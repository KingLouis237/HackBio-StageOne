from typing import Dict, List, Tuple, Literal

# --- Genetic code (RNA codons: U not T). Using ALL CAPS to signal constants. ---
CODON_TABLE: Dict[str, str] = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}

ValidDNAPolicy = Literal["strict", "ignore", "replace"]
NoStartPolicy  = Literal["return_empty", "raise"]
NoStopPolicy   = Literal["truncate", "return_full", "raise"]

def _validate_dna(seq: str,
                  policy: ValidDNAPolicy = "strict",
                  replace_with: str = "A") -> str:
    """
    Validate/normalize DNA: uppercase & handle invalid symbols per policy.
    - strict  : raise ValueError if anything not in A/C/G/T
    - ignore  : drop invalid characters silently
    - replace : replace invalid characters with `replace_with` (default 'A')
    """
    seq = seq.upper()
    valid = set("ACGT")
    out = []
    for i, c in enumerate(seq):
        if c in valid:
            out.append(c)
        else:
            if policy == "strict":
                raise ValueError(f"Invalid DNA character '{c}' at pos {i}. Allowed: A/C/G/T.")
            elif policy == "ignore":
                continue
            elif policy == "replace":
                out.append(replace_with)
    return "".join(out)

def dna_to_rna(dna: str) -> str:
    """Transcribe 5'→3' DNA coding strand to RNA (T→U)."""
    return dna.replace("T", "U")

def reverse_complement(dna: str) -> str:
    """Reverse-complement DNA."""
    comp = str.maketrans({"A":"T","T":"A","C":"G","G":"C"})
    return dna.translate(comp)[::-1]

def translate_rna(
    rna: str,
    start_at_first_aug: bool = True,
    no_start_policy: NoStartPolicy = "return_empty",
    no_stop_policy: NoStopPolicy = "truncate",
) -> str:
    """
    RNA → protein.
    - start_at_first_aug=True: begin at first AUG (M). Otherwise start at index 0.
    - no_start_policy: what to do if no AUG is found.
        'return_empty' | 'raise'
    - no_stop_policy: when no stop codon is encountered:
        'truncate' (translate to end), 'return_full' (same as truncate), or 'raise'
    """
    if not rna:
        return ""

    start = 0
    if start_at_first_aug:
        start = rna.find("AUG")
        if start == -1:
            if no_start_policy == "return_empty":
                return ""
            raise ValueError("No AUG start codon found in RNA.")

    protein: List[str] = []
    i = start
    while i <= len(rna) - 3:
        codon = rna[i:i+3]
        aa = CODON_TABLE.get(codon, "?")  # '?' flags unknown triplets
        if aa == "*":  # stop codon
            break
        protein.append(aa)
        i += 3

    # Handle lack of stop codon based on policy
    if i > len(rna) - 3:
        if no_stop_policy == "raise":
            raise ValueError("No stop codon encountered before RNA end.")
        # 'truncate' / 'return_full' -> just return what we have

    return "".join(protein)

def translate_dna(
    dna: str,
    *,
    start_at_first_atg: bool = True,
    invalid_policy: ValidDNAPolicy = "strict",
    no_start_policy: NoStartPolicy = "return_empty",
    no_stop_policy: NoStopPolicy = "truncate",
) -> str:
    """
    DNA → (validate/normalize) → RNA → protein.
    - invalid_policy: how to treat non-ACGT symbols in DNA
    - start_at_first_atg: True = enforce canonical start (ATG/AUG)
    """
    dna = _validate_dna(dna, policy=invalid_policy)
    rna = dna_to_rna(dna)
    return translate_rna(
        rna,
        start_at_first_aug=start_at_first_atg,
        no_start_policy=no_start_policy,
        no_stop_policy=no_stop_policy,
    )




from typing import Literal

UnequalPolicy = Literal["raise", "pad"]

def hamming_distance(a: str, b: str, *, policy: UnequalPolicy = "raise",
                     pad_char: str = "_", case_insensitive: bool = True) -> int:
    """
    Hamming distance between two strings.
    - If lengths differ:
        - policy='raise' -> raise ValueError (as grader suggested)
        - policy='pad'   -> pad the shorter with `pad_char`
    - case_insensitive=True -> compare lowercased versions
    """
    if case_insensitive:
        a, b = a.lower(), b.lower()

    if len(a) != len(b):
        if policy == "raise":
            raise ValueError(f"Strings must have equal length (got {len(a)} vs {len(b)}).")
        elif policy == "pad":
            if len(a) < len(b):
                a = a + pad_char * (len(b) - len(a))
            else:
                b = b + pad_char * (len(a) - len(b))

    return sum(1 for x, y in zip(a, b) if x != y)




slack = "Ashu Louis"
x_handle = "KingLouis237"
try:
    print("Hamming (strict):", hamming_distance(slack, x_handle, policy="raise"))
except ValueError as e:
    print("Strict mode error:", e)

print("Hamming (padded):", hamming_distance(slack, x_handle, policy="pad"))
