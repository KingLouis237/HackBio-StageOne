# Clean, beginner-friendly DNA→protein translator with extras
from typing import Dict, List, Tuple

# 1) Standard codon table (RNA codons: U instead of T)
CODON_TABLE: Dict[str, str] = {
    # Phenylalanine/Leucine
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    # Serine
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    # Tyrosine / Stop
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    # Cysteine / Stop / Tryptophan
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    # Leucine
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    # Proline
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    # Histidine / Glutamine
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Arginine
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    # Isoleucine / Methionine (START)
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    # Threonine
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    # Asparagine / Lysine
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    # Serine / Arginine
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    # Valine
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    # Alanine
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    # Aspartic acid / Glutamic acid
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    # Glycine
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}

def clean_sequence(seq: str, alphabet: str = "DNA") -> str:
    """
    Keep only valid letters for DNA or RNA and uppercase the sequence.
    - DNA valid: A C G T
    - RNA valid: A C G U
    """
    seq = seq.upper()
    valid = set("ACGT") if alphabet == "DNA" else set("ACGU")
    return "".join([c for c in seq if c in valid])


def dna_to_rna(dna: str) -> str:
    """Transcribe DNA (5' to 3') to mRNA by replacing T with U."""
    return dna.upper().replace("T", "U")


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA string."""
    comp = str.maketrans({"A":"T","T":"A","C":"G","G":"C"})
    return dna.upper().translate(comp)[::-1]


def translate_rna(rna: str, start_at_first_aug: bool = True) -> str:
    """
    Translate an RNA sequence to amino acids.
    - If start_at_first_aug is True, translation begins at the first AUG.
    - Stops at stop codon (*) or end of sequence.
    """
    rna = clean_sequence(rna, alphabet="RNA")
    if not rna:
        return ""
    # Optionally find first AUG
    start = 0
    if start_at_first_aug:
        start = rna.find("AUG")
        if start == -1:
            return ""  # no start codon
    protein = []
    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i+3]
        aa = CODON_TABLE.get(codon, "?")
        if aa == "*":  # stop
            break
        if len(codon) == 3:
            protein.append(aa)
    return "".join(protein)


def translate_dna(dna: str, start_at_first_atg: bool = True) -> str:
    """
    Convenience wrapper: DNA to mRNA to protein.
    """
    dna = clean_sequence(dna, alphabet="DNA")
    rna = dna_to_rna(dna)
    return translate_rna(rna, start_at_first_aug=start_at_first_atg)


def translate_in_frame(dna: str, frame: int = 0, strand: str = "+") -> str:
    """
    Translate a specific reading frame without enforcing a start codon.
    - frame: 0, 1, or 2
    - strand: '+' or '-' (reverse-complement)
    """
    dna = clean_sequence(dna, alphabet="DNA")
    if strand == "-":
        dna = reverse_complement(dna)
    dna = dna[frame:]
    rna = dna_to_rna(dna)
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        aa = CODON_TABLE.get(codon, "?")
        protein.append(aa if aa != "*" else "*")
    return "".join(protein)


def find_orfs(dna: str, min_len: int = 0) -> List[Tuple[str, int, int, str]]:
    """
    Find simple ORFs on both strands. Returns a list of tuples:
    (protein, start_nt, end_nt, strand)
    - start_nt/end_nt are 0-based coordinates on the + strand reference.
    - min_len filters proteins shorter than min_len aa.
    """
    dna = clean_sequence(dna, alphabet="DNA")
    results: List[Tuple[str, int, int, str]] = []
    for strand in ["+", "-"]:
        seq = dna if strand == "+" else reverse_complement(dna)
        rna = dna_to_rna(seq)
        i = 0
        while i <= len(rna) - 3:
            codon = rna[i:i+3]
            if codon == "AUG":
                # extend until stop
                j = i
                peptides = []
                while j <= len(rna) - 3:
                    c = rna[j:j+3]
                    aa = CODON_TABLE.get(c, "?")
                    if aa == "*" or aa is None:
                        j += 3
                        break
                    peptides.append(aa)
                    j += 3
                protein = "".join(peptides)
                if len(protein) >= min_len:
                    # Map back to + strand coordinates
                    if strand == "+":
                        start_nt = i
                        end_nt = j
                    else:
                        # For reverse strand, map to + coordinates
                        start_nt = len(dna) - j
                        end_nt = len(dna) - i
                    results.append((protein, start_nt, end_nt, strand))
                i = j  # skip past this ORF
            else:
                i += 3 if codon in CODON_TABLE else 1
    return results


# DEMOS

# Classic test sequence (should begin with 'M')
dna_test = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
print("translate_dna():", translate_dna(dna_test))

# Show all six-frame translations (simple, without start/stop handling)
print("\nSix-frame translations (frames 0/1/2 on + and -):")
for strand in ["+", "-"]:
    for frame in [0, 1, 2]:
        prot = translate_in_frame(dna_test, frame=frame, strand=strand)
        print(f"strand {strand}, frame {frame}: {prot[:60]}")

# Find ORFs of length >= 5 aa
print("\nORFs (>=5 aa) on both strands:")
for prot, s, e, st in find_orfs(dna_test, min_len=5):
    print(f"{st}-strand {s}-{e} -> {prot}")




#HAMMING DISTANCE
def hamming_distance(str1, str2):
    """
    hERE I Calculate the Hamming distance between two strings.
    If they are not the same length, pad the shorter one with '_'.
    """
    # Step 1: make both lowercase so case doesn’t affect comparison
    str1 = str1.lower()
    str2 = str2.lower()

    # Step 2: find their lengths
    len1 = len(str1)
    len2 = len(str2)

    # Step 3: if lengths differ, pad the shorter one
    if len1 < len2:
        str1 = str1 + "_" * (len2 - len1)
    elif len2 < len1:
        str2 = str2 + "_" * (len1 - len2)

    # Step 4: set a counter to count the differences
    distance = 0

    # Step 5: loop through both strings at once
    for a, b in zip(str1, str2):
        # zip() pairs letters together, like:
        # "abc" and "xyz" to ('a','x'), ('b','y'), ('c','z')
        if a != b:
            distance += 1  # add 1 if letters are different

    # Step 6: return the final count
    return distance


slack_name = "Ashu Louis"
twitter_handle = "KingLouis237"

dist = hamming_distance(slack_name, twitter_handle)
print(f"Hamming distance between '{slack_name}' and '{twitter_handle}' is:", dist)

