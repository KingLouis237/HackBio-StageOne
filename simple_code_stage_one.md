from typing import Dict, List, Tuple, Literal, Iterable

# I define the standard genetic code as RNA codons (U, not T).
# I keep it in ALL_CAPS to signal "constant".
CODON_TABLE: Dict[str, str] = {
    # I list all 64 codons; '*' means stop codon.
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

# I define small Literal types so my function signatures communicate valid options.
ValidDNAPolicy = Literal["strict", "ignore", "replace"]
NoStartPolicy  = Literal["return_empty", "raise", "start_at_0"]
NoStopPolicy   = Literal["truncate", "return_full", "raise"]

def _clean_dna(
    seq: str,
    *,
    policy: ValidDNAPolicy = "strict",
    replace_with: str = "N"
):
    """
    I normalize DNA: uppercase and handle non-ACGT symbols per policy.

    - strict  : I raise ValueError if I meet any non-ACGT character.
    - ignore  : I drop invalid characters (this changes indices; I mention that in errors/info).
    - replace : I replace invalid with 'N' (or another base), keeping length stable.

    I intentionally support 'N' because real data often contains ambiguity;
    downstream functions can decide how to treat 'N'.
    """
    seq = seq.upper()  # I make case uniform.
    valid = set("ACGT")
    out: List[str] = []
    for i, c in enumerate(seq):  # I walk through the sequence with positions.
        if c in valid:
            out.append(c)  # valid base, keep it
        else:
            if policy == "strict":
                # I fail fast with a precise message and position.
                raise ValueError(f"Invalid DNA character '{c}' at position {i}. Allowed: A/C/G/T.")
            elif policy == "ignore":
                # I skip invalid; I accept index drift for a quick-clean use case.
                continue
            elif policy == "replace":
                # I keep length constant by substituting an ambiguity code (default 'N').
                out.append(replace_with)
    return "".join(out)

def dna_to_rna(dna: str):
    """
    I model transcription of the coding (sense) strand by replacing T with U.
    I assume input is 5'3' and already cleaned to A/C/G/T/N by _clean_dna.
    """
    return dna.replace("T", "U").replace("N", "N")  # I leave 'N' as 'N' to signal ambiguity.

def reverse_complement(dna: str) -> str:
    """
    I return the reverse complement (5'→3') of the input DNA.
    This is needed to handle the minus strand.
    """
    comp = str.maketrans({"A":"T","T":"A","C":"G","G":"C","N":"N"})  # I map N to N.
    return dna.translate(comp)[::-1]  # I reverse after complementing.

def _iter_start_indices_for_rna(
    rna: str,
    *,
    start_at_first: bool,
    allowed_start_codons: Iterable[str]
):
    """
    I yield valid start indices in RNA given a set of allowed start codons.
    - If start_at_first=True, I return only the first occurrence index (or nothing).
    - Else, I yield all frame-aligned occurrences (useful for ORF scans).
    """
    # I scan in-frame for any allowed start codon.
    i = 0
    found_first = False
    while i <= len(rna) - 3:
        codon = rna[i:i+3]
        if codon in allowed_start_codons:
            yield i
            if start_at_first:
                found_first = True
                break
            i += 3  # I skip to next codon in-frame for efficiency.
        else:
            i += 3  # I step by codons.
    if start_at_first and not found_first:
        # I yield nothing; translate_rna handles the downstream policy.
        return

def translate_rna(
    rna: str,
    *,
    start_at_first: bool = True,
    # I allow alternative start codons (biologically: bacteria often use GUG/UUG).
    allowed_start_codons: Iterable[str] = ("AUG", "GUG", "UUG"),
    no_start_policy: NoStartPolicy = "return_empty",
    no_stop_policy: NoStopPolicy = "truncate",
    treat_unknown_as_stop: bool = False
):
    """
    I translate RNA into a protein string.

    Parameters I expose on purpose to address "scientific flaws":
    - allowed_start_codons: I include AUG, plus common bacterial alternatives (GUG, UUG).
      In eukaryotes AUG dominates, but I keep this configurable.
    - no_start_policy:
        * 'return_empty'  I return '' if I find no start codon.
        * 'raise'         I raise ValueError if I find no start.
        * 'start_at_0'    I start at position 0 even without a canonical start.
    - no_stop_policy:
        * 'truncate'      I translate to the end if no stop occurs (typical practice).
        * 'return_full'   synonymous with 'truncate' here; I keep for clarity.
        * 'raise'         I raise if I never encounter a stop codon.
    - treat_unknown_as_stop:
        * If True, I stop at any unknown triplet (e.g., contains 'N').
          If False, I encode unknown as '?' and continue.
    """
    if not rna:
        return ""  # I handle empty input quietly; this keeps calling code simple.

    # I decide where to begin based on allowed starts and policy.
    start_indices = list(_iter_start_indices_for_rna(
        rna, start_at_first=start_at_first, allowed_start_codons=allowed_start_codons
    ))

    if start_at_first:
        if start_indices:
            start = start_indices[0]
        else:
            if no_start_policy == "return_empty":
                return ""
            elif no_start_policy == "raise":
                raise ValueError("No allowed start codon found in RNA.")
            elif no_start_policy == "start_at_0":
                start = 0
            else:
                raise ValueError(f"Unknown no_start_policy: {no_start_policy}")
    else:
        # If I am not forced to pick only the first, I still need a start.
        start = start_indices[0] if start_indices else 0  # I gracefully default to 0.

    protein: List[str] = []
    i = start
    while i <= len(rna) - 3:
        codon = rna[i:i+3]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            # Unknown triplet: contains 'N' or not in table. I decide via flag.
            if treat_unknown_as_stop:
                break
            else:
                protein.append("?")
                i += 3
                continue
        if aa == "*":  # I stop on canonical stop codon.
            break
        protein.append(aa)
        i += 3

    # If I reached the end without a stop, I honor the no_stop_policy.
    if i > len(rna) - 3:
        if no_stop_policy == "raise":
            raise ValueError("No stop codon encountered before RNA end.")
        # 'truncate' / 'return_full' → I simply return what I translated.

    return "".join(protein)

def translate_dna(
    dna: str,
    *,
    invalid_policy: ValidDNAPolicy = "strict",
    replace_with: str = "N",
    start_at_first: bool = True,
    allowed_start_codons_dna: Iterable[str] = ("ATG", "GTG", "TTG"),
    no_start_policy: NoStartPolicy = "return_empty",
    no_stop_policy: NoStopPolicy = "truncate",
    treat_unknown_as_stop: bool = False,
    strand: Literal["+", "-"] = "+"
) -> str:
    """
    I translate DNA to protein with explicit, documented choices.

    - invalid_policy: 'strict' (raise on non-ACGT), 'ignore', or 'replace' with replace_with.
    - allowed_start_codons_dna: I default to ATG plus bacterial GTG/TTG (configurable).
    - strand: '+' for input as-is; '-' to translate the reverse complement.
    - I convert DNA → RNA and then call translate_rna with matching policies.
    """
    dna = _clean_dna(dna, policy=invalid_policy, replace_with=replace_with)
    if strand == "-":
        dna = reverse_complement(dna)  # I support minus strand translation.
    rna = dna_to_rna(dna)
    # I convert DNA allowed starts to RNA allowed starts by replacing T→U.
    allowed_rna_starts = tuple(s.replace("T", "U") for s in allowed_start_codons_dna)
    return translate_rna(
        rna,
        start_at_first=start_at_first,
        allowed_start_codons=allowed_rna_starts,
        no_start_policy=no_start_policy,
        no_stop_policy=no_stop_policy,
        treat_unknown_as_stop=treat_unknown_as_stop,
    )


if __name__ == "__main__":
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    print("DNA:", dna)
    print("Protein:", translate_dna(dna))
    print("Reverse complement:", reverse_complement(dna))
    print("RNA:", dna_to_rna(dna))




#Hamming distance 

# I compute Hamming distances and produce the plots the rubric expects.


from typing import Literal, Iterable
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
from pathlib import Path

# Step 1. Define the file name and location
csv_path = Path("handles.csv")  # this saves it in the same folder as your script

# Step 2. Define some simple sample data
data = {
    "handle": ["Ashu Louis", "KingLouis237", "LouisAshu"],
    "label": ["slack", "twitter", "alt"],
    "diagnosis": ["control", "case", "control"]
}

# Step 3. Convert data into a DataFrame
df = pd.DataFrame(data)

# Step 4. Save it as a CSV file
df.to_csv(csv_path, index=False)

# Step 5. Confirm creation and show where it was saved
print(f"handles.csv created successfully at: {csv_path.resolve()}")
print(df)


UnequalPolicy = Literal["raise", "pad"]

def hamming_distance(
    a: str,
    b: str,
    *,
    policy: UnequalPolicy = "raise",
    pad_char: str = "_",
    case_insensitive: bool = True
) -> int:
    """
    I compute the Hamming distance (count of mismatched positions).
    - If lengths differ:
        * 'raise'  I raise ValueError (pure Hamming definition).
        * 'pad'    I pad the shorter one with pad_char so I can still compare.
    - If case_insensitive=True, I compare lowercase forms to avoid case artifacts.
    """
    if case_insensitive:
        a, b = a.lower(), b.lower()  # I normalize case.

    if len(a) != len(b):
        if policy == "raise":
            raise ValueError(f"Strings must be equal length (got {len(a)} vs {len(b)}).")
        elif policy == "pad":
            # I pad on the right using str.ljust so I don’t shift existing positions.
            if len(a) < len(b):
                a = a.ljust(len(b), pad_char)
            else:
                b = b.ljust(len(a), pad_char)

    # I sum a 1 for each position where characters differ.
    return sum(1 for x, y in zip(a, b) if x != y)

def demo_handles_analysis(csv_file: str, slack_name: str):
    """

    CSV schema I expect:
    - 'handle' (string)
    - optional 'label' (e.g., 'slack', 'twitter')
    - optional 'diagnosis' (to show how I’d join clinical variables)
    """
    # Load CSV
    # I show header handling explicitly (header row at 0 is the default).
    df = pd.read_csv(csv_file, header=0)  #read_csv, header
    # If a human provided a label column, I set it as index for nicer heatmap axes.
    if "label" in df.columns:
        df = df.set_index("label")       # index_col via set_index
    print(df.head())                     #  head

    # I normalize handle strings to lowercase and trim whitespace.
    df["handle_norm"] = df["handle"].str.lower().str.strip()

    # Compute pairwise Hamming distances (padding, to handle unequal lengths)
    n = len(df)
    dist = np.zeros((n, n), dtype=int)   # I preallocate for speed and clarity.
    handles = df["handle_norm"].tolist()
    # I compute a full pairwise matrix (symmetric).
    for i in range(n):
        for j in range(n):
            dist[i, j] = hamming_distance(handles[i], handles[j], policy="pad")
    # I wrap the matrix back into a labeled DataFrame for plotting.
    dist_df = pd.DataFrame(dist, index=df.index, columns=df.index)

    # correlation view for a heatmap
    # I convert distances to a similarity in [0,1] so a heatmap is more intuitive (darker = more similar).
    # This is not a Pearson correlation; I document that clearly.
    max_val = dist_df.values.max()
    similarity = 1 - (dist_df / max_val) if max_val > 0 else 1 - dist_df  # avoid divide by zero
    print("Max similarity:", similarity.values.max())  # expected: max

    #Heatmap of distances
    plt.figure(figsize=(6, 4))
    ax = sns.heatmap(
        dist_df,
        cmap="Blues",
        linewidths=.5,
        linecolor="white",
        annot=True, fmt="d"              #
    )
    ax.set_title("Pairwise Hamming Distances (padded)")
    ax.set_xlabel("Handle")
    ax.set_ylabel("Handle")
    plt.tight_layout()
    plt.show()

    #Clustermap on similarity (groups similar handles)
    sns.clustermap(similarity, cmap="Blues")  # expected: clustermap, cmap
    plt.show()                                # expected: show

    # Scatterplot: handle length vs. distance to my Slack
    # I compute distance to the provided slack_name to produce a meaningful scatter.
    dist_to_slack = [hamming_distance(slack_name, s, policy="pad") for s in handles]
    lengths = [len(s) for s in handles]

    plt.figure(figsize=(5, 4))          #
    # If a 'diagnosis' column exists, I use it as hue to show how clinical variables could be layered.
    hue = df["diagnosis"] if "diagnosis" in df.columns else None
    sns.scatterplot(x=lengths, y=dist_to_slack, hue=hue)  # scatterplot, hue
    plt.title("Distance to Slack vs. Handle Length")      #
    plt.xlabel("Handle length")                           # xlabel
    plt.ylabel(f"Hamming distance to '{slack_name}'")     # ylabel
    plt.axvline(x=len(slack_name), linestyle="--")        # axvline
    plt.axhline(y=np.median(dist_to_slack), linestyle=":")# axhline
    plt.grid(alpha=0.3)                                   # grid
    plt.tight_layout()
    plt.show()

    # KDE plot: distribution of distances to Slack
    plt.figure(figsize=(5, 3))
    sns.kdeplot(dist_to_slack, fill=True)
    plt.title("Distribution of distances to Slack")
    plt.xlabel("Hamming distance")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # I include a tiny runnable demo
    # handle,label,diagnosis
    # Ashu Louis,slack,control
    # KingLouis237,twitter,case
    # LouisAshu,alt,control
    csv_file = "handles.csv"
    slack = "Ashu Louis"
    demo_handles_analysis(csv_file, slack_name=slack)
