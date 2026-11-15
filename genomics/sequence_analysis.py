
# ==========================================
# OMICSBRIDGE
# DNA SEQUENCE ANALYSIS ENGINE
# ==========================================

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def analyze_sequence(sequence):

    # --------------------------------------
    # 1. CLEAN DNA SEQUENCE
    # --------------------------------------

    sequence = "".join(sequence.split()).upper()

    # --------------------------------------
    # 2. VALIDATE DNA SEQUENCE
    # --------------------------------------

    # Full standard DNA IUPAC alphabet
    allowed_bases = set("ATGCRYSWKMBDHVN")

    if not sequence:
        raise ValueError("The DNA sequence is empty.")

    invalid_bases = set(sequence) - allowed_bases

    if invalid_bases:
        invalid_list = ", ".join(sorted(invalid_bases))

        raise ValueError(
            f"Invalid DNA characters detected: {invalid_list}. "
            "Use standard DNA IUPAC codes: "
            "A, T, G, C, R, Y, S, W, K, M, B, D, H, V and N."
        )

    # Create Biopython sequence object
    dna = Seq(sequence)

    # --------------------------------------
    # 3. NUCLEOTIDE COMPOSITION
    # --------------------------------------

    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    n_count = sequence.count("N")

    # --------------------------------------
    # 4. AMBIGUOUS NUCLEOTIDE COMPOSITION
    # --------------------------------------

    ambiguity_codes = {
        "R": "A or G",
        "Y": "C or T",
        "S": "G or C",
        "W": "A or T",
        "K": "G or T",
        "M": "A or C",
        "B": "C, G or T",
        "D": "A, G or T",
        "H": "A, C or T",
        "V": "A, C or G",
        "N": "Any nucleotide"
    }

    ambiguous_bases = {
        code: sequence.count(code)
        for code in ambiguity_codes
    }

    ambiguous_count = sum(ambiguous_bases.values())

    # --------------------------------------
    # 5. GC CONTENT
    # --------------------------------------

    # Ambiguous positions are excluded from
    # the GC calculation.
    gc_content = gc_fraction(
        dna,
        ambiguous="remove"
    ) * 100

    # --------------------------------------
    # 6. REVERSE COMPLEMENT
    # --------------------------------------

    reverse_complement = str(
        dna.reverse_complement()
    )

    # --------------------------------------
    # 7. RETURN ANALYSIS RESULTS
    # --------------------------------------

    results = {

        "sequence_length": len(dna),

        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,
        "N": n_count,

        "ambiguous_count": ambiguous_count,

        "ambiguous_bases": ambiguous_bases,

        "GC_content": round(gc_content, 2),

        "reverse_complement": reverse_complement

    }

    return results
