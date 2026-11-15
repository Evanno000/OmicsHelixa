
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

    allowed_bases = set("ATGCN")

    if not sequence:
        raise ValueError("The DNA sequence is empty.")

    if not set(sequence).issubset(allowed_bases):
        raise ValueError(
            "Invalid DNA sequence. Use only A, T, G, C and N."
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
    # 4. GC CONTENT
    # --------------------------------------

    gc_content = gc_fraction(dna, ambiguous="remove") * 100

    # --------------------------------------
    # 5. REVERSE COMPLEMENT
    # --------------------------------------

    reverse_complement = str(dna.reverse_complement())

    # --------------------------------------
    # 6. RETURN ANALYSIS RESULTS
    # --------------------------------------

    results = {

        "sequence_length": len(dna),

        "A": a_count,
        "T": t_count,
        "G": g_count,
        "C": c_count,
        "N": n_count,

        "GC_content": round(gc_content, 2),

        "reverse_complement": reverse_complement

    }

    return results
