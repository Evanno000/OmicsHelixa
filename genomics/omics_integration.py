
from Bio.Seq import Seq


# ==========================================
# DNA VALIDATION
# ==========================================

def clean_coding_dna(sequence):

    if not sequence:
        raise ValueError("Please provide a DNA sequence.")

    sequence = "".join(sequence.split()).upper()

    allowed_bases = set("ATGCN")

    invalid_bases = set(sequence) - allowed_bases

    if invalid_bases:
        raise ValueError(
            "Invalid DNA characters: "
            + ", ".join(sorted(invalid_bases))
            + ". Use only A, T, G, C and N."
        )

    return sequence


# ==========================================
# CENTRAL DOGMA ANALYSIS
# ==========================================

def analyze_central_dogma(sequence):

    dna = clean_coding_dna(sequence)

    # --------------------------------------
    # 1. DNA COMPOSITION
    # --------------------------------------

    dna_length = len(dna)

    a_count = dna.count("A")
    t_count = dna.count("T")
    g_count = dna.count("G")
    c_count = dna.count("C")
    n_count = dna.count("N")

    known_bases = dna_length - n_count

    gc_content = (
        ((g_count + c_count) / known_bases) * 100
        if known_bases > 0
        else 0
    )

    # --------------------------------------
    # 2. TRANSCRIPTION
    # --------------------------------------

    # Input is the DNA coding strand.
    # Therefore, transcription replaces T with U.

    mrna = dna.replace("T", "U")

    # --------------------------------------
    # 3. TRANSLATION
    # --------------------------------------

    complete_length = len(dna) - (len(dna) % 3)

    complete_coding_sequence = dna[:complete_length]

    if complete_length >= 3:

        protein = str(
            Seq(complete_coding_sequence).translate(
                to_stop=False
            )
        )

    else:
        protein = ""

    # --------------------------------------
    # 4. PROTEIN COMPOSITION
    # --------------------------------------

    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    amino_acid_counts = {
        amino_acid: protein.count(amino_acid)
        for amino_acid in amino_acids
    }

    amino_acid_counts["X"] = protein.count("X")

    stop_count = protein.count("*")

    # --------------------------------------
    # 5. SUMMARY
    # --------------------------------------

    results = {

        "dna_sequence": dna,
        "dna_length": dna_length,

        "dna_A": a_count,
        "dna_T": t_count,
        "dna_G": g_count,
        "dna_C": c_count,
        "dna_N": n_count,

        "gc_content": round(gc_content, 2),

        "mrna_sequence": mrna,
        "mrna_length": len(mrna),

        "protein_sequence": protein,
        "protein_length": len(protein),

        "amino_acid_counts": amino_acid_counts,

        "stop_codons": stop_count,

        "incomplete_bases": len(dna) - complete_length,

        "complete_codons": complete_length // 3
    }

    return results
