
from Bio.Seq import Seq


# ==========================================
# DNA VALIDATION
# ==========================================

def clean_coding_dna(sequence):

    if not sequence:
        raise ValueError("Please provide a DNA sequence.")

    sequence = "".join(sequence.split()).upper()

    # Full standard DNA IUPAC alphabet
    allowed_bases = set("ATGCRYSWKMBDHVN")

    invalid_bases = set(sequence) - allowed_bases

    if invalid_bases:
        raise ValueError(
            "Invalid DNA characters: "
            + ", ".join(sorted(invalid_bases))
            + ". Use standard DNA IUPAC codes: "
            "A, T, G, C, R, Y, S, W, K, M, B, D, H, V and N."
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

    # Count all ambiguous nucleotide codes.
    ambiguity_codes = "RYSWKMBDHVN"

    ambiguous_bases = {

        code: dna.count(code)

        for code in ambiguity_codes

    }

    ambiguous_count = sum(
        ambiguous_bases.values()
    )

    # Count positions with definite nucleotide identities.
    known_bases = (
        a_count
        + t_count
        + g_count
        + c_count
    )

    # GC content excludes ambiguous positions.
    gc_content = (

        ((g_count + c_count) / known_bases) * 100

        if known_bases > 0

        else 0

    )

    # --------------------------------------
    # 2. TRANSCRIPTION
    # --------------------------------------

    # Input represents the DNA coding strand.
    # Therefore, transcription replaces T with U.
    # IUPAC ambiguity codes are preserved.

    mrna = dna.replace("T", "U")

    # --------------------------------------
    # 3. TRANSLATION
    # --------------------------------------

    complete_length = len(dna) - (
        len(dna) % 3
    )

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

    # X represents an unknown or ambiguous amino acid.
    amino_acid_counts["X"] = protein.count("X")

    # Stop codons are not amino acid residues.
    stop_count = protein.count("*")

    # Number of translated amino acid residues,
    # excluding stop codon positions.
    amino_acid_residue_count = (
        len(protein) - stop_count
    )

    # --------------------------------------
    # 5. SUMMARY
    # --------------------------------------

    results = {

        # DNA information

        "dna_sequence": dna,

        "dna_length": dna_length,

        "dna_A": a_count,
        "dna_T": t_count,
        "dna_G": g_count,
        "dna_C": c_count,
        "dna_N": n_count,

        "ambiguous_count": ambiguous_count,

        "ambiguous_bases": ambiguous_bases,

        "known_bases": known_bases,

        "gc_content": round(gc_content, 2),

        # RNA information

        "mrna_sequence": mrna,

        "mrna_length": len(mrna),

        # Protein information

        "protein_sequence": protein,

        # Retained for compatibility with existing templates.
        # This includes stop codon symbols.
        "protein_length": len(protein),

        "amino_acid_residue_count": amino_acid_residue_count,

        "amino_acid_counts": amino_acid_counts,

        "stop_codons": stop_count,

        # Translation information

        "incomplete_bases": dna_length - complete_length,

        "complete_codons": complete_length // 3

    }

    return results
