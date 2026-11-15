
from Bio.Seq import Seq


# ==========================================
# RNA SEQUENCE VALIDATION
# ==========================================

def clean_rna_sequence(sequence):

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError("The RNA sequence is empty.")

    allowed_bases = set("AUGCN")

    if not set(sequence).issubset(allowed_bases):
        raise ValueError(
            "Invalid RNA sequence. Use only A, U, G, C and N."
        )

    return sequence


# ==========================================
# RNA SEQUENCE ANALYSIS
# ==========================================

def analyze_rna_sequence(sequence):

    sequence = clean_rna_sequence(sequence)

    a_count = sequence.count("A")
    u_count = sequence.count("U")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    n_count = sequence.count("N")

    gc_bases = g_count + c_count
    known_bases = len(sequence) - n_count

    gc_content = (
        (gc_bases / known_bases) * 100
        if known_bases > 0 else 0
    )

    # RNA base-pairing:
    # A pairs with U
    # U pairs with A
    # G pairs with C
    # C pairs with G
    # N remains N

    complement_table = str.maketrans(
        "AUGCN",
        "UACGN"
    )

    complementary_sequence = sequence.translate(
        complement_table
    )

    reverse_complement = complementary_sequence[::-1]

    return {
        "sequence_length": len(sequence),
        "A": a_count,
        "U": u_count,
        "G": g_count,
        "C": c_count,
        "N": n_count,
        "GC_content": round(gc_content, 2),
        "complementary_sequence": complementary_sequence,
        "reverse_complement": reverse_complement
    }


# ==========================================
# RNA TO DNA CONVERSION
# ==========================================

def convert_rna_to_dna(sequence):

    sequence = clean_rna_sequence(sequence)

    return str(Seq(sequence).back_transcribe())


# ==========================================
# RNA TRANSLATION
# ==========================================

def translate_rna(sequence):

    sequence = clean_rna_sequence(sequence)

    if len(sequence) % 3 != 0:
        raise ValueError(
            "RNA sequence length must be divisible by 3 "
            "for complete codon translation."
        )

    return str(Seq(sequence).translate())
