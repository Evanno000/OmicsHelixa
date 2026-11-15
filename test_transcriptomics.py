
from genomics.transcriptomics_analysis import (
    analyze_rna_sequence,
    convert_rna_to_dna,
    translate_rna
)


rna_sequence = "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGA"


print("================================")
print("OMICSBRIDGE TRANSCRIPTOMICS TEST")
print("================================")


print("\nRNA SEQUENCE ANALYSIS")

results = analyze_rna_sequence(rna_sequence)

for key, value in results.items():

    print(f"{key}: {value}")


print("\nRNA TO DNA CONVERSION")

print(convert_rna_to_dna(rna_sequence))


print("\nRNA TRANSLATION")

print(translate_rna(rna_sequence))


print("\nTEST COMPLETED")
