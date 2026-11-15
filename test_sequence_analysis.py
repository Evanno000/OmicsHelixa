
from genomics.sequence_analysis import analyze_sequence

# Test DNA sequence
dna = "ATGCGCAT"

# Analyse sequence
results = analyze_sequence(dna)

print("OMICSBRIDGE DNA ANALYSIS")
print("------------------------")

for key, value in results.items():
    print(f"{key}: {value}")
