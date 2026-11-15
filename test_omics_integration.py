
from genomics.omics_integration import analyze_central_dogma


print("===================================")
print("OMICSBRIDGE INTEGRATION TEST")
print("===================================")

dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA"

results = analyze_central_dogma(dna)

print("\nDNA SEQUENCE:")
print(results["dna_sequence"])

print("\nDNA LENGTH:")
print(results["dna_length"])

print("\nGC CONTENT:")
print(results["gc_content"])

print("\nTRANSCRIBED mRNA:")
print(results["mrna_sequence"])

print("\nTRANSLATED PROTEIN:")
print(results["protein_sequence"])

print("\nAMINO ACID COMPOSITION:")
print(results["amino_acid_counts"])

print("\nSTOP CODONS:")
print(results["stop_codons"])

print("\n===================================")
print("TEST COMPLETED")
print("===================================")
