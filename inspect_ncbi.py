
from Bio import Entrez, SeqIO
from collections import Counter

Entrez.email = "your-email@example.com"

accession = "CP002686.1"

print("OMICSBRIDGE - NCBI ANNOTATION SUMMARY")
print("=" * 45)

handle = Entrez.efetch(
    db="nuccore",
    id=accession,
    rettype="gb",
    retmode="text"
)

record = SeqIO.read(handle, "genbank")
handle.close()

print("Accession:", record.id)
print("Organism:", record.annotations.get("organism"))
print("Sequence length:", len(record.seq))

feature_counts = Counter(
    feature.type for feature in record.features
)

print("\nFEATURE COUNTS")
for feature_type, count in sorted(feature_counts.items()):
    print(f"{feature_type}: {count}")

cds_features = [
    feature for feature in record.features
    if feature.type == "CDS"
]

mrna_features = [
    feature for feature in record.features
    if feature.type == "mRNA"
]

print("\nCDS SUMMARY")
print("Total CDS features:", len(cds_features))

with_translation = sum(
    1 for feature in cds_features
    if "translation" in feature.qualifiers
)

print("CDS with translations:", with_translation)

print("\nmRNA SUMMARY")
print("Total mRNA features:", len(mrna_features))

print("\nFIRST 5 CDS EXAMPLES")
for feature in cds_features[:5]:
    qualifiers = feature.qualifiers

    print(
        "Gene:",
        qualifiers.get("gene", qualifiers.get("locus_tag", ["Unknown"]))[0]
    )
    print(
        "Product:",
        qualifiers.get("product", ["Unknown"])[0]
    )
    print(
        "Has translation:",
        "translation" in qualifiers
    )
    print("-" * 30)

print("\nInspection completed.")

