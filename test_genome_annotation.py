
from genomics.ncbi_annotation import annotate_genome

try:

    results = annotate_genome("CP002686.1")

    print("\nGENOME ANNOTATION TEST")
    print("=" * 45)

    print("Accession:", results["accession"])
    print("Organism:", results["organism"])
    print("Record length:", results["record_length"])

    print("Gene features:", results["gene_feature_count"])
    print("mRNA features:", results["mrna_feature_count"])
    print("Annotated CDS:", results["cds_count"])
    print("Unique gene identifiers:", results["unique_gene_count"])

    print("\nFIRST FIVE CDS RECORDS")
    print("=" * 45)

    for cds in results["cds_records"][:5]:

        print("\nGene:", cds["gene"])
        print("Locus tag:", cds["locus_tag"])
        print("Product:", cds["product"])
        print("Protein ID:", cds["protein_id"])
        print("Coordinates:", cds["start"], "-", cds["end"])
        print("Strand:", cds["strand"])
        print("Protein length:", cds["protein_length"], "aa")

except Exception as error:

    print("TEST FAILED:", error)
