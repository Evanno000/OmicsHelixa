
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from io import StringIO
import os


# ==========================================
# PROTEIN ALIGNMENT CONFIGURATION
# ==========================================

aligner = PairwiseAligner()

aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5


# ==========================================
# COMPARE TWO PROTEIN SEQUENCES
# ==========================================

def compare_proteins(query, reference):

    if not query or not reference:
        return None

    query = query.upper()
    reference = reference.upper()

    alignment = aligner.align(query, reference)[0]

    indices = alignment.indices

    aligned_positions = 0
    identical_positions = 0

    query_aligned_positions = 0
    reference_aligned_positions = 0

    for query_position, reference_position in zip(
        indices[0],
        indices[1]
    ):

        # Count residues from each sequence included
        # in the local alignment.
        if query_position >= 0:
            query_aligned_positions += 1

        if reference_position >= 0:
            reference_aligned_positions += 1

        # Identity requires a residue from both sequences.
        if query_position < 0 or reference_position < 0:
            continue

        aligned_positions += 1

        if query[query_position] == reference[reference_position]:
            identical_positions += 1

    if aligned_positions == 0:
        return None

    identity = (
        identical_positions / aligned_positions
    ) * 100

    query_coverage = (
        query_aligned_positions / len(query)
    ) * 100

    reference_coverage = (
        reference_aligned_positions / len(reference)
    ) * 100

    return {
        "identity": round(identity, 2),
        "query_coverage": round(query_coverage, 2),
        "reference_coverage": round(reference_coverage, 2),
        "alignment_score": round(float(alignment.score), 2),
        "aligned_positions": aligned_positions
    }


# ==========================================
# FETCH GENBANK RECORD FROM NCBI
# ==========================================

def fetch_ncbi_record(accession):

    ncbi_email = os.environ.get(
        "NCBI_EMAIL",
        ""
    ).strip()

    if not ncbi_email:

        raise ValueError(
            "NCBI email is not configured. "
            "Set the NCBI_EMAIL environment variable."
        )

    Entrez.email = ncbi_email

    if not accession.strip():

        raise ValueError(
            "Please provide an NCBI accession number."
        )

    try:

        handle = Entrez.efetch(
            db="nuccore",
            id=accession.strip(),
            rettype="gb",
            retmode="text"
        )

        try:
            record_text = handle.read()

        finally:
            handle.close()

        if not record_text.strip():

            raise ValueError(
                "NCBI returned an empty record."
            )

        record = SeqIO.read(
            StringIO(record_text),
            "genbank"
        )

        return record

    except Exception as error:

        raise ValueError(
            f"Unable to retrieve the NCBI record: {error}"
        )


# ==========================================
# EXTRACT ANNOTATED PROTEINS
# ==========================================

def extract_annotated_proteins(record):

    proteins = []

    for feature in record.features:

        if feature.type != "CDS":
            continue

        qualifiers = feature.qualifiers

        translation = qualifiers.get(
            "translation",
            [""]
        )[0]

        if not translation:
            continue

        proteins.append({

            "gene": qualifiers.get(
                "gene",
                ["Unknown"]
            )[0],

            "product": qualifiers.get(
                "product",
                ["Unannotated protein"]
            )[0],

            "protein_id": qualifiers.get(
                "protein_id",
                ["Not provided"]
            )[0],

            "translation": translation,

            "location": str(feature.location)

        })

    return proteins


# ==========================================
# ANNOTATE PREDICTED ORFS
# ==========================================

def annotate_orfs(accession, predicted_orfs):

    record = fetch_ncbi_record(accession)

    annotated_proteins = extract_annotated_proteins(
        record
    )

    if not annotated_proteins:

        raise ValueError(
            "The NCBI record contains no annotated CDS translations."
        )

    matches = []

    for orf in predicted_orfs:

        query = orf.get(
            "protein_sequence",
            ""
        )

        if not query:
            continue

        best_match = None

        for reference in annotated_proteins:

            comparison = compare_proteins(
                query,
                reference["translation"]
            )

            if comparison is None:
                continue

            candidate = {
                **reference,
                **comparison
            }

            # Rank candidates using both query and
            # reference coverage, then identity and score.
            candidate_rank = (
                min(
                    candidate["query_coverage"],
                    candidate["reference_coverage"]
                ),
                candidate["identity"],
                candidate["alignment_score"]
            )

            if best_match is None:

                best_match = candidate

            else:

                best_rank = (
                    min(
                        best_match["query_coverage"],
                        best_match["reference_coverage"]
                    ),
                    best_match["identity"],
                    best_match["alignment_score"]
                )

                if candidate_rank > best_rank:
                    best_match = candidate

        if best_match:

            matches.append({

                "orf_number": orf.get("number"),
                "orf_length": orf.get("length"),
                "orf_strand": orf.get("strand"),
                "orf_frame": orf.get("frame"),
                "match": best_match

            })

    # Sort ORFs by the same combined quality criteria.
    matches.sort(
        key=lambda item: (
            min(
                item["match"]["query_coverage"],
                item["match"]["reference_coverage"]
            ),
            item["match"]["identity"],
            item["match"]["alignment_score"]
        ),
        reverse=True
    )

    return {

        "accession": record.id,

        "description": record.description,

        "organism": record.annotations.get(
            "organism",
            "Unknown"
        ),

        "record_length": len(record.seq),

        "annotated_protein_count": len(annotated_proteins),

        "matches": matches

    }
