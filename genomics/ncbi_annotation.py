
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from io import StringIO
from urllib.error import HTTPError, URLError
from http.client import IncompleteRead, RemoteDisconnected
from functools import lru_cache
import os
import re
import time


# ============================================================
# NCBI DOWNLOAD CONFIGURATION
# ============================================================

NCBI_MAX_ATTEMPTS = 3
NCBI_RETRY_DELAY = 5
NCBI_TIMEOUT = 180


# ============================================================
# PROTEIN ALIGNMENT CONFIGURATION
# ============================================================

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5


# ============================================================
# PROTEIN COMPARISON SCREENING THRESHOLDS
# ============================================================

MIN_IDENTITY = 30.0
MIN_QUERY_COVERAGE = 50.0
MIN_REFERENCE_COVERAGE = 20.0
MIN_ALIGNED_POSITIONS = 20


# ============================================================
# NCBI ACCESSION VALIDATION
# ============================================================

def validate_accession(accession):

    if not accession or not accession.strip():
        raise ValueError(
            "Please provide an NCBI accession number."
        )

    accession = accession.strip()

    pattern = r"^[A-Za-z]{1,4}_?\d{5,9}(?:\.\d+)?$"

    if not re.fullmatch(pattern, accession):
        raise ValueError(
            "Invalid accession format. Please enter a valid "
            "NCBI nucleotide accession, for example "
            "U91966.1, CP002686.1 or NC_003074.8."
        )

    return accession


# ============================================================
# FETCH GENBANK RECORD FROM NCBI WITH RETRIES
# ============================================================

def fetch_ncbi_record(accession):

    accession = validate_accession(accession)

    ncbi_email = os.environ.get(
        "NCBI_EMAIL", ""
    ).strip()

    if not ncbi_email:
        raise ValueError(
            "NCBI email is not configured. "
            "Please configure the NCBI_EMAIL environment variable."
        )

    Entrez.email = ncbi_email

    # Configure Biopython's own retry behaviour.
    Entrez.max_tries = 3
    Entrez.sleep_between_tries = 10

    last_error = None

    retryable_errors = (
        IncompleteRead,
        RemoteDisconnected,
        ConnectionResetError,
        TimeoutError,
        URLError,
        OSError
    )

    for attempt in range(1, NCBI_MAX_ATTEMPTS + 1):

        handle = None

        try:

            print(
                f"Connecting to NCBI: {accession} "
                f"(attempt {attempt}/{NCBI_MAX_ATTEMPTS})",
                flush=True
            )

            handle = Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gb",
                retmode="text",
                timeout=NCBI_TIMEOUT
            )

            print(
                "NCBI response received. "
                "Downloading GenBank record...",
                flush=True
            )

            record_text = handle.read()

            if not record_text or not record_text.strip():

                raise ValueError(
                    "NCBI returned an empty record. "
                    "Please check the accession number."
                )

            print(
                f"Downloaded {len(record_text):,} characters. "
                "Parsing GenBank record...",
                flush=True
            )

            try:

                record = SeqIO.read(
                    StringIO(record_text),
                    "genbank"
                )

            except (ValueError, TypeError) as error:

                raise ValueError(
                    f"GenBank parsing failed: {error}"
                ) from error

            print(
                f"Successfully retrieved {record.id}. "
                f"Sequence length: {len(record.seq):,} bp. "
                f"Features: {len(record.features):,}.",
                flush=True
            )

            return record

        except HTTPError as error:

            if error.code == 400:

                raise ValueError(
                    "NCBI could not process this accession. "
                    "Please check that it is a valid nucleotide accession."
                ) from error

            elif error.code == 404:

                raise ValueError(
                    "The requested NCBI record could not be found."
                ) from error

            elif error.code == 429:

                last_error = error

                print(
                    "NCBI rate limit reached. "
                    "Waiting before retrying...",
                    flush=True
                )

            elif error.code >= 500:

                last_error = error

                print(
                    f"NCBI server error {error.code}. "
                    "Waiting before retrying...",
                    flush=True
                )

            else:

                raise ValueError(
                    f"NCBI returned an HTTP error ({error.code}). "
                    "Please try again later."
                ) from error

        except retryable_errors as error:

            last_error = error

            print(
                f"NCBI download interrupted: "
                f"{type(error).__name__}: {error}",
                flush=True
            )

        except ValueError:
            raise

        except Exception as error:

            raise ValueError(
                f"NCBI retrieval error: "
                f"{type(error).__name__}: {error}"
            ) from error

        finally:

            if handle is not None:

                try:
                    handle.close()
                except Exception:
                    pass

        if attempt < NCBI_MAX_ATTEMPTS:

            delay = NCBI_RETRY_DELAY * attempt

            print(
                f"Retrying NCBI download in {delay} seconds...",
                flush=True
            )

            time.sleep(delay)

    raise ValueError(
        f"NCBI could not provide the complete record "
        f"after {NCBI_MAX_ATTEMPTS} attempts. "
        f"Last error: {type(last_error).__name__}: {last_error}. "
        "Please try again later or check your internet connection."
    )


# ============================================================
# EXTRACT ANNOTATED PROTEINS
# ============================================================

def extract_annotated_proteins(record):

    proteins = []

    for feature in record.features:

        if feature.type != "CDS":
            continue

        qualifiers = feature.qualifiers

        translation = qualifiers.get(
            "translation", [""]
        )[0]

        translation = re.sub(
            r"\s+", "", translation
        ).upper()

        translation = translation.rstrip("*")

        if not translation:
            continue

        proteins.append({

            "gene": qualifiers.get(
                "gene", ["Unknown"]
            )[0],

            "product": qualifiers.get(
                "product", ["Unannotated protein"]
            )[0],

            "protein_id": qualifiers.get(
                "protein_id", ["Not provided"]
            )[0],

            "translation": translation,

            "location": str(feature.location)

        })

    return proteins


# ============================================================
# COMPARE TWO PROTEIN SEQUENCES
# ============================================================

def compare_proteins(query, reference):

    if not query or not reference:
        return None

    query = re.sub(r"\s+", "", query).upper()
    reference = re.sub(r"\s+", "", reference).upper()

    if not query or not reference:
        return None

    try:

        alignment = aligner.align(
            query,
            reference
        )[0]

    except Exception:
        return None

    indices = alignment.indices

    aligned_positions = 0
    identical_positions = 0

    query_aligned_positions = 0
    reference_aligned_positions = 0

    for query_position, reference_position in zip(
        indices[0],
        indices[1]
    ):

        if query_position >= 0:
            query_aligned_positions += 1

        if reference_position >= 0:
            reference_aligned_positions += 1

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

        "alignment_score": round(
            float(alignment.score), 2
        ),

        "aligned_positions": aligned_positions

    }


# ============================================================
# ALIGNMENT SCREENING
# ============================================================

def passes_alignment_thresholds(comparison):

    if comparison is None:
        return False

    return (

        comparison["identity"] >= MIN_IDENTITY

        and comparison["query_coverage"] >= MIN_QUERY_COVERAGE

        and comparison["reference_coverage"] >= MIN_REFERENCE_COVERAGE

        and comparison["aligned_positions"] >= MIN_ALIGNED_POSITIONS

    )


def get_alignment_rejection_reasons(comparison):

    reasons = []

    if comparison["identity"] < MIN_IDENTITY:

        reasons.append(
            f"Identity below {MIN_IDENTITY}%"
        )

    if comparison["query_coverage"] < MIN_QUERY_COVERAGE:

        reasons.append(
            f"Query coverage below {MIN_QUERY_COVERAGE}%"
        )

    if comparison["reference_coverage"] < MIN_REFERENCE_COVERAGE:

        reasons.append(
            f"Reference coverage below {MIN_REFERENCE_COVERAGE}%"
        )

    if comparison["aligned_positions"] < MIN_ALIGNED_POSITIONS:

        reasons.append(
            f"Fewer than {MIN_ALIGNED_POSITIONS} aligned amino-acid positions"
        )

    return reasons


# ============================================================
# RANK PROTEIN ALIGNMENTS
# ============================================================

def alignment_rank(candidate):

    return (

        min(
            candidate["query_coverage"],
            candidate["reference_coverage"]
        ),

        candidate["identity"],

        candidate["alignment_score"],

        candidate["aligned_positions"]

    )


# ============================================================
# ANNOTATE PREDICTED ORFS
# ============================================================

def annotate_orfs(accession, predicted_orfs):

    """
    Compare predicted ORF protein sequences against
    annotated CDS translations in an NCBI GenBank record.

    Candidate matches and filtered alignments are returned
    separately. These are computational results, not
    experimentally confirmed gene annotations.
    """

    if not isinstance(predicted_orfs, list):

        raise ValueError(
            "Predicted ORFs must be provided as a list. "
            "Run ORF prediction before protein comparison."
        )

    if not predicted_orfs:

        raise ValueError(
            "No predicted ORFs were provided. "
            "Run ORF prediction first."
        )

    for index, orf in enumerate(predicted_orfs):

        if not isinstance(orf, dict):

            raise ValueError(
                f"Invalid predicted ORF at position {index + 1}. "
                "Each ORF must be a dictionary containing "
                "its sequence information."
            )

        protein_sequence = orf.get(
            "protein_sequence", ""
        )

        if (
            protein_sequence is not None
            and not isinstance(protein_sequence, str)
        ):

            raise ValueError(
                f"Invalid protein sequence in ORF "
                f"at position {index + 1}. "
                "The protein sequence must be text."
            )

    record = fetch_ncbi_record(accession)

    annotated_proteins = extract_annotated_proteins(
        record
    )

    if not annotated_proteins:

        raise ValueError(
            "The NCBI record was retrieved successfully, "
            "but it contains no CDS features with annotated "
            "protein translations. Try a record with annotated "
            "coding sequences."
        )

    matches = []
    filtered_matches = []

    orfs_without_alignments = 0

    print(
        f"Comparing {len(predicted_orfs)} predicted ORFs "
        f"against {len(annotated_proteins)} annotated proteins...",
        flush=True
    )

    for orf in predicted_orfs:

        query = orf.get(
            "protein_sequence", ""
        )

        if not query:
            orfs_without_alignments += 1
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

            if (
                best_match is None
                or alignment_rank(candidate)
                > alignment_rank(best_match)
            ):

                best_match = candidate

        if best_match is None:

            orfs_without_alignments += 1
            continue

        result = {

            "orf_number": orf.get("number"),

            "orf_length": orf.get("length"),

            "orf_strand": orf.get("strand"),

            "orf_frame": orf.get("frame"),

            "orf_start": orf.get("start"),

            "orf_end": orf.get("end"),

            "gene": best_match["gene"],

            "product": best_match["product"],

            "protein_id": best_match["protein_id"],

            "identity": best_match["identity"],

            "query_coverage": best_match["query_coverage"],

            "reference_coverage": best_match["reference_coverage"],

            "alignment_score": best_match["alignment_score"],

            "aligned_positions": best_match["aligned_positions"],

            "reference_location": best_match["location"]

        }

        if passes_alignment_thresholds(best_match):

            result["match_status"] = "Candidate protein match"
            result["rejection_reasons"] = []

            matches.append(result)

        else:

            result["match_status"] = "Below screening thresholds"

            result["rejection_reasons"] = (
                get_alignment_rejection_reasons(best_match)
            )

            filtered_matches.append(result)

    matches.sort(
        key=lambda item: (
            min(
                item["query_coverage"],
                item["reference_coverage"]
            ),
            item["identity"],
            item["alignment_score"],
            item["aligned_positions"]
        ),
        reverse=True
    )

    filtered_matches.sort(
        key=lambda item: (
            min(
                item["query_coverage"],
                item["reference_coverage"]
            ),
            item["identity"],
            item["alignment_score"],
            item["aligned_positions"]
        ),
        reverse=True
    )

    print(
        f"Protein comparison complete. "
        f"{len(matches)} candidate matches passed screening. "
        f"{len(filtered_matches)} alignments were filtered. "
        f"{orfs_without_alignments} ORFs had no usable alignment.",
        flush=True
    )

    return {

        "accession": record.id,

        "description": record.description,

        "organism": record.annotations.get(
            "organism", "Unknown"
        ),

        "record_length": len(record.seq),

        "predicted_orf_count": len(predicted_orfs),

        "annotated_protein_count": len(annotated_proteins),

        "matches": matches,

        "candidate_match_count": len(matches),

        "filtered_matches": filtered_matches,

        "filtered_match_count": len(filtered_matches),

        "orfs_without_alignments": orfs_without_alignments,

        "comparison_thresholds": {

            "minimum_identity": MIN_IDENTITY,

            "minimum_query_coverage": MIN_QUERY_COVERAGE,

            "minimum_reference_coverage": MIN_REFERENCE_COVERAGE,

            "minimum_aligned_positions": MIN_ALIGNED_POSITIONS

        }

    }


# ============================================================
# ANNOTATED GENOME EXPLORER
# ============================================================

@lru_cache(maxsize=2)
def annotate_genome(accession):

    """
    Retrieve a GenBank record and extract existing
    annotated CDS features.

    Successful results are cached in memory.
    """

    accession = validate_accession(accession)

    record = fetch_ncbi_record(accession)

    cds_records = []

    gene_features = 0
    mrna_features = 0

    unique_gene_ids = set()

    print(
        "Extracting genome annotations...",
        flush=True
    )

    for feature in record.features:

        if feature.type == "gene":
            gene_features += 1

        elif feature.type == "mRNA":
            mrna_features += 1

        if feature.type != "CDS":
            continue

        qualifiers = feature.qualifiers

        gene = qualifiers.get(
            "gene", [""]
        )[0]

        locus_tag = qualifiers.get(
            "locus_tag", [""]
        )[0]

        product = qualifiers.get(
            "product", ["Unspecified protein product"]
        )[0]

        protein_id = qualifiers.get(
            "protein_id", ["Not available"]
        )[0]

        translation = qualifiers.get(
            "translation", [""]
        )[0]

        if not translation:
            continue

        gene_identifier = locus_tag or gene

        if gene_identifier:
            unique_gene_ids.add(gene_identifier)

        start = int(feature.location.start) + 1
        end = int(feature.location.end)

        strand_value = feature.location.strand

        if strand_value == 1:
            strand = "+"

        elif strand_value == -1:
            strand = "-"

        else:
            strand = "Unknown"

        extracted_sequence = feature.extract(
            record.seq
        )

        cds_records.append({

            "gene": gene or "Unknown",

            "locus_tag": locus_tag or "Not available",

            "product": product,

            "protein_id": protein_id,

            "start": start,

            "end": end,

            "strand": strand,

            "location": str(feature.location),

            "cds_length": len(extracted_sequence),

            "protein_length": len(translation),

            "protein_sequence": translation

        })

    if not cds_records:

        raise ValueError(
            "The NCBI record was retrieved, but no CDS "
            "features with annotated protein translations "
            "were found. Try an accession with annotated "
            "coding sequences, such as CP002686.1."
        )

    cds_records.sort(
        key=lambda item: (
            item["start"],
            item["end"]
        )
    )

    print(
        f"Genome annotation extraction complete: "
        f"{len(cds_records):,} CDS features.",
        flush=True
    )

    return {

        "accession": record.id,

        "description": record.description,

        "organism": record.annotations.get(
            "organism", "Unknown organism"
        ),

        "record_length": len(record.seq),

        "gene_feature_count": gene_features,

        "mrna_feature_count": mrna_features,

        "cds_count": len(cds_records),

        "unique_gene_count": len(unique_gene_ids),

        "cds_records": cds_records

    }
