
from flask import Flask, render_template, request, send_file, send_from_directory
from werkzeug.exceptions import RequestEntityTooLarge
from Bio.Seq import Seq
from io import BytesIO
import math
import re

from genomics.sequence_analysis import analyze_sequence
from genomics.ncbi_annotation import annotate_orfs, annotate_genome
from genomics.omics_integration import analyze_central_dogma

# ============================================================
# APPLICATION CONFIGURATION
# ============================================================

app = Flask(__name__)

MAX_UPLOAD_SIZE = 5 * 1024 * 1024
MAX_UPLOAD_ERROR = (
    "The uploaded file is too large. "
    "The maximum allowed size is 5 MB."
)

app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_SIZE

PER_PAGE = 50

DNA_BASES = set("ATGCRYSWKMBDHVN")
RNA_BASES = set("AUGCRYSWKMBDHVN")
PROTEIN_BASES = set("ACDEFGHIKLMNPQRSTVWYBZJXUO*")

# ============================================================
# SEQUENCE VALIDATION
# ============================================================

def remove_fasta_formatting(sequence):

    # Detect separate FASTA records before removing headers.
    fasta_headers = re.findall(
        r"(?m)^\s*>",
        sequence
    )

    if len(fasta_headers) > 1:

        raise ValueError(
            "Multiple FASTA records detected. "
            "Please submit one sequence at a time."
        )

    sequence = re.sub(
        r"^>.*$",
        "",
        sequence,
        flags=re.MULTILINE
    )

    return re.sub(r"\s+", "", sequence).upper()


def clean_dna_sequence(sequence):

    sequence = remove_fasta_formatting(sequence)

    if not sequence:
        raise ValueError("Please enter a DNA sequence.")

    invalid = set(sequence) - DNA_BASES

    if invalid:
        raise ValueError(
            "Invalid DNA characters detected: "
            + ", ".join(sorted(invalid))
            + ". Use standard DNA IUPAC codes: "
            "A, T, G, C, R, Y, S, W, K, M, B, D, H, V and N."
        )

    return sequence


def clean_rna_sequence(sequence):

    sequence = remove_fasta_formatting(sequence)

    if not sequence:
        raise ValueError("Please enter an RNA sequence.")

    invalid = set(sequence) - RNA_BASES

    if invalid:
        raise ValueError(
            "Invalid RNA characters detected: "
            + ", ".join(sorted(invalid))
            + ". Use standard RNA IUPAC codes: "
            "A, U, G, C, R, Y, S, W, K, M, B, D, H, V and N."
        )

    return sequence


def clean_protein_sequence(sequence):

    sequence = remove_fasta_formatting(sequence)

    if not sequence:
        raise ValueError("Please provide a protein sequence.")

    invalid = set(sequence) - PROTEIN_BASES

    if invalid:
        raise ValueError(
            "Invalid protein characters detected: "
            + ", ".join(sorted(invalid))
            + ". Use standard amino acid codes, "
            "including X for unknown residues and * for stop codons."
        )

    return sequence


def clean_dna_or_rna_sequence(sequence):

    sequence = remove_fasta_formatting(sequence)

    if not sequence:
        raise ValueError(
            "Please provide a DNA or RNA sequence."
        )

    if "T" in sequence and "U" in sequence:
        raise ValueError(
            "The sequence contains both T and U. "
            "Please provide either DNA or RNA."
        )

    if "U" in sequence:
        return clean_rna_sequence(sequence)

    return clean_dna_sequence(sequence)


# ============================================================
# UPLOAD SIZE HANDLING AND INPUT READING
# ============================================================

def read_dna_input():

    try:

        sequence = request.form.get(
            "sequence", ""
        ).strip()

        uploaded_file = request.files.get("file")

    except RequestEntityTooLarge:

        raise ValueError(MAX_UPLOAD_ERROR)

    if uploaded_file and uploaded_file.filename:

        try:

            sequence = uploaded_file.read().decode("utf-8")

        except UnicodeDecodeError:

            raise ValueError(
                "The uploaded file is not valid UTF-8 text."
            )

        except RequestEntityTooLarge:

            raise ValueError(MAX_UPLOAD_ERROR)

    return clean_dna_sequence(sequence)


def read_rna_input():

    try:

        sequence = request.form.get(
            "sequence", ""
        ).strip()

        uploaded_file = request.files.get("file")

    except RequestEntityTooLarge:

        raise ValueError(MAX_UPLOAD_ERROR)

    if uploaded_file and uploaded_file.filename:

        try:

            sequence = uploaded_file.read().decode("utf-8")

        except UnicodeDecodeError:

            raise ValueError(
                "The uploaded file is not valid UTF-8 text."
            )

        except RequestEntityTooLarge:

            raise ValueError(MAX_UPLOAD_ERROR)

    return clean_rna_sequence(sequence)


# ============================================================
# DNA COMPOSITION
# ============================================================

def calculate_dna_results(sequence):

    return analyze_sequence(sequence)


# ============================================================
# SIX-FRAME TRANSLATION
# ============================================================

def translate_six_frames(sequence):

    dna = Seq(sequence)
    translations = []

    for frame in range(3):

        frame_sequence = str(dna[frame:])

        usable_length = len(frame_sequence) - (
            len(frame_sequence) % 3
        )

        protein = str(
            Seq(frame_sequence[:usable_length]).translate()
        )

        translations.append({
            "strand": "+",
            "frame": frame + 1,
            "protein_sequence": protein
        })

    reverse = str(dna.reverse_complement())

    for frame in range(3):

        frame_sequence = reverse[frame:]

        usable_length = len(frame_sequence) - (
            len(frame_sequence) % 3
        )

        protein = str(
            Seq(frame_sequence[:usable_length]).translate()
        )

        translations.append({
            "strand": "-",
            "frame": frame + 1,
            "protein_sequence": protein
        })

    return translations


# ============================================================
# ORF OVERLAP AND NESTING CLASSIFICATION
# ============================================================

def classify_orf_relationships(orfs):

    """
    Classify relationships between predicted ORFs.

    Classifications:

    No overlap:
        No other ORF intersects this ORF.

    Nested ORF:
        This ORF is completely contained within another
        same-strand ORF.

    Contains nested ORF:
        This ORF completely contains another same-strand ORF.

    Overlapping ORF:
        ORFs intersect, but neither completely contains
        the other on the same strand.

    Coordinates are one-based and inclusive.
    """

    for orf in orfs:

        orf["nested_in"] = []
        orf["nested_orfs"] = []
        orf["overlapping_orfs"] = []
        orf["overlap_status"] = "No overlap"

    for index, orf in enumerate(orfs):

        for other_index, other in enumerate(orfs):

            if index == other_index:
                continue

            overlap_start = max(
                orf["start"],
                other["start"]
            )

            overlap_end = min(
                orf["end"],
                other["end"]
            )

            if overlap_start > overlap_end:
                continue

            same_strand = (
                orf["strand"] == other["strand"]
            )

            if (
                same_strand
                and other["start"] <= orf["start"]
                and other["end"] >= orf["end"]
            ):

                orf["nested_in"].append(
                    other["number"]
                )

            elif (
                same_strand
                and orf["start"] <= other["start"]
                and orf["end"] >= other["end"]
            ):

                orf["nested_orfs"].append(
                    other["number"]
                )

            else:

                orf["overlapping_orfs"].append(
                    other["number"]
                )

    for orf in orfs:

        if orf["nested_in"]:

            orf["overlap_status"] = "Nested ORF"

        elif orf["nested_orfs"]:

            orf["overlap_status"] = "Contains nested ORF"

        elif orf["overlapping_orfs"]:

            orf["overlap_status"] = "Overlapping ORF"

        else:

            orf["overlap_status"] = "No overlap"

    return orfs


# ============================================================
# ORF PREDICTION
# ============================================================

def find_orfs(sequence, min_length=30):

    """
    Identify candidate ORFs in all six reading frames.

    Rules:
    - ATG is the accepted start codon.
    - TAA, TAG and TGA are accepted stop codons.
    - Minimum length is measured in amino acids,
      excluding the terminal stop codon.
    - Multiple start codons before a stop codon
      can produce nested ORF candidates.
    - An ambiguous codon interrupts candidate ORFs.
    - Coordinates are one-based and inclusive.
    - Nested and overlapping relationships are classified
      separately after ORF prediction.

    This is a basic computational ORF predictor,
    not a gene-expression or gene-function detector.
    """

    if not isinstance(min_length, int) or min_length < 1:

        raise ValueError(
            "Minimum ORF length must be a positive integer."
        )

    sequence = clean_dna_sequence(sequence)

    stop_codons = {"TAA", "TAG", "TGA"}

    canonical_codons = {"A", "T", "G", "C"}

    all_orfs = []

    sequence_length = len(sequence)

    for strand in ["+", "-"]:

        if strand == "+":

            working_sequence = sequence

        else:

            working_sequence = str(
                Seq(sequence).reverse_complement()
            )

        for frame in range(3):

            start_positions = []

            for position in range(
                frame,
                len(working_sequence) - 2,
                3
            ):

                codon = working_sequence[
                    position:position + 3
                ]

                if not set(codon).issubset(canonical_codons):

                    start_positions = []
                    continue

                if codon == "ATG":

                    start_positions.append(position)

                elif codon in stop_codons:

                    for start_position in start_positions:

                        protein_length = (
                            position - start_position
                        ) // 3

                        if protein_length < min_length:
                            continue

                        coding_dna = working_sequence[
                            start_position:position
                        ]

                        protein = str(
                            Seq(coding_dna).translate()
                        )

                        if strand == "+":

                            genomic_start = (
                                start_position + 1
                            )

                            genomic_end = (
                                position + 3
                            )

                        else:

                            genomic_start = (
                                sequence_length - position - 2
                            )

                            genomic_end = (
                                sequence_length - start_position
                            )

                        all_orfs.append({

                            "start": genomic_start,
                            "end": genomic_end,

                            "strand": strand,
                            "frame": frame + 1,

                            "length": protein_length,

                            "dna_sequence": coding_dna,

                            "protein_sequence": protein,

                            "protein_length": protein_length,

                            "overlap_status": "No overlap",

                            "nested_in": [],
                            "nested_orfs": [],
                            "overlapping_orfs": []

                        })

                    start_positions = []

    all_orfs.sort(
        key=lambda item: (
            item["start"],
            item["end"],
            item["strand"]
        )
    )

    for number, orf in enumerate(all_orfs, start=1):

        orf["number"] = number
        orf["orf_number"] = number

    all_orfs = classify_orf_relationships(
        all_orfs
    )

    return all_orfs


# ============================================================
# SIX-FRAME ORF TRACKS
# ============================================================

def build_orf_tracks(orfs, sequence_length):

    tracks = []

    number_of_ticks = 5

    ticks = sorted(set(
        max(
            1,
            round(
                1 + (sequence_length - 1) * i / number_of_ticks
            )
        )
        for i in range(number_of_ticks + 1)
    ))

    for strand in ["+", "-"]:

        for frame in range(1, 4):

            frame_orfs = [

                dict(orf)

                for orf in orfs

                if (
                    orf["strand"] == strand
                    and orf["frame"] == frame
                )

            ]

            frame_orfs.sort(
                key=lambda item: (
                    item["start"],
                    item["end"]
                )
            )

            lane_ends = []

            for orf in frame_orfs:

                assigned_lane = None

                for lane_number, lane_end in enumerate(lane_ends):

                    if orf["start"] > lane_end:

                        assigned_lane = lane_number
                        break

                if assigned_lane is None:

                    assigned_lane = len(lane_ends)

                    lane_ends.append(
                        orf["end"]
                    )

                else:

                    lane_ends[assigned_lane] = (
                        orf["end"]
                    )

                orf["lane"] = assigned_lane

            tracks.append({

                "strand": strand,
                "frame": frame,
                "orfs": frame_orfs,
                "lane_count": max(1, len(lane_ends)),
                "ticks": ticks

            })

    return tracks


# ============================================================
# RNA COMPOSITION AND TRANSFORMATIONS
# ============================================================

def calculate_rna_results(sequence):

    counts = {

        "A": sequence.count("A"),
        "U": sequence.count("U"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
        "N": sequence.count("N")

    }

    length = len(sequence)

    unambiguous_count = (
        counts["A"]
        + counts["U"]
        + counts["G"]
        + counts["C"]
    )

    gc_content = (

        (
            (counts["G"] + counts["C"])
            / unambiguous_count
        ) * 100

        if unambiguous_count else 0

    )

    ambiguity_codes = "RYSWKMBDHVN"

    ambiguous_bases = {

        code: sequence.count(code)

        for code in ambiguity_codes

    }

    ambiguous_count = sum(
        ambiguous_bases.values()
    )

    return {

        "sequence_length": length,

        "A": counts["A"],
        "U": counts["U"],
        "G": counts["G"],
        "C": counts["C"],
        "N": counts["N"],

        "ambiguous_count": ambiguous_count,

        "ambiguous_bases": ambiguous_bases,

        "GC_content": round(gc_content, 2),

        "complementary_sequence": str(
            Seq(sequence).complement_rna()
        ),

        "reverse_complement": str(
            Seq(sequence).reverse_complement_rna()
        )

    }


def translate_rna_complete(sequence):

    usable_length = len(sequence) - (
        len(sequence) % 3
    )

    if usable_length == 0:
        return ""

    return str(
        Seq(sequence[:usable_length]).translate()
    )


# ============================================================
# FASTA DOWNLOAD
# ============================================================

def fasta_download(sequence, filename, description):

    fasta_text = f">{description}\n{sequence}\n"

    file_buffer = BytesIO(
        fasta_text.encode("utf-8")
    )

    file_buffer.seek(0)

    return send_file(
        file_buffer,
        mimetype="text/plain",
        as_attachment=True,
        download_name=filename
    )


# ============================================================
# HOME PAGE
# ============================================================

@app.route("/")
def home():

    return render_template("index.html")


# ============================================================
# SEO FILES: ROBOTS.TXT AND SITEMAP.XML
# ============================================================

@app.route("/robots.txt")
def robots_txt():

    return send_from_directory(
        app.root_path,
        "robots.txt",
        mimetype="text/plain"
    )


@app.route("/sitemap.xml")
def sitemap_xml():

    return send_from_directory(
        app.root_path,
        "sitemap.xml",
        mimetype="application/xml"
    )


# ============================================================
# GENOMICS WORKSPACE
# ============================================================

@app.route("/genomics", methods=["GET", "POST"])
def genomics():

    if request.method == "GET":

        return render_template("genomics.html")

    sequence = ""

    try:

        sequence = read_dna_input()

        results = calculate_dna_results(sequence)

        predicted_orfs = find_orfs(sequence)

        translations = translate_six_frames(sequence)

        orf_tracks = build_orf_tracks(
            predicted_orfs,
            len(sequence)
        )

        return render_template(
            "genomics.html",
            results=results,
            sequence_input=sequence,
            sequence_length=len(sequence),
            predicted_orfs=predicted_orfs,
            translations=translations,
            orf_tracks=orf_tracks,
            error=None
        )

    except Exception as error:

        return render_template(
            "genomics.html",
            sequence_input=sequence,
            error=str(error)
        )


# ============================================================
# NCBI GENOME EXPLORER
# ============================================================

def prepare_genome_explorer(
    accession,
    search="",
    page=1
):

    genome_data = annotate_genome(accession)

    all_records = genome_data.get(
        "cds_records", []
    )

    if search:

        search_lower = search.lower()

        filtered_records = [

            record for record in all_records

            if search_lower in " ".join([

                str(record.get("gene", "")),
                str(record.get("locus_tag", "")),
                str(record.get("product", "")),
                str(record.get("protein_id", ""))

            ]).lower()

        ]

    else:

        filtered_records = all_records

    total_records = len(filtered_records)

    total_pages = max(
        1,
        math.ceil(total_records / PER_PAGE)
    )

    page = max(
        1,
        min(page, total_pages)
    )

    start_index = (page - 1) * PER_PAGE

    displayed_records = filtered_records[
        start_index:start_index + PER_PAGE
    ]

    safe_records = []

    for record in displayed_records:

        safe_record = dict(record)

        safe_record.pop(
            "protein_sequence",
            None
        )

        safe_records.append(safe_record)

    genome_summary = dict(genome_data)

    genome_summary.pop(
        "cds_records",
        None
    )

    return {

        "genome_summary": genome_summary,
        "cds_records": safe_records,
        "search_term": search,
        "total_records": total_records,
        "page": page,
        "total_pages": total_pages

    }


# ============================================================
# NCBI ANNOTATION WORKSPACE
# ============================================================

@app.route("/annotation", methods=["GET", "POST"])
def annotation():

    if request.method == "GET":

        analysis_mode = request.args.get(
            "analysis_mode",
            "protein_comparison"
        )

        if analysis_mode == "genome_explorer":

            accession = request.args.get(
                "accession", ""
            ).strip()

            search = request.args.get(
                "search", ""
            ).strip()

            try:

                page = max(
                    1,
                    int(request.args.get("page", 1))
                )

            except ValueError:

                page = 1

            if not accession:

                return render_template(
                    "annotation.html",
                    analysis_mode="genome_explorer",
                    accession=accession,
                    search_term=search,
                    error="Please enter an NCBI accession number."
                )

            try:

                explorer_data = prepare_genome_explorer(
                    accession,
                    search,
                    page
                )

                return render_template(
                    "annotation.html",
                    analysis_mode="genome_explorer",
                    accession=accession,
                    error=None,
                    **explorer_data
                )

            except Exception as error:

                return render_template(
                    "annotation.html",
                    analysis_mode="genome_explorer",
                    accession=accession,
                    search_term=search,
                    error=str(error)
                )

        return render_template(
            "annotation.html",
            analysis_mode="protein_comparison"
        )

    accession = request.form.get(
        "accession", ""
    ).strip()

    sequence_input = ""

    try:

        sequence = read_dna_input()

        sequence_input = sequence

        predicted_orfs = find_orfs(sequence)

        results = annotate_orfs(
            accession,
            predicted_orfs
        )

        return render_template(
            "annotation.html",
            analysis_mode="protein_comparison",
            accession=accession,
            sequence_input=sequence_input,
            analysis=calculate_dna_results(sequence),
            predicted_orfs=predicted_orfs,
            results=results,
            matches=results.get("matches", []),
            best_match=results.get("matches", [None])[0]
            if results.get("matches") else None,
            error=None
        )

    except Exception as error:

        return render_template(
            "annotation.html",
            analysis_mode="protein_comparison",
            accession=accession,
            sequence_input=sequence_input,
            error=str(error)
        )


# ============================================================
# GENOMICS DOWNLOAD ROUTES
# ============================================================

@app.route("/download-reverse-complement", methods=["POST"])
def download_reverse_complement():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    result = str(
        Seq(sequence).reverse_complement()
    )

    return fasta_download(
        result,
        "reverse_complement.fasta",
        "DNA_reverse_complement"
    )


@app.route("/download-protein-fasta", methods=["POST"])
def download_protein_fasta():

    protein = request.form.get(
        "protein_sequence", ""
    ).strip()

    if protein:

        protein = clean_protein_sequence(protein)

    else:

        sequence = clean_dna_sequence(
            request.form.get("sequence", "")
        )

        usable_length = len(sequence) - len(sequence) % 3

        protein = str(
            Seq(sequence[:usable_length]).translate()
        )

    return fasta_download(
        protein,
        "translated_protein.fasta",
        "translated_protein"
    )


@app.route("/download-longest-orf", methods=["POST"])
def download_longest_orf():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    orfs = find_orfs(sequence)

    if not orfs:

        return "No qualifying ORFs were found.", 400

    longest_orf = max(
        orfs,
        key=lambda item: item["length"]
    )

    return fasta_download(
        longest_orf["protein_sequence"],
        "longest_orf.fasta",
        "longest_ORF"
    )


# ============================================================
# TRANSCRIPTOMICS WORKSPACE
# ============================================================

@app.route("/transcriptomics", methods=["GET", "POST"])
def transcriptomics():

    if request.method == "GET":

        return render_template("transcriptomics.html")

    sequence_input = ""

    try:

        rna_sequence = read_rna_input()

        sequence_input = rna_sequence

        results = calculate_rna_results(
            rna_sequence
        )

        dna_sequence = rna_sequence.replace(
            "U", "T"
        )

        if len(rna_sequence) % 3 == 0:

            protein_sequence = translate_rna_complete(
                rna_sequence
            )

        else:

            protein_sequence = ""

        return render_template(
            "transcriptomics.html",
            results=results,
            sequence_input=sequence_input,
            dna_sequence=dna_sequence,
            protein_sequence=protein_sequence,
            error=None
        )

    except Exception as error:

        return render_template(
            "transcriptomics.html",
            sequence_input=sequence_input,
            error=str(error)
        )


# ============================================================
# RNA DOWNLOAD ROUTES
# ============================================================

@app.route("/download-rna-complement", methods=["POST"])
def download_rna_complement():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    return fasta_download(
        str(Seq(sequence).complement_rna()),
        "rna_complement.fasta",
        "RNA_complement"
    )


@app.route("/download-rna-reverse-complement", methods=["POST"])
def download_rna_reverse_complement():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    return fasta_download(
        str(Seq(sequence).reverse_complement_rna()),
        "rna_reverse_complement.fasta",
        "RNA_reverse_complement"
    )


@app.route("/download-rna-dna", methods=["POST"])
def download_rna_dna():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    return fasta_download(
        sequence.replace("U", "T"),
        "rna_dna_equivalent.fasta",
        "DNA_equivalent"
    )


@app.route("/download-rna-protein", methods=["POST"])
def download_rna_protein():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    if len(sequence) % 3 != 0:

        return (
            "Protein translation is unavailable because "
            "the RNA sequence length is not divisible by three.",
            400
        )

    protein = translate_rna_complete(sequence)

    return fasta_download(
        protein,
        "rna_translated_protein.fasta",
        "translated_protein"
    )


# ============================================================
# OMICS INTEGRATION
# ============================================================

@app.route("/integration", methods=["GET", "POST"])
def integration():

    if request.method == "GET":

        return render_template("integration.html")

    sequence = ""

    try:

        sequence = read_dna_input()

        results = analyze_central_dogma(sequence)

        return render_template(
            "integration.html",
            results=results,
            sequence=sequence,
            error=None
        )

    except Exception as error:

        return render_template(
            "integration.html",
            sequence=sequence,
            error=str(error)
        )


# ============================================================
# INTEGRATION DOWNLOAD ROUTES
# ============================================================

@app.route("/download-integration-dna", methods=["POST"])
def download_integration_dna():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    return fasta_download(
        sequence,
        "integration_dna.fasta",
        "DNA_sequence"
    )


@app.route("/download-integration-rna", methods=["POST"])
def download_integration_rna():

    sequence = clean_dna_or_rna_sequence(
        request.form.get("sequence", "")
    )

    if "U" in sequence:

        mrna_sequence = sequence

    else:

        mrna_sequence = sequence.replace(
            "T", "U"
        )

    return fasta_download(
        mrna_sequence,
        "integration_mrna.fasta",
        "mRNA_sequence"
    )


@app.route("/download-integration-protein", methods=["POST"])
def download_integration_protein():

    protein_input = request.form.get(
        "protein_sequence", ""
    ).strip()

    if protein_input:

        protein = clean_protein_sequence(
            protein_input
        )

    else:

        legacy_input = request.form.get(
            "sequence", ""
        ).strip()

        if not legacy_input:

            raise ValueError(
                "No protein or nucleotide sequence was provided."
            )

        cleaned_legacy = remove_fasta_formatting(
            legacy_input
        )

        if set(cleaned_legacy).issubset(DNA_BASES | RNA_BASES):

            sequence = clean_dna_or_rna_sequence(
                cleaned_legacy
            )

            usable_length = len(sequence) - (
                len(sequence) % 3
            )

            if usable_length == 0:

                return (
                    "The sequence does not contain a complete codon.",
                    400
                )

            protein = str(
                Seq(sequence[:usable_length]).translate()
            )

        else:

            protein = clean_protein_sequence(
                cleaned_legacy
            )

    return fasta_download(
        protein,
        "integration_protein.fasta",
        "translated_protein"
    )


# ============================================================
# FILE SIZE ERROR HANDLER
# ============================================================

@app.errorhandler(RequestEntityTooLarge)
@app.errorhandler(413)
def file_too_large(error):

    return (
        MAX_UPLOAD_ERROR,
        413
    )


# ============================================================
# APPLICATION START
# ============================================================

if __name__ == "__main__":

    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True,
        use_reloader=False
    )
