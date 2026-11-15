
from flask import (
    Flask,
    render_template,
    request,
    send_file,
    abort
)

from Bio import SeqIO
from Bio.Seq import Seq

from io import StringIO, BytesIO
import heapq
import math

from genomics.sequence_analysis import analyze_sequence
from genomics.ncbi_annotation import (
    annotate_orfs,
    annotate_genome
)

from genomics.transcriptomics_analysis import (
    analyze_rna_sequence,
    convert_rna_to_dna,
    translate_rna
)

from genomics.omics_integration import analyze_central_dogma


# ============================================================
# APPLICATION CONFIGURATION
# ============================================================

app = Flask(__name__)

app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024

CDS_PER_PAGE = 50


# ============================================================
# DNA INPUT VALIDATION
# ============================================================

def clean_dna_sequence(raw_sequence):

    if not raw_sequence:
        raise ValueError("Please provide a DNA sequence.")

    lines = raw_sequence.splitlines()

    sequence = "".join(
        line.strip()
        for line in lines
        if not line.strip().startswith(">")
    )

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError("The DNA sequence is empty.")

    invalid = sorted(set(sequence) - set("ATGCN"))

    if invalid:
        raise ValueError(
            "Invalid DNA characters detected: "
            + ", ".join(invalid)
            + ". Use only A, T, G, C and N."
        )

    return sequence


def read_dna_input(pasted_sequence, uploaded_file):

    if uploaded_file and uploaded_file.filename:

        content = uploaded_file.read().decode(
            "utf-8",
            errors="replace"
        )

        if not content.strip():
            raise ValueError("The uploaded DNA file is empty.")

        try:
            records = list(
                SeqIO.parse(StringIO(content), "fasta")
            )

            if records:
                return clean_dna_sequence(str(records[0].seq))

        except Exception:
            pass

        return clean_dna_sequence(content)

    return clean_dna_sequence(pasted_sequence)


# ============================================================
# RNA INPUT VALIDATION
# ============================================================

def clean_rna_sequence(raw_sequence):

    if not raw_sequence:
        raise ValueError("Please provide an RNA sequence.")

    lines = raw_sequence.splitlines()

    sequence = "".join(
        line.strip()
        for line in lines
        if not line.strip().startswith(">")
    )

    sequence = "".join(sequence.split()).upper()

    if not sequence:
        raise ValueError("The RNA sequence is empty.")

    invalid = sorted(set(sequence) - set("AUGCN"))

    if invalid:
        raise ValueError(
            "Invalid RNA characters detected: "
            + ", ".join(invalid)
            + ". Use only A, U, G, C and N."
        )

    return sequence


def read_rna_input(pasted_sequence, uploaded_file):

    if uploaded_file and uploaded_file.filename:

        content = uploaded_file.read().decode(
            "utf-8",
            errors="replace"
        )

        if not content.strip():
            raise ValueError("The uploaded RNA file is empty.")

        try:
            records = list(
                SeqIO.parse(StringIO(content), "fasta")
            )

            if records:
                return clean_rna_sequence(str(records[0].seq))

        except Exception:
            pass

        return clean_rna_sequence(content)

    return clean_rna_sequence(pasted_sequence)


# ============================================================
# SIX-FRAME TRANSLATION
# ============================================================

def translate_six_frames(sequence):

    dna = Seq(sequence)

    reverse_complement = dna.reverse_complement()

    frames = {}

    for strand_name, strand_sequence in (
        ("+", dna),
        ("-", reverse_complement)
    ):

        for frame in range(3):

            translated = str(
                strand_sequence[frame:].translate(
                    to_stop=False
                )
            )

            frames[f"{strand_name}{frame + 1}"] = translated

    return frames


# ============================================================
# OPEN READING FRAME PREDICTION
# ============================================================

def find_orfs(sequence, minimum_length=30):

    sequence = sequence.upper()

    stop_codons = {"TAA", "TAG", "TGA"}

    orfs = []

    sequence_length = len(sequence)

    strands = [
        ("+", sequence),
        ("-", str(Seq(sequence).reverse_complement()))
    ]

    for strand, strand_sequence in strands:

        for frame in range(3):

            starts = []

            for position in range(
                frame,
                len(strand_sequence) - 2,
                3
            ):

                codon = strand_sequence[
                    position:position + 3
                ]

                if codon == "ATG":
                    starts.append(position)

                elif codon in stop_codons:

                    for start_position in starts:

                        protein_sequence = str(
                            Seq(
                                strand_sequence[
                                    start_position:position
                                ]
                            ).translate()
                        )

                        if len(protein_sequence) < minimum_length:
                            continue

                        if strand == "+":

                            start_coordinate = start_position + 1
                            end_coordinate = position + 3

                        else:

                            start_coordinate = (
                                sequence_length - (position + 3) + 1
                            )

                            end_coordinate = (
                                sequence_length - start_position
                            )

                        orfs.append({

                            "start": min(
                                start_coordinate,
                                end_coordinate
                            ),

                            "end": max(
                                start_coordinate,
                                end_coordinate
                            ),

                            "strand": strand,

                            "frame": frame + 1,

                            "length": len(
                                protein_sequence
                            ),

                            "protein_sequence": protein_sequence,

                            "sequence": strand_sequence[
                                start_position:position + 3
                            ]

                        })

                    starts = []

    orfs.sort(
        key=lambda item: (
            -item["length"],
            item["start"]
        )
    )

    for number, orf in enumerate(orfs, start=1):
        orf["orf_number"] = number

    # Determine overlap relationships.

    for orf in orfs:

        orf["overlaps"] = []

        for other in orfs:

            if orf["orf_number"] == other["orf_number"]:
                continue

            if (
                orf["start"] <= other["end"]
                and other["start"] <= orf["end"]
            ):

                orf["overlaps"].append(
                    other["orf_number"]
                )

    return orfs


# ============================================================
# ORF VISUALIZATION TRACKS
# ============================================================

def build_orf_tracks(predicted_orfs, sequence_length):

    track_names = [
        "+1", "+2", "+3",
        "-1", "-2", "-3"
    ]

    tracks = []

    for track_name in track_names:

        strand = track_name[0]
        frame = int(track_name[1])

        matching_orfs = [
            orf for orf in predicted_orfs
            if (
                orf["strand"] == strand
                and orf["frame"] == frame
            )
        ]

        matching_orfs.sort(
            key=lambda item: item["start"]
        )

        lanes = []
        lane_end_positions = []

        for orf in matching_orfs:

            assigned_lane = None

            for lane_index, lane_end in enumerate(
                lane_end_positions
            ):

                if orf["start"] > lane_end:

                    assigned_lane = lane_index
                    break

            if assigned_lane is None:

                assigned_lane = len(lanes)

                lanes.append([])
                lane_end_positions.append(0)

            lanes[assigned_lane].append(orf)

            lane_end_positions[assigned_lane] = orf["end"]

        tracks.append({

            "name": track_name,

            "strand": strand,

            "frame": frame,

            "lanes": lanes

        })

    return tracks


# ============================================================
# FASTA DOWNLOAD HELPER
# ============================================================

def fasta_download(filename, header, sequence):

    fasta_text = f">{header}\n"

    for index in range(0, len(sequence), 70):
        fasta_text += sequence[index:index + 70] + "\n"

    buffer = BytesIO(
        fasta_text.encode("utf-8")
    )

    buffer.seek(0)

    return send_file(
        buffer,
        mimetype="text/plain",
        as_attachment=True,
        download_name=filename
    )


# ============================================================
# GENOME EXPLORER PREPARATION
# ============================================================

def prepare_genome_explorer(
    accession,
    search="",
    page=1
):

    accession = accession.strip().upper()

    search = search.strip()

    page = max(1, int(page))

    # annotate_genome uses an LRU cache.
    # Moving between pages for the same accession
    # reuses the cached result while Flask remains running.

    genome_data = annotate_genome(accession)

    all_records = genome_data.get(
        "cds_records",
        []
    )

    if search:

        search_lower = search.lower()

        filtered_records = []

        for record in all_records:

            searchable_fields = [

                record.get("gene", ""),

                record.get("locus_tag", ""),

                record.get("product", ""),

                record.get("protein_id", "")

            ]

            searchable_text = " ".join(
                str(value)
                for value in searchable_fields
            ).lower()

            if search_lower in searchable_text:
                filtered_records.append(record)

    else:

        filtered_records = all_records

    total_records = len(filtered_records)

    total_pages = max(
        1,
        math.ceil(total_records / CDS_PER_PAGE)
    )

    page = min(page, total_pages)

    start_index = (
        page - 1
    ) * CDS_PER_PAGE

    end_index = start_index + CDS_PER_PAGE

    page_records = filtered_records[
        start_index:end_index
    ]

    # Do not send all protein sequences to the browser.
    # Only send the records required for this page.

    display_records = []

    for record in page_records:

        display_records.append({

            key: value

            for key, value in record.items()

            if key != "protein_sequence"

        })

    summary = {
        key: value
        for key, value in genome_data.items()
        if key != "cds_records"
    }

    return {

        "genome_summary": summary,

        "cds_records": display_records,

        "search_term": search,

        "total_records": total_records,

        "page": page,

        "total_pages": total_pages

    }


# ============================================================
# HOME
# ============================================================

@app.route("/")
def home():

    return render_template("index.html")


# ============================================================
# GENOMICS WORKSPACE
# ============================================================

@app.route("/genomics", methods=["GET", "POST"])
def genomics():

    if request.method == "GET":

        return render_template(
            "genomics.html"
        )

    sequence_input = request.form.get(
        "sequence",
        ""
    )

    uploaded_file = request.files.get(
        "file"
    )

    try:

        sequence = read_dna_input(
            sequence_input,
            uploaded_file
        )

        analysis = analyze_sequence(sequence)

        predicted_orfs = find_orfs(sequence)

        translations = translate_six_frames(sequence)

        orf_tracks = build_orf_tracks(
            predicted_orfs,
            len(sequence)
        )

        reverse_complement = str(
            Seq(sequence).reverse_complement()
        )

        return render_template(

            "genomics.html",

            analysis=analysis,

            sequence=sequence,

            sequence_input=sequence_input,

            predicted_orfs=predicted_orfs,

            predicted_orf_count=len(predicted_orfs),

            translations=translations,

            orf_tracks=orf_tracks,

            sequence_length=len(sequence),

            reverse_complement=reverse_complement

        )

    except Exception as error:

        return render_template(
            "genomics.html",
            error=str(error),
            sequence_input=sequence_input
        )


# ============================================================
# NCBI ANNOTATION WORKSPACE
# ============================================================

@app.route("/annotation", methods=["GET", "POST"])
def annotation():

    # --------------------------------------------------------
    # GET REQUEST: GENOME EXPLORER
    # --------------------------------------------------------

    if request.method == "GET":

        analysis_mode = request.args.get(
            "analysis_mode",
            "protein_comparison"
        )

        if analysis_mode == "genome_explorer":

            accession = request.args.get(
                "accession",
                ""
            )

            search = request.args.get(
                "search",
                ""
            )

            try:

                page = int(
                    request.args.get("page", 1)
                )

                if page < 1:
                    page = 1

            except (ValueError, TypeError):

                page = 1

            if not accession:

                return render_template(
                    "annotation.html",
                    analysis_mode="genome_explorer",
                    accession="",
                    genome_summary=None,
                    cds_records=[],
                    search_term=search,
                    total_records=0,
                    page=1,
                    total_pages=1
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
                    **explorer_data
                )

            except Exception as error:

                return render_template(
                    "annotation.html",
                    analysis_mode="genome_explorer",
                    accession=accession,
                    error=str(error),
                    genome_summary=None,
                    cds_records=[],
                    search_term=search,
                    total_records=0,
                    page=1,
                    total_pages=1
                )

        return render_template(
            "annotation.html",
            analysis_mode="protein_comparison",
            accession="U91966.1"
        )

    # --------------------------------------------------------
    # POST REQUEST: PROTEIN COMPARISON
    # --------------------------------------------------------

    analysis_mode = request.form.get(
        "analysis_mode",
        "protein_comparison"
    )

    if analysis_mode != "protein_comparison":

        return render_template(
            "annotation.html",
            error="Unsupported annotation analysis mode."
        )

    accession = request.form.get(
        "accession",
        "U91966.1"
    ).strip()

    sequence_input = request.form.get(
        "sequence",
        ""
    )

    uploaded_file = request.files.get(
        "file"
    )

    try:

        sequence = read_dna_input(
            sequence_input,
            uploaded_file
        )

        sequence_analysis = analyze_sequence(
            sequence
        )

        predicted_orfs = find_orfs(
            sequence
        )

        orf_tracks = build_orf_tracks(
            predicted_orfs,
            len(sequence)
        )

        # Fetch reference annotations and compare proteins.

        results = annotate_orfs(
            accession,
            predicted_orfs
        )

        # Attach the predicted ORF details to each match.

        orf_lookup = {
            orf["orf_number"]: orf
            for orf in predicted_orfs
        }

        for match in results.get("matches", []):

            orf_number = match.get(
                "orf_number"
            )

            matching_orf = orf_lookup.get(
                orf_number,
                {}
            )

            match["start"] = matching_orf.get(
                "start"
            )

            match["end"] = matching_orf.get(
                "end"
            )

            match["strand"] = matching_orf.get(
                "strand"
            )

            match["frame"] = matching_orf.get(
                "frame"
            )

            match["orf_length"] = matching_orf.get(
                "length"
            )

            match["overlaps"] = matching_orf.get(
                "overlaps",
                []
            )

        best_orf = None

        if results.get("matches"):

            best_orf = max(
                results["matches"],
                key=lambda item: (
                    min(
                        item.get("query_coverage", 0),
                        item.get("reference_coverage", 0)
                    ),
                    item.get("identity", 0),
                    item.get("score", 0)
                )
            )

        return render_template(

            "annotation.html",

            analysis_mode="protein_comparison",

            accession=accession,

            sequence_input=sequence_input,

            sequence=sequence,

            sequence_analysis=sequence_analysis,

            predicted_orfs=predicted_orfs,

            predicted_orf_count=len(predicted_orfs),

            orf_tracks=orf_tracks,

            sequence_length=len(sequence),

            results=results,

            best_orf_number=(
                best_orf.get("orf_number")
                if best_orf
                else None
            )

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
        "reverse_complement.fasta",
        "Reverse_Complement",
        result
    )


@app.route("/download-protein-fasta", methods=["POST"])
def download_protein_fasta():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    protein = str(
        Seq(sequence).translate(
            to_stop=False
        )
    )

    return fasta_download(
        "translated_protein.fasta",
        "Translated_Protein",
        protein
    )


@app.route("/download-longest-orf", methods=["POST"])
def download_longest_orf():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    predicted_orfs = find_orfs(sequence)

    if not predicted_orfs:

        return "No ORFs were detected.", 400

    longest_orf = predicted_orfs[0]

    return fasta_download(
        "longest_orf.fasta",
        f"ORF_{longest_orf['orf_number']}",
        longest_orf["protein_sequence"]
    )


# ============================================================
# TRANSCRIPTOMICS WORKSPACE
# ============================================================

@app.route("/transcriptomics", methods=["GET", "POST"])
def transcriptomics():

    if request.method == "GET":

        return render_template(
            "transcriptomics.html"
        )

    sequence_input = request.form.get(
        "sequence",
        ""
    )

    uploaded_file = request.files.get(
        "file"
    )

    try:

        sequence = read_rna_input(
            sequence_input,
            uploaded_file
        )

        analysis = analyze_rna_sequence(
            sequence
        )

        dna_equivalent = convert_rna_to_dna(
            sequence
        )

        protein = translate_rna(
            sequence
        )

        complement = str(
            Seq(sequence).complement_rna()
        )

        reverse_complement = str(
            Seq(sequence).reverse_complement_rna()
        )

        return render_template(

            "transcriptomics.html",

            sequence=sequence,

            sequence_input=sequence_input,

            analysis=analysis,

            dna_equivalent=dna_equivalent,

            protein=protein,

            complement=complement,

            reverse_complement=reverse_complement

        )

    except Exception as error:

        return render_template(
            "transcriptomics.html",
            error=str(error),
            sequence_input=sequence_input
        )


# ============================================================
# TRANSCRIPTOMICS DOWNLOAD ROUTES
# ============================================================

@app.route("/download-rna-complement", methods=["POST"])
def download_rna_complement():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    complement = str(
        Seq(sequence).complement_rna()
    )

    return fasta_download(
        "rna_complement.fasta",
        "RNA_Complement",
        complement
    )


@app.route("/download-rna-reverse-complement", methods=["POST"])
def download_rna_reverse_complement():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    reverse_complement = str(
        Seq(sequence).reverse_complement_rna()
    )

    return fasta_download(
        "rna_reverse_complement.fasta",
        "RNA_Reverse_Complement",
        reverse_complement
    )


@app.route("/download-rna-dna", methods=["POST"])
def download_rna_dna():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    dna = convert_rna_to_dna(
        sequence
    )

    return fasta_download(
        "rna_dna_equivalent.fasta",
        "DNA_Equivalent",
        dna
    )


@app.route("/download-rna-protein", methods=["POST"])
def download_rna_protein():

    sequence = clean_rna_sequence(
        request.form.get("sequence", "")
    )

    protein = translate_rna(
        sequence
    )

    return fasta_download(
        "rna_protein.fasta",
        "Translated_Protein",
        protein
    )


# ============================================================
# OMICS INTEGRATION WORKSPACE
# ============================================================

@app.route("/integration", methods=["GET", "POST"])
def integration():

    if request.method == "GET":

        return render_template(
            "integration.html"
        )

    sequence_input = request.form.get(
        "sequence",
        ""
    )

    uploaded_file = request.files.get(
        "file"
    )

    try:

        sequence = read_dna_input(
            sequence_input,
            uploaded_file
        )

        results = analyze_central_dogma(
            sequence
        )

        return render_template(

            "integration.html",

            sequence=sequence,

            sequence_input=sequence_input,

            results=results

        )

    except Exception as error:

        return render_template(
            "integration.html",
            error=str(error),
            sequence_input=sequence_input
        )


# ============================================================
# OMICS INTEGRATION DOWNLOAD ROUTES
# ============================================================

@app.route("/download-integration-dna", methods=["POST"])
def download_integration_dna():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    return fasta_download(
        "integration_dna.fasta",
        "DNA_Sequence",
        sequence
    )


@app.route("/download-integration-rna", methods=["POST"])
def download_integration_rna():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    rna = sequence.replace("T", "U")

    return fasta_download(
        "integration_mrna.fasta",
        "mRNA_Sequence",
        rna
    )


@app.route("/download-integration-protein", methods=["POST"])
def download_integration_protein():

    sequence = clean_dna_sequence(
        request.form.get("sequence", "")
    )

    protein = str(
        Seq(sequence).translate(
            to_stop=False
        )
    )

    return fasta_download(
        "integration_protein.fasta",
        "Translated_Protein",
        protein
    )


# ============================================================
# FILE UPLOAD ERROR HANDLING
# ============================================================

@app.errorhandler(413)
def request_entity_too_large(error):

    return (
        "The uploaded file is too large. "
        "The maximum permitted upload size is 5 MB.",
        413
    )


# ============================================================
# APPLICATION ENTRY POINT
# ============================================================

if __name__ == "__main__":

    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True,
        use_reloader=False
    )
