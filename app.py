
from flask import (
    Flask,
    render_template,
    request,
    send_file
)

from Bio import SeqIO
from Bio.Seq import Seq

from io import StringIO, BytesIO

from genomics.sequence_analysis import analyze_sequence
from genomics.ncbi_annotation import annotate_orfs

from genomics.transcriptomics_analysis import (
    analyze_rna_sequence,
    convert_rna_to_dna,
    translate_rna
)

from genomics.omics_integration import analyze_central_dogma


app = Flask(__name__)

app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024


# ==========================================
# CLEAN AND VALIDATE DNA SEQUENCE
# ==========================================

def clean_dna_sequence(raw_sequence):

    if not raw_sequence:
        raise ValueError("Please provide a DNA sequence.")

    raw_sequence = raw_sequence.strip()

    lines = [
        line.strip()
        for line in raw_sequence.splitlines()
        if line.strip() and not line.strip().startswith(">")
    ]

    sequence = "".join(lines).upper()
    sequence = "".join(sequence.split())

    if not sequence:
        raise ValueError("No DNA sequence was found.")

    allowed_bases = set("ATGCN")
    invalid_bases = set(sequence) - allowed_bases

    if invalid_bases:
        raise ValueError(
            "Invalid DNA characters detected: "
            + ", ".join(sorted(invalid_bases))
            + ". Use only A, T, G, C and N."
        )

    return sequence


# ==========================================
# READ DNA FROM TEXT OR FASTA FILE
# ==========================================

def read_dna_input(pasted_sequence, uploaded_file):

    if uploaded_file and uploaded_file.filename:

        try:
            file_content = uploaded_file.read().decode("utf-8")

        except UnicodeDecodeError:
            raise ValueError(
                "The uploaded file must be a text-based FASTA file."
            )

        if file_content.lstrip().startswith(">"):

            records = list(
                SeqIO.parse(StringIO(file_content), "fasta")
            )

            if not records:
                raise ValueError(
                    "No valid FASTA sequence was found."
                )

            raw_sequence = str(records[0].seq)

        else:
            raw_sequence = file_content

    else:
        raw_sequence = pasted_sequence

    return clean_dna_sequence(raw_sequence)


# ==========================================
# CLEAN AND VALIDATE RNA SEQUENCE
# ==========================================

def clean_rna_sequence(raw_sequence):

    if not raw_sequence:
        raise ValueError("Please provide an RNA sequence.")

    raw_sequence = raw_sequence.strip()

    lines = [
        line.strip()
        for line in raw_sequence.splitlines()
        if line.strip() and not line.strip().startswith(">")
    ]

    sequence = "".join(lines).upper()
    sequence = "".join(sequence.split())

    if not sequence:
        raise ValueError("No RNA sequence was found.")

    allowed_bases = set("AUGCN")
    invalid_bases = set(sequence) - allowed_bases

    if invalid_bases:
        raise ValueError(
            "Invalid RNA characters detected: "
            + ", ".join(sorted(invalid_bases))
            + ". Use only A, U, G, C and N."
        )

    return sequence


# ==========================================
# READ RNA FROM TEXT OR FASTA FILE
# ==========================================

def read_rna_input(pasted_sequence, uploaded_file):

    if uploaded_file and uploaded_file.filename:

        try:
            file_content = uploaded_file.read().decode("utf-8")

        except UnicodeDecodeError:
            raise ValueError(
                "The uploaded file must be a text-based FASTA file."
            )

        if file_content.lstrip().startswith(">"):

            records = list(
                SeqIO.parse(StringIO(file_content), "fasta")
            )

            if not records:
                raise ValueError(
                    "No valid FASTA sequence was found."
                )

            raw_sequence = str(records[0].seq)

        else:
            raw_sequence = file_content

    else:
        raw_sequence = pasted_sequence

    return clean_rna_sequence(raw_sequence)


# ==========================================
# SIX-FRAME TRANSLATION
# ==========================================

def translate_six_frames(sequence):

    dna = Seq(sequence)
    reverse_dna = dna.reverse_complement()

    translations = []

    strands = [
        ("+", dna),
        ("-", reverse_dna)
    ]

    for strand, strand_sequence in strands:

        for frame in range(3):

            protein = strand_sequence[frame:].translate(
                to_stop=False
            )

            translations.append({
                "strand": strand,
                "frame": frame + 1,
                "protein_sequence": str(protein)
            })

    return translations


# ==========================================
# IMPROVED ORF DETECTION AND OVERLAP ANALYSIS
# ==========================================

def find_orfs(sequence, minimum_length=30):

    sequence = clean_dna_sequence(sequence)

    dna = Seq(sequence)
    reverse_dna = dna.reverse_complement()

    stop_codons = {"TAA", "TAG", "TGA"}

    orfs = []
    sequence_length = len(sequence)

    strands = [
        ("+", str(dna)),
        ("-", str(reverse_dna))
    ]

    # ======================================
    # STEP 1: PREDICT ORFs ON BOTH STRANDS
    # ======================================

    for strand, strand_sequence in strands:

        for frame in range(3):

            for start in range(
                frame,
                len(strand_sequence) - 2,
                3
            ):

                if strand_sequence[start:start + 3] != "ATG":
                    continue

                for stop in range(
                    start + 3,
                    len(strand_sequence) - 2,
                    3
                ):

                    codon = strand_sequence[stop:stop + 3]

                    if codon in stop_codons:

                        dna_sequence = strand_sequence[
                            start:stop + 3
                        ]

                        protein_sequence = str(
                            Seq(dna_sequence).translate(
                                to_stop=True
                            )
                        )

                        protein_length = len(protein_sequence)

                        if protein_length >= minimum_length:

                            if strand == "+":

                                start_position = start + 1
                                end_position = stop + 3

                            else:

                                start_position = (
                                    sequence_length - (stop + 3) + 1
                                )

                                end_position = (
                                    sequence_length - start
                                )

                            orfs.append({

                                "strand": strand,
                                "frame": frame + 1,
                                "start": start_position,
                                "end": end_position,
                                "length": protein_length,
                                "dna_sequence": dna_sequence,
                                "protein_sequence": protein_sequence

                            })

                        # Stop at the first in-frame stop codon.
                        break

    # ======================================
    # STEP 2: SORT AND NUMBER ORFs
    # ======================================

    orfs.sort(
        key=lambda item: item["length"],
        reverse=True
    )

    for index, orf in enumerate(orfs):

        orf["number"] = index + 1

        orf["overlaps_with"] = []
        orf["nested_in"] = []
        orf["contains"] = []

        orf["overlap_status"] = "No overlap detected"

    # ======================================
    # STEP 3: COMPARE ALL GENOMIC INTERVALS
    # ======================================

    for index, first in enumerate(orfs):

        for second in orfs[index + 1:]:

            overlap = (
                first["start"] <= second["end"]
                and second["start"] <= first["end"]
            )

            if not overlap:
                continue

            first["overlaps_with"].append(
                second["number"]
            )

            second["overlaps_with"].append(
                first["number"]
            )

            # ==================================
            # STEP 4: IDENTIFY STRICT CONTAINMENT
            # ==================================

            first_contains_second = (
                first["start"] <= second["start"]
                and first["end"] >= second["end"]
            )

            second_contains_first = (
                second["start"] <= first["start"]
                and second["end"] >= first["end"]
            )

            if (
                first_contains_second
                and not second_contains_first
            ):

                first["contains"].append(
                    second["number"]
                )

                second["nested_in"].append(
                    first["number"]
                )

            elif (
                second_contains_first
                and not first_contains_second
            ):

                second["contains"].append(
                    first["number"]
                )

                first["nested_in"].append(
                    second["number"]
                )

            # Identical intervals are recorded as overlaps,
            # but not as strict containment.

    # ======================================
    # STEP 5: CLASSIFY OVERLAP STATUS
    # ======================================

    for orf in orfs:

        if orf["nested_in"]:

            orf["overlap_status"] = "Nested ORF"

        elif orf["contains"]:

            orf["overlap_status"] = "Contains nested ORF"

        elif orf["overlaps_with"]:

            orf["overlap_status"] = "Overlapping ORF"

        else:

            orf["overlap_status"] = "No overlap detected"

        orf["overlaps_with"].sort()
        orf["nested_in"].sort()
        orf["contains"].sort()

    return orfs


# ==========================================
# BUILD SIX ORF VISUALIZATION TRACKS
# ==========================================

def build_orf_tracks(predicted_orfs, sequence_length):

    tracks = []

    if sequence_length <= 0:
        return tracks

    ticks = sorted(set([
        1,
        max(1, round(sequence_length * 0.25)),
        max(1, round(sequence_length * 0.50)),
        max(1, round(sequence_length * 0.75)),
        sequence_length
    ]))

    for strand in ["+", "-"]:

        for frame in range(1, 4):

            frame_orfs = [
                dict(orf)
                for orf in predicted_orfs
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

                lane = None

                for index, last_end in enumerate(lane_ends):

                    if orf["start"] > last_end:
                        lane = index
                        break

                if lane is None:

                    lane = len(lane_ends)
                    lane_ends.append(orf["end"])

                else:

                    lane_ends[lane] = orf["end"]

                orf["lane"] = lane

            tracks.append({

                "strand": strand,
                "frame": frame,
                "orfs": frame_orfs,
                "lane_count": max(1, len(lane_ends)),
                "ticks": ticks

            })

    return tracks


# ==========================================
# CREATE FASTA DOWNLOAD
# ==========================================

def fasta_download(filename, header, sequence):

    if not sequence:
        return "No sequence is available to download.", 400

    fasta_content = f">{header}\n{sequence}\n"

    file_data = BytesIO(
        fasta_content.encode("utf-8")
    )

    return send_file(
        file_data,
        mimetype="text/plain",
        as_attachment=True,
        download_name=filename
    )


# ==========================================
# HOME PAGE
# ==========================================

@app.route("/")
def home():

    return render_template("index.html")


# ==========================================
# GENOMICS WORKSPACE
# ==========================================

@app.route("/genomics", methods=["GET", "POST"])
def genomics():

    results = None
    error = None
    sequence_input = ""
    translations = []
    predicted_orfs = []
    orf_tracks = []
    sequence_length = 0

    if request.method == "POST":

        try:

            uploaded_file = (
                request.files.get("file")
                or request.files.get("fasta_file")
            )

            sequence = read_dna_input(
                request.form.get("sequence", ""),
                uploaded_file
            )

            results = analyze_sequence(sequence)

            sequence_input = sequence
            sequence_length = len(sequence)

            translations = translate_six_frames(sequence)

            predicted_orfs = find_orfs(sequence)

            orf_tracks = build_orf_tracks(
                predicted_orfs,
                sequence_length
            )

        except Exception as e:

            error = str(e)

    return render_template(
        "genomics.html",
        results=results,
        error=error,
        sequence_input=sequence_input,
        translations=translations,
        predicted_orfs=predicted_orfs,
        orf_tracks=orf_tracks,
        sequence_length=sequence_length
    )


# ==========================================
# NCBI PROTEIN ANNOTATION
# ==========================================

@app.route("/annotation", methods=["GET", "POST"])
def annotation():

    results = None
    error = None
    accession = "U91966.1"
    predicted_orfs = []
    predicted_orf_count = 0
    sequence_input = ""
    orf_tracks = []
    sequence_length = 0
    best_orf_number = None

    if request.method == "POST":

        try:

            accession = request.form.get(
                "accession",
                "U91966.1"
            ).strip()

            uploaded_file = (
                request.files.get("file")
                or request.files.get("fasta_file")
            )

            sequence = read_dna_input(
                request.form.get("sequence", ""),
                uploaded_file
            )

            sequence_input = sequence
            sequence_length = len(sequence)

            analyze_sequence(sequence)

            # Predict ORFs and identify genomic overlaps.

            predicted_orfs = find_orfs(sequence)

            predicted_orf_count = len(predicted_orfs)

            if not predicted_orfs:

                raise ValueError(
                    "No ORFs of at least 30 amino acids were detected."
                )

            # Build the six genomic visualization tracks.

            orf_tracks = build_orf_tracks(
                predicted_orfs,
                sequence_length
            )

            print("Predicted ORFs:", predicted_orf_count)

            # Retrieve NCBI annotations and compare proteins.

            results = annotate_orfs(
                accession,
                predicted_orfs
            )

            results["predicted_orf_count"] = predicted_orf_count

            # Attach coordinates and overlap data to each match.

            for match_item in results["matches"]:

                orf_number = match_item.get("orf_number")

                matching_orf = next(
                    (
                        orf for orf in predicted_orfs
                        if orf["number"] == orf_number
                    ),
                    None
                )

                if matching_orf:

                    match_item["orf_start"] = matching_orf["start"]
                    match_item["orf_end"] = matching_orf["end"]

                    match_item["overlap_status"] = (
                        matching_orf["overlap_status"]
                    )

                    match_item["overlaps_with"] = (
                        matching_orf["overlaps_with"]
                    )

                    match_item["nested_in"] = (
                        matching_orf["nested_in"]
                    )

                    match_item["contains"] = (
                        matching_orf["contains"]
                    )

                    match_item["orf_dna_sequence"] = (
                        matching_orf["dna_sequence"]
                    )

                    match_item["orf_protein_sequence"] = (
                        matching_orf["protein_sequence"]
                    )

            # Select the best match using the same
            # ranking principle as the interpretation dashboard.

            if results["matches"]:

                best_match = max(
                    results["matches"],
                    key=lambda item: (
                        min(
                            item["match"]["query_coverage"],
                            item["match"]["reference_coverage"]
                        ),
                        item["match"]["identity"],
                        item["match"]["alignment_score"]
                    )
                )

                best_orf_number = best_match["orf_number"]

            print("NCBI annotation completed.")

            print(
                "Annotated matches:",
                len(results["matches"])
            )

            print(
                "Best matching ORF:",
                best_orf_number
            )

        except Exception as e:

            error = str(e)
            results = None

            print("ANNOTATION ERROR:", error)

    return render_template(
        "annotation.html",
        results=results,
        error=error,
        accession=accession,
        sequence_input=sequence_input,
        predicted_orfs=predicted_orfs,
        predicted_orf_count=predicted_orf_count,
        orf_tracks=orf_tracks,
        sequence_length=sequence_length,
        best_orf_number=best_orf_number
    )


# ==========================================
# DOWNLOAD DNA REVERSE COMPLEMENT
# ==========================================

@app.route("/download-reverse-complement", methods=["POST"])
def download_reverse_complement():

    try:

        sequence = clean_dna_sequence(
            request.form.get("sequence", "")
        )

        reverse_complement = str(
            Seq(sequence).reverse_complement()
        )

        return fasta_download(
            "reverse_complement.fasta",
            "OmicsBridge_Reverse_Complement",
            reverse_complement
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# DOWNLOAD PROTEIN FASTA
# ==========================================

@app.route("/download-protein-fasta", methods=["POST"])
def download_protein_fasta():

    protein_sequence = request.form.get(
        "protein_sequence",
        request.form.get("sequence", "")
    )

    return fasta_download(
        "protein_sequence.fasta",
        "OmicsBridge_Protein",
        protein_sequence
    )


# ==========================================
# DOWNLOAD LONGEST ORF
# ==========================================

@app.route("/download-longest-orf", methods=["POST"])
def download_longest_orf():

    try:

        sequence = clean_dna_sequence(
            request.form.get("sequence", "")
        )

        predicted_orfs = find_orfs(sequence)

        if not predicted_orfs:
            return "No ORFs were detected.", 404

        longest_orf = predicted_orfs[0]

        return fasta_download(
            "longest_orf.fasta",
            f"OmicsBridge_ORF_{longest_orf['number']}",
            longest_orf["protein_sequence"]
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# TRANSCRIPTOMICS WORKSPACE
# ==========================================

@app.route("/transcriptomics", methods=["GET", "POST"])
def transcriptomics():

    results = None
    error = None
    sequence_input = ""
    dna_sequence = ""
    protein_sequence = ""

    if request.method == "POST":

        print("TRANSCRIPTOMICS POST REQUEST RECEIVED")

        try:

            uploaded_file = (
                request.files.get("file")
                or request.files.get("fasta_file")
            )

            sequence = read_rna_input(
                request.form.get("sequence", ""),
                uploaded_file
            )

            results = analyze_rna_sequence(sequence)

            dna_sequence = convert_rna_to_dna(sequence)

            sequence_input = sequence

            if len(sequence) % 3 == 0:
                protein_sequence = translate_rna(sequence)

        except Exception as e:

            error = str(e)

            print("TRANSCRIPTOMICS ERROR:", error)

    return render_template(
        "transcriptomics.html",
        results=results,
        error=error,
        sequence_input=sequence_input,
        dna_sequence=dna_sequence,
        protein_sequence=protein_sequence
    )


# ==========================================
# TRANSCRIPTOMICS DOWNLOADS
# ==========================================

@app.route("/download-rna-complement", methods=["POST"])
def download_rna_complement():

    try:

        sequence = clean_rna_sequence(
            request.form.get("sequence", "")
        )

        results = analyze_rna_sequence(sequence)

        return fasta_download(
            "rna_complement.fasta",
            "OmicsBridge_RNA_Complement",
            results["complementary_sequence"]
        )

    except Exception as e:

        return str(e), 400


@app.route("/download-rna-reverse-complement", methods=["POST"])
def download_rna_reverse_complement():

    try:

        sequence = clean_rna_sequence(
            request.form.get("sequence", "")
        )

        results = analyze_rna_sequence(sequence)

        return fasta_download(
            "rna_reverse_complement.fasta",
            "OmicsBridge_RNA_Reverse_Complement",
            results["reverse_complement"]
        )

    except Exception as e:

        return str(e), 400


@app.route("/download-rna-dna", methods=["POST"])
def download_rna_dna():

    try:

        sequence = clean_rna_sequence(
            request.form.get("sequence", "")
        )

        dna_sequence = convert_rna_to_dna(sequence)

        return fasta_download(
            "rna_to_dna.fasta",
            "OmicsBridge_DNA_From_RNA",
            dna_sequence
        )

    except Exception as e:

        return str(e), 400


@app.route("/download-rna-protein", methods=["POST"])
def download_rna_protein():

    try:

        sequence = clean_rna_sequence(
            request.form.get("sequence", "")
        )

        protein_sequence = translate_rna(sequence)

        return fasta_download(
            "rna_protein.fasta",
            "OmicsBridge_Translated_Protein",
            protein_sequence
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# OMICS INTEGRATION
# ==========================================

@app.route("/integration", methods=["GET", "POST"])
def integration():

    results = None
    error = None
    sequence_input = ""

    if request.method == "POST":

        try:

            uploaded_file = (
                request.files.get("file")
                or request.files.get("fasta_file")
            )

            sequence = read_dna_input(
                request.form.get("sequence", ""),
                uploaded_file
            )

            results = analyze_central_dogma(sequence)

            sequence_input = sequence

            print("OMICS INTEGRATION ANALYSIS COMPLETED.")

        except Exception as e:

            error = str(e)

            print("OMICS INTEGRATION ERROR:", error)

    return render_template(
        "integration.html",
        results=results,
        error=error,
        sequence_input=sequence_input
    )


# ==========================================
# INTEGRATION DNA FASTA DOWNLOAD
# ==========================================

@app.route("/download-integration-dna", methods=["POST"])
def download_integration_dna():

    try:

        sequence = clean_dna_sequence(
            request.form.get("sequence", "")
        )

        return fasta_download(
            "omicsbridge_dna.fasta",
            "OmicsBridge_DNA_coding_strand",
            sequence
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# INTEGRATION mRNA FASTA DOWNLOAD
# ==========================================

@app.route("/download-integration-rna", methods=["POST"])
def download_integration_rna():

    try:

        sequence = clean_rna_sequence(
            request.form.get("sequence", "")
        )

        return fasta_download(
            "omicsbridge_mRNA.fasta",
            "OmicsBridge_generated_mRNA",
            sequence
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# INTEGRATION PROTEIN FASTA DOWNLOAD
# ==========================================

@app.route("/download-integration-protein", methods=["POST"])
def download_integration_protein():

    try:

        sequence = request.form.get(
            "sequence",
            ""
        ).strip().upper()

        sequence = "".join(sequence.split())

        if not sequence:

            raise ValueError(
                "No protein sequence was provided."
            )

        allowed_amino_acids = set(
            "ACDEFGHIKLMNPQRSTVWYX*"
        )

        invalid_characters = set(sequence) - allowed_amino_acids

        if invalid_characters:

            raise ValueError(
                "Invalid protein characters detected: "
                + ", ".join(sorted(invalid_characters))
            )

        return fasta_download(
            "omicsbridge_protein.fasta",
            "OmicsBridge_translated_protein",
            sequence
        )

    except Exception as e:

        return str(e), 400


# ==========================================
# FILE SIZE ERROR HANDLER
# ==========================================

@app.errorhandler(413)
def request_entity_too_large(error):

    return (
        "The uploaded file exceeds the 5 MB limit.",
        413
    )


# ==========================================
# RUN APPLICATION
# ==========================================

if __name__ == "__main__":

    print("\n====================================")
    print("          OMICSBRIDGE ROUTES")
    print("====================================")

    for route in app.url_map.iter_rules():

        print(
            route.rule,
            sorted(route.methods)
        )

    print("====================================")
    print("Starting OmicsBridge...")
    print("====================================\n")

    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True,
        use_reloader=False
    )
