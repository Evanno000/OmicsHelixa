# OmicsHelixa

### Connecting Genomes, Transcriptomes and Biological Insights

**OmicsHelixa** is an integrated bioinformatics platform developed to explore biological sequence data through computational analysis.

The platform brings together genomics, transcriptomics, sequence translation, and NCBI annotation exploration in one web application.

Its initial biological focus is *Plasmodium falciparum*, the parasite responsible for malaria.

---

## About the Author

**Created by Evans Kioko**

Evans Kioko is a Bachelor of Science in Biochemistry student at Kisii University with a strong interest in bioinformatics, computational biology, genomics, and systems biology.

OmicsHelixa was developed as a learning and portfolio project to combine biochemical knowledge with programming and computational approaches to biological data analysis.

The project reflects an interest in applying computational tools to understand biological systems and explore the relationships between genomic, transcriptomic, and proteomic information.

**Developer:** Evans Kioko  
**Institution:** Kisii University  
**Academic Programme:** BSc Biochemistry  
**Project:** OmicsHelixa  
**Field:** Bioinformatics and Computational Biology

---

## Features

### 1. Genomics

- DNA sequence analysis
- Nucleotide composition
- GC content calculation
- Reverse complement generation
- Six-frame translation
- FASTA file input and output

### 2. Transcriptomics

- RNA sequence analysis
- RNA nucleotide composition
- GC content calculation
- Complement and reverse complement generation
- Translation into protein sequences

### 3. Integrated Omics Analysis

- Connects DNA, mRNA, and protein sequence analysis
- Generates corresponding molecular sequences
- Summarises sequence statistics
- Explores relationships between the central dogma components

### 4. NCBI Annotation Explorer

- Retrieves nucleotide records using NCBI accession numbers
- Explores annotated coding sequences (CDS)
- Examines annotated protein information
- Supports exploratory ORF and protein comparison
- Provides a genome explorer with paginated annotation records
- Supports searches using available annotation identifiers and descriptions

---

## Technology Stack

- Python
- Flask
- Biopython
- HTML
- CSS
- JavaScript
- NCBI Entrez databases

---

## Requirements

- Python 3.10 or later
- pip
- Internet connection for NCBI data retrieval

---

## Installation

### 1. Clone the repository

```bash
git clone YOUR_GITHUB_REPOSITORY_URL
```

Replace `YOUR_GITHUB_REPOSITORY_URL` with the actual repository URL after publishing the project.

### 2. Enter the project directory

```bash
cd OmicsHelixa
```

### 3. Create a virtual environment

```bash
python -m venv venv
```

### 4. Activate the virtual environment

**Windows PowerShell:**

```powershell
.\venv\Scripts\Activate.ps1
```

### 5. Install dependencies

```bash
python -m pip install -r requirements.txt
```

### 6. Start the application

```bash
python app_test.py
```

Open the following address in your browser:

```text
http://127.0.0.1:5000/
```

---

## Project Structure

```text
OmicsHelixa/
│
├── app.py
├── app_test.py
├── requirements.txt
├── .gitignore
│
├── genomics/
│
├── templates/
│
├── static/
│
└── test_files/
```

---

## Scientific Scope and Limitations

OmicsHelixa is an educational and exploratory bioinformatics platform.

Its computational outputs are intended for learning, sequence exploration, and preliminary analysis. They should not be interpreted as clinically validated results or definitive biological conclusions.

Protein comparison thresholds used in exploratory screening are implementation-specific and should not be treated as universal standards for establishing homology or protein function.

NCBI annotations and sequence records are provided by external databases. Their accuracy and completeness depend on the underlying source records.

---

## Acknowledgements

The project uses Biopython for biological sequence processing and NCBI resources for nucleotide and annotation data retrieval.

- [Biopython](https://biopython.org/)
- [NCBI](https://www.ncbi.nlm.nih.gov/)
- [NCBI Entrez](https://www.ncbi.nlm.nih.gov/books/NBK25501/)

---

## Project Status

OmicsHelixa is under active development. Its current focus is completing functional testing, improving documentation, and preparing the application for publication.

---

## Licence

A licence will be specified when the project licence file is added.