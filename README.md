# Bioinformatics Python Toolkit 🧬

A collection of Python scripts demonstrating core bioinformatics workflows, structural biology parsing, and sequence analysis. This toolkit leverages **Biopython** and **NumPy** to automate data retrieval, perform sequence alignments, and analyze biological data.

## 🚀 What's Inside (Core Functionalities)

This repository contains standalone scripts that perform the following tasks:

### 1. Structural Biology & Sequence Alignment
* **SMCRA Parsing & Sequence Extraction:** Downloads and parses mmCIF/PDB files to navigate the Structure-Model-Chain-Residue-Atom (SMCRA) architecture. Extracts 1D amino acid sequences directly from 3D structures.
* **Biophysical Calculations:** Uses spatial coordinates to calculate metrics like the 3D Center of Mass and identifies highly flexible regions using atomic B-factors.
* **Custom Pairwise Alignment:** Performs local sequence alignments (Smith-Waterman style) with customizable match, mismatch, and gap penalty parameters. 

### 2. Database Fetching & BLAST Automation
* **NCBI Entrez Integration:** Interacts with the Entrez API to search and fetch nucleotide or protein sequences in batches, handling pagination for large datasets.
* **Remote BLAST pipeline:** Automates remote BLAST queries via NCBIWWW, parses the resulting XML outputs, filters hits based on custom E-value thresholds, and saves significant hits into ready-to-use FASTA files.
* **BLAST XML Analysis:** Reads BLAST XML outputs to extract pairwise alignments (HSPs) and generates summary statistics such as average identity, maximum score, and best E-values.

### 3. File Parsing & Data Filtering
* **Advanced FASTA Processing:** Reads FASTA files to calculate summary statistics (max/min/avg sequence lengths), computes GC content, and filters sequences based on user-defined length thresholds.
* **GenBank Metadata Extraction:** Parses GenBank (.gb/.gbk) files to extract metadata, summarize feature counts, and specifically isolate Coding Sequences (CDS) and protein IDs.

## 🛠️ Prerequisites & Installation

To run these scripts locally, you need Python 3.x installed along with the following libraries:

```bash
pip install biopython numpy
