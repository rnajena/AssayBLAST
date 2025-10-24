<img src="https://raw.github.com/rnajena/AssayBLAST/main/logo/assayblast_logo_low.png" alt="logo" width="400">

# AssayBLAST - for *in silico* analysis of PCR oligos

[![build status](https://github.com/rnajena/AssayBLAST/workflows/tests/badge.svg)](https://github.com/rnajena/AssayBLAST/actions)
[![codecov](https://codecov.io/gh/rnajena/AssayBLAST/branch/main/graph/badge.svg)](https://codecov.io/gh/rnajena/AssayBLAST)
[![pypi version](https://img.shields.io/pypi/v/assay_blast.svg)](https://pypi.python.org/pypi/assay_blast)
[![python version](https://img.shields.io/pypi/pyversions/assay_blast.svg)](https://python.org)

This tool provides in silico predictions of microarray hybridization results, calculating the expected binding interactions between query DNA sequences (primers and probes) and a genome database. Based on the BLAST hits and their mismatch numbers, the functionality of the primers and probes can be estimated. The tool also checks the strand specificity of the primers and probes.

The tool has been evaluated and was published at MDPI applied biosciences (https://doi.org/10.3390/applbiosci4020018)


## Key Features
- BLAST Database Creation: Automatically generates a custom BLAST database from user-provided genome files.
- BLAST Search Execution: Runs BLAST search to evaluate how primers or probes interact with both DNA strands.
- Mismatch Analysis: Identifies and filters matches with acceptable levels of mismatches, simulating how primers perform in the presence of mutations or single nucleotide polymorphisms (SNPs).
- Off-Target Detection: Helps identify off-target binding sites to ensure primer/probe specificity.
- Support for Multiple Genomes: Capable of performing searches across multiple genome sequences, ideal for comparative genomics or multi-strain analyses.
- Results in Multiple Formats: Outputs detailed results in TSV, and text formats for easy review and further analysis.

## Requirements
- Python 3.11+
- rnajena-sugar
- NCBI BLAST+

## Installation

Install this repository alongside a BLAST installation:

```bash
pip install https://github.com/rnajena/AssayBLAST/archive/refs/heads/main.zip
```

One way to install BLAST and AssayBLAST with conda:

```bash
conda create -c bioconda -n assay_blast python=="3.13" blast
conda activate assay_blast
pip install https://github.com/rnajena/AssayBLAST/archive/refs/heads/main.zip
```


## Usage
Version 2.0 comes with two separate scripts: one for running BLAST and another for analyzing the results.

### Command-Line Arguments

```bash
assay_blast <genome_files_glob_pattern> -q <query_file.fasta> -o <BLAST output file> [options]
assay_analyze <BLAST output file>
```
#### Required Arguments for assay_blast
- `genomes`: Glob pattern for the genome FASTA or GenBank files.
- `-q, --queries`: Path to the query FASTA file containing primers or probes.

For a description of optional arguments please run `assay_blast -h`.
The `--num-threads` parameter can be used to specify the number of threads used by the BLAST search.

#### Required Arguments for assay_analyze

- `fname`: BLAST output file from `assay_blast`

For a description of optional arguments please run `assay_analyze -h`.

#### Example

```bash
assay_test -d .  # Download the two example files
assay_blast example_database.fasta -q example_queries.fasta --mismatch 2
assay_analyze blast_results.tsv --mismatch 2
```

These commands:
- Uses the FASTA files `example_database.fasta` to build the BLAST database.
- Runs the BLAST search using the primers/probes in `example_queries.fasta`.
- Allows a maximum of 2 mismatches in alignments.


## Outputs
- **BLAST TSV Results**: Detailed output of BLAST hits including scores, E-values, mismatches.
- **BLAST Mismatching Alignments**: Overview of mismatching alignments
- **TSV Overview**: Summary table containing mismatch counts and growth analysis for each probe/primer.
- **TSV Details**: Detailed table with all found probes/primers, mismatch counts, distances and positions

## How It Works
1. **BLAST Database Creation**: The tool generates a BLAST database from the user-provided genome files.
2. **Forward and Reverse Complement BLAST**: Runs two separate BLAST searches: one for the detection of the sequences and a second one to output alignments of the matches.
3. **Mismatch Filtering and Analysis**: Alignments are filtered based on user-defined mismatch thresholds. Results include mismatch counts and binding positions.
4. **Result Generation**: Outputs the results in various formats including TSV (for easy parsing and analysis), and text (for alignments).

## Run tests

Run the tests with `assay_test`.

## Contributing
Feel free to contribute to this project by submitting pull requests, reporting issues, or suggesting improvements.

## License
This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact
For any questions, issues, or suggestions, please reach out via [GitHub Issues](https://github.com/rnajena/AssayBLAST/issues).
