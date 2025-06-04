# assayBLAST - for *in silico* analysis of PCR oligos

This tool provides in silico predictions of microarray hybridization results, calculating the expected binding interactions between query DNA sequences (primers and probes) and a genome database. Based on the BLAST hits and their mismatch numbers, the functionality of the primers and probes can be estimated. The tool also checks the strand specificity of the primers and probes.

The tool has been evaluated and was published at MDPI applied biosciences (https://doi.org/10.3390/applbiosci4020018)


## Key Features
- BLAST Database Creation: Automatically generates a custom BLAST database from user-provided genome files.
- BLAST Search Execution: Runs BLAST search to evaluate how primers or probes interact with both DNA strands.
- Mismatch Analysis: Identifies and filters matches with acceptable levels of mismatches, simulating how primers perform in the presence of mutations or single nucleotide polymorphisms (SNPs).
- Off-Target Detection: Helps identify off-target binding sites to ensure primer/probe specificity.
- Support for Multiple Genomes: Capable of performing searches across multiple genome sequences, ideal for comparative genomics or multi-strain analyses.
- Results in Multiple Formats: Outputs detailed results in TSV, and text formats for easy review and further analysis.

## Python Requirements
- Python 3.11+
- rnajena-sugar
- NCBI BLAST+

## Usage
Version 2.0 comes with two separate scripts: one for running BLAST and another for analyzing the results.
To run the tool, download the assay_blast.py and assay_analyze.py files and run it with python as explained below.

### Command-Line Arguments

```bash
python assay_blast.py <genome_files_glob_pattern> -q <query_file.fasta> -o <BLAST output file> [options]
python assay_analyze.py <BLAST output file>
```
#### Required Arguments for assay_blast.py
- `genomes`: Glob pattern for the genome FASTA or GenBank files.
- `-q, --queries`: Path to the query FASTA file containing primers or probes.

For a description of optional arguments please run `python assay_blast.py -h`.

#### Required Arguments for assay_analyze.py

- `fname`: BLAST output file from `assay_blast.py`

For a description of optional arguments please run `python assay_analyze.py -h`.

#### Example

```bash
python assay_blast.py example_database.fasta -q example_queries.fasta --max_mismatches 2
python assay_analyze.py blast_results.tsv --max_mismatches 2
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
2. **Forward and Reverse Complement BLAST**: Runs two separate BLAST searches: one for the forward strand and one for the reverse complement of the query sequences.
3. **Mismatch Filtering and Analysis**: Alignments are filtered based on user-defined mismatch thresholds. Results include mismatch counts and binding positions.
4. **Result Generation**: Outputs the results in various formats including TSV (for easy parsing and analysis), and text (for alignments).

## Run tests

Run the tests with `python test_assay.py`.

## Contributing
Feel free to contribute to this project by submitting pull requests, reporting issues, or suggesting improvements.

## License
This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact
For any questions, issues, or suggestions, please reach out via [GitHub Issues](https://github.com/mcollatz/AssayBLAST/issues).
