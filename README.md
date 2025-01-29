# assayBLAST - for *in silico* analysis of PCR oligos

This tool provides in silico predictions of microarray hybridization results, calculating the expected binding interactions between query DNA sequences (primers and probes) and a genome database. Based on the BLAST hits and their mismatch numbers, the functionality of the primers and probes can be estimated. The tool also checks the strand specificity of the primers and probes.



## Key Features
- BLAST Database Creation: Automatically generates a custom BLAST database from user-provided genome files.
- BLAST Search Execution: Runs forward and reverse complement BLAST searches to evaluate how primers or probes interact with both DNA strands.
- Mismatch Analysis: Identifies and filters matches with acceptable levels of mismatches, simulating how primers perform in the presence of mutations or single nucleotide polymorphisms (SNPs).
- Melting Temperature (Tm) Calculation: Computes the melting temperature of sequences using SantaLucia's thermodynamic nearest-neighbor model.
- Off-Target Detection: Helps identify off-target binding sites to ensure primer/probe specificity.
- Support for Multiple Genomes: Capable of performing searches across multiple genome sequences, ideal for comparative genomics or multi-strain analyses.
- Results in Multiple Formats: Outputs detailed results in XML, TSV, and text formats for easy review and further analysis.

## Python Requirements
- Python 3.7+
- Biopython
- NCBI BLAST+

## Usage
To run the tool, download the assayBLAST.py file and run it with python as explained below.
### Command-Line Arguments

```bash
python assayBLAST.py -g <"genome_files_glob_pattern.fasta"> -q <query_file.fasta> [options]
```
#### Required Arguments
- `-g, --genome`: Glob pattern for the genome FASTA files in "".
- `-q, --queries`: Path to the query FASTA file containing primers or probes.

#### Optional Arguments
- `-d, --db_name`: Name of the BLAST database (default: `genome_db`).
- `-o, --output`: Output file name for BLAST results (default: `blast_results.xml`).
- `-c, --tsv_output`: Output file for mismatch matrix in TSV format (default: `blast_results.tsv`).
- `-a, --alignments_output`: Output file for alignments with mismatches (default: `alignments.txt`).
- `-mh, --multi_hits_output`: Output file for multi-hit BLAST results (default: `multi_hits.txt`).
- `-db, --db_dir`: Directory to store BLAST database files (default: `blast_db`).
- `-m, --max_mismatches`: Maximum number of allowed mismatches in BLAST alignments (default: `4`).
- `-cc, --concatenate`: Concatenate input sequences into one (default: `False`).
- `-k, --keep_blast_db`: Keep the previously created BLAST database (default: `False`).

#### Example

```bash
python assayBLAST.py --genome "example_database.fasta" --queries example_queries.fasta --max_mismatches 2
```
This command:
- Uses the FASTA files `example_database.fasta` to build the BLAST database.
- Runs the BLAST search using the primers/probes in `example_queries.fasta`.
- Allows a maximum of 2 mismatches in alignments.
- The --concatenate parameter is not used, thus every entry in the `example_database.fasta` is a possible target in the output files.


![grafik](https://github.com/user-attachments/assets/0965f69c-6bf7-4155-9ca3-73287fd1c162)
	primer_lukF_11b_forward	primer_lukF_11b_revcomp	probe_lukF_10_forward	probe_lukF_10_revcomp	primer_entN_51_forward	primer_entN_51_revcomp	probe_entN_11_forward	probe_entN_11_revcomp
Query_Sequences	TGTCTGCCAGCTAAGAAG	CTTCTTAGCTGGCAGACA	AGCTTCCACCCAACATATGGTAATGA	TCATTACCATATGTTGGGTGGAAGCT	CTCTTCATCTAATTGATTTCCA	TGGAAATCAATTAGATGAAGAG	TCATGCTTATACGGAGGAGTTACGA	TCGTAACTCCTCCGTATAAGCATGA
SantaLucia_Tm	51.48	51.48	59.71	59.71	49.19	49.19	58.69	58.69
mm=0	2	1	1	2	0	1	1	0
mm=1	0	0	0	0	0	0	0	0
mm=2	0	0	0	0	0	0	0	0
example_database.fasta|CP102956.fna		0 (pos: 1942449-1942466)	0 (pos: 1942420-1942445)					
example_database.fasta|CP102957.fna	0 (pos: 2204112-2204129)			0 (pos: 2204133-2204158)				
example_database.fasta|CP102958.fna	0 (pos: 2344146-2344163)			0 (pos: 2344167-2344192)		0 (pos: 184189-184210)	0 (pos: 184157-184181)	
![grafik](https://github.com/user-attachments/assets/ace37924-becc-457c-b958-c08dbf42b288)

The Figure shows the output of the blast_result.tsv of the example data. The header of the table contains the names of the provided query oligos.
Every oligo is presented once as used for the forward BLAST search (_forward) and once with its reverse complementary sequence as used for the reverse BLAST search (_revcomp).
Below the oligo names are the Query sequences followd by thier respective Tm values calculated using the Santa Lucia formula.
Followed by the sum of BLAST hits with the respective mismatchcount (mm=0, mm=1 ... mm=--max_mismatches)
Below are the hits against every entry in the provided target sequences (example_database.fasta). If the query sequence was found in the target sequence it is listed with the number of mismatches followed by the 
start and end position within the target sequence. If there is no entry, there were no hits of the query sequence against the respective target sequence. For example the first privided oligo sequence primer_lukF_11b_forward

## Outputs
- **BLAST XML Results**: Detailed output of BLAST alignments including scores, E-values, mismatches, and alignments.
- **TSV Output**: Summary table containing mismatch counts and melting temperatures for each query-target alignment.
- **Alignments**: Text file containing all alignments with mismatches, useful for further review or filtering.
- **Multi-Hits**: Text file listing all targets that have multiple hits from the same query sequence.

## How It Works
1. **BLAST Database Creation**: The tool generates a BLAST database from the user-provided genome files.
2. **Forward and Reverse Complement BLAST**: Runs two separate BLAST searches: one for the forward strand and one for the reverse complement of the query sequences.
3. **Mismatch Filtering and Analysis**: Alignments are filtered based on user-defined mismatch thresholds. Results include mismatch counts, binding positions, and melting temperatures.
4. **Result Generation**: Outputs the results in various formats including XML (for detailed BLAST results), TSV (for easy parsing and analysis), and text (for alignments and multi-hits).

## Contributing
Feel free to contribute to this project by submitting pull requests, reporting issues, or suggesting improvements.

## License
This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact
For any questions, issues, or suggestions, please reach out via [GitHub Issues](https://github.com/mcollatz/AssayBLAST/issues).



















