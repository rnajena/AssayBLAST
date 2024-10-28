# AssayBLAST - for in silico analysis of PCR oligos

This tool provides in silico predictions of microarray hybridization results, calculating the expected binding interactions between query DNA sequences (primers and probes) and a genome database. Using mismatch counts and probe abundance, it identifies likely positive or negative hybridization results based on a threshold interpretation, which is especially useful for microarray design and diagnostic applications.



## Key Features
- BLAST Database Creation: Automatically generates a custom BLAST database from user-provided genome files.
- BLAST Search Execution: Runs forward and reverse complement BLAST searches to evaluate how primers or probes interact with both DNA strands.
- Mismatch Analysis: Identifies and filters matches with acceptable levels of mismatches, simulating how primers perform in the presence of mutations or single nucleotide polymorphisms (SNPs).
- Melting Temperature (Tm) Calculation: Computes the melting temperature of sequences using SantaLucia's thermodynamic nearest-neighbor model.
- Off-Target Detection: Helps identify off-target binding sites to ensure primer/probe specificity.
- Support for Multiple Genomes: Capable of performing searches across multiple genome sequences, ideal for comparative genomics or multi-strain analyses.
- Results in Multiple Formats: Outputs detailed results in XML, TSV, and text formats for easy review and further analysis.

---

## Requirements
- Python 3.7+
- Biopython
- NCBI BLAST+

## Usage
Download and run the AssayBLAST.py
