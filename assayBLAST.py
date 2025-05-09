#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
import argparse
import os
import glob
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import MeltingTemp as mt

    # Create BLAST database
    print("Creating BLAST database...")
    db_file = create_blast_db(args.genome, args.db_name, args.db_dir, concatenate=args.concatenate)
    print("BLAST database created successfully.")
    
    with open(combined_fasta, "w") as out_handle:
        if concatenate:
            # Combine all sequences from each FASTA file into one sequence under a single header
            for db_file in db_files:
                super_contig = ""
                for record in SeqIO.parse(db_file, "fasta"):
                    super_contig += str(record.seq)
                out_handle.write(f">{os.path.basename(db_file)}\n{super_contig}\n")
        else:
            # Combine all sequences from each FASTA file without concatenation
            for db_file in db_files:
                for record in SeqIO.parse(db_file, "fasta"):
                    new_header = f"{os.path.basename(db_file).replace(' ','_')}|{record.id}" # also replaces white space from the file name
                    out_handle.write(f">{new_header}\n{record.seq}\n")

    cmd = f"makeblastdb -in {combined_fasta} -dbtype nucl -out {db_path}"
    os.system(cmd)
    return db_path

# Read the query sequences and create the reverse complement
def generate_reverse_complement_query(query_file, output_file):
    query_sequences = {}
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(query_file, "fasta"):
            query_sequences[f"{record.id}_forward"] = record.seq
            reverse_complement_seq = record.seq.reverse_complement()
            query_sequences[f"{record.id}_revcomp"] = reverse_complement_seq
            out_handle.write(f">{record.id}_revcomp\n{reverse_complement_seq}\n")
    return query_sequences

# Runs the BLAST search with the forward or the reverse complementary sequences respectively
def run_blast(query_file, db_name, output_file, db_dir, reverse_complement=False):
    db_path = os.path.join(db_dir, db_name)
    query_sequences = {}
    if reverse_complement:
        query_sequences = generate_reverse_complement_query(query_file, "reverse_complement_queries.fasta")
        query_file = "reverse_complement_queries.fasta"
    blast_cline = NcbiblastnCommandline(query=query_file, db=db_path, out=output_file, outfmt=5,
                                         dust='no', word_size=7, gapopen=10, gapextend=6,
                                         evalue=100000, reward=5, penalty=-4, strand='plus',
                                         max_target_seqs=50000)
    print(f"Running BLAST search for {'reverse complement ' if reverse_complement else ''}queries...")
    stdout, stderr = blast_cline()
    print(f"BLAST search for {'reverse complement ' if reverse_complement else ''}queries completed.")
    return query_sequences, stdout, stderr

# Analyze and filter BLAST results
def parse_blast_results(xml_file, max_mismatches):
    results = {}
    alignments = []
    all_hits = []
    mismatch_counts = {i: {} for i in range(max_mismatches + 1)}
    with open(xml_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            query_id = blast_record.query.split()[0]
            query_length = blast_record.query_length
            results[query_id] = {}
            for alignment in blast_record.alignments:
                target_id = alignment.title.split()[1]
                for hsp in alignment.hsps:
                    mismatches = hsp.align_length - hsp.identities
                    if mismatches <= max_mismatches and hsp.align_length == query_length:
                        if target_id not in results[query_id]:
                            results[query_id][target_id] = []
                        results[query_id][target_id].append((mismatches, hsp.sbjct_start, hsp.sbjct_end, hsp))
                        if mismatches > 0:  # Only consider alignments with mismatches
                            alignments.append((query_id, target_id, hsp))
                        if target_id not in mismatch_counts[mismatches]:
                            mismatch_counts[mismatches][target_id] = 0
                        mismatch_counts[mismatches][target_id] += 1
                        # Add to all_hits for additional output
                        all_hits.append((query_id, target_id, hsp))
    return results, alignments, mismatch_counts, all_hits

def combine_results(original_results, revcomp_results):
    combined_results = {}

    # Iterate through both original and reverse complement results
    for query_id in original_results.keys():
        combined_results[query_id + "_forward"] = original_results[query_id]
        combined_results[query_id + "_revcomp"] = revcomp_results.get(query_id, {})

    # Add any additional reverse complement results not present in original
    for query_id in revcomp_results.keys():
        if query_id not in original_results:
            combined_results[query_id] = revcomp_results[query_id]

    return combined_results

# Writes all hits that contain at least one mismatch to a file
def write_alignments(alignments, output_file):
    with open(output_file, "w") as out_handle:
        for query_id, target_id, hsp in alignments:
            out_handle.write(f">{query_id} vs {target_id}\n")
            out_handle.write(f"Score: {hsp.score}\n")
            out_handle.write(f"E-value: {hsp.expect}\n")
            out_handle.write(f"Identities: {hsp.identities}/{hsp.align_length}\n")
            out_handle.write(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n\n")

# Creates a file with only hits that occur more than once per target sequence
def write_multi_hits(all_hits, output_file, max_mismatches):
    # Filter hits to include only those queries with multiple hits against the same target
    filtered_hits = {}
    for query_id, target_id, hsp in all_hits:
        mismatches = hsp.align_length - hsp.identities
        if mismatches <= max_mismatches:
            if query_id not in filtered_hits:
                filtered_hits[query_id] = {}
            if target_id not in filtered_hits[query_id]:
                filtered_hits[query_id][target_id] = []
            filtered_hits[query_id][target_id].append(hsp)

    # Write filtered hits to file
    with open(output_file, "w") as out_handle:
        for query_id, targets in filtered_hits.items():
            for target_id, hsps in targets.items():
                if len(hsps) > 1:  # Only include if there are multiple hits against the same target
                    for hsp in hsps:
                        mismatches = hsp.align_length - hsp.identities
                        out_handle.write(f">{query_id} vs {target_id}\n")
                        out_handle.write(f"Score: {hsp.score}\n")
                        out_handle.write(f"E-value: {hsp.expect}\n")
                        out_handle.write(f"Identities: {hsp.identities}/{hsp.align_length}\n")
                        out_handle.write(f"Position: {hsp.sbjct_start}-{hsp.sbjct_end}\n")
                        out_handle.write(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}\n\n")


def main():
    parser = argparse.ArgumentParser(description="Perform BLAST search with genome FASTA file and query FASTA file.")
    parser.add_argument("-g", "--genome", help="Glob pattern for genome FASTA files")
    parser.add_argument("-q", "--queries", help="Path to query FASTA file")
    parser.add_argument("-d", "--db_name", default="genome_db", help="Name for BLAST database (default: genome_db)")
    parser.add_argument("-o", "--output", default="blast_results.xml", help="Output file name for BLAST results (default: blast_results.xml)")
    parser.add_argument("-c", "--tsv_output", default="blast_results.tsv", help="Output file name for mismatches matrix in TSV format (default: blast_results.tsv)")
    parser.add_argument("-a", "--alignments_output", default="alignments.txt", help="Output file name for alignments with mismatches (default: alignments.txt)")
    parser.add_argument("-mh", "--multi_hits_output", default="multi_hits.txt", help="Output file name for multi BLAST hits (default: multi_hits.txt)")
    parser.add_argument("-db", "--db_dir", default="blast_db", help="Directory to store BLAST database files (default: blast_db)")
    parser.add_argument("-m", "--max_mismatches", type=int, default=4, help="Maximum number of accepted mismatches (default: 4)")
    parser.add_argument("-cc", "--concatenate", action="store_true", help="Concatenate input sequences into one (default: False)")
    args = parser.parse_args()

    # Create BLAST database if it doesn't exist
    if not os.path.exists(os.path.join(args.db_dir, args.db_name + ".nhr")):
        print("Creating BLAST database...")
        db_file = create_blast_db(args.genome, args.db_name, args.db_dir, concatenate=args.concatenate)
        print("BLAST database created successfully.")

    # Run BLAST search for original queries
    query_sequences, stdout, stderr = run_blast(args.queries, args.db_name, args.output, args.db_dir)
    if stderr:
        print("An error occurred during BLAST search for original queries:")
        print(stderr)
    else:
        print("BLAST search for original queries completed successfully.")

    # Run BLAST search for reverse complement queries
    query_sequences, stdout, stderr = run_blast(args.queries, args.db_name, "reverse_complement_" + args.output, args.db_dir, reverse_complement=True)
    if stderr:
        print("An error occurred during BLAST search for reverse complement queries:")
        print(stderr)
    else:
        print("BLAST search for reverse complement queries completed successfully.")

    # Parse BLAST results for original queries
    print("Parsing BLAST results for original queries...")
    original_results, original_alignments, original_mismatch_counts, original_all_hits = parse_blast_results(args.output, args.max_mismatches)

    # Parse BLAST results for reverse complement queries
    print("Parsing BLAST results for reverse complement queries...")
    revcomp_results, revcomp_alignments, revcomp_mismatch_counts, revcomp_all_hits = parse_blast_results("reverse_complement_" + args.output, args.max_mismatches)

    # Combine results
    print("Combining results...")
    combined_results = combine_results(original_results, revcomp_results)

    # Combine mismatch counts
    combined_mismatch_counts = {i: {} for i in range(args.max_mismatches + 1)}
    for i in range(args.max_mismatches + 1):
        for target_id, count in original_mismatch_counts[i].items():
            combined_mismatch_counts[i][target_id] = count
        for target_id, count in revcomp_mismatch_counts[i].items():
            if target_id in combined_mismatch_counts[i]:
                combined_mismatch_counts[i][target_id] += count
            else:
                combined_mismatch_counts[i][target_id] = count

    # Calculate mismatch counts for each query
    query_mismatch_counts = {query_id: {i: 0 for i in range(args.max_mismatches + 1)} for query_id in combined_results.keys()}
    for query_id in combined_results.keys():
        for target_id, mismatches_list in combined_results[query_id].items():
            for mismatches, _, _, _ in mismatches_list:
                if mismatches in query_mismatch_counts[query_id]:
                    query_mismatch_counts[query_id][mismatches] += 1
                else:
                    query_mismatch_counts[query_id][mismatches] = 1

    # Write combined results to table
    print("Writing combined results to table...")
    with open(args.tsv_output, "w") as out_table:
        # Write header
        out_table.write("\t" + "\t".join(combined_results.keys()) + "\n")
        # Write query sequences
        query_sequences = [str(query_sequences.get(query_id, "")) for query_id in combined_results.keys()]
        out_table.write("Query_Sequences\t" + "\t".join(query_sequences) + "\n")
        # Write SantaLucia Tm values
        tm_values = [str(round(mt.Tm_NN(sequence, Na=100), 2)) for sequence in query_sequences]
        out_table.write("SantaLucia_Tm\t" + "\t".join(tm_values) + "\n")
        # Write mismatch counts
        for mm in range(args.max_mismatches + 1):
            mm_key = f"mm={mm}"
            mm_counts = [str(query_mismatch_counts[query_id].get(mm, 0)) for query_id in combined_results.keys()]
            out_table.write(mm_key + "\t" + "\t".join(mm_counts) + "\n")
        # Write data rows
        for target_id in sorted(set(target for targets in combined_results.values() for target in targets.keys())):
            row_values = []
            for query_id in combined_results.keys():
                mismatches_list = combined_results[query_id].get(target_id, [])
                if not mismatches_list:
                    row_values.append("")
                else:
                    row_values.append("; ".join([f"{mismatches} (pos: {start}-{end})" for mismatches, start, end, _ in mismatches_list]))
            out_table.write(target_id + "\t" + "\t".join(row_values) + "\n")

    print(f"Combined results table written to {args.tsv_output}")

    # Write alignments with mismatches to file
    all_alignments = original_alignments + revcomp_alignments
    print("Writing alignments with mismatches to file...")
    write_alignments(all_alignments, args.alignments_output)
    print(f"Alignments with mismatches written to {args.alignments_output}")
    
    # Write multi hits to file
    all_hits = original_all_hits + revcomp_all_hits
    print("Writing multi BLAST hits to file...")
    write_multi_hits(all_hits, args.multi_hits_output, args.max_mismatches)
    print(f"Multi BLAST hits written to {args.multi_hits_output}")


if __name__ == "__main__":
    main()


