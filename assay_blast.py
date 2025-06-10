#!/usr/bin/env python
# (C) 2025, Tom Eulenfeld, Max Collatz, MIT license
"""
Start BLAST primer/probe assay by running BLAST
"""
_epilog = """
This script needs rnajena-sugar.
Install with
    `pip install --no-deps rnajena-sugar`.
"""

import argparse
from functools import reduce
from operator import add
import os
from pathlib import Path
import re
import time
from warnings import warn


__version__ = '2.0'


class ParseError(Exception):
    pass


OUTFMT7 = '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen'
REWARD = 5
PENALTY = -4
BLAST = (
    'blastn '
    '-query {query} -db {db} '
    '-word_size 7 -gapopen 1000 -gapextend 1000 -reward {reward} -penalty {penalty} '
    '-perc_identity {perc:.2f} -qcov_hsp_perc 100 '
    '-evalue {evalue} -max_target_seqs 100000 '
    '-num_threads {num_threads} '
    "-outfmt '{outfmt}' -dust no -out {out}"
    )


_CACHE = {}

def _read(fname):
    if fname in _CACHE:
        return _CACHE[fname]
    try:
        from sugar import read
    except ImportError:
        import sys
        sys.exit('Script needs rnajena-sugar. Install with:\npip install rnajena-sugar')
    seqs = _CACHE[fname] = read(fname)
    return seqs


def _read_and_add_id(fname, add_id=False):
    """Read file and optionally prepend seqid with the file stem"""
    seqs = _read(fname)
    if add_id:
        stem = Path(fname).stem
        for seq in seqs:
            seq.id = stem + '--' + seq.id
    return seqs


def _evalue_factor(n, mismatch):
    """Calculate factor for evalue from max mismatch and query length"""
    score = REWARD * (n - mismatch) + PENALTY * mismatch
    bitscore = 0.2764 * score + 2.474
    return n / 2 ** bitscore


def _get_blast_evalue(query, combined, mismatch):
    """Calculate BLAST evalue and perc"""
    len_queries = [len(seq) for seq in _read(query)]
    N = min(len_queries)
    Nmax = max(len_queries)
    M = sum(len(seq) for seq in _read(combined))
    perc = 100 * max(0, N - mismatch - 0.5) / N
    evalue = 1.5 * M * max(_evalue_factor(N, mismatch), _evalue_factor(Nmax, mismatch))  # use 1.5x the expected upper bound of evalue
    if evalue > 0.1:
        evalue = round(evalue, 2)
    return {'perc': perc, 'evalue': evalue}


def _adapt_outfmt7(fname, call, mismatch, query_ids, source_ids, **blast_config):
    """Add comments to BLAST outfmt 7 file"""
    intro = (f'# File was written with BLAST outfmt=7, assay_blast v{__version__}. '
                'More information at the end of the file.\n')
    extro = (
        f'#\n# assay_blast: Calculated evalue {blast_config["evalue"]} from option {mismatch=}\n'
        f'# BLAST call: {call}\n'
        f'# queryids: {" ".join(query_ids)}\n'
        f'# sourceids: {" ".join(source_ids)}\n'
        )
    with open(fname) as f:
        data = f.read()
    with open(fname, 'w') as f:
        f.write(intro + data + extro)


def _filter_outfmt0(fname):
    """Filter BLAST outfmt 0 file, only select alignments with mismatches"""
    with open(fname) as f:
        data = f.read()
    regex = r'(?sm)Query=.*?$|>[^>]*?Identities = [\d/]+ \((?!100)\d+%\).*?(?=Lambda|>)'
    data2 = re.findall(regex, data)
    with open(fname, 'w') as f:
        f.write(f'# File was written with BLAST outfmt=0 and filtered by assay_blast v{__version__}.\n'
                '# Only alignments with mismatches are kept.\n')
        f.write('\n'.join(data2) + '\n')


def run_blast(query, genomes, out, db=None, super_contig=False, mismatch=2, num_threads=1, keep_db=False):
    """
    Run blast for primer/probe detection, see CLI help for description of arguments
    """
    if not db and len(genomes) == 0:
        raise ParseError('Please specify either db or genomes argument')
    if db:
        db = Path(db)
    else:
        db = Path(genomes[0])
        db = (db.with_name(db.stem + '_database') / db.stem)
    combined = db.with_suffix('.combined.fasta')
    if not keep_db or not Path(str(db) + '.nsq').exists():
        if len(genomes) == 0:
            raise ParseError('To create the BLAST database, please specify genome files')
        db.parent.mkdir(exist_ok=True)
        srcs = [_read_and_add_id(fname, add_id=super_contig) for fname in genomes]
        seqs = _CACHE[combined] = reduce(add, srcs)
        seqs.write(combined)
        del seqs
        call = f'makeblastdb -in {combined} -dbtype nucl -out {db}'
        print(call)
        os.system(call)
    else:
        print(f'Use BLAST database at {db}')
    blast_config = _get_blast_evalue(query, combined, mismatch)
    blast_config.update(query=query, db=db, reward=REWARD, penalty=PENALTY, num_threads=num_threads)
    out = Path(out)
    out2 = out.with_name(out.stem + '_mismatching_alignments.txt')
    t1 = time.time()
    call = BLAST.format(out=out, outfmt=OUTFMT7, **blast_config)
    print(call)
    os.system(call)
    call = BLAST.format(out=out2, outfmt=0, **blast_config)
    print(call)
    os.system(call)
    t2 = time.time()
    print(f'blast time used: {t2-t1}s')
    source_ids = [seq.id.split('--')[0] for seq in _read(combined)]
    _adapt_outfmt7(out, call=call, mismatch=mismatch, query_ids=_read(query).ids, source_ids=source_ids, **blast_config)
    _filter_outfmt0(out2)
    print()
    print(f'BLAST file created at {out}.')
    print('The BLAST file can be used as input for the assay_analyze script.')
    print(f'Mismatching alignments file created at at {out2}.')


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog=_epilog)
    parser.add_argument('genomes', nargs='*', help='Files with genomes (FASTA, GenBank, ...)')
    parser.add_argument('-q', '--query', help='File with primers and probes (FASTA)', required=True)
    parser.add_argument('-o', '--out', help='BLAST output file', default='blast_results.tsv')
    parser.add_argument('-n', '--num-threads', type=int, default=1, help='Number of threads used for BLAST')
    parser.add_argument('--super-contig', action='store_true', help='Treat sequences in each genome file as a super-contig and not as individual sequences')
    parser.add_argument('--db', help='BLAST DB file, by default derived from first genome file')
    parser.add_argument('--keep-db', action='store_true', help='Overwrite database if it exists')
    parser.add_argument('--mismatch', type=int, default=2, help='Maximum allowed mismatches, defaults to 2, used to calculate the evalue passed to BLAST')
    args = vars(parser.parse_args())
    try:
        run_blast(**args)
    except ParseError as ex:
        parser.error(ex)


if __name__ == "__main__":
    main()
