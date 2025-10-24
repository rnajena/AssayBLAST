#!/usr/bin/env python
# (C) 2025, Tom Eulenfeld, Max Collatz, MIT license
"""
Start BLAST primer/probe assay by running BLAST
"""
import argparse
from functools import reduce
from operator import add
import os
from pathlib import Path
import re
import time
import warnings
from warnings import warn


__version__ = '2.0'


def _formatwarning(message, category, filename, lineno, file=None, line=None):
    return f'{filename}:{lineno}: {category.__name__}: {message}\n'

warnings.formatwarning = _formatwarning


class ParseError(Exception):
    pass


OUTFMT7 = '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen'
REWARD = 5
PENALTY = -4
MAX_TARGET_SEQS = 1_000_000_000
BLAST = (
    'blastn '
    '-query {query} -db {db} '
    '-word_size 7 -gapopen 1000 -gapextend 1000 -reward {reward} -penalty {penalty} '
    '-perc_identity {perc:.2f} -qcov_hsp_perc {perc:.2f} '
    '-evalue {evalue} -max_target_seqs {max_target_seqs} '
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
    # The filtering by identity and query coverage percentage should not be necessary,
    # because we already filter with E-values, but it does not hurt either
    perc = 100 * max(0, N - mismatch - 0.5) / N
    # the first term should be the largest, we are pedantic here
    # use 1.5x the expected upper bound of evalue
    evalue = 1.5 * M * max(_evalue_factor(N, mismatch), _evalue_factor(N - mismatch, 0),
                           _evalue_factor(Nmax, mismatch), _evalue_factor(Nmax - mismatch, 0))
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


def _check_duplicates(seqs):
    nts = {}
    for seq in seqs:
        if seq.data in nts:
            warn(f'Detected duplicate: {seq.id} and {nts[seq.data]}')
        elif seq.rc().data in nts:
            warn(f'Detected duplicate: {seq.id} is the reverse complement of {nts[seq.data]}')
        else:
            nts[seq.data] = seq.id
            for seq2 in seqs:
                if seq2 != seq:
                    if seq2.data in seq.data:
                        warn(f'{seq2.id} is part of {seq.id}')
                    elif seq2.rc().data in seq.data:
                        warn(f'The reverse complement of {seq2.id} is part of {seq.id}')


def run_blast(query, genomes, out, db=None, filename_as_id=False, mismatch=2, num_threads=1,
              keep_db=False, mismatch_alignments=False):
    """
    Run blast for primer/probe detection, see CLI help for description of arguments
    """
    _check_duplicates(_read(query))
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
        srcs = [_read_and_add_id(fname, add_id=filename_as_id) for fname in genomes]
        seqs = _CACHE[combined] = reduce(add, srcs)
        seqs.write(combined)
        del seqs
        call = f'makeblastdb -in {combined} -dbtype nucl -out {db}'
        print(call)
        os.system(call)
    else:
        print(f'Use BLAST database at {db}')
    blast_config = _get_blast_evalue(query, combined, mismatch)
    blast_config.update(query=query, db=db, reward=REWARD, penalty=PENALTY,
                        num_threads=num_threads, max_target_seqs=MAX_TARGET_SEQS)
    out = Path(out)
    out2 = out.with_name(out.stem + '_mismatching_alignments.txt')
    t1 = time.time()
    if not mismatch_alignments:
        call = BLAST.format(out=out, outfmt=OUTFMT7, **blast_config)
        print(call)
        os.system(call)
        source_ids = [seq.id.split('--')[0] for seq in _read(combined)]
        source_ids = list(dict.fromkeys(source_ids))  # remove duplicates and keep order
        _adapt_outfmt7(out, call=call, mismatch=mismatch, query_ids=_read(query).ids, source_ids=source_ids, **blast_config)
    else:
        call = BLAST.format(out=out2, outfmt=0, **blast_config)
        print(call)
        os.system(call)
        _filter_outfmt0(out2)
    t2 = time.time()
    print(f'blast time used: {t2-t1}s\n')
    if not mismatch_alignments:
        print(f'BLAST file created at {out}.')
        print('The BLAST file can be used as input for the assay_analyze script.')
    else:
        print(f'Mismatching alignments file created at at {out2}.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('genomes', nargs='*', help='Files with genomes (FASTA, GenBank, ...)')
    parser.add_argument('-q', '--query', help='File with primers and probes (FASTA), required argument', required=True)
    parser.add_argument('-o', '--out', help='BLAST output file, default is blast_results.tsv', default='blast_results.tsv')
    parser.add_argument('-n', '--num-threads', type=int, default=1, help='Number of threads used for BLAST, defaults to 1')
    parser.add_argument('--filename-as-id', action='store_true', help=
                        'If set, treat all sequences in a genome file as belonging to the same source (e.g. genome or organism). '
                        'In this case, the source ID is derived from the file name. '
                        'If not set, the ID is derived from the FASTA header and every sequence is treated separately.')
    parser.add_argument('--db', help='BLAST DB prefix, by default derived from first genome file')
    parser.add_argument('--keep-db', action='store_true', help='Keep the BLAST database if it exists (the database only needs to be re-created for different genomes input)')
    msg = ('Maximum allowed mismatches, defaults to 2, used to calculate the e-value passed to BLAST. '
           'Note, that the BLAST output file may include hits with a higher mismatch.')
    parser.add_argument('--mismatch', type=int, default=2, help=msg)
    msg = ('Run BLAST with a different --outfmt parameter to create a file with mismatch alignments instead of the usual results file. '
           'Consider using the --keep-db option. The file specified with the --out option will automatically receive an appropriate file suffix.')
    parser.add_argument('--mismatch-alignments', action='store_true', help=msg)
    args = vars(parser.parse_args())
    try:
        run_blast(**args)
    except ParseError as ex:
        parser.error(ex)


if __name__ == "__main__":
    main()
