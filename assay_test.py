
from contextlib import contextmanager
from pathlib import Path
from subprocess import run
import tempfile
import urllib.request

try:
    import pandas as pd
except ImportError:
    print('pandas not found, please install for full tests')
    pd = None
from sugar import read



@contextmanager
def _tmpdir(out=None):
    if out:
        tmpdir = Path(out)
        tmpdir.mkdir(exist_ok=True)
        yield tmpdir
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)


_URL = 'https://raw.githubusercontent.com/rnajena/AssayBLAST/refs/heads/main/'


def _download(fname):
    if not fname.exists():
        url = _URL + fname.name
        print(f'Download file {url}')
        urllib.request.urlretrieve(url, fname)
    return fname


def _call(cmd):
    print(cmd)
    run(cmd.split())
    print()


def _read_result(out, type_):
    return pd.read_csv(str(out).replace('.blastn', f'_assay_{type_}.tsv'), sep='\t')

def _count_linear(df):
    return df.apply(lambda x: x.str.startswith('lin')).values.sum()


def test_assay(out=None, testit=True):
    with _tmpdir(out) as tmp:
        query = _download(tmp / 'example_queries.fasta')
        genomes = _download(tmp / 'example_database.fasta')
        if testit:
            out = tmp / 'probes.blastn'
            db = tmp / 'db/db.db'
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db}')
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db} --mismatch-alignments --keep-db')
            _call(f'assay_analyze {out}')
            _call(f'assay_analyze {out} --zero-based-numbering -o {tmp / "probes_assay_0based"}')
            _call(f'assay_analyze {out} --only-primer -o {tmp / "primer_assay"}')
            if pd:
                df1 = _read_result(out, 'overview')
                num_linear = _count_linear(df1)
                df2 = _read_result(out, 'details')
                assert len(df2) ==  4
            print()
            # test filename_as_id parameters
            out = tmp / 'probes_filename_as_id.blastn'
            db = tmp / 'db/db.db'
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db} --filename-as-id')
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db} --filename-as-id  --mismatch-alignments --keep-db')
            _call(f'assay_analyze {out}')
            _call(f'assay_analyze {out} --only-primer -o {tmp / "primer_filename_as_id"}')
            if pd:
                df = _read_result(out, 'overview')
                num_linear = _count_linear(df)
                assert num_linear == 2
                df = _read_result(out, 'details')
                assert len(df) ==  4
            print()
            # test multiple genome files
            seqs = read(genomes)
            genomes1 = tmp / 'genome1.fasta'
            genomes2 = tmp / 'genome2.fasta'
            seqs[:1].write(genomes1)
            seqs[1:].write(genomes2)
            out = tmp / 'probes_multiple_genomes.blastn'
            _call(f'assay_blast {genomes1} {genomes2} -q {query} -o {out} --db {db}')
            _call(f'assay_analyze {out}')
            if pd:
                df = _read_result(out, 'overview')
                num_linear = _count_linear(df)
                assert num_linear == 4
                # results are the same independent of whether we use multiple genome files
                # or a single combined file
                assert df.values.tolist() == df1.values.tolist()
                df = _read_result(out, 'details')
                assert len(df) ==  4
                assert df.values.tolist() == df2.values.tolist()
            print()
            out = tmp / 'probes_multiple_genomes_filename_as_id.blastn'
            _call(f'assay_blast {genomes1} {genomes2} -q {query} -o {out} --db {db} --filename-as-id')
            _call(f'assay_analyze {out}')
            if pd:
                df = _read_result(out, 'overview')
                num_linear = _count_linear(df)
                assert num_linear == 3
                df = _read_result(out, 'details')
                assert len(df) ==  4
            print('Tests run successful.')


def main():
    import argparse
    p = argparse.ArgumentParser('assay_tests')
    p.add_argument('out', help='output directory (by default a temporary directory)', nargs='?')
    p.add_argument('-d', '--download', help='only download test files, do not run tests', action='store_true')
    args = p.parse_args()
    if args.download and not args.out:
        p.error('[out] argument needed for --download option')
    test_assay(args.out, not args.download)


if __name__ == '__main__':
    main()
