
from contextlib import contextmanager
from pathlib import Path
from subprocess import check_output
import tempfile
import urllib.request


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
    check_output(cmd.split())


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
            print()
            out = tmp / 'probes_filename_as_id.blastn'
            db = tmp / 'db/db.db'
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db} --filename-as-id')
            _call(f'assay_blast {genomes} -q {query} -o {out} --db {db} --filename-as-id  --mismatch-alignments --keep-db')
            _call(f'assay_analyze {out}')
            _call(f'assay_analyze {out} --only-primer -o {tmp / "primer_filename_as_id"}')
            print()
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
