
from contextlib import contextmanager
import os
from pathlib import Path
import sys
import tempfile


@contextmanager
def _tmpdir():
    if '-h' in sys.argv:
        print('Run tests for assay_blast and assay_analyze.\nOnly accepted argument is a tmpdir.')
        sys.exit()
    if len(sys.argv) == 2:
        tmpdir = Path(sys.argv[1])
        tmpdir.mkdir(exist_ok=True)
        yield tmpdir
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)


def test_assay():
    example = Path(__file__).parent
    query = example / 'example_queries.fasta'
    genomes = example / 'example_database.fasta'
    with _tmpdir() as tmp:
        out = tmp / 'probes.blastn'
        db = tmp / 'db.db'
        call = f'assay_blast.py {genomes} -q {query} -o {out} --db {db}'
        print(call)
        os.system(call)
        os.system(call)
        call = f'assay_analyze.py {out}'
        print(call)
        os.system(call)
        call = f'assay_analyze.py {out} --only-primer -o {tmp / "primer"}'
        print(call)
        os.system(call)
    print()
    print('Tests run successful.')


if __name__ == '__main__':
    test_assay()
