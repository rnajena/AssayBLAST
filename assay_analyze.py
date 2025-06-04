#!/usr/bin/env python
# (C) 2025, Tom Eulenfeld, Max Collatz, MIT license
"""
Read BLAST/GFF file and report probes and primers.

Find primer-probe-primer triplets on +,+,- or +,-,- strand (exponential growth)
or probe-primer pairs (linear growth).
You can also search for primer pairs only.
For BLAST files, the script expects the string "primer" or "probe"
to be part of the source name (e.g., the sseqid field).
We also recommend inclusion of the qlen field as created with the assay_blast.py script.
For GFF files, the type must be primer or probe and mismatch has to be defined.
"""
_epilog = """
This script needs rnajena-sugar.
Install with
    `pip install --no-deps rnajena-sugar`.
"""

import argparse
from pathlib import Path
from warnings import warn


__version__ = '2.0'


# import sys
# def _debug(type, value, tb):
#     if hasattr(sys, 'ps1') or not sys.stderr.isatty():
#     # we are in interactive mode or we don't have a tty-like
#     # device, so we call the default hook
#         sys.__excepthook__(type, value, tb)
#     else:
#         import traceback, pdb
#         # we are NOT in interactive mode, print the exception...
#         traceback.print_exception(type, value, tb)
#         print
#         # ...then start the debugger in post-mortem mode.
#         # pdb.pm() # deprecated
#         pdb.post_mortem(tb) # more "modern"
# sys.excepthook = _debug


def _find_probe_and_primers(superfts, distance=250, fast=False):
    """Find probes and corresponding primers"""
    r = []
    probes_done = set()
    for fts in superfts.groupby('seqid').values():
        for i, probe in enumerate(fts):
            if probe.type != 'probe' or fast and probe.meta.name in probes_done:
                continue
            primer1 = []
            for ft in fts[i-1::-1]:  # look for primers before probe on the + strand
                if ft.type != 'primer' or ft.loc.strand != '+':
                    continue
                if probe.distance(ft) > distance:
                    break
                if probe.overlaps(ft):
                    warn(f'Detected overlap between {probe.meta.name} and {ft.meta.name}')
                primer1.append(ft)
            primer2 = []
            for ft in fts[i+1:]:   # look for primers after probe on the - strand
                if ft.type != 'primer' or ft.loc.strand != '-':
                    continue
                if probe.distance(ft) > distance:
                    break
                if probe.overlaps(ft):
                    warn(f'Detected overlap between {probe.meta.name} and {ft.meta.name}')
                primer2.append(ft)
            # for linear amplification remove probe primer pairs located on the same strand
            if len(primer1) == 0 and len(primer2) > 0:
                primer2 = [p for p in primer2 if p.loc.strand != probe.loc.strand]
            elif len(primer1) > 0 and len(primer2) == 0:
                primer1 = [p for p in primer1 if p.loc.strand != probe.loc.strand]
            # depending on amplification add probe-primer singlets/doublets/triplets
            if len(primer1) + len(primer2) == 0:
                r.append((probe.meta.name, (probe, )))
            elif len(primer1) == 0:
                for ft in primer2:
                    r.append((probe.meta.name, (probe, ft)))
            elif len(primer2) == 0:
                for ft in primer1:
                    r.append((probe.meta.name, (ft, probe)))
            else:
                r2 = []
                for ft in primer1:
                    for ft2 in primer2:
                        r2.append((probe.meta.name, (ft, probe, ft2)))
                        if fast and ft.meta.mismatch + probe.meta.mismatch + ft2.probe.mismatch == 0:
                            probes_done.add(probe.meta.name)
                r.extend(r2)
    return r


def _find_primers(fts, distance=250):
    """Find primer pairs"""
    r = []
    done = set()
    for i, ft in enumerate(fts):
        found_pair = False
        if ft.loc in done:
            # primer on - strand already used in pair
            continue
        if ft.loc.strand == '-':
            # add not used single primer on - strand
            r.append((ft.meta.name, (ft,)))
        else:
            for ft2 in fts[i+1:]:
                if ft.loc.range == ft2.loc.range:
                    done.add(ft2.loc)
                    continue
                if ft2.loc.strand == '+':
                    continue
                if ft2.distance(ft) > distance:
                    break
                found_pair = True
                done.add(ft2.loc)
                name = ft.meta.name + ',' + ft2.meta.name
                r.append((name, (ft, ft2)))
            if not found_pair:
                r.append((ft.meta.name, (ft, )))
    return r


def _reassign_primer_header(results):
    """
    Assign linear primer to other primer names
    """
    primer_pairs = set(name for r in results.values() for name, _ in r if ',' in name)
    new_names = {p: name for name in primer_pairs for p in name.split(',')}
    for rs in results.values():
        for i, r in enumerate(rs):
            name, pps = r
            if ',' not in name and name in new_names:
                rs[i] = (new_names[name], pps)


def _sortkey(result):
    """Key to sort by goodness"""
    name, r = result
    return name, -len(r), sum(ft.meta.mismatch for ft in r), r[0].distance(r[-1])


def _groupby(ft):
    return ft.seqid.split('--')[0]


def find_probe_primer(allfts, mismatch=2, only_primer=False, distance=150):
    """Find primer-probe-primer triplets or primer pairs"""
    allfts = allfts.select(type_in=('primer', 'probe'), mismatch_lt=mismatch).sort()
    results = {}
    for superseqid, fts in allfts.groupby(_groupby).items():
        if only_primer:
            fts = fts.select(type='primer')
            r = _find_primers(fts, distance=distance)
        else:
            r = _find_probe_and_primers(fts, distance=distance, fast=False)
        results[superseqid] = sorted(r, key=_sortkey)
    if only_primer:
        _reassign_primer_header(results)
    return results


NUM2STR_PROBE = {0: 'none', 1: 'noprimer', 2: 'lin', 3: 'exp'}
NUM2STR_PRIMER = {0: 'none', 1: 'lin', 2: 'exp'}


def output_assay_overview(results, out, only_primer=False):
    """Write the assay overview file"""
    num2str = NUM2STR_PRIMER if only_primer else NUM2STR_PROBE
    lines = []
    all_probes = sorted(set(proben for res in results.values() for proben, _ in res))
    lines.append('seqid\t' + '\t'.join(all_probes))
    for superseqid, r2 in results.items():
        probes = {}
        for proben, r in r2:
            if proben not in probes:
                if only_primer:
                    primerstr = ' '.join(f'{"F" if p.loc.strand == "+" else "R"}{p.meta.mismatch}{p.loc.strand}' for p in r)
                else:
                    primerstr = ' '.join(f'{"P" if p.type == "probe" else "F" if i == 0 else "R"}{p.meta.mismatch}{p.loc.strand}' for i, p in enumerate(r))
                growth = num2str[len(r)]
                # mms = ''.join(str(p.meta.mismatch) for p in r)
                # primerstr = ''.join('P' if p.type == 'probe' else '>' if p.loc.strand == '+' else '<' for p in r)
                probes[proben] = f'{growth} {primerstr}'
        line = f'{superseqid}\t' + '\t'.join(probes.get(proben, NUM2STR_PROBE[0]) for proben in all_probes)
        lines.append(line)
    with open(out, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def output_assay_details(results, out, only_primer=False, verbose=False):
    """Write the assay details file"""
    num2str = NUM2STR_PRIMER if only_primer else NUM2STR_PROBE
    lines = []
    for superseqid, res in results.items():
        max_growth = 0
        for proben, pps in res:
            max_growth = max(max_growth, len(pps))
            growth = num2str[len(pps)]
            names = ','.join(p.meta.name + p.loc.strand for p in pps)
            mms = ','.join(str(p.meta.mismatch) for p in pps)
            dists = ','.join(str(p1.distance(p2)) for p1, p2 in zip(pps[1:], pps[:-1]))
            pos = ','.join(f'{p.loc.start}:{p.loc.stop}:{p.loc.strand}' for p in pps)
            contig = pps[0].seqid.split('--')[-1]
            contig = '' if contig == superseqid else contig + '_'
            lines.append(f'{superseqid}\t{proben}\t{growth}\t{names}\tM{mms}\tD{dists}\t{contig}{pos}\n')
        if verbose:
            print(f'{superseqid}\t{num2str[max_growth]}')
    with open(out, 'w') as f:
        f.write(''.join(lines))


def find_probe_primer_cli(fname, out=None, only_primer=False, **kw):
    """Read BLAST file or GFF file and find+report probes/primers"""
    try:
        from sugar import read_fts
    except ImportError:
        import sys
        sys.exit('Script needs rnajena-sugar. Install with:\npip install rnajena-sugar')
    fts = read_fts(fname)
    print(f'Successfully parsed BLAST/GFF file at {fname}.')
    fmt = fts[0].meta._fmt
    if fmt == 'blast':
        for ft in fts:
            ft.meta.mismatch = ft.meta._blast.mismatch
            if 'qlen' in ft.meta._blast:
                if len(ft) != ft.meta._blast.qlen:
                    warn('Ignore BLAST hits shorter than query length')
                    continue
            else:
                warn('No qlen row found, cannot check if BLAST hit spans full query length')
            if 'probe' in ft.meta.name:
                ft.type = 'probe'
            elif 'primer' in ft.meta.name:
                ft.type = 'primer'
            else:
                warn(f'Ignore BLAST hits which are neither primer nor probe')
    elif fmt == 'gff':
        for ft in fts:
            if 'mismatch' in ft.meta._gff:
                ft.meta.mismatch = int(ft.meta._gff.mismatch)
            else:
                warn(f'No mismatch information found in GFF file, use mismatch 0')
                ft.meta.mismatch = 0
    results = find_probe_primer(fts, only_primer=only_primer, **kw)
    if out is None:
        fname = Path(fname)
        out = fname.with_name(fname.stem + '_assay' + '_only_primer' * only_primer)
    out = Path(out)
    out1 = out.with_name(out.stem + '_overview.tsv')
    out2 = out.with_name(out.stem + '_details.tsv')
    output_assay_overview(results, out1, only_primer=only_primer)
    print(f'Assay overview file created at {out1}.')
    output_assay_details(results, out2, only_primer=only_primer)
    print(f'Assay details file created at {out2}.')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog=_epilog)
    parser.add_argument('fname', help='BLAST outfile outfmt 7, 6 or 10 or GFF file, e.g. created with assay_blast.py')
    parser.add_argument("--mismatch", type=int, default=2, help='Maximum allowed mismatches (default: 2)')
    parser.add_argument("--distance", type=int, default=250, help='Distance threshold (bp) between adjacent oligos (default: 250)')
    parser.add_argument('--only-primer', action='store_true', help='Find primer pairs instead of primer-probe-primer triplets')
    parser.add_argument('-o', '--out', help='Prefix of output files, by default derived from fname')

    args = vars(parser.parse_args())
    find_probe_primer_cli(**args)


if __name__ == "__main__":
    main()
