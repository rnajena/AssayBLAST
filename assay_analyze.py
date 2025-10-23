#!/usr/bin/env python
# (C) 2025, Tom Eulenfeld, Max Collatz, MIT license
"""
Read BLAST/GFF file and report probes and primers.

Find primer-probe-primer triplets on +,+,- or +,-,- strand (exponential growth)
or probe-primer pairs (linear growth).
You can also search for primer pairs only.
For BLAST files, the script expects the string "primer" or "probe"
to be part of the source name (e.g., the sseqid field).
We also recommend inclusion of the qlen field as created with the assay_blast script.
For GFF files, the type must be primer or probe and mismatch has to be defined.
"""

import argparse
from copy import deepcopy
from pathlib import Path
from warnings import warn


from assay_blast import __version__


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
            if i > 0:
                for ft in fts[i-1::-1]:  # look for primers before probe on the + strand
                    if ft.type == 'primer':
                        if probe.distance(ft) > distance:
                            break
                        if probe.overlaps(ft):
                            warn(f'Detected overlap between {probe.meta.name} and {ft.meta.name}')
                        ft = deepcopy(ft)
                        ft.meta.strand_check = ft.loc.strand == '+'
                        primer1.append(ft)
            primer2 = []
            for ft in fts[i+1:]:   # look for primers after probe on the - strand
                if ft.type == 'primer':
                    if probe.distance(ft) > distance:
                        break
                    if probe.overlaps(ft):
                        warn(f'Detected overlap between {probe.meta.name} and {ft.meta.name}')
                    ft = deepcopy(ft)
                    ft.meta.strand_check = ft.loc.strand == '-'
                    primer2.append(ft)
            # depending on amplification add probe-primer singlets/doublets/triplets
            if len(primer1) + len(primer2) == 0:
                r.append((probe.meta.name, (probe, )))
            elif len(primer1) == 0:
                for ft in primer2:
                    # for linear amplification correct strand check for probe primer pairs located on the same strand
                    if probe.loc.strand == ft.loc.strand:
                        ft.meta.strand_check = False
                    r.append((probe.meta.name, (probe, ft)))
            elif len(primer2) == 0:
                for ft in primer1:
                    # for linear amplification correct strand check for probe primer pairs located on the same strand
                    if probe.loc.strand == ft.loc.strand:
                        ft.meta.strand_check = False
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
                if ft2.distance(ft) > distance:
                    break
                if ft2.loc.strand == '+':
                    # second primer should be on - strand
                    ft2 = deepcopy(ft2)
                    ft2.meta.strand_check = False
                else:
                    done.add(ft2.loc)
                found_pair = True
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


def _strand_check(r):
    return all(getattr(ft.meta, 'strand_check', True) for ft in r)

def _growth(r):
    return sum(getattr(ft.meta, 'strand_check', True) for ft in r)

def _sum_mismatch(r):
    return sum(ft.meta.mismatch for ft in r)

def _sortkey(result):
    """Key to sort by goodness"""
    name, r = result
    return name, -_growth(r), _sum_mismatch(r), r[0].distance(r[-1])


def _groupby(ft):
    return ft.seqid.split('--')[0]


def find_probe_primer(allfts, mismatch=2, only_primer=False, distance=150):
    """Find primer-probe-primer triplets or primer pairs"""
    allfts = allfts.select(type_in=('primer', 'probe'), mismatch_le=mismatch).sort()
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


def output_assay_overview(results, out, query_ids=None, source_ids=None, only_primer=False):
    """Write the assay overview file"""
    num2str = NUM2STR_PRIMER if only_primer else NUM2STR_PROBE
    lines = []
    all_probes = sorted(query_ids or set(proben for res in results.values() for proben, _ in res))
    lines.append('Genome\t' + '\t'.join(all_probes))
    for superseqid in (source_ids or results):
        probes = {}
        if superseqid in results:
            for proben, r in results[superseqid]:
                if proben not in probes:
                    if only_primer:
                        primerstr = ' '.join(f'{"F" if p.loc.strand == "+" else "R"}{p.meta.mismatch}{p.loc.strand}' for p in r)
                    else:
                        primerstr = ' '.join(f'{"P" if p.type == "probe" else "F" if i == 0 else "R"}{p.meta.mismatch}{p.loc.strand}' for i, p in enumerate(r))
                    # growth = num2str[_growth(r)]
                    # mms = ''.join(str(p.meta.mismatch) for p in r)
                    # primerstr = ''.join('P' if p.type == 'probe' else '>' if p.loc.strand == '+' else '<' for p in r)
                    probes[proben] = f'{num2str[_growth(r)]} {primerstr}'
        line = f'{superseqid}\t' + '\t'.join(probes.get(proben, NUM2STR_PROBE[0]) for proben in all_probes)
        lines.append(line)
    with open(out, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def output_assay_details(results, out, source_ids=None, only_primer=False, verbose=False, zero_based_numbering=False):
    """Write the assay details file"""
    num2str = NUM2STR_PRIMER if only_primer else NUM2STR_PROBE
    header = (f'Genome\t{"Primer" if only_primer else "Probe"}\tAmplification\t'
              'Pairing and Strand direction ((+)/(-))\tStrand direction\tStrand check\tMismatch\tDistance\tLocation\n')
    lines = [header]
    for superseqid in (source_ids or results):
        max_growth = 0
        if superseqid in results:
            for proben, pps in results[superseqid]:
                max_growth = max(max_growth, len(pps))
                growth = num2str[_growth(pps)]
                names = ','.join(f'{p.meta.name} ({p.loc.strand})' for p in pps)
                strands = ''.join(p.loc.strand for p in pps)
                strand_check = 'pass' if _strand_check(pps) else 'FAIL'
                mms = ','.join(str(p.meta.mismatch) + '*' * p.meta.warn_end for p in pps)
                dists = ','.join(str(p1.distance(p2)) for p1, p2 in zip(pps[1:], pps[:-1]))
                pos = ','.join(f'{p.loc.start+(not zero_based_numbering)}:{p.loc.stop}:{p.loc.strand}' for p in pps)
                contig = pps[0].seqid.split('--')[-1]
                contig = '' if contig == superseqid else contig + ','
                lines.append(f'{superseqid}\t{proben}\t{growth}\t{names}\t{strands}\t{strand_check}\tM{mms}\tD{dists}\t{contig}{pos}\n')
        else:
            lines.append(f'{superseqid}\t\t{num2str[0]}\t\t\t\tM\tD\t\n')
        if verbose:
            print(f'{superseqid}\t{num2str[max_growth]}')
    with open(out, 'w') as f:
        f.write(''.join(lines))


def find_probe_primer_cli(fname, out=None, only_primer=False, zero_based_numbering=False, **kw):
    """Read BLAST file or GFF file and find+report probes/primers"""
    try:
        from sugar import read, read_fts
    except ImportError:
        import sys
        sys.exit('Script needs rnajena-sugar. Install with:\npip install rnajena-sugar')
    fts = read_fts(fname, comments=(comments:=[]))
    query_ids = ([line.removeprefix(pre).split() for line in comments if line.startswith(pre := '# queryids:')] + [None])[0]
    if query_ids:
        query_ids = sorted({qid for qid in query_ids if 'probe' in qid.lower()})
    source_ids = ([line.removeprefix(pre).split() for line in comments if line.startswith(pre := '# sourceids:')] + [None])[0]

    fmt = fts[0].meta._fmt
    print(f'Successfully parsed {fmt.upper()} file at {fname}.')
    for ft in fts:
        if ft.type is not None:
            ft.type = ft.type.lower()
        if ft.type not in ('probe', 'primer'):
            if 'probe' in ft.meta.get('name'):
                ft.type = 'probe'
            elif 'primer' in ft.meta.get('name'):
                ft.type = 'primer'
            else:
                msg = f'At least one {fmt.upper()} hit is missing keyword "primer" or "probe"'
                catchy_no_kw_warning_msg = f'\n\n{len(msg) * "#"}\n{msg}\n{len(msg) * "#"}\n'
                warn(catchy_no_kw_warning_msg)
        try:
            ft.meta.mismatch = ft.meta[f'_{fmt}'].mismatch
        except KeyError:
            warn(f'No mismatch information found in {fmt.upper()} file, use mismatch 0')
            ft.meta.mismatch = 0
        ft.meta.warn_end = False
        try:
            qlen = ft.meta[f'_{fmt}'].qlen
            qstart = ft.meta[f'_{fmt}'].qstart
            qend = ft.meta[f'_{fmt}'].qend
        except KeyError:
            warn(f'No qstart, qend or qlen information found in {fmt.upper()} file, cannot check for gaps at query ends')
        else:
            if len(ft) > qlen:
                raise ValueError('Detected a gap, please contact developers')
            elif len(ft) < qlen:
                ft.meta.mismatch += qlen - len(ft)
                # warn for mismatches at 3' primer ends or 5' probe ends
                warn_end = ft.type == 'probe' and qstart > 1 or ft.type == 'primer' and qend < qlen
                if warn_end:
                    ft.meta.warn_end = True
                    which_end = 5 if ft.type == 'probe' else 3
                    msg = (f'Detected mismatch at {which_end}\' end of {ft.type}. '
                           'Will be marked with "*" in the details table.')
                    warn(msg)
    results = find_probe_primer(fts, only_primer=only_primer, **kw)
    if out is None:
        fname = Path(fname)
        out = fname.with_name(fname.stem + '_assay' + '_only_primer' * only_primer)
    out = Path(out)
    out1 = out.with_name(out.stem + '_overview.tsv')
    out2 = out.with_name(out.stem + '_details.tsv')
    output_assay_overview(results, out1, query_ids=None if only_primer else query_ids, source_ids=source_ids, only_primer=only_primer)
    print(f'Assay overview file created at {out1}.')
    output_assay_details(results, out2, source_ids=source_ids, only_primer=only_primer, zero_based_numbering=zero_based_numbering)
    print(f'Assay details file created at {out2}.')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fname', help='BLAST outfile outfmt 7, 6 or 10 or GFF file, e.g. created with assay_blast')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument("--mismatch", type=int, default=2, help='Maximum allowed mismatches (default: 2)')
    parser.add_argument("--distance", type=int, default=250, help='Distance threshold (bp) between adjacent oligos (default: 250)')
    parser.add_argument('--only-primer', action='store_true', help='Find primer pairs instead of primer-probe-primer triplets')
    parser.add_argument('-o', '--out', help='Prefix of output files, by default derived from fname')
    msg = 'By default the locations in the output are given by one-based numbering. Switch to zero-based numbering.'
    parser.add_argument('--zero-based-numbering', action='store_true', help=msg)

    args = vars(parser.parse_args())
    find_probe_primer_cli(**args)


if __name__ == "__main__":
    main()
