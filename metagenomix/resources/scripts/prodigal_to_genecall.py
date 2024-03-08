# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse
from skbio.io import read


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', nargs=1, required=True, help='input gff3 file')
    parser.add_argument(
        '-f', nargs=1, required=True, help='input fasta file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output file')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_r(r):
    if int(r):
        r = '1'
    else:
        r = '0'
    return r


def get_direction(pl):
    if pl == '+':
        di = 'f'
    else:
        di = 'r'
    return di


def reformat(gff, fasta, filou):
    d = load_seqs(fasta)
    head = 'gene_callers_id\tcontig\tstart\tstop\tdirection'
    head += '\tpartial\tcall_type\tsource\tversion\taa_sequence\n'
    with open(filou, 'w') as o, open(gff) as f:
        o.write(head)
        i = 0
        for line in f:
            if line[0] == '#':
                continue
            i += 1
            na, tv, cd, st, en, sc, pl, ph, at = line.strip().split('\t')
            di = get_direction(pl)
            t, v = tv.split('_', 1)
            me = at.split(';')[0].split('_')[-1]
            pa = get_r(at.split(';')[1].split('=')[-1])
            o.write('\t'.join([
                i, na, str(int(st)-1), en, di, pa, '1',
                t.lower(), v, d['%s_%s' % (na, me)]
            ]))


def load_seqs(fasta):
    fasta_d = {}
    for entry in read(fasta, format="fasta"):
        fasta_d[entry.metadata['id']] = str(entry)
    return fasta_d


if __name__ == '__main__':
    args = get_args()
    reformat(args['g'][0], args['f'][0], args['o'][0])
