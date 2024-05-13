# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', nargs=1, required=True, help='path to the fasta files')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output fasta file (mobmess input)')
    parser.add_argument(
        '-c', nargs=1, required=True, help='plasmid boolean (mobmess complete)')
    parser.add_argument(
        '-p', nargs="?", required=False, help='File with circular contig names')
    parser.add_argument(
        '-n', nargs="?", required=False, help='File with target contig names')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_circles(circs_fp=''):
    circles_per_group = {}
    if circs_fp:
        circs_pd = pd.read_csv(circs_fp)
        circles_per_group = circs_pd.groupby('group').apply(
            lambda x: x['contig'].tolist()).to_dict()
    return circles_per_group


def get_contigs(names_fp=''):
    contigs_per_group = {}
    if names_fp:
        names_pd = pd.read_csv(names_fp)
        contigs_per_group = names_pd.groupby('group').apply(
            lambda x: x['contig'].tolist()).to_dict()
    return contigs_per_group


def write_seq(circs_fp, circles, seq, grp, o2):
    if circs_fp:
        if seq in circles.get(grp, {}):
            o2.write('%s__%s\t1\n' % (grp, seq))
        else:
            o2.write('%s__%s\t0\n' % (grp, seq))
    else:
        o2.write('%s__%s\t1\n' % (grp, seq))


def make_complete(
        filin: str,
        fastou: str,
        filou: str,
        circs_fp: str,
        names_fp: str
):
    circles, contigs = get_circles(circs_fp), get_contigs(names_fp)
    fastas_pd = pd.read_csv(filin)
    with open(fastou, 'w') as o1, open(filou, 'w') as o2:
        for grp, grp_pd in fastas_pd.groupby('group'):
            fas = grp_pd['filepath'].item()
            with gzip.open(fas) as f:
                for line_ in f:
                    line = line_.decode()
                    if line[0] == ">":
                        write = False
                        seq = line[1:].strip().split()[0]
                        if names_fp:
                            if seq in contigs.get(grp, {}):
                                write = True
                                write_seq(circs_fp, circles, seq, grp, o2)
                        else:
                            write = True
                            write_seq(circs_fp, circles, seq, grp, o2)
                    if write:
                        o1.write(line)


if __name__ == '__main__':
    args = get_args()
    make_complete(
        args['f'][0],
        args['o'][0],
        args['c'][0],
        args['p'],
        args['n'])
