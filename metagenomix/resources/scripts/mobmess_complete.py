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
        '-n', nargs="?", required=False, help='File with target contigs names')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_contigs_per_group(names_fp=''):
    contigs_per_group = {}
    if names_fp:
        names_pd = pd.read_table(names_fp)
        contigs_per_group = names_pd.groupby('sample_group').apply(
            lambda x: x['contig'].tolist()).to_dict()
    return contigs_per_group


def make_complete(filin, fastou, filou, names_fp):
    contigs_per_group = get_contigs_per_group(names_fp)
    fastas_pd = pd.read_csv(filin)
    with open(fastou, 'w') as o1, open(filou, 'w') as o2:
        for group, group_pd in fastas_pd.groupby('group'):
            fasta = group_pd['filepath'].item()
            with gzip.open(fasta) as f:
                for line_ in f:
                    line = line_.decode()
                    if line[0] == ">":
                        write = False
                        contig = line[1:].strip().split()[0]
                        if names_fp:
                            if contig in contigs_per_group.get(group, {}):
                                write = True
                                o2.write('%s\t1\n' % contig)
                        else:
                            write = True
                            o2.write('%s\t1\n' % contig)
                    if write:
                        o1.write(line)


if __name__ == '__main__':
    args = get_args()
    filin = args['f'][0]
    fastou = args['o'][0]
    filou = args['c'][0]
    names_fp = args['n']
    make_complete(filin, fastou, filou, names_fp)
