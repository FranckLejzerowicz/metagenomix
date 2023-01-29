# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pysam
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True, help='input table file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output file')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def assembling_counts(tab):
    dat = []
    for vals in tab.values:
        alignment = pysam.AlignmentFile(vals[1], "rb")
        with open(vals[0]) as f:
            for line in f:
                if line[0] == '>':
                    contig = line[1:].split()[0]
                    reads = alignment.count(contig)
                    dat.append((vals[2:] + [contig, reads]))
    dat = pd.DataFrame(dat, columns=['tech', 'sample', 'ali', 'target',
                                     'prev', 'mode', 'contig', 'reads'])
    return dat


def count_reads(filin, filou):
    tab = pd.read_table(filin)
    mode = tab['mode'].unique().tolist()[0]
    if mode == 'assembling':
        assembling_counts(tab).write(filou, index=False, sep='\t')


if __name__ == '__main__':
    args = get_args()
    count_reads(args['i'][0], args['o'][0])
