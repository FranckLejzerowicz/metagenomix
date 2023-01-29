# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pysam
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True, help='input table file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output file')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def assembling_counts(tab, filou):
    dat = []
    for vals in tab:
        alignment = pysam.AlignmentFile(vals[1], "rb")
        with open(vals[0]) as f:
            for line in f:
                if line[0] == '>':
                    contig = line[1:].split()[0]
                    reads = alignment.count(contig)
                    dat.append((vals[2:] + [contig, str(reads)]))
    with open(filou, 'w') as o:
        o.write('tech\tsample\tali\ttarget\tprev\tmode\tcontig\treads\n')
        for row in dat:
            o.write('%s\n' % '\t'.join(row))
    return dat


def count_reads(filin, filou):
    tab = [x.strip() for x in open(filin).readlines()]
    if tab[-1][-1] == 'assembling':
        assembling_counts(tab, filou)


if __name__ == '__main__':
    args = get_args()
    count_reads(args['i'][0], args['o'][0])
