#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs=1, help='Input contig sequences in fasta')
    parser.add_argument('-s', nargs=1, help='Contig names to subset')
    parser.add_argument('-o', nargs=1, help='Output contig sequences in fasta')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def write_seq(seq, seq_id, s, o):
    if seq and seq_id in s:
        o.write('>%s\n%s\n' % (seq_id, seq))


def subset_contigs(filin, filou, names):
    s = set([x.strip() for x in open(names).readlines()])
    with open(filou, 'w') as o:
        if filin.endswith('gz'):
            func = gzip.open(filin, 'rt')
        else:
            func = open(filin)
        with func as f:
            seq, seq_id = '', ''
            for line in f:
                if line[0] == '>':
                    write_seq(seq, seq_id, s, o)
                    seq_id = line[1:].strip().split()[0]
                else:
                    seq += line.strip()
        write_seq(seq, seq_id, s, o)


if __name__ == '__main__':
    args = get_args()
    subset_contigs(args['i'][0], args['o'][0], args['s'][0])
