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
    parser.add_argument('-i', nargs=1, help='Input fasta file')
    parser.add_argument('-o', nargs=1, help='Output table')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def write_seq(seq, seq_id, o):
    if seq:
        o.write('%s\t%s\n' % (seq_id, len(seq)))


def measure_length(filin, filou):
    with open(filou, 'w') as o:
        if filin.endswith('gz'):
            func = gzip.open(filin, 'rt')
        else:
            func = open(filin)
        with func as f:
            seq, seq_id = '', ''
            for line in f:
                if line[0] == '>':
                    write_seq(seq, seq_id, o)
                    seq_id = line[1:].strip().split()[0]
                    seq = ''
                else:
                    seq += line.strip()
        write_seq(seq, seq_id, o)


if __name__ == '__main__':
    args = get_args()
    measure_length(args['i'][0], args['o'][0])
