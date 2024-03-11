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
    parser.add_argument('-o', nargs=1, help='Output fasta file')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def reformat(filin, filou):
    with open(filou, 'w') as o:
        if filin.endswith('gz'):
            func = gzip.open(filin, 'rt')
        else:
            func = open(filin)
        with func as f:
            for line in f:
                if line[0] == '>':
                    o.write('%s\n' % line.strip().split()[0])
                else:
                    o.write(line)


if __name__ == '__main__':
    args = get_args()
    reformat(args['i'][0], args['o'][0])
