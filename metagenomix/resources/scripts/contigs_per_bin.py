#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob
import argparse
from os.path import basename, splitext


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs=1, help='Output of a binning tool')
    parser.add_argument('-b', nargs=1, help='Binning tool')
    parser.add_argument('-o', nargs=1, help='Output file')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_contigs_per_bin(filin, filou, b):
    if b == 'metawrap_refine':
        with open(filou, 'w') as o:
            for bin_fp in glob.glob('%s/*.fa' % filin):
                bin = splitext(basename(bin_fp))[0]
                with open(bin_fp) as f:
                    for line in f:
                        if line[0] == '>':
                            o.write('%s\t%s\n' % (line.strip()[1:], bin))
    else:
        sys.exit('[contigs_per_bin.py] Binning tool "%s" not supported yet' % b)


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    b = args['b'][0]
    get_contigs_per_bin(filin, filou, b)