#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
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
        '-i', nargs=1, required=True, help='fasta file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output fasta file name')
    parser.add_argument(
        '-t', nargs=1, required=True, help='length threshold in bases')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def do_filtering(filin, filou, thresh):
    kept = 0
    removed = 0
    with open(filou, 'w') as o:
        for entry in read(filin, format='fasta'):
            if len(str(entry)) < thresh:
                removed += 1
                continue
            else:
                kept += 1
                o.write('>%s\n%s\n' % (entry.metadata['id'], str(entry)))
    print('Filter on length %s:' % thresh)
    print(' - removed:', removed)
    print(' - kept:', kept)


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    thresh = int(args['t'][0])
    do_filtering(filin, filou, thresh)