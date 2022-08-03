#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True, help='fasta file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output fasta file name')
    parser.add_argument(
        '--n', action='store_true', default=False,
        help='instead add a numeric id')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def do_edit(filin, filou, N):
    idx = 0
    with open(filou, 'w') as o, open(filin) as f:
        for line in f:
            if line[0] == '>':
                idx += 1
                if N:
                    o.write('%s%s\n' % (line.strip().split()[0], idx))
                else:
                    o.write(line.replace(' ', ''))
            else:
                o.write(line)


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    N = args['n']
    do_edit(filin, filou, N)
