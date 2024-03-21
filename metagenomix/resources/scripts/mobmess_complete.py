# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', nargs=1, required=True, help='input fasta file')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output complete file for mobmess')
    parser.add_argument(
        '-n', nargs="*", required=False, help='File with target contigs names')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def make_complete(filin, filou, names):
    contigs = set()
    if names:
        contigs.update([
            x.split()[0].strip() for x in open(names[0]).readlines()])
    with open(filou, 'w') as o, open(filin) as f:
        for line in f:
            if line[0] == ">":
                contig = line[1:].strip().split()[0]
                if contigs:
                    if contig in contigs:
                        o.write('%s\t1\n' % contig)
                else:
                    o.write('%s\t1\n' % contig)


if __name__ == '__main__':
    args = get_args()
    filin = args['f'][0]
    filou = args['o'][0]
    names = args['n']
    make_complete(filin, filou, names)
