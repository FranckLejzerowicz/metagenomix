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
from os.path import splitext
from numpy import array_split


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', nargs=1, required=True, help='fasta file')
    parser.add_argument(
        '-c', nargs=1, type=int, required=False, help='number of chunks')
    parser.add_argument(
        '--no', action='store_true', default=False, help='no chunk, just bed')
    parser.add_argument(
        '--array', action='store_true', default=False, help='keep array number')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def make_bed(fasta, chunks, bed_only, keep_array):
    if bed_only:
        one_bed(fasta, keep_array)
    else:
        multi_beds(fasta, chunks)


def one_bed(fasta, keep_array):
    if keep_array:
        bed_fpo = "%s.bed.%s" % tuple(splitext(fasta))
    else:
        bed_fpo = "%s.bed" % splitext(fasta)[0]
    with open(bed_fpo, 'w') as o:
        for e in read(fasta, format="fasta"):
            o.write("%s\t0\t%s\n" % (e.metadata['id'], len(e)))


def multi_beds(fasta, chunks):
    seqs = {e.metadata['id']: len(e) for e in read(fasta, format="fasta")}
    seqs_chunks = array_split(seqs, chunks)
    for sdx, seqs_chunk in enumerate(seqs_chunks):
        bed_fpo = "%s_chunk%s.bed" % (splitext(fasta)[0], sdx)
        with open(bed_fpo, 'w') as o:
            for seq in seqs_chunk:
                o.write("%s\t0\t%s\n" % (seq, seqs[seq]))


if __name__ == '__main__':
    args = get_args()
    fasta = args['i'][0]
    chunks = args['c'][0]
    bed_only = args['no']
    keep_array = args['array']
    make_bed(fasta, chunks, bed_only, keep_array)
