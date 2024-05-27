# ----------------------------------------------------------------------------
# Copyright (c) 2023, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', nargs=1, required=True, help='path to the fasta files')
    parser.add_argument(
        '-o', nargs=1, required=True, help='output fasta file (mobmess input)')
    parser.add_argument(
        '-c', nargs=1, required=True, help='plasmid boolean (mobmess complete)')
    parser.add_argument(
        '-s', nargs=1, required=False, type=int, help='Minimum contig length')
    parser.add_argument(
        '-p', nargs="?", required=False, help='File with circular contig names')
    parser.add_argument(
        '-n', nargs="?", required=False, help='File with target contig names')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def get_circles(circs_fp=''):
    circles = {}
    if circs_fp:
        circs_pd = pd.read_csv(circs_fp)
        for group, group_pd in circs_pd.groupby('group'):
            group_fp = group_pd['filepath'].item()
            circles[group] = [x.strip() for x in open(group_fp).readlines()]
    return circles


def get_contigs(names_fp=''):
    contigs = {}
    if names_fp:
        names_pd = pd.read_table(names_fp)
        contigs = names_pd.groupby('sample_name').apply(
            lambda x: x['contig'].tolist()).to_dict()
    return contigs


def write_seq(circ_fp, circs, name, group, seq, o1, o2):
    new_name = '%s__%s' % (group, name)
    if circ_fp:
        if name in circs.get(group, {}):
            o2.write('%s\t1\n' % new_name)
        else:
            o2.write('%s\t0\n' % new_name)
    else:
        o2.write('%s\t1\n' % new_name)
    o1.write('%s\n%s\n' % (new_name, seq))


def check_circ(contigs, name, group, name_fp):
    if name_fp:
        if name in contigs.get(group, {}):
            return True
        else:
            return False
    else:
        return True


def prep(
        filin: str,
        fastou: str,
        filou: str,
        minlen: int,
        circ_fp: str,
        name_fp: str
):
    circs = get_circles(circ_fp)
    contigs = get_contigs(name_fp)
    fastas_pd = pd.read_csv(filin)
    with open(fastou, 'w') as o1, open(filou, 'w') as o2:
        for group, group_pd in fastas_pd.groupby('group'):
            fasta = group_pd['filepath'].item()
            seq = ''
            with gzip.open(fasta) as f:
                for line_ in f:
                    line = line_.decode().strip()
                    if line[0] == ">":
                        write = False
                        name = line[1:].split()[0]
                        if not (seq and len(seq) >= minlen):
                            continue
                        write = check_circ(contigs, name, group, name_fp)
                        if write:
                            write_seq(circ_fp, circs, name, group, seq, o1, o2)
                        seq = ''
                    else:
                        seq += line
            if write:
                write_seq(circ_fp, circs, name, group, seq, o1, o2)


if __name__ == '__main__':
    args = get_args()
    prep(
        args['f'][0],
        args['o'][0],
        args['c'][0],
        args['s'][0],
        args['p'],
        args['n'])
