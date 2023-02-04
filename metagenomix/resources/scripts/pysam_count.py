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


def assembling_counts(tab):
    dat = []
    for vdx, vals in enumerate(tab):
        if not vdx:
            continue
        alignment = pysam.AlignmentFile(vals[1], "rb")
        with open(vals[0]) as f:
            for line in f:
                if line[0] == '>':
                    contig = line[1:].split()[0]
                    reads = alignment.count(contig)
                    dat.append((vals[2:] + [contig, str(reads)]))
    return dat


def prodigal_counts(tab):
    dat = []
    for vals in tab:
        alignment = pysam.AlignmentFile(vals[1], "rb")
        with open(vals[0]) as f:
            for line in f:
                if line[0] == '>':
                    contig_split = line[1:].split()
                    gene = contig_split[0]
                    contig = gene.rsplit('_')[0]
                    start = contig_split[2]
                    end = contig_split[4]
                    reads = alignment.count(contig, int(start), int(end))
                    dat.append((vals[2:] + [gene, str(reads)]))
    return dat


def binning_counts(tab):
    pass
    # dat = []
    # for vals in tab:
    #     alignment = pysam.AlignmentFile(vals[1], "rb")
    #     with open(vals[0]) as f:
    #         for line in f:
    #             if line[0] == '>':
    #                 contig = line[1:].split()[0]
    #                 reads = alignment.count(contig)
    #                 dat.append((vals[2:] + [contig, str(reads)]))
    # with open(filou, 'w') as o:
    #     o.write('tech\tsample\tali\ttarget\tprev\tmode\tcontig\treads\n')
    #     for row in dat:
    #         o.write('%s\n' % '\t'.join(row))
    # return dat


def write_output(dat, unit, filou):
    head = 'tech\tsample\tali\tcoassembly\tgroup\t'
    head += 'target\tprev\tmode\t%s\treads\n' % unit
    with open(filou, 'w') as o:
        o.write(head)
        for row in dat:
            o.write('%s\n' % '\t'.join(row))


def count_reads(filin, filou):
    tab = [x.strip().split('\t') for x in open(filin).readlines()]
    if tab[-1][-1] == 'assembling':
        unit = 'contig'
        dat = assembling_counts(tab)
    elif tab[-1][-1] == 'prodigal':
        unit = 'gene'
        dat = prodigal_counts(tab)
    elif tab[-1][-1] == 'binning':
        unit = 'bin'
        dat = binning_counts(tab)
    else:
        print(skdfkjb)
    write_output(dat, unit, filou)


if __name__ == '__main__':
    args = get_args()
    count_reads(args['i'][0], args['o'][0])
