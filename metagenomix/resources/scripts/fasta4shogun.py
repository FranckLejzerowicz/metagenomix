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
    parser.add_argument('-i', nargs=1, required=True, help='fasta file')
    parser.add_argument('-o', nargs=1, required=True, help='output fasta file')
    parser.add_argument('-s', nargs=1, required=True, help='sample name')
    parser.add_argument('-r', nargs='*', choices=['1', '2'], help='orientation')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def write_formatted(filin, filou, sample, orient):
    with open(filou, 'w') as o, open(filin) as f:
        idx = 0
        for line in f:
            if line[0] == '>':
                idx += 1
                o.write('>%s_%s%s %s' % (sample, idx, orient, line[1:]))
            else:
                o.write(line)


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    sample = args['s'][0]
    orient = ''
    if 'r' in args:
        orient = '/%s' % args['r'][0]
    write_formatted(filin, filou, sample, orient)
