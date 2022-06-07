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
    parser.add_argument('-i', nargs=1, required=True, help='sam file')
    parser.add_argument('-o', nargs=1, required=True, help='output sam file')
    parser.add_argument('-s', nargs=1, required=True, help='sample name')
    parse = parser.parse_args()
    args = vars(parse)
    return args


def write_formatted(filin, filou, sample):
    with open(filou, 'w') as o, open(filin) as f:
        idx = 0
        for line in f:
            if line[0] == '@':
                continue
            idx += 1
            ls = line.split()
            o.write('%s_%s_%s\t%s\n' % (sample, idx, ls[0], ls[-1]))


if __name__ == '__main__':
    args = get_args()
    filin = args['i'][0]
    filou = args['o'][0]
    sample = args['s'][0]
    write_formatted(filin, filou, sample)
