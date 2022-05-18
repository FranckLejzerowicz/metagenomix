# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs=1, required=True, help='Upper triangular')
    parse = parser.parse_args()
    arguments = vars(parse)
    return arguments


def get_samples(mat_fp: str) -> list:
    """Get the sample names which are in the
    first line of the upper triangular matrix.

    Parameters
    ----------
    mat_fp : str
        Input path to the upper triangular matrix.

    Returns
    -------
    first_line : list
        Sample names.
    """
    with open(mat_fp) as f:
        for ldx, line in enumerate(f):
            first_line = line.strip().split(';')[1:]
            return first_line


def get_matrix_in_pandas(mat_fp: str) -> pd.DataFrame:
    """Read the upper triangular matrix into pandas.

    Parameters
    ----------
    mat_fp : str
        Input path to the upper triangular matrix.

    Returns
    -------
    mat_pd : pd.DataFrame
        Upper triangular matrix table.
    """
    print('reading matrix in pandas...')
    mat_pd = pd.read_csv(
        mat_fp,
        skiprows=1,
        header=None,
        index_col=0,
        sep=';',
        low_memory=False,
        dtype=str
    )
    return mat_pd


def write_matrix_from_pandas(mat_pd: pd.DataFrame, mat_fp: str, samples: list):
    """Write the matrix as symmetric table.

    Parameters
    ----------
    mat_pd : pd.DataFrame
        Upper triangular matrix table.
    mat_fp : str
        Input path to the upper triangular matrix.
    samples : list
        Sample names.
    """
    mat_o = mat_fp.replace('.csv', '_sym.tsv')
    with open(mat_o, 'w') as o:
        o.write('\t%s\n' % '\t'.join(samples))
        for idx, row in enumerate(mat_pd.values):
            for jdx in range(len(row)):
                if jdx < idx:
                    row[jdx] = mat_pd.iloc[jdx, idx]
            o.write('%s\t%s\n' % (samples[idx], '\t'.join(row)))


if __name__ == '__main__':
    args = get_args()
    fp = args['i'][0]
    samples = get_samples(fp)
    mat_pd = get_matrix_in_pandas(fp)
    write_matrix_from_pandas(mat_pd, fp, samples)
