# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd


def read_metadata(
        meta_fp: str
) -> pd.DataFrame:
    """Read the metadata passed by the user.

    Parameters
    ----------
    meta_fp : str
        Metadata file path.

    Returns
    -------
    meta : pd.DataFrame
        Metadata table.
    """
    first_col = get_first_column(meta_fp)
    meta = pd.read_csv(meta_fp, dtype=str, sep='\t', low_memory=False)
    meta.rename(columns={first_col: 'sample_name'}, inplace=True)
    return meta


def get_first_column(
        meta_fp: str
) -> str:
    """Get the first column of the metadata file.

    Parameters
    ----------
    meta_fp : str
        Metadata file path.

    Returns
    -------
    first_column : str
        Name of the first column of the metadata input file.
    """
    line = ''
    with open(meta_fp) as f:
        for line in f:
            break
    if line.split():
        first_column = line.split()[0]
        return first_column
    else:
        raise IOError('File "%s" is empty' % meta_fp)


def get_cur_type(
        cur_input: str
) -> str:
    """Get the string "nucl" or "prot" depending
    on the type of data in the fasta file.

    Parameters
    ----------
    cur_input : str

    Returns
    -------
    typ : str
        Type of sequence data.
    """
    len_set = []
    with open(cur_input) as f:
        for ldx, line in enumerate(f):
            if line.startswith('>'):
                continue
            len_set.append(len(set(line.strip())))
            if ldx > 15:
                break
    if np.mean(len_set) < 6:
        typ = 'nucl'
    else:
        typ = 'prot'
    return typ

