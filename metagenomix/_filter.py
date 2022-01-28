# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def filter_table(preval: str, abund: str,
                 tsv_pd: pd.DataFrame) -> pd.DataFrame:
    preval = float(preval)
    abund = float(abund)
    if preval + abund == 0:
        return tsv_pd
    tsv_filt_pd = tsv_pd.copy()
    # get the min number of samples based on prevalence percent
    if preval < 1:
        n_perc = tsv_pd.shape[1] * preval
    else:
        n_perc = preval
    # abundance filter in terms of min reads counts
    tsv_pd_perc = tsv_filt_pd.copy()
    if abund < 1:
        tsv_pd_perc = tsv_filt_pd / tsv_filt_pd.sum()
    tsv_filt_pd = tsv_filt_pd.loc[(tsv_pd_perc > abund).sum(1) > n_perc, :]
    tsv_filt_pd = tsv_filt_pd.loc[tsv_filt_pd.sum(1) > 0,
                                  tsv_filt_pd.sum(0) > 0]
    return tsv_filt_pd


def write_filtered_tsv(tsv_out: str, tsv_pd: pd.DataFrame) -> None:
    tsv_sams_col = tsv_pd.reset_index().columns[0]
    tsv_pd = tsv_pd.reset_index().rename(
        columns={tsv_sams_col: 'Feature ID'}).set_index('Feature ID')
    tsv_pd.reset_index().to_csv(tsv_out, index=False, sep='\t')


def get_unique_filterings(songbird_filtering):
    unique_filterings = {}
    for filt_name, dat_d in songbird_filtering[''].items():
        for dat, (preval, abund) in dat_d.items():
            unique_filterings.setdefault(dat, set()).add(
                (filt_name, preval, abund))
    return unique_filterings
