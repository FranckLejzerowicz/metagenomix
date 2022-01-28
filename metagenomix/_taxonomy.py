# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from os.path import isfile

from prep_songbird._io_utils import (
    get_taxonomy_classifier, get_analysis_folder, parse_g2lineage)
from prep_songbird._cmds import (
    write_seqs_fasta, write_taxonomy_sklearn, run_export, run_import)


def edit_split_taxonomy(
        ranks: dict, split_taxa_pd: pd.DataFrame) -> pd.DataFrame:
    if len(ranks) == split_taxa_pd.shape[1]:
        split_taxa_pd = split_taxa_pd.rename(columns=ranks)
    else:
        alpha = 'ABCDEFGHIJKLMNOPQRST'
        cols = [alpha[x] for x in range(split_taxa_pd.shape[1])]
        split_taxa_pd = pd.DataFrame([
            pd.Series([
                '%s__%s' % (cols[idx], str(x).replace(' ', '_'))
                for idx, x in enumerate(row) if str(x) != 'nan']
            ) for row in split_taxa_pd.values],
            columns=cols,
            index=split_taxa_pd.index.tolist()
        )
    return split_taxa_pd


def get_ranks_from_split_taxonomy(split_taxa_pd, col):
    rank = [x.split('__')[0] for x in split_taxa_pd[col]
            if str(x) not in ['nan', 'None', 'Unassigned']]
    if len(set(rank)) == 1:
        return rank[0]
    else:
        rank = [x.split('_')[0] for x in split_taxa_pd[col]
                if str(x) not in ['nan', 'None', 'Unassigned']]
        if len(set(rank)) == 1:
            return rank[0]
    return ''


def parse_split_taxonomy(split_taxa_pd: pd.DataFrame) -> dict:
    torm = []
    ranks = {}
    # not_collapsable = False
    # parse each column/levels to determine which is
    # (i) to be removed or (ii) labelled with consistent rank (e.g. "k__"))
    for col in split_taxa_pd.columns:
        # if some levels of the split taxonomy are the feature IDs themselves:
        # remove these levels (as it would be a copy of the non-collapsed table)
        col_features = split_taxa_pd[col].tolist()
        if split_taxa_pd.index.tolist() == col_features:
            torm.append(col)
        else:
            rank = get_ranks_from_split_taxonomy(split_taxa_pd, col)
            if rank and rank != col:
                ranks[col] = rank
    # remove columns to be removed
    if torm:
        split_taxa_pd.drop(columns=torm, inplace=True)
    return ranks


def get_split_taxonomy(taxa, taxo_sep=';'):
    split_lens = set()
    split_taxa = []
    for tdx, taxon in enumerate(taxa):
        if str(taxon) == 'nan':
            split_lens.add(1)
            split_taxa.append(pd.Series(['Unassigned']))
        else:
            taxon_split = [x.strip() for x in str(taxon).split(taxo_sep)
                           if len(x.strip()) and not x.startswith('x__')]
            split_lens.add(len(taxon_split))
            split_taxa.append(pd.Series(taxon_split))
    # if the parsed and split taxonomies have  very variable number of fields
    # or very long split results  it is terminated here as a taxonomy that
    # does not make sense
    if len(split_lens) > 15 or max(split_lens) > 15:
        return pd.DataFrame([[x] for x in taxa], columns=['not_really_taxon'])
    # build a dataframe from the spit taxonomy
    split_taxa_pd = pd.DataFrame(split_taxa)
    # add a column label
    ALPHA = 'ABCDEFGHIJKLMNOPQRST'
    split_taxa_pd.columns = ['Taxolevel_%s' % (ALPHA[idx])
                             for idx in range(split_taxa_pd.shape[1])]
    return split_taxa_pd


def get_tax_tables(tax_fp: str) -> tuple:
    # read taxonomy with features as index, format and collect features IDs list
    tax_pd = pd.read_csv(tax_fp, header=0, sep='\t', dtype=str)
    tax_pd.rename(columns={tax_pd.columns[0]: 'Feature ID'}, inplace=True)
    # perform taxonomic split on the Taxon list and give the Feature as index
    split_taxa_pd = get_split_taxonomy(tax_pd.Taxon.tolist())
    split_taxa_pd.index = tax_pd['Feature ID'].tolist()
    return tax_pd, split_taxa_pd


def get_taxa_edit(taxo):
    taxo_edit = []
    for tax in taxo:
        if str(tax) == 'nan':
            taxo_edit.append(tax)
        elif not tax.strip('_'):
            taxo_edit.append(tax)
        else:
            taxo_edit.append(tax.replace(',', '_'))
    return taxo_edit


def get_edit_taxonomy_command(data):
    cmd = ''
    out_pd = pd.read_csv(data.tax[2], dtype=str, sep='\t')
    taxo = out_pd['Taxon'].tolist()
    taxo_edit = get_taxa_edit(taxo)
    if taxo != taxo_edit:
        out_pd['Taxon'] = taxo_edit
        out_pd.to_csv(data.tax[2], index=False, sep='\t')
        cmd = run_import(
            data.tax[2], data.tax[1], 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_wol(force: bool, tsv_pd: pd.DataFrame, out_qza: str,
                     out_tsv: str, cur_datasets_features: dict) -> str:
    cmd = ''
    if force or not isfile(out_qza):
        g2lineage = parse_g2lineage()
        rev_cur_datasets_features = dict(
            (y, x) for x, y in cur_datasets_features.items())
        if not isfile(out_tsv):
            with open(out_tsv, 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in tsv_pd.index:
                    if rev_cur_datasets_features[feat] in g2lineage:
                        o.write('%s\t%s\n' % (
                            feat, g2lineage[rev_cur_datasets_features[feat]]))
                    else:
                        o.write('%s\t%s\n' % (feat, feat.replace('|', '; ')))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def run_taxonomy_amplicon(
        dat: str, i_datasets_folder: str, force: bool, tsv_pd: pd.DataFrame,
        out_qza: str, out_tsv: str, i_classifier: str) -> str:
    cmd = ''
    if isfile(out_tsv) and not isfile(out_qza):
        cmd += run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    else:
        ref_classifier_qza = get_taxonomy_classifier(i_classifier)
        odir_seqs = get_analysis_folder(i_datasets_folder, 'seqs/%s' % dat)
        out_fp_seqs_rad = '%s/seq_%s' % (odir_seqs, dat)
        out_fp_seqs_fasta = '%s.fasta' % out_fp_seqs_rad
        out_fp_seqs_qza = '%s.qza' % out_fp_seqs_rad
        if force or not isfile(out_fp_seqs_qza):
            cmd += write_seqs_fasta(out_fp_seqs_fasta, out_fp_seqs_qza, tsv_pd)
        if force or not isfile(out_qza):
            cmd += write_taxonomy_sklearn(out_qza, out_fp_seqs_qza,
                                          ref_classifier_qza)
            cmd += run_export(out_qza, out_tsv, '')
    return cmd


def run_taxonomy_others(force: bool, tsv_pd: pd.DataFrame,
                        out_qza: str, out_tsv: str) -> str:
    cmd = ''
    if force or not isfile(out_qza):
        if not isfile(out_tsv):
            with open(out_tsv, 'w') as o:
                o.write('Feature ID\tTaxon\n')
                for feat in tsv_pd.index:
                    o.write('%s\t%s\n' % (feat, feat))
        cmd = run_import(out_tsv, out_qza, 'FeatureData[Taxonomy]')
    return cmd


def get_taxonomy_command(dat, config, data):
    cmd = ''
    if data.tax[0] == 'wol':
        cmd = run_taxonomy_wol(
            config.force, data.data[0], data.tax[1],
            data.tax[2], data.features)
    elif data.tax[0] in ['amplicon', 'sklearn']:
        if config.i_classifier:
            cmd = run_taxonomy_amplicon(
                dat, config.i_datasets_folder, config.force, data.data[0],
                data.tax[1], data.tax[2], config.i_classifier)
        else:
            print('No classifier passed for 16S data\nExiting...')
    else:
        cmd = run_taxonomy_others(
            config.force, data.data[0], data.tax[1], data.tax[2])
    return cmd


