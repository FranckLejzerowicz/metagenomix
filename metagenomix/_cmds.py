# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import sys
import pandas as pd
from os.path import isfile, splitext


def run_import(input_path: str, output_path: str, typ: str) -> str:
    """
    Return the import qiime2 command.

    :param input_path: input file path.
    :param output_path: output file path.
    :param typ: qiime2 type.
    :return: command to qiime2.
    """
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not input_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(input_path)[0]
            cmd += 'biom convert \\\n'
            cmd += '  -i %s \\\n' % input_path
            cmd += '  -o %s \\\n' % cur_biom
            cmd += '  --table-type="OTU table" \\\n'
            cmd += '  --to-hdf5\n\n'
            cmd += 'qiime tools import \\\n'
            cmd += '  --input-path %s \\\n' % cur_biom
            cmd += '  --output-path %s \\\n' % output_path
            cmd += '  --type "FeatureTable[Frequency]"\n'
        else:
            cmd += 'qiime tools import \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s \\\n' % output_path
            cmd += '  --type "FeatureTable[Frequency]"\n'
    else:
        cmd += 'qiime tools import \\\n'
        cmd += '  --input-path %s \\\n' % input_path
        cmd += '  --output-path %s \\\n' % output_path
        cmd += '  --type "%s"\n' % typ
    return cmd


def run_export(input_path: str, output_path: str, typ: str) -> str:
    cmd = ''
    if typ.startswith("FeatureTable"):
        if not output_path.endswith('biom'):
            cur_biom = '%s.biom' % splitext(output_path)[0]
            cmd += 'qiime tools export \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s\n' % splitext(output_path)[0]
            cmd += 'mv %s/*.biom %s\n' % (splitext(output_path)[0], cur_biom)
            cmd += 'biom convert'
            cmd += '  -i %s \\\n' % cur_biom
            cmd += '  -o %s.tmp \\\n' % output_path
            cmd += '  --to-tsv\n\n'
            cmd += 'tail -n +2 %s.tmp > %s\n\n' % (output_path, output_path)
            cmd += 'rm -rf %s %s.tmp\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'qiime tools export \\\n'
            cmd += '  --input-path %s \\\n' % input_path
            cmd += '  --output-path %s\n' % splitext(output_path)[0]
            cmd += 'mv %s/*.biom %s\n' % (splitext(input_path)[0], output_path)
            cmd += 'rm -rf %s\n' % splitext(input_path)[0]
    else:
        cmd += 'qiime tools export \\\n'
        cmd += '  --input-path %s \\\n' % input_path
        cmd += '  --output-path %s\n' % splitext(output_path)[0]
        if 'songbird' in typ:
            cmd += 'mv %s/index.html %s\n' % (splitext(output_path)[0], output_path)
        else:
            cmd += 'mv %s/*.tsv %s\n' % (splitext(output_path)[0], output_path)
        cmd += 'rm -rf %s\n' % splitext(output_path)[0]
    return cmd


def write_seqs_fasta(out_fp_seqs_fasta: str, out_fp_seqs_qza: str,
                     tsv_pd: pd.DataFrame, tsv_fp: str = '') -> str:
    """
    Write the fasta sequences.

    :param out_fp_seqs_fasta: output sequences fasta file name.
    :param out_fp_seqs_qza: output sequences qiime2 Artefact file name.
    :param tsv_pd: table which feature names are sequences.
    :param cur_sh: writing file handle.
    """
    with open(out_fp_seqs_fasta, 'w') as fas_o:
        for seq in tsv_pd.index:
            fas_o.write('>%s\n%s\n' % (seq.strip(), seq.strip()))
    cmd = '# Write features as fasta file:\n'
    cmd += '#  - Features from: %s\n' % tsv_fp
    cmd += '# Snippet:\n'
    cmd += '# ```:\n'
    cmd += "# with open(fasta_out, 'w') as o:\n"
    cmd += "#     for seq in tsv_pd.index:\n"
    cmd += "#         o.write('>%s\\n%s\\n' % (seq.strip(), seq.strip()))\n"
    cmd += '# ```:\n'
    cmd += run_import(
        out_fp_seqs_fasta, out_fp_seqs_qza, 'FeatureData[Sequence]')
    return cmd


def write_taxonomy_sklearn(out_qza: str, out_fp_seqs_qza: str,
                           ref_classifier_qza: str) -> str:
    """
    Classify reads by taxon using a fitted classifier.
    https://docs.qiime2.org/2020.2/plugins/available/feature-classifier/classify-sklearn

    :param out_qza: Taxonomic classifications.
    :param out_fp_seqs_qza: The features sequences that should be classified.
    :param ref_classifier_qza: Taxonomic classifier.
    """
    cmd = 'qiime feature-classifier classify-sklearn \\\n'
    cmd += '--i-reads %s \\\n' % out_fp_seqs_qza
    cmd += '--i-classifier %s \\\n' % ref_classifier_qza
    cmd += '--p-n-jobs %s \\\n' % '4'
    cmd += '--o-classification %s\n' % out_qza
    return cmd


def get_case(
        case_vals: list,
        case_var: str,
        form: str = None) -> str:
    if len(case_vals):
        case = '%s_%s' % (case_var, '-'.join(
            [x.replace('<', 'below').replace('>', 'above') for x in case_vals]))
    else:
        case = case_var
    if form:
        case = '%s_%s' % (case, form)
    case = re.sub('[ /()]', '', case.replace('__', '_'))
    return case


def get_new_meta_pd(meta_pd: pd.DataFrame, case: str,
                    case_var: str, case_vals: list) -> pd.DataFrame:
    if 'ALL' in case:
        new_meta_pd = meta_pd.copy()
    elif len([x for x in case_vals if x[0] == '>' or x[0] == '<']):
        new_meta_pd = meta_pd.copy()
        for case_val in case_vals:
            if case_val[0] == '>':
                new_meta_pd = new_meta_pd[
                    new_meta_pd[case_var].astype(float) >= float(case_val[1:])
                ].copy()
            elif case_val[0] == '<':
                new_meta_pd = new_meta_pd[
                    new_meta_pd[case_var].astype(float) <= float(case_val[1:])
                ].copy()
    else:
        new_meta_pd = meta_pd[meta_pd[case_var].isin(case_vals)].copy()
    return new_meta_pd


def filter_feature_table(qza: str, new_qza: str, meta: str) -> str:
    cmd = '\nqiime feature-table filter-samples \\\n'
    cmd += '--i-table %s \\\n' % qza
    cmd += '--m-metadata-file %s \\\n' % meta
    cmd += '--o-filtered-table %s\n' % new_qza
    return cmd


def songbird_cmd(
        qza, new_qza, new_meta, nsams, params,
        formula, bformula, out_paths, force) -> tuple:

    fcmd = ''
    if force or not isfile(new_qza):
        fcmd += filter_feature_table(qza, new_qza, new_meta)

    batches = int(params['batches'])

    cmd = ''
    if force or not isfile(out_paths['diff_qza']):
        cmd = '# model\n'
        cmd += '\nqiime songbird multinomial \\\n'
        cmd += ' --i-table %s \\\n' % new_qza
        cmd += ' --m-metadata-file %s \\\n' % new_meta
        cmd += ' --p-formula "%s" \\\n' % formula
        cmd += ' --p-epochs %s \\\n' % params['epochs']
        if batches > 0.8 * nsams:
            cmd += ' --p-batch-size %s \\\n' % str(int(nsams * 0.8))
        else:
            cmd += ' --p-batch-size %s \\\n' % params['batches']
        cmd += ' --p-differential-prior %s \\\n' % params['diff_priors']
        cmd += ' --p-learning-rate %s \\\n' % params['learns']
        cmd += ' --p-min-sample-count %s \\\n' % params['thresh_samples']
        cmd += ' --p-min-feature-count %s \\\n' % params['thresh_feats']
        if 'examples' in params:
            cmd += ' --p-num-random-test-examples %s \\\n' % params['examples']
        else:
            cmd += ' --p-training-column %s \\\n' % params['train']
        cmd += ' --p-summary-interval %s \\\n' % params['summary_interval']
        cmd += ' --o-differentials %s \\\n' % out_paths['diff_qza']
        cmd += ' --o-regression-stats %s \\\n' % out_paths['stat']
        cmd += ' --o-regression-biplot %s\n' % out_paths['plot']

    if force or not isfile(out_paths['diff']):
        cmd += run_export(out_paths['diff_qza'], out_paths['diff'], '')

    stat_tsv = '%s.txt' % splitext(out_paths['stat'])[0]
    if force or not isfile(stat_tsv):
        cmd += run_export(out_paths['stat'], stat_tsv, '')

    bcmd = ''
    if force or len(out_paths['bdiff_qza']) and not isfile(out_paths['bstat']):
        bcmd += '\nqiime songbird multinomial \\\n'
        bcmd += ' --i-table %s \\\n' % new_qza
        bcmd += ' --m-metadata-file %s \\\n' % new_meta
        bcmd += ' --p-formula "%s" \\\n' % bformula
        bcmd += ' --p-epochs %s \\\n' % params['epochs']
        if batches > 0.8 * nsams:
            bcmd += ' --p-batch-size %s \\\n' % str(int(nsams * 0.8))
        else:
            bcmd += ' --p-batch-size %s \\\n' % params['batches']
        bcmd += ' --p-differential-prior %s \\\n' % params['diff_priors']
        bcmd += ' --p-learning-rate %s \\\n' % params['learns']
        bcmd += ' --p-min-sample-count %s \\\n' % params['thresh_samples']
        bcmd += ' --p-min-feature-count %s \\\n' % params['thresh_feats']
        if 'examples' in params:
            bcmd += ' --p-num-random-test-examples %s \\\n' % params['examples']
        else:
            bcmd += ' --p-training-column %s \\\n' % params['train']
        bcmd += ' --p-summary-interval %s \\\n' % params['summary_interval']
        bcmd += ' --o-differentials %s \\\n' % out_paths['bdiff_qza']
        bcmd += ' --o-regression-stats %s \\\n' % out_paths['bstat']
        bcmd += ' --o-regression-biplot %s\n' % out_paths['bplot']

    if force or not isfile(out_paths['tens']):
        cmd += '\nqiime songbird summarize-paired \\\n'
        cmd += ' --i-regression-stats %s \\\n' % out_paths['stat']
        cmd += ' --i-baseline-stats %s \\\n' % out_paths['bstat']
        cmd += ' --o-visualization %s\n' % out_paths['tens']

    if force or not isfile(out_paths['html']):
        cmd += run_export(out_paths['tens'], out_paths['html'], 'songbird')

    return cmd, fcmd, bcmd


def caller(self, namespace):
    """Calls as function the part of the software name that follows the first
    underscore.
    For example, software name "search_diamond" would call function`diamond()`.
    """
    func = self.soft.name.rsplit('_', 1)[1]
    module = sys.modules[namespace]
    if hasattr(module, func) and callable(getattr(module, func)):
        module_call = getattr(module, func)
        module_call(self)
