# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from os.path import isdir, isfile
from metagenomix._io_utils import mkdr

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def get_simka_input(dir_path, inputs) -> str:
    """Write the file containing the paths to each sample.

    Parameters
    ----------
    dir_path : str
        Path to the output directory for the simka analyses.
    inputs : dict
        Lists of fastq files per sample name.

    Returns
    -------
    sim_out : str
        Path to the file to write.
    """
    sim_out = '%s/samples_files.txt' % dir_path
    mkdr(sim_out, True)
    with open(sim_out, 'w') as o:
        for sam in inputs:
            fs = inputs[sam]
            if 'after_None' not in dir_path:
                fs = ['%sq' % x[:-1] if 'notCombined_' in x else x for x in fs]
            o.write('%s: %s\n' % (sam, '; '.join(fs)))
    return sim_out


def simka_cmd(soft, smin: bool, sim_in: str, out_dir: str,
              k: str, n: int, tmp_dir: str) -> str:
    """

    Parameters
    ----------
    soft : pipeline.Soft
        Class instance for the current software.
    smin : bool
        Whether to use SimkaMin (True) or base Simka (False).
    sim_in : str
        Input file containing to the fastq file paths for Simka.
    out_dir : str
        Output folder for Simka.
    k : str
        Length of the k-mer.
    n : int
        Number of sequences to inject in the Simka analysis.
    tmp_dir : str
        Temporary, scratch directory.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = ''
    if isdir('%s/simkamin' % out_dir):
        cmd = 'rm -rf %s/simkamin\n' % out_dir
    if smin:
        if isfile('%s/mat_abundance_braycurtis.csv' % out_dir):
            return ''
        elif isfile('%s/mat_abundance_braycurtis.csv.gz' % out_dir):
            return ''
        cmd += simka_min_cmd(soft, sim_in, out_dir, k, str(n))
    else:
        cmd += simka_base_cmd(soft, sim_in, out_dir, k, str(n), tmp_dir)
    return cmd


def simka_min_cmd(soft, sim_in: str, out_dir: str, k: str, n: str) -> str:
    """Write the Simka command for the SimkaMin algorithm.

    Parameters
    ----------
    soft : pipeline.Soft
        Class instance for the current software.
    sim_in : str
        Input file containing to the fastq file paths for Simka.
    out_dir : str
        Output folder for Simka.
    k : str
        Length of the k-mer.
    n : str
        Number of sequences to inject in the Simka analysis.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = 'python %s/simkaMin/simkaMin.py' % soft.params['path']
    cmd += ' -bin %s/bin/simkaMinCore' % soft.params['path']
    cmd += ' -in %s -out %s' % (sim_in, out_dir)
    cmd += ' -kmer-size %s' % k
    cmd += ' -max-reads %s' % n
    cmd += ' -nb-kmers 50000'
    cmd += ' -max-memory 180000'
    cmd += ' -min-read-size 100'
    cmd += ' -filter'
    cmd += ' -nb-cores %s\n' % soft.params['cpus']
    cmd += 'rm -rf %s/simkamin' % out_dir
    return cmd


def simka_base_cmd(
        soft, sim_in: str, out_dir: str, k: str, n: str, tmpdir: str) -> str:
    """Write the Simka command for the Simka base algorithm.

    Parameters
    ----------
    soft : pipeline.Soft
        Class instance of for the current software.
    sim_in : str
        Input file containing to the fastq file paths for Simka.
    out_dir : str
        Output folder for Simka.
    k : str
        Length of the k-mer.
    n : str
        Number of sequences to inject in the Simka analysis.
    tmpdir : str
        Temporary, scratch directory.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = '%s/build/bin/simka' % soft.path
    cmd += ' -in %s' % sim_in
    cmd += ' -out %s' % out_dir
    cmd += ' -out-tmp %s' % tmpdir
    cmd += ' -abundance-min 5'
    cmd += ' -kmer-size %s' % k
    cmd += ' -max-reads %s' % n
    cmd += ' -data-info'
    cmd += ' -simple-dist'
    cmd += ' -nb-cores %s' % soft.params['cpus']
    cmd += ' -max-count %s' % soft.params['cpus']
    cmd += ' -max-memory %s' % str((int(soft.params['mem_num'][0])*1000)-1000)
    return cmd


def simka_pcoa_cmd(mat: str, meta_fpo: str) -> str:
    """Write the Simka command for the pcoa based on the distance matrices.

    Parameters
    ----------
    mat : str
        Path to the input matrix.
    meta_fpo: str
        Path to the output matrix.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = ''
    if 'sym' in mat:
        return cmd
    if mat.endswith('gz'):
        mat_fp = mat.replace('.csv.gz', '.csv')
        if not isfile(mat_fp):
            cmd += 'gunzip %s\n' % mat
    else:
        mat_fp = mat

    mat_o = mat_fp.replace('.csv', '_sym.tsv')
    if not isfile(mat_o):
        cmd += '%s/symmetrize_simka_matrix.py' % RESOURCES
        cmd += ' -i %s\n' % mat_fp

    mat_dm = mat_o.replace('.tsv', '_dm.qza')
    if not isfile(mat_dm):
        cmd += 'qiime tools import'
        cmd += ' --input-path %s' % mat_o
        cmd += ' --output-path %s' % mat_dm
        cmd += ' --type DistanceMatrix\n'

    mat_pcoa = mat_dm.replace('.qza', '_pcoa.qza')
    if not isfile(mat_pcoa):
        cmd += 'qiime diversity pcoa'
        cmd += ' --i-distance-matrix %s' % mat_dm
        cmd += ' --o-pcoa %s\n' % mat_pcoa

    mat_pcoa_dir = mat_pcoa.replace('.qza', '')
    mat_pcoa_txt = mat_pcoa.replace('.qza', '.txt')
    if not isfile(mat_pcoa_txt):
        cmd += 'qiime tools export'
        cmd += ' --input-path %s' % mat_pcoa
        cmd += ' --output-path %s\n' % mat_pcoa_dir
        cmd += 'mv %s/ordination.txt %s\n' % (mat_pcoa_dir, mat_pcoa_txt)
        cmd += 'rm -rf %s\n' % mat_pcoa_dir

    mat_emp = mat_pcoa.replace('.qza', '_emp.qzv')
    if not isfile(mat_emp):
        cmd += 'qiime emperor plot'
        cmd += ' --i-pcoa %s' % mat_pcoa
        cmd += ' --m-metadata-file %s' % meta_fpo
        cmd += ' --o-visualization %s' % mat_emp
    return cmd
