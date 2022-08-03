# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import pkg_resources
from os.path import isdir

from metagenomix._io_utils import io_update, to_do
from metagenomix.parameters import tech_params

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def get_simka_input(
        self,
        tech: str
) -> tuple:
    """Write the file containing the paths to each sample.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
        .inputs : dict
            Input files
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    -------
    cmd : str
        Command to create the inputs file
    out : str
        Path to the file to write
    """
    out = '%s/%s/samples_files.txt' % (self.dir, tech)
    cmd = ''
    for sdx, sam in enumerate(self.inputs):
        fs = self.inputs[sam][(tech, sam)]
        io_update(self, i_f=fs, key=tech)
        if not fs:
            continue
        if sdx:
            cmd += 'echo -e "%s: %s\\n" >> %s\n' % (sam, '; '.join(fs), out)
        else:
            cmd += 'echo -e "%s: %s\\n" > %s\n' % (sam, '; '.join(fs), out)
    if cmd:
        cmd += 'envsubst < %s > %s.tmp\n' % (out, out)
        cmd += 'mv %s.tmp %s\n' % (out, out)
    return cmd, out


def simka_cmd(self, params: dict, sim_in: str, out_dir: str,
              k: int, n: int) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
        .config
            Configurations
    params : dict
        Parameters for the current technology
    sim_in : str
        Input file containing to the fastq file paths for Simka
    out_dir : str
        Output folder for Simka
    k : int
        Length of the k-mer
    n : int
        Number of sequences to inject in the Simka analysis

    Returns
    -------
    cmd : str
        Simka command line
    """
    cmd = ''
    if isdir('%s/simkamin' % out_dir):
        cmd = 'rm -rf %s/simkamin\n' % out_dir
    if params['simkaMin']:
        if not self.config.force:
            if not to_do('%s/mat_abundance_braycurtis.csv' % out_dir):
                return ''
            elif not to_do('%s/mat_abundance_braycurtis.csv.gz' % out_dir):
                return ''
        cmd += simka_min_cmd(params, sim_in, out_dir, k, str(n))
    else:
        cmd += simka_base_cmd(params, sim_in, out_dir, k, str(n))
    return cmd


def simka_min_cmd(params: dict, sim_in: str, out_dir: str,
                  k: int, n: str) -> str:
    """Write the Simka command for the SimkaMin algorithm.

    Parameters
    ----------
    params : dict
        Simka parameters
    sim_in : str
        Input file containing to the fastq file paths for Simka.
    out_dir : str
        Output folder for Simka.
    k : int
        Length of the k-mer.
    n : str
        Number of sequences to inject in the Simka analysis.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = 'python %s/simkaMin/simkaMin.py' % params['path']
    cmd += ' -bin %s/bin/simkaMinCore' % params['path']
    cmd += ' -in %s' % sim_in
    cmd += ' -out %s' % out_dir
    cmd += ' -max-reads %s' % n
    cmd += ' -kmer-size %s' % str(k)
    cmd += ' -nb-kmers %s' % params['nb_kmers']
    cmd += ' -min-read-size %s\n' % params['min_read_size']
    if int(params['mem']) == 1:
        mem = '1'
    else:
        mem = str((params['mem'] - 1))
    if params['mem_dim'].lower()[0] == 'g':
        mem += '000'
    cmd += ' -filter'
    cmd += ' -max-memory %s' % mem
    cmd += ' -nb-cores %s\n' % params['cpus']
    cmd += 'rm -rf %s/simkamin\n' % out_dir
    return cmd


def simka_base_cmd(
        params: dict, sim_in: str, out_dir: str, k: int, n: str) -> str:
    """Write the Simka command for the Simka base algorithm.

    Parameters
    ----------
    params : dict
        Simka parameters
    sim_in : str
        Input file containing to the fastq file paths for Simka.
    out_dir : str
        Output folder for Simka.
    k : int
        Length of the k-mer.
    n : str
        Number of sequences to inject in the Simka analysis.

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = '%s/bin/simka' % params['path']
    cmd += ' -in %s.tmp' % sim_in
    cmd += ' -out %s' % out_dir
    cmd += ' -out-tmp %s_tmp' % out_dir
    cmd += ' -abundance-min 5'
    cmd += ' -kmer-size %s' % int(k)
    cmd += ' -max-reads %s' % n
    cmd += ' -data-info'
    cmd += ' -simple-dist'
    cmd += ' -nb-cores %s' % params['cpus']
    cmd += ' -max-count %s' % params['cpus']
    cmd += ' -max-memory %s\n' % str((int(params['mem'][0])*1000)-1000)
    cmd += 'rm -rf %s_tmp\n' % out_dir
    cmd += 'rm %s.tmp\n' % sim_in
    return cmd


def simka_pcoa_cmd(mat: str, config) -> str:
    """Write the Simka command for the pcoa based on the distance matrices.

    Parameters
    ----------
    mat : str
        Path to the input matrix.
    config
        Configurations

    Returns
    -------
    cmd : str
        Simka command line.
    """
    cmd = ''
    if not config.force and 'sym' in mat:
        return cmd
    if mat.endswith('gz'):
        mat_fp = mat.replace('.csv.gz', '.csv')
        if to_do(mat_fp):
            cmd += 'gunzip %s\n' % mat
    else:
        mat_fp = mat

    mat_o = mat_fp.replace('.csv', '_sym.tsv')
    if config.force or to_do(mat_o):
        cmd += 'python %s/symmetrize_simka_matrix.py' % RESOURCES
        cmd += ' -i %s\n' % mat_fp

    mat_dm = mat_o.replace('.tsv', '_dm.qza')
    if config.force or to_do(mat_dm):
        cmd += 'qiime tools import'
        cmd += ' --input-path %s' % mat_o
        cmd += ' --output-path %s' % mat_dm
        cmd += ' --type DistanceMatrix\n'

    mat_pcoa = mat_dm.replace('.qza', '_pcoa.qza')
    if config.force or to_do(mat_pcoa):
        cmd += 'qiime diversity pcoa'
        cmd += ' --i-distance-matrix %s' % mat_dm
        cmd += ' --o-pcoa %s\n' % mat_pcoa

    mat_pcoa_dir = mat_pcoa.replace('.qza', '')
    mat_pcoa_txt = mat_pcoa.replace('.qza', '.txt')
    if config.force or to_do(mat_pcoa_txt):
        cmd += 'qiime tools export'
        cmd += ' --input-path %s' % mat_pcoa
        cmd += ' --output-path %s\n' % mat_pcoa_dir
        cmd += 'mv %s/ordination.txt %s\n' % (mat_pcoa_dir, mat_pcoa_txt)
        cmd += 'rm -rf %s\n' % mat_pcoa_dir

    mat_emp = mat_pcoa.replace('.qza', '_emp.qzv')
    if config.force or to_do(mat_emp):
        cmd += 'qiime emperor plot'
        cmd += ' --i-pcoa %s' % mat_pcoa
        cmd += ' --m-metadata-file %s' % config.meta_fp
        cmd += ' --o-visualization %s' % mat_emp
    return cmd


def simka(self) -> None:
    """Create command lines for Simka.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    for tech in self.config.techs:
        params = tech_params(self, tech)
        input_cmd, input_file = get_simka_input(self, tech)
        for k in map(int, params['kmer']):
            for n in map(int, params['log_reads']):
                out_d = '%s/%s/k%s/n%s' % (self.dir, tech, k, n)
                cmd = simka_cmd(self, params, input_file, out_d, k, n)
                if cmd:
                    cmd = input_cmd + cmd
                    self.outputs['dirs'].append(out_d)
                    self.outputs['cmds'].setdefault(tech, []).append(cmd)
                    io_update(self, o_d=out_d, key=tech)
                for mdx, mat in enumerate(glob.glob('%s/mat_*.csv*' % out_d)):
                    cmd = simka_pcoa_cmd(mat, self.config)
                    if cmd:
                        self.outputs['cmds'].setdefault(tech, []).append(cmd)
                        io_update(self, o_d=out_d, key=tech)
