# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from os.path import isdir

from metagenomix._io_utils import io_update, to_do, status_update
from metagenomix.core.parameters import tech_params

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
    inputs : list
        Path(s) to the input fastq files
    out : str
        Path to the file to write
    """
    cmd = ''
    inputs = []
    out = '%s/%s/samples_files.txt' % (self.dir, tech)
    for sdx, sam in enumerate(self.inputs):
        fs = self.inputs[sam][(tech, sam)]
        inputs.extend(fs)
        if not fs:
            continue
        if sdx:
            cmd += 'echo -e "%s: %s\\n" >> %s\n' % (sam, '; '.join(fs), out)
        else:
            cmd += 'echo -e "%s: %s\\n" > %s\n' % (sam, '; '.join(fs), out)
    if cmd:
        cmd += 'envsubst < %s > %s.tmp\n' % (out, out)
        cmd += 'mv %s.tmp %s\n' % (out, out)
    return cmd, inputs, out


def simka_cmd(
        self,
        params: dict,
        sim_in: str,
        out_dir: str,
        k: int,
        n: int
) -> str:
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


def simka_min_cmd(
        params: dict,
        sim_in: str,
        out_dir: str,
        k: int,
        n: str
) -> str:
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
    cmd += ' -min-read-size %s' % params['min_read_size']
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
        params: dict,
        sim_in: str,
        out_dir: str,
        k: int,
        n: str
) -> str:
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


def simka_pcoa_cmd(
        self,
        mat: str
) -> str:
    """Write the Simka command for the pcoa/tsne based on the distance matrices.

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    mat : str
        Distance matrix

    Returns
    -------
    cmd : str
        Simka command line
    """
    sym_cmd = ''
    mat_o = '%s_sym.tsv' % mat.split('.csv')[0]
    if self.config.force or to_do(mat_o):
        if mat.endswith('gz'):
            mat_fp = mat.replace('.csv.gz', '.csv')
            if to_do(mat_fp):
                sym_cmd += 'gunzip %s\n' % mat
        else:
            mat_fp = mat
        sym_cmd += '\n%s/symmetrize_simka_matrix.py' % RESOURCES
        sym_cmd += ' -i %s\n' % mat_fp

    imp_cmd = ''
    mat_dm = mat_o.replace('.tsv', '_dm.qza')
    if self.config.force or to_do(mat_dm):
        imp_cmd += '\nqiime softwares import'
        imp_cmd += ' --input-path %s' % mat_o
        imp_cmd += ' --output-path %s' % mat_dm
        imp_cmd += ' --type DistanceMatrix\n'

    ordi_cmd = ''
    exp_cmd = ''
    emp_cmd = ''
    for ordi in ['pcoa', 'tsne']:
        ordi_fp = mat_dm.replace('.qza', '_%s.qza' % ordi)
        if self.config.force or to_do(ordi_fp):
            ordi_cmd += '\nqiime diversity %s' % ordi
            ordi_cmd += ' --i-distance-matrix %s' % mat_dm
            ordi_cmd += ' --o-%s %s' % (ordi, ordi_fp)
            if ordi == 'tsne':
                ordi_cmd += ' --p-perplexity 25'
                ordi_cmd += ' --p-early-exaggeration 10'
                ordi_cmd += ' --p-learning-rate 200\n'

        ordi_dir = ordi_fp.replace('.qza', '')
        ordi_txt = ordi_fp.replace('.qza', '.txt')
        if self.config.force or to_do(ordi_txt):
            exp_cmd += '\nqiime softwares export'
            exp_cmd += ' --input-path %s' % ordi_fp
            exp_cmd += ' --output-path %s\n' % ordi_dir
            exp_cmd += 'mv %s/ordination.txt %s\n' % (ordi_dir, ordi_txt)
            exp_cmd += 'rm -rf %s\n' % ordi_dir

        emp_fp = ordi_fp.replace('.qza', '_emp.qzv')
        if self.config.force or to_do(emp_fp):
            emp_cmd += '\nqiime emperor plot'
            emp_cmd += ' --i-pcoa %s' % ordi_fp
            emp_cmd += ' --m-metadata-file %s' % self.config.meta_fp
            emp_cmd += ' --o-visualization %s\n' % emp_fp

    cmd = ''
    if emp_cmd:
        cmd = sym_cmd + imp_cmd + ordi_cmd + exp_cmd + emp_cmd
    return cmd


def simka(self) -> None:
    """Simka is a de novo comparative metagenomics tool. Simka represents
    each dataset as a k-mer spectrum and compute several classical ecological
    distances between them.

    References
    ----------
    Benoit, Gaëtan, et al. "Multiple comparative metagenomics using multiset
    k-mer counting." PeerJ Computer Science 2 (2016): e94.

    Notes
    -----
    GitHub  : https://github.com/GATB/simka
    Docs    : https://gatb.inria.fr/software/simka
    Paper   : https://doi.org/10.7717/peerj-cs.94

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
    """
    for tech in self.config.techs:
        params = tech_params(self, tech)
        input_cmd, input_fastqs, input_file = get_simka_input(self, tech)
        status_update(self, tech, input_fastqs)
        cmds = ''
        for k in map(int, params['kmer']):
            for n in map(int, params['log_reads']):
                out_dir = '%s/%s/k%s/n%s' % (self.dir, tech, k, n)
                cmd = simka_cmd(self, params, input_file, out_dir, k, n)
                if cmd:
                    cmds += cmd
                    self.outputs['dirs'].append(out_dir)
                    io_update(self, i_f=input_fastqs, o_d=out_dir, key=tech)
                else:
                    io_update(self, i_d=out_dir, key=tech)

                for met in ['abundance_braycurtis', 'presenceAbsence_jaccard']:
                    mat = '%s/mat_%s.csv.gz' % (out_dir, met)
                    cmd = simka_pcoa_cmd(self, mat)
                    if cmd:
                        cmds += cmd
                        self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                        io_update(self, o_d=out_dir, key=tech)
        if cmds:
            cmd = input_cmd + cmd
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)
            self.soft.add_status(tech, 'all samples', 1)
        else:
            self.soft.add_status(tech, 'all samples', 0)
