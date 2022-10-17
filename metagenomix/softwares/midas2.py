# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import isdir, isfile
from metagenomix.core.parameters import tech_params
from metagenomix._io_utils import (caller, status_update, io_update,
                                   to_do, tech_specificity)


def get_midas2_cmd(
        self,
        fastqs: list,
        focus_dir: str,
        db_name: str,
        db_path: str,
        species_list: str,
        params: dict,
        run_step: str
) -> str:
    """Collect the relevant MIDAS2 command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastqs : list
        Path to the input files
    focus_dir : str
        Path to the output folder
    db_name : str
        MIDAS Database name
    db_path : str
        Path to local MIDAS Database
    species_list : str
        path to file containing the list of species IDs to focus on
    params : dict
        Run parameters
    run_step : str
        Name of the MIDAS2 module

    Returns
    -------
    cmd : str
        MIDAS2 command
    """
    cmd = 'midas2 run_%s' % run_step
    cmd += ' --sample_name %s' % self.sam_pool
    if self.config.force:
        cmd += ' --force'
    cmd += ' -1 %s' % fastqs[0]
    if len(fastqs) == 2:
        cmd += ' -2 %s' % fastqs[1]

    if run_step in ['genes', 'snps']:
        for param in ['prebuilt_bowtie2_indexes', 'prebuilt_bowtie2_species']:
            if params[param]:
                cmd += ' --%s %s' % (param, params[param])

        if params['aln_interleaved']:
            cmd += ' --aln_interleaved'

        if species_list:
            cmd += ' --species_list %s' % params['species_list']

        for param in ['select_by', 'select_threshold', 'aln_speed',
                      'aln_mode', 'fragment_length', 'aln_readq']:
            cmd += ' --%s %s' % (param, params[param])

        if not params['aln_mapid']:
            cmd += ' --aln_mapid 0.94'
        else:
            cmd += ' --aln_mapid %s' % params['aln_mapid']

        if params['aln_mapq']:
            cmd += ' --aln_mapq %s' % params['aln_mapq']
        else:
            if run_step == 'genes':
                cmd += ' --aln_mapq 2'
            if run_step in 'snps':
                cmd += ' --aln_mapq 10'

        if params['chunk_size']:
            cmd += ' --chunk_size %s' % params['chunk_size']
        else:
            if run_step == 'genes':
                cmd += ' --chunk_size 50000'
            if run_step in 'snps':
                cmd += ' --chunk_size 1000000'

        if run_step == 'genes':
            cmd += ' --read_depth %s' % params['read_depth']
        else:
            for boolean in ['paired_only', 'ignore_ambiguous',
                            'advanced', 'analysis_ready']:
                if params[boolean]:
                    cmd += ' --%s' % boolean
            for param in ['aln_baseq', 'aln_trim', 'site_depth', 'snp_maf']:
                cmd += ' --%s %s' % (param, params[param])

    else:
        cmd += ' --word_size %s' % params['word_size']
        cmd += ' --marker_reads %s' % params['marker_reads']
        cmd += ' --marker_covered %s' % params['marker_covered']
        if params['aln_mapid']:
            cmd += ' --aln_mapid %s' % params['aln_mapid']

    for param in ['aln_cov', 'max_reads']:
        if params[param]:
            cmd += ' --%s %s' % (param, params[param])

    cmd += ' --midasdb_name %s' % db_name
    cmd += ' --midasdb_dir %s' % db_path
    cmd += ' --num_cores %s' % params['cpus']
    cmd += ' %s' % focus_dir

    return cmd


def get_midas2(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        params: dict,
        run_step: str
) -> None:
    """Get the MIDAS2 inputs and outputs for the current tech and
    sample and collect the commands and data structures.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for checkm tetra
        .sam_pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    run_step : str
        Name of the MIDAS2 step
    """
    species_lists = {'all_species': ''}
    if run_step == 'species':
        out = '%s/species/species_profile.tsv' % self.sam_pool
    else:
        if params['species_list']:
            species_lists = params['species_list']
        if run_step == 'genes':
            out = '%s/genes/genes_summary.tsv' % self.sam_pool
        else:
            out = '%s/snps/snps_summary.tsv' % self.sam_pool

    to_dos = status_update(self, tech, fastqs)

    db_folder = self.databases.paths['midas2']
    for db_name in self.soft.params['databases']:
        if self.config.dev:
            db_path = '%s/%s' % (db_folder, db_name)
        else:
            db_path = '%s/%s' % (db_folder, db_name)
            if not isdir(db_path):
                sys.exit('[midas2] Path to database "%s not found' % db_path)

        for focus, species_list in species_lists.items():
            focus_dir = out_dir + '/' + db_name + '/' + focus
            out_fp = focus_dir + '/' + out
            if self.config.force or to_do(out_fp):
                cmd = get_midas2_cmd(self, fastqs, focus_dir, db_name, db_path,
                                     species_list, params, run_step)
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                    i_f = fastqs
                    if isfile(species_list):
                        i_f += [species_list]
                    io_update(self, i_f=i_f, o_d=focus_dir, key=(tech,))


def midas2_snps(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        params: dict
) -> None:
    """Predict single-nucleotide-polymorphism.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder
        .sam_pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, fastqs, out_dir, params, 'snps')


def midas2_genes(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        params: dict
) -> None:
    """Metagenomic pan-genome profiling.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder
        .sam_pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, fastqs, out_dir, params, 'genes')


def midas2_species(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        params: dict
) -> None:
    """Estimate species abundance profile for given sample.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder
        .sam_pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, fastqs, out_dir, params, 'species')


def midas2_(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        params: dict
) -> None:
    """Perform all the key steps of the MIDAS2 analysis.

    Notes
    -----
    This is done if the pipeline specifies "cutadapt midas2" and not specific
    modules, such as "cutadapt midas2_species" or "cutadapt midas2_snps".

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for checkm tetra
        .sam_pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    for run_step in ['species', 'snps', 'genes']:
        get_midas2(self, tech, fastqs, out_dir, params, run_step)


def midas2(self) -> None:
    """Metagenomic Intra-Species Diversity Analysis System 2 (MIDAS2) is an
    integrated pipeline for profiling single nucleotide variants (SNVs) and
    gene copy number variants (CNVs) in shotgun metagenomic reads. MIDAS2
    implements the same analyses as the original MIDAS, but re-engineered to
    addresses the computational challenges presented by increasingly large
    reference genome databases.

    References
    ----------
    Zhao, C., Dimitrov, B., Goldman, M., Nayfach, S. and Pollard, K.S.,
    2022. MIDAS2: Metagenomic Intra-species Diversity Analysis System. bioRxiv.

    Notes
    -----
    GitHub  : https://github.com/czbiohub/MIDAS2
    Paper   : https://www.biorxiv.org/content/10.1101/2022.06.16.496510v2
    Docs    : https://midas2.readthedocs.io/en/latest/index.html

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    if set(self.sam_pool) == {''}:
        sys.exit('[midas2] Run on reads not MAGs (not "%s")' % self.soft.prev)
    elif self.sam_pool in self.pools:
        sys.exit('[midas2] Run on reads (not after "%s")' % self.soft.prev)

    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam):
            continue
        out_dir = self.dir + '/' + tech
        module_call = caller(self, __name__)
        params = tech_params(self, tech)
        module_call(self, tech, fastqs, out_dir, params)
