# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import dirname, isdir
from metagenomix.core.parameters import tech_params
from metagenomix._io_utils import (caller, status_update, io_update,
                                   to_do, tech_specificity)


def get_merge_cmd(
        self,
        tech: str,
        step: str,
        db: str,
        spc_list: str,
        out_dir: str,
        sam_fp: str
) -> str:
    """Collect the command line for merging MIDAS2 analysis results per sample.

    Parameters
    ----------
    self : Commands class instance
    tech : str
        Technology
    step : str
        Name of the MIDAS2 module
    db : str
        MIDAS Database name
    spc_list : str
        path to file containing the list of species IDs to focus on
    out_dir : str
        Path to the output folder
    sam_fp : dtr
        Path to samples input folders per tech, database and species focus

    Returns
    -------
    cmd : str
        Command line for merging the results of a midas run_* step
    """
    params = tech_params(self, tech)
    cmd = step.replace('_', ' merge_')
    if step == 'midas2_species':
        cmd += ' --min_cov %s' % params['min_cov']
    else:
        cmd += ' --num_cores %s' % params['cpus']
        if 'uhgg' in db:
            cmd += ' --midasdb_name uhgg'
        elif 'gtdb' in db:
            cmd += ' --midasdb_name gtdb'
        elif 'newdb' in db:
            cmd += ' --midasdb_name newdb'
        elif 's3db' in db:
            cmd += ' --midasdb_name s3db'
        cmd += ' --midasdb_dir %s' % self.databases.paths[db]

        if spc_list:
            cmd += ' --species_list %s' % spc_list

        if params['genome_depth']:
            cmd += ' --genome_depth %s' % params['genome_depth']
        elif step == 'midas2_genes':
            cmd += ' --genome_depth 1.0'
        else:
            cmd += ' --genome_depth 5.0'

        if params['sample_counts']:
            cmd += ' --sample_counts %s' % int(params['sample_counts'])
        elif step == 'midas2_genes':
            cmd += ' --sample_counts 1'
        else:
            cmd += ' --sample_counts 2'

        if step == 'midas2_genes':
            cmd += ' --min_copy %s' % params['min_copy']
            cmd += ' --cluster_pid %s' % params['cluster_pid']
        else:
            for param in [
                'site_depth', 'site_ratio', 'site_prev', 'snv_type',
                'snp_pooled_method', 'snp_maf', 'snp_type', 'locus_type'
            ]:
                if param == 'snp_type':
                    cmd += ' --snp_type {%s}' % ','.join(params[param])
                else:
                    cmd += ' --%s %s' % (param, params[param])
            for boolean in ['chunk_size', 'advanced', 'robust_chunk']:
                if params[boolean]:
                    cmd += ' --%s' % boolean
    cmd += ' --force'
    cmd += ' --samples_list %s' % sam_fp
    cmd += ' %s\n' % out_dir
    return cmd


def get_sams_list(
        self,
        tech: str,
        step: str,
        sam_fp: str,
        sams_folders: list,
        db: str,
        focus: str
) -> tuple:
    """Collect the command to make a samples list input file.

    Parameters
    ----------
    self : Commands class instance
    tech : str
        Technology
    step : str
        Name of the MIDAS2 module
    sam_fp : str
        Path to the file to contain the samples and their input folders
    sams_folders : list
        Samples and their input folders
    db : str
        MIDAS Database name
    focus : str
        Name of the species IDs list to focus on

    Returns
    -------
    cmd : str
        Command to create the samples list input file to MIDAS2 merge commands
    i_d : list
        Paths to the per-sample outputs of the midas2 run_* commands
    to_dos : list
        Paths to input file that need to be created (or whether they are links)
    """
    cmd = ''
    i_d, to_dos = [], []
    for sdx, (sam, folder) in enumerate(sams_folders):

        folder_sam = '%s/%s' % (folder, sam)
        i_d.append(folder_sam)

        if step == 'midas2_species':
            out = '%s/species/species_profile.tsv' % folder_sam
        elif step == 'midas2_genes':
            out = '%s/genes/genes_summary.tsv' % folder_sam
        else:
            out = '%s/snps/snps_summary.tsv' % folder_sam
        to_dos.extend(status_update(self, tech, [out], group=db, genome=focus))

        if not sdx:
            cmd += 'echo -e "sample_name\\tmidas_outdir" > %s\n' % sam_fp
        cmd += 'echo -e "%s\\t%s" >> %s\n' % (sam, folder, sam_fp)
    return cmd, i_d, to_dos


def get_samples(
        self,
        sams_inputs: dict,
        tech_step_samples: dict
) -> None:
    """Make a dict data structure per merging, with samples to merge as values.

    Parameters
    ----------
    self : Commands class instance
    sams_inputs : dict
        Output data structure from the MIDAS analysis steps (midas2 run_*)
    tech_step_samples : dict
    """
    for sam, step_inputs in sams_inputs.items():
        for step in ['midas2_species', 'midas2_genes', 'midas2_snps']:
            if step not in step_inputs:
                continue
            self.outputs['outs'][step] = {}
            inputs = step_inputs[step]
            for (tech, db, focus, spc_list), folder in inputs.items():
                if (tech, step) not in tech_step_samples:
                    tech_step_samples[(tech, step)] = {}
                tech_step_samples[(tech, step)].setdefault(
                    (db, focus, spc_list), []).append([sam, folder])


def get_merge_outs(
        out_dir: str,
        step: str
) -> str:
    """Get the name of the expected output files for MIDAS2 merge commands.

    Parameters
    ----------
    out_dir : str
        Path to the output folder
    step : str
        Name of the MIDAS2 step

    Returns
    -------
    out : str
        Path to the current output
    """
    if step == 'midas2_species':
        out = '%s/species/species_prevalence.tsv' % out_dir
    elif step == 'midas2_genes':
        out = '%s/genes/genes_summary.tsv' % out_dir
    else:
        out = '%s/snps/snps_summary.tsv' % out_dir
    return out


def merge(self):
    """Collect commands to perform all necessary MIDAS2 merging, including:
        - midas2 merge_species
            merge MIDAS species abundance results across metagenomic samples
        - midas2 merge_genes
            metagenomic pan-genome profiling
        - midas2 merge_snps
            pooled-samples SNPs calling

    Parameters
    ----------
    self : Commands class instance
    """
    tech_step_samples = {}
    for step, inputs in self.softs.items():
        if step.startswith('midas2') and step != self.soft.name:
            get_samples(self, inputs.outputs, tech_step_samples)

    for (tech, step), samples in tech_step_samples.items():
        for (db, focus, spc_list), sams_dirs in samples.items():
            out_dir = '/'.join([self.dir, tech, db, focus])
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'][step][(tech, db, focus, spc_list)] = out_dir
            sam_fp = '%s/sample_list.txt' % out_dir
            cmd_sam, i_d, to_dos = get_sams_list(
                self, tech, step, sam_fp, sams_dirs, db, focus)
            out = get_merge_outs(out_dir, step)
            if self.config.force or to_do(out):
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    cmd = get_merge_cmd(
                        self, tech, step, db, spc_list, out_dir, sam_fp)
                    c = cmd_sam + '\n' + cmd
                    self.outputs['cmds'].setdefault((tech, step), []).append(c)
                    io_update(self, i_d=i_d, o_d=out_dir, key=(tech, step))
                self.soft.add_status(tech, step, 1)
            else:
                self.soft.add_status(tech, step, 0)


def get_midas2_cmd(
        self,
        fastqs: list,
        focus_dir: str,
        db: str,
        db_path: str,
        spc_list: str,
        params: dict,
        step: str
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
    db : str
        MIDAS Database name
    db_path : str
        Path to local MIDAS Database
    spc_list : str
        path to file containing the list of species IDs to focus on
    params : dict
        Run parameters
    step : str
        Name of the MIDAS2 module

    Returns
    -------
    cmd : str
        MIDAS2 command
    """
    cmd = step.replace('_', ' run_')
    cmd += ' --sample_name %s' % self.sam_pool
    cmd += ' --force'
    cmd += ' -1 %s' % fastqs[0]
    if len(fastqs) == 2:
        cmd += ' -2 %s' % fastqs[1]
    if step in ['midas2_genes', 'midas2_snps']:
        cmd += ' --debug'
        for param in ['prebuilt_bowtie2_indexes', 'prebuilt_bowtie2_species']:
            if params[param]:
                cmd += ' --%s %s' % (param, params[param])

        if params['aln_interleaved']:
            cmd += ' --aln_interleaved'

        if spc_list:
            cmd += ' --species_list %s' % spc_list

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
            if step == 'midas2_genes':
                cmd += ' --aln_mapq 2'
            if step in 'midas2_snps':
                cmd += ' --aln_mapq 10'

        if params['chunk_size']:
            cmd += ' --chunk_size %s' % params['chunk_size']
        else:
            if step == 'midas2_genes':
                cmd += ' --chunk_size 50000'
            if step in 'midas2_snps':
                cmd += ' --chunk_size 1000000'

        if step == 'midas2_genes':
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

    if 'uhgg' in db:
        cmd += ' --midasdb_name uhgg'
    elif 'gtdb' in db:
        cmd += ' --midasdb_name gtdb'
    elif 'newdb' in db:
        cmd += ' --midasdb_name newdb'
    elif 's3db' in db:
        cmd += ' --midasdb_name s3db'
    cmd += ' --midasdb_dir %s' % db_path
    cmd += ' --num_cores %s' % params['cpus']
    cmd += ' %s\n' % focus_dir

    if step in ['midas2_genes', 'midas2_snps']:
        cmd += 'rm -rf %s/%s/bt2_indexes\n' % (focus_dir, self.sam_pool)
        cmd += 'rm %s/%s/temp/*/*.bam*\n' % (focus_dir, self.sam_pool)
        cmd += 'rm %s/%s/temp/*/*/*.lz4\n' % (focus_dir, self.sam_pool)
    return cmd


def get_species_lists(
        self,
        params,
        run_step
) -> tuple:
    """

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
    params : dict
        Run parameters
    run_step : str
        Name of the MIDAS2 step

    Returns
    -------
    species_lists : dict
        Path to species IDs list per species list name
    out : str
        End of the path to the current output
    """
    species_lists = {'all_species': ''}
    if run_step == 'midas2_species':
        out = '%s/species/species_profile.tsv' % self.sam_pool
    else:
        species_lists = params['species_list']
        if run_step == 'midas2_genes':
            out = '%s/genes/genes_summary.tsv' % self.sam_pool
        else:
            out = '%s/snps/snps_summary.tsv' % self.sam_pool
    return species_lists, out


def get_midas2(
        self,
        tech: str,
        sam: str,
        fastqs: list,
        out_dir: str,
        params: dict,
        step: str
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
    sam : str
        Sample name
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    step : str
        Name of the MIDAS2 step
    """
    if step not in self.outputs['outs']:
        self.outputs['outs'][step] = {}

    species_lists, out = get_species_lists(self, params, step)
    to_dos = status_update(self, tech, fastqs)
    for db in self.soft.params['databases']:

        db_path = self.databases.paths[db]
        if not self.config.dev and not isdir(db_path):
            sys.exit('[midas2] Path to database "%s not found' % db_path)

        for focus, spc_list in species_lists.items():
            focus_dir = out_dir + '/' + db + '/' + focus
            self.outputs['dirs'].append(focus_dir)
            if sam in self.soft.params['skip_samples']:
                continue
            self.outputs['outs'][step][(tech, db, focus, spc_list)] = focus_dir
            out_fp = focus_dir + '/' + out
            spc_dir = focus_dir + '/%s/species' % self.sam_pool
            if self.config.force or to_do(out_fp):
                cmd = get_midas2_cmd(self, fastqs, focus_dir, db, db_path,
                                     spc_list, params, step)
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                    i_f, i_d = fastqs, []
                    if not to_do(spc_list):
                        i_f += [spc_list]
                    if not to_do(spc_dir + '/species_profile.tsv'):
                        i_d.append(spc_dir)
                    io_update(self, i_f=i_f, i_d=i_d, o_d=focus_dir, key=(tech,))
                self.soft.add_status(tech, self.sam_pool, 1)
            else:
                self.soft.add_status(tech, self.sam_pool, 0)


def snps(
        self,
        tech: str,
        sam: str,
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
    sam : str
        Sample name
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, sam, fastqs, out_dir, params, 'midas2_snps')


def genes(
        self,
        tech: str,
        sam: str,
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
    sam : str
        Sample name
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, sam, fastqs, out_dir, params, 'midas2_genes')


def species(
        self,
        tech: str,
        sam: str,
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
    sam : str
        Sample name
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    get_midas2(self, tech, sam, fastqs, out_dir, params, 'midas2_species')


def midas2_(
        self,
        tech: str,
        sam: str,
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
    sam : str
        Sample name
    fastqs : list
        Path(s) to the input files
    out_dir : str
        Path to the output folder
    params : dict
        Run parameters
    """
    for step in ['midas2_species', 'midas2_genes', 'midas2_snps']:
        get_midas2(self, tech, sam, fastqs, out_dir, params, step)


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

    if self.sam_pool == '':
        merge(self)
    else:
        if self.soft.prev.startswith('midas2'):
            if 'midas2' in self.softs:
                reads = self.softs[self.softs['midas2'].prev].outputs
            elif 'midas2_species' in self.softs:
                reads = self.softs[self.softs['midas2_species'].prev].outputs
            else:
                sys.exit('[%s] No reads found used for "midas2[_species]"' %
                         self.soft.name)
        else:
            reads = self.inputs

        for (tech, sam), fastqs in reads[self.sam_pool].items():
            if tech_specificity(self, fastqs, tech, sam):
                continue
            out_dir = self.dir + '/' + tech
            module_call = caller(self, __name__)
            params = tech_params(self, tech)
            module_call(self, tech, sam, fastqs, out_dir, params)
