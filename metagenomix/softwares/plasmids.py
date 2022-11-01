# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from os.path import dirname
from metagenomix._inputs import (sample_inputs, group_inputs,
                                 genome_key, genome_out_dir)
from metagenomix._io_utils import io_update, to_do, status_update
from metagenomix.core.parameters import tech_params


def plasforest_cmd(
        self,
        in_fp: str,
        out_fp: str,
        out_dir: str
) -> str:
    """Collect plasforest command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    in_fp : str
        Path to the input file
    out_fp : str
        Path to the output file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        plasforest command
    """
    binary = self.soft.params['binary']
    cmd = 'cd %s\n' % out_dir
    cmd += 'cp %s/plasforest.sav %s/.\n' % (dirname(binary), out_dir)
    cmd += 'cp %s/*.fasta* %s/.\n' % (dirname(binary), out_dir)
    cmd += 'python3 %s' % binary
    cmd += ' -i %s' % in_fp
    cmd += ' -o %s' % out_fp
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    cmd += '\nrm %s/plasforest.sav\n' % out_dir
    cmd += 'rm %s/*.fasta*\n' % out_dir
    return cmd


def get_plasforest(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample or co-assembly name
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        out_fp = '%s/plasmids.csv' % out_dir
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_fp)
        to_dos = status_update(self, tech, [fasta[0]], group=sam_group,
                               genome=genome)

        if self.config.force or to_do(out_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = plasforest_cmd(self, fasta[0], out_fp, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta[0], i_d=out_dir, o_f=out_fp, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def plasmidfinder_cmd(
        self,
        fasta: str,
        out_dir: str,
        key: tuple
) -> str:
    """Collect PlasmidFinder command.

    Parameters
    ----------
    self : Commands class instance
        .pool : str
            Pool name.
        .soft.params
            PlasmidFinder parameters
    fasta : str
        Path to the input file
    out_dir : str
        Path to the output folder for the current sample/MAG
    key : str
        Technology and/or co-assembly pool group name

    Returns
    -------
    cmd : str
        PlasmidFinder command
    """
    tmp_dir = '$TMPDIR/plasmidfinder_%s' % '_'.join(key)
    cmd = 'mkdir -p %s\n' % tmp_dir
    cmd += '%s' % self.soft.params['binary']
    if len(fasta) == 2:
        cmd += ' --infile %s' % ' '.join(fasta)
    else:
        cmd += ' --infile %s' % fasta[0]
    cmd += ' --outputPath %s' % out_dir
    cmd += ' --tmp_dir %s' % tmp_dir
    cmd += ' --methodPath %s' % self.soft.params['methodPath']
    cmd += ' --databasePath %s' % self.soft.params['databasePath']
    if 'databases' in self.soft.params:
        cmd += ' --databases %s' % self.soft.params['databases']
    cmd += ' --mincov %s' % (float(self.soft.params['mincov']) / 100)
    cmd += ' --threshold %s' % (float(self.soft.params['threshold']) / 100)
    if self.soft.params['extented_output']:
        cmd += ' --extented_output'
    cmd += '\nrm -rf %s\n' % tmp_dir
    return cmd


def get_plasmidfinder(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample or co-assembly name
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)
        to_dos = status_update(
            self, tech, [fasta[0]], group=sam_group, genome=genome)

        json_fp = '%s/data.json' % out_dir
        if self.config.force or to_do(json_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = plasmidfinder_cmd(self, fasta, out_dir, key)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, i_d=out_dir, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def dispatch(self) -> None:
    """Classify assembly contigs or genomes/MAGs as plasmids or not.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
    """
    __plasmid_tool__ = getattr(sys.modules[__name__], 'get_%s' % self.soft.name)

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            __plasmid_tool__(self, tech, fastas, group)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            __plasmid_tool__(self, tech, fastas, mags)

    else:
        if self.soft.name == 'plasmidfinder':
            tech_fastas = sample_inputs(self, raw=True)
        else:
            tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            __plasmid_tool__(self, tech, fastas, self.sam_pool)


def plasmidfinder(self) -> None:
    """PlasmidFinder is a tool for the identification and typing of Plasmid
    Replicons in Whole-Genome Sequencing (WGS).

    References
    ----------
    Carattoli, Alessandra, and Henrik Hasman. "PlasmidFinder and in silico
    pMLST: identification and typing of plasmid replicons in whole-genome
    sequencing (WGS)." Horizontal gene transfer. Humana, New York, NY,
    2020. 285-294.

    Notes
    -----
    BitBucket   : https://bitbucket.org/genomicepidemiology/plasmidfinder
    Paper       : https://doi.org/10.1007/978-1-4939-9877-7_20

    Parameters
    ----------
    self : Commands class instance
        .soft.name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for PlasmidFinder
        .soft.prev : str
            Previous software in the pipeline
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
        .sam_pool
            Sample of co-assembly group name
    """
    dispatch(self)


def plasforest(self) -> None:
    """A random forest classifier of contigs to identify contigs of plasmid
    origin in contig and scaffold genomes.

    References
    ----------
    Pradier, Léa, et al. "PlasForest: a homology-based random forest
    classifier for plasmid detection in genomic datasets." BMC bioinformatics
    22.1 (2021): 1-17.

    Notes
    -----
    GitHub  : https://github.com/leaemiliepradier/PlasForest
    Paper   : https://doi.org/10.1186/s12859-021-04270-w

    Parameters
    ----------
    self : Commands class instance
        .soft.name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for PlasForest
        .soft.prev : str
            Previous software in the pipeline
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
        .sam_pool
            Sample of co-assembly group name
    """
    dispatch(self)


def mob_typer_cmd(
        self,
        tech: str,
        fasta: str,
        typer_out: str,
        mge_report: str,
        key: tuple
) -> str:
    """Collect the mob_typer command line.

    Parameters
    ----------
    self
    tech : str
    fasta : str
    typer_out : str
    mge_report : str
    key : tuple

    Returns
    -------
    cmd : str
        mob_typer command line
    """
    params = tech_params(self, tech)
    tmp_dir = '$TMPDIR/mob_typer_%s' % '_'.join(key)
    cmd = 'mob_typer'
    cmd += ' --infile %s' % fasta
    cmd += ' --out_file %s' % typer_out
    cmd += ' --mge_report_file %s' % mge_report
    cmd += ' --analysis_dir %s' % tmp_dir
    cmd += ' --num_threads %s' % params['cpus']
    cmd += ' --sample_id %s_%s_%s' % (tech, self.sam_pool, key[-1])
    cmd += ' --force'
    for boolean in ['multi', 'keep_tmp']:
        if params[boolean]:
            cmd += ' --%s' % boolean

    for param in [
        'min_rep_evalue', 'min_mob_evalue', 'min_con_evalue', 'min_rpp_evalue',
        'min_rep_ident', 'min_mob_ident', 'min_con_ident', 'min_rpp_ident',
        'min_rep_cov', 'min_mob_cov', 'min_con_cov', 'min_rpp_cov',
        'min_length', 'min_overlap', 'primary_cluster_dist'
    ]:
        if params[param]:
            cmd += ' --%s %s' % (param, params[param])

    for path in [
        'plasmid_mash_db', 'plasmid_meta', 'plasmid_db_type',
        'plasmid_replicons', 'repetitive_mask', 'plasmid_mob',
        'plasmid_mpf', 'plasmid_orit', 'database_directory'
    ]:
        if params[path]:
            cmd += ' --%s %s' % (path, params[path])
    return cmd


def mob_recon_cmd(
        self,
        tech: str,
        fasta: str,
        out_dir: str,
        key: tuple
) -> str:
    """Collect the mob_recon command line.

    Parameters
    ----------
    self
    tech : str
    fasta : str
    out_dir : str
    mge_report : str
    key : tuple

    Returns
    -------
    cmd : str
        mob_recon command line
    """
    params = tech_params(self, tech)
    cmd = 'mob_recon'
    cmd += ' --infile %s' % fasta
    cmd += ' --outdir %s' % out_dir
    cmd += ' --num_threads %s' % params['cpus']

    if params['prefix']:
        cmd += ' --prefix %s' % params['prefix']
    else:
        cmd += ' --prefix %s_%s_%s' % (tech, self.sam_pool, key[-1])

    cmd += ' --force'
    for boolean in ['multi', 'keep_tmp', 'run_overhang']:
        if params[boolean]:
            cmd += ' --%s' % boolean

    if self.soft.prev == 'unicycler':
        cmd += ' --unicycler_contigs'

    for param in [
        'min_rep_evalue', 'min_mob_evalue', 'min_con_evalue', 'min_rpp_evalue',
        'min_rep_ident', 'min_mob_ident', 'min_con_ident', 'min_rpp_ident',
        'min_rep_cov', 'min_mob_cov', 'min_con_cov', 'min_rpp_cov',
        'min_length', 'min_overlap', 'primary_cluster_dist',
        'secondary_cluster_dist', 'mash_genome_neighbor_threshold',
        'max_contig_size', 'max_plasmid_size',
    ]:
        if params[param]:
            cmd += ' --%s %s' % (param, params[param])

    for path in [
        'plasmid_db', 'plasmid_mash_db', 'plasmid_meta', 'plasmid_db_type',
        'plasmid_replicons', 'repetitive_mask', 'plasmid_mob', 'plasmid_mpf',
        'plasmid_orit', 'database_directory',
        # Mask frequent plasmid sequences that are not in the organism
        'filter_db',
        # Mash check of close genomes to only blast to these chromosomes.
        # (nuanced/automatic to filter in specific genomic contexts)
        'genome_filter_db_prefix'
    ]:
        if params[path]:
            cmd += ' --%s %s' % (path, params[path])
    return cmd


def mob_cluster_cmd(
        self,
        tech: str,
        typer_report: str,
        out_dir: str
) -> str:
    """Collect the mob_cluster command line.

    Parameters
    ----------
    self
    tech : str
    typer_report : str
    out_dir : str

    Returns
    -------
    cmd : str
        mob_cluster command line
    """
    params = tech_params(self, tech)
    cmd = 'mob_cluster'
    if params['mode']:
        cmd += ' --mode %s' % params['mode']
    cmd += ' --infile %s' % params['new_plasmids']
    if params['mob_typer_file']:
        cmd += ' --mob_typer_file %s' % params['mob_typer_file']
    else:
        cmd += ' --mob_typer_file %s' % typer_report
    if params['taxonomy']:
        cmd += ' --taxonomy %s' % params['taxonomy']
    cmd += ' --outdir %s' % out_dir
    if params['ref_cluster_file']:
        cmd += ' --ref_cluster_file %s' % params['ref_cluster_file']
    if params['ref_fasta_file']:
        cmd += ' --ref_fasta_file %s' % params['ref_fasta_file']
    cmd += ' --num_threads %s' % params['cpus']
    cmd += ' --primary_cluster_dist %s' % params['primary_cluster_dist']
    cmd += ' --secondary_cluster_dist %s' % params['secondary_cluster_dist']
    return cmd


def mobsuite_cmds(
        self,
        tech: str,
        fasta: str,
        out_dir: str,
        typer_out: str,
        key: tuple
) -> str:
    """

    Parameters
    ----------
    self
    tech : str
    fasta : str
    out_dir : str
    typer_out : str
    key : tuple

    Returns
    -------
    cmd : str
    """
    cmd = ''
    typer_report = '%s/typer_mge_report.txt' % out_dir
    if to_do(typer_report):
        cmd += mob_typer_cmd(self, tech, fasta, typer_out, typer_report, key)

    mge_report = '%s/mge.report.txt' % out_dir
    if to_do(mge_report):
        cmd += mob_recon_cmd(self, tech, fasta, out_dir, key)

    params = tech_params(self, tech)
    if params['cluster'] or params['new_plasmids']:
        cmd += mob_cluster_cmd(self, tech, typer_report, out_dir)
    return cmd


def get_mobsuite(
        self,
        tech: str,
        contigs_dict: dict,
        group: str
) -> None:
    """

    Parameters
    ----------
    self
    tech : str
    contigs_dict : dict
    group : str
    """
    for genome, contigs in contigs_dict.items():

        out_dir = genome_out_dir(self, tech, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        fasta = contigs[0]
        to_dos = status_update(self, tech, [fasta], group=group)

        typer_out = '%s/mobtyper_results.txt' % out_dir
        if self.config.force or to_do(typer_out):
            key = (tech, group)
            cmds = mobsuite_cmds(self, tech, fasta, out_dir, typer_out, key)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmds)
            io_update(self, i_f=fasta, i_d=out_dir, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def mobsuite(self) -> None:
    """MOB-suite: Software softwares for clustering, reconstruction and typing of
    plasmids from draft assemblies.

    Plasmids are mobile genetic elements (MGEs), which allow for rapid
    evolution and adaption of bacteria to new niches through horizontal
    transmission of novel traits to different genetic backgrounds. The
    MOB-suite is designed to be a modular set of softwares for the typing and
    reconstruction of plasmid sequences from WGS assemblies.

    The MOB-suite depends on a series of databases which are too large to be
    hosted in git-hub. They can be downloaded or updated by running mob_init
    or if running any of the softwares for the first time, the databases will
    download and initialize automatically if you do not specify an alternate
    database location. However, they are quite large so the first run will
    take a long time depending on your connection and speed of your computer.

    References
    ----------
    Robertson, James, and John HE Nash. "MOB-suite: software softwares for
    clustering, reconstruction and typing of plasmids from draft assemblies."
    Microbial genomics 4.8 (2018).

    Notes
    -----
    GitHub  : https://github.com/phac-nml/mob-suite
    Paper   : https://doi.org/10.1099/mgen.0.000206

    Parameters
    ----------
    self
    """
    assemblers = self.config.tools['assembling']
    if self.soft.prev not in assemblers:
        sys.exit('[mobsuite] Only after assembly (%s)' % ''.join(assemblers))

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            contigs_dict = group_inputs(self, inputs)
            get_mobsuite(self, tech, contigs_dict, group)

