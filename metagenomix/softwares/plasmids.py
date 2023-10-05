# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import pkg_resources
from os.path import basename, dirname, splitext

from metagenomix._inputs import (
    sample_inputs, group_inputs, genome_key, get_contigs_from_path,
    genome_out_dir, get_plasmids)
from metagenomix._io_utils import caller, io_update, to_do, status_update
from metagenomix.core.parameters import tech_params

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def plasforest_cmd(
        self,
        fasta: str,
        out_dir: str
) -> str:
    """Collect plasforest command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fasta : str
        Path to the input file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        plasforest command
    """
    binary = self.soft.params['binary']
    cmd = 'cd %s\n' % out_dir
    cmd_rm = ''
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += '\ngunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        fasta = fasta.rstrip('.gz')
        cmd_rm += 'rm %s\n' % fasta

    cmd += 'cp %s/plasforest.sav %s/.\n' % (dirname(binary), out_dir)
    cmd += 'cp %s/*.fasta* %s/.\n' % (dirname(binary), out_dir)
    cmd += 'python3 %s' % binary
    cmd += ' -i %s' % fasta
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    cmd += ' -o %s/plasmids.csv\n' % out_dir

    cmd += 'rm %s/plasforest.sav\n' % out_dir
    cmd += 'rm %s/*.fasta*\n' % out_dir
    cmd += 'gzip -q %s/plasmids.csv\n' % out_dir
    cmd += cmd_rm
    return cmd


def get_plasforest(
        self,
        tech: str,
        sam_group: str,
        fastas: dict,
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample or co-assembly name
    fastas : dict
        Paths to the input fasta files per genome/MAG
    """
    for genome, fastas in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        out_fp = '%s/plasmids.csv.gz' % out_dir
        self.outputs['outs'][(tech, sam_group)][genome] = out_dir
        fasta = fastas[0]
        to_dos = status_update(self, tech, [fasta], group=sam_group,
                               genome=genome)

        if self.config.force or to_do(out_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = plasforest_cmd(self, fasta, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, i_d=out_dir, o_f=out_fp, key=key)
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
    cmd_rm = ''
    cmd = 'mkdir -p %s\n' % tmp_dir
    if len(fasta) == 2:
        infile = ' '.join(fasta)
    elif fasta[0].endswith('.fa.gz') or fasta[0].endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta[0], fasta[0].rstrip('.gz'))
        infile = fasta[0].rstrip('.gz')
        cmd_rm += 'rm %s\n' % infile
    else:
        infile = fasta[0]
    cmd += 'export PATH=$PATH:%s\n' % dirname(self.soft.params['methodPath'])
    if self.soft.params['binary']:
        cmd += '%s' % self.soft.params['binary']
    else:
        cmd += 'plasmidfinder.py'
    cmd += ' --infile %s' % infile
    cmd += ' --tmp_dir %s' % tmp_dir
    cmd += ' --databasePath %s' % self.soft.params['databasePath']
    if 'databases' in self.soft.params:
        cmd += ' --databases %s' % self.soft.params['databases']
    cmd += ' --mincov %s' % (float(self.soft.params['mincov']) / 100)
    cmd += ' --threshold %s' % (float(self.soft.params['threshold']) / 100)
    if self.soft.params['extented_output']:
        cmd += ' --extented_output'
    cmd += ' --outputPath %s\n' % out_dir

    cmd += 'rm -rf %s\n' % tmp_dir
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    cmd += cmd_rm
    return cmd


def get_plasmidfinder(
        self,
        tech: str,
        sam_group: str,
        fastas: dict,
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample or co-assembly name
    fastas : dict
        Paths to the input fasta files per genome/MAG
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][(tech, sam_group)][genome] = out_dir
        to_dos = status_update(
            self, tech, [fasta[0]], group=sam_group, genome=genome)

        json_fp = '%s/data.json.gz' % out_dir
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


def platon_cmd(
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
            platon parameters
    fasta : str
        Path to the input file
    out_dir : str
        Path to the output folder for the current sample/MAG
    key : str
        Technology and/or co-assembly pool group name

    Returns
    -------
    cmd : str
        platon command
    """
    tmp_dir = '$TMPDIR/platon_%s' % '_'.join(key)
    cmd_rm = ''
    cmd = 'mkdir -p %s\n' % tmp_dir
    if len(fasta) == 2:
        infile = ' '.join(fasta)
    elif fasta[0].endswith('.fa.gz') or fasta[0].endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta[0], fasta[0].rstrip('.gz'))
        infile = fasta[0].rstrip('.gz')
        cmd_rm += 'rm %s\n' % infile
    else:
        infile = fasta[0]

    cmd += 'platon'
    cmd += ' --prefix output'
    cmd += ' --db %s' % self.databases.paths['platon']
    for boolean in ['characterize', 'meta']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    cmd += ' --mode %s' % self.soft.params['mode']
    cmd += ' --output %s' % out_dir
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' %s\n' % infile

    cmd += 'rm -rf %s\n' % tmp_dir
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    cmd += cmd_rm
    return cmd


def get_platon(
        self,
        tech: str,
        sam_group: str,
        fastas: dict,
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample or co-assembly name
    fastas : dict
        Paths to the input fasta files per genome/MAG
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][(tech, sam_group)][genome] = out_dir
        to_dos = status_update(
            self, tech, [fasta[0]], group=sam_group, genome=genome)

        out_fp = '%s/output.tsv.gz' % out_dir
        if self.config.force or to_do(out_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = platon_cmd(self, fasta, out_dir, key)
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
            self.outputs['outs'][(tech, group)] = {}
            fastas = group_inputs(self, inputs)
            __plasmid_tool__(self, tech, group, fastas)
    else:
        if self.soft.name == 'plasmidfinder':
            tech_fastas = sample_inputs(self, raw=True)
        else:
            tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            self.outputs['outs'][(tech, self.sam_pool)] = {}
            __plasmid_tool__(self, tech, self.sam_pool, fastas)


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


def platon(self):
    """Platon detects plasmid-borne contigs within bacterial draft (meta)
    genomes assemblies. Therefore, Platon analyzes the distribution bias of
    protein-coding gene families among chromosomes and plasmids.
    This analysis is complemented by comprehensive contig characterizations
    followed by heuristic filters.

    Notes
    -----
    GitHub  : https://github.com/oschwengers/platon
    Paper   : https://doi.org/10.1099/mgen.0.000398

    Parameters
    ----------
    self
    """
    dispatch(self)


def mob_typer_cmd(
        self,
        tech: str,
        fasta: str,
        typer_dir: str,
        key: tuple
) -> str:
    """Collect the mob_typer command line.

    Parameters
    ----------
    self
    tech : str
    fasta : str
    typer_dir : str
    key : tuple

    Returns
    -------
    cmd : str
        mob_typer command line
    """
    params = tech_params(self, tech)
    tmp_dir = '$TMPDIR/mob_typer_%s' % '_'.join(key)

    cmd, cmd_rm = '', ''
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        fasta = fasta.rstrip('.gz')
        cmd_rm += 'rm %s\n' % fasta

    cmd_typer = '\nmob_typer'
    cmd_typer += ' --infile %s' % fasta
    cmd_typer += ' --analysis_dir %s' % tmp_dir
    cmd_typer += ' --num_threads %s' % params['cpus']
    cmd_typer += ' --force'
    for boolean in ['multi', 'keep_tmp']:
        if params[boolean]:
            cmd_typer += ' --%s' % boolean

    for param in [
        'min_rep_evalue', 'min_mob_evalue', 'min_con_evalue', 'min_rpp_evalue',
        'min_rep_ident', 'min_mob_ident', 'min_con_ident', 'min_rpp_ident',
        'min_rep_cov', 'min_mob_cov', 'min_con_cov', 'min_rpp_cov',
        'min_length', 'min_overlap', 'primary_cluster_dist'
    ]:
        if params[param]:
            cmd_typer += ' --%s %s' % (param, params[param])

    for path in [
        'plasmid_mash_db', 'plasmid_meta', 'plasmid_db_type',
        'plasmid_replicons', 'repetitive_mask', 'plasmid_mob',
        'plasmid_mpf', 'plasmid_orit', 'database_directory'
    ]:
        if params[path]:
            cmd_typer += ' --%s %s' % (path, params[path])

    cmd_typer += ' --mge_report_file %s/mge.report.txt' % typer_dir
    cmd_typer += ' --out_file %s/mobtyper_results.txt\n' % typer_dir

    cmd += 'if [ -s %s ]; then\n' % fasta
    cmd += cmd_typer
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % typer_dir
    cmd += 'fi\n'
    cmd += cmd_rm
    return cmd


def mob_recon_cmd(
        self,
        tech: str,
        fasta: str,
        recon_dir: str,
) -> str:
    """Collect the mob_recon command line.

    Parameters
    ----------
    self
    tech : str
    fasta : str
    recon_dir : str

    Returns
    -------
    cmd : str
        mob_recon command line
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        fasta = fasta.rstrip('.gz')
        cmd_rm += 'rm %s\n' % fasta

    cmd += '\nmob_recon'
    cmd += ' --infile %s' % fasta
    cmd += ' --num_threads %s' % params['cpus']

    if params['prefix']:
        cmd += ' --prefix %s' % params['prefix']

    cmd += ' --force'
    for boolean in ['keep_tmp', 'run_overhang']:
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
    cmd += ' --outdir %s\n' % recon_dir

    cmd += 'for i in %s/*; do gzip -q $i; done\n' % recon_dir
    cmd += cmd_rm
    return cmd


def mob_cluster_cmd(
        self,
        tech: str,
        cluster_dir: str
) -> str:
    """Collect the mob_cluster command line.

    Parameters
    ----------
    self
    tech : str
    cluster_dir : str

    Returns
    -------
    cmd : str
        mob_cluster command line
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''

    typer_report = '%s/../typer/mge.report.txt.gz' % cluster_dir
    if params['mob_typer_file']:
        typer_report = params['mob_typer_file']

    if typer_report.endswith('.gz'):
        cmd += 'gunzip -c %s > %s' % (typer_report, typer_report.rstrip('.gz'))
        typer_report = typer_report.rstrip('.gz')
        cmd_rm += 'rm %s\n' % typer_report.rstrip('.gz')

    cmd += 'mob_cluster'
    if params['mode']:
        cmd += ' --mode %s' % params['mode']

    cmd += ' --infile %s' % params['new_plasmids']
    cmd += ' --mob_typer_file %s' % typer_report
    if params['taxonomy']:
        cmd += ' --taxonomy %s' % params['taxonomy']

    cmd += ' --outdir %s' % cluster_dir
    if params['ref_cluster_file']:
        cmd += ' --ref_cluster_file %s' % params['ref_cluster_file']
    if params['ref_fasta_file']:
        cmd += ' --ref_fasta_file %s' % params['ref_fasta_file']

    cmd += ' --num_threads %s' % params['cpus']
    cmd += ' --primary_cluster_dist %s' % params['primary_cluster_dist']
    cmd += ' --secondary_cluster_dist %s\n' % params['secondary_cluster_dist']

    cmd += 'for i in %s/*; do gzip -q $i; done\n' % cluster_dir
    cmd += cmd_rm
    return cmd


def mobsuite_cmds(
        self,
        tech: str,
        fasta: str,
        out_dir: str,
        key: tuple
) -> str:
    """

    Parameters
    ----------
    self
    tech : str
    fasta : str
    out_dir : str
    key : tuple

    Returns
    -------
    cmd : str
    """
    cmd = ''
    typer_dir = '%s/typer' % out_dir
    self.outputs['dirs'].append(typer_dir)
    if to_do('%s/mge.report.txt.gz' % typer_dir):
        cmd += mob_typer_cmd(self, tech, fasta, typer_dir, key)

    if self.config.tools[self.soft.prev] != 'annotation (plasmid)':
        recon_dir = '%s/recon' % out_dir
        self.outputs['dirs'].append(recon_dir)
        if to_do('%s/mge.report.txt.gz' % recon_dir):
            cmd += mob_recon_cmd(self, tech, fasta, recon_dir)

        params = tech_params(self, tech)
        if params['cluster'] or params['new_plasmids']:
            cluster_dir = '%s/cluster' % out_dir
            self.outputs['dirs'].append(cluster_dir)
            cmd += mob_cluster_cmd(self, tech, cluster_dir)

    return cmd


def get_mobsuite(
        self,
        tech: str,
        group: str,
        inputs_dict: dict,
        contigs: str
) -> None:
    """

    Parameters
    ----------
    self
    tech : str
    group : str
    inputs_dict : dict
    contigs: str
        Empty of not run after a plasmid-detection tool
    """
    for genome, fastas in inputs_dict.items():

        key = (tech, group)
        out_dir = genome_out_dir(self, tech, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        if contigs:
            plasmids, fasta, cmds, rms = get_plasmids(self, fastas, contigs)
            i_f = [plasmids, contigs]
        else:
            fasta, cmds, rms = fastas[0], '', ''
            i_f = [fasta]
        to_dos = status_update(self, tech, i_f, group=group)

        typer_out = '%s/typer/mobtyper_results.txt.gz' % out_dir
        if self.config.force or to_do(typer_out):
            cmds += mobsuite_cmds(self, tech, fasta, out_dir, key)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmds)
            io_update(self, i_f=i_f, i_d=out_dir, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def mobsuite(self) -> None:
    """MOB-suite: Software softwares for clustering, reconstruction and
    typing of plasmids from draft assemblies.

    Plasmids are mobile genetic elements (MGEs), which allow for rapid
    evolution and adaption of bacteria to new niches through horizontal
    transmission of novel traits to different genetic backgrounds. The
    MOB-suite is designed to be a modular set of softwares for the typing
    and reconstruction of plasmid sequences from WGS assemblies.

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
    previous = self.config.tools[self.soft.prev]
    if not (self.soft.prev in assemblers or previous == 'annotation (plasmid)'):
        sys.exit('[mobsuite] Only after assembly (%s)' % ', '.join(assemblers))

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            inputs_dict = group_inputs(self, inputs)
            contigs = ''
            if previous == 'annotation (plasmid)':
                contigs = get_contigs_from_path(self, tech, group)
            get_mobsuite(self, tech, group, inputs_dict, contigs)


def oritfinder_cmd(
        gbk: str,
) -> str:
    """Get oriTfinder command.

    Parameters
    ----------
    gbk : str
        Path to the input gbk file

    Returns
    -------
    cmd : str
        oriTfinder command
    """
    cmd = 'echo "%s" > list.txt\n' % gbk.rstrip('.gz')
    cmd += 'perl oriTfinder_local.pl list.txt'
    return cmd


def get_oritfinder(
        self,
        tech: str,
        fastas: dict,
        group: str
):
    """

    Parameters
    ----------
    self
    tech: str
    fastas : dict
    group: str
    """
    for genome, folder in fastas.items():
        out_dir = genome_out_dir(self, tech, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        gbk = '%s/gene.coords.gbk.gz' % folder[0]
        to_dos = status_update(self, tech, [gbk], group=group)

        cmd = ''
        out = '%s/gene.coords_summary.txt.gz' % out_dir
        key = (tech, group)
        if self.config.force or to_do(out):
            cmd += oritfinder_cmd(gbk)

        if cmd:
            path = self.soft.params['path']
            full_cmd = 'cd %s\n' % out_dir
            full_cmd += 'mkdir gbk\n'
            full_cmd += 'cp %s gbk/.\n' % gbk
            full_cmd += 'gunzip -c %s > %s\n' % (gbk, gbk.rstrip('.gz'))
            full_cmd += 'cp -r %s/data .\n' % path
            full_cmd += 'cp -r %s/tools .\n' % path
            full_cmd += 'cp -r %s/scripts .\n' % path
            full_cmd += 'cp %s/oriTfinder_local.pl .\n' % path
            full_cmd += cmd
            full_cmd += 'rm %s\n' % gbk.rstrip('.gz')
            full_cmd += 'rm -rf data\n'
            full_cmd += 'rm -rf tools\n'
            full_cmd += 'rm -rf scripts\n'
            full_cmd += 'rm oriTfinder_local.pl\n'
            full_cmd += 'for i in %s/*; do gzip -q %s; done\n' % out_dir
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(full_cmd)
            io_update(self, i_f=gbk, i_d=out_dir, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def oritfinder(self):
    """Identification of origin of transfers in DNA sequences of bacterial
    mobile genetic elements.

    References
    ----------
    Li, X., Xie, Y., Liu, M., Tai, C., Sun, J., Deng, Z. and Ou, H.Y.,
    2018. oriTfinder: a web-based tool for the identification of origin of
    transfers in DNA sequences of bacterial mobile genetic elements. Nucleic
    acids research, 46(W1), pp.W229-W234.

    Notes
    -----
    Docs    : http://bioinfo-mml.sjtu.edu.cn/oriTfinder
    Paper   : https://doi.org/10.1093/nar/gky352

    Parameters
    ----------
    self
    """
    name = self.soft.name
    if self.sam_pool in self.pools:
        if self.soft.prev != 'prodigal':
            sys.exit('[%s] Runs on protein data (plass, prodigal...)' % name)
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_oritfinder(self, tech, fastas, group)
    else:
        prev = self.soft.prev
        if prev not in ['plass', 'prodigal']:
            sys.exit('[%s] Runs on protein data (plass, prodigal...)' % name)
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_oritfinder(self, tech, fastas, self.sam_pool)


def get_genomad(
        self,
        tech: str,
        fastas: list,
        group: str
) -> None:
    """Get the geNomad command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Databases
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample name or group for the current co-assembly
    """
    for genome, fasta_ in fastas.items():

        key = genome_key(tech, group, genome)
        out = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out)
        self.outputs['outs'].setdefault((tech, group), []).append(out)

        fasta = fasta_[0]
        if self.soft.prev == 'prodigal':
            fasta = '%s/nucleotide.sequences.fasta.gz' % fasta_[0]

        print(fasta)
        to_dos = status_update(self, tech, [fasta], group=group, genome=genome)

        fpo = '%s/Results_Integron_Finder_mysequences/mysequences.summary' % out
        if self.config.force or to_do(fpo):
            # cmd = hmms_cmd + integronfinder_cmd(self, tech, fasta, out)
            cmd = genomad_cmd(self, tech, fasta, out)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            # io_update(self, i_f=[fasta, hmms_sh], o_d=out, key=key)
            io_update(self, i_f=fasta, o_d=out, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def genomad_end_to_end_cmd(
        self,
        tech: str,
        fasta: str,
        out_dir: str,
        prefix: str
):
    """

    Parameters
    ----------
    self
    tech : str
    fasta : str
    out_dir : str
    prefix : str

    Returns
    -------
    cmd : str
        genomad end-to-end command
    """
    params = tech_params(self, tech)
    cmd = 'genomad end-to-end'
    cmd += ' --threads %s' % params['cpus']
    if params['presets']:
        cmd += ' --%s' % params['presets']
    else:
        for param in [
            'min_score',
            'max_fdr',
            'min_plasmid_marker_enrichment',
            'min_virus_marker_enrichment',
            'min_plasmid_hallmarks',
            'min_plasmid_hallmarks_short_seqs',
            'min_virus_hallmarks',
            'min_virus_hallmarks_short_seqs',
            'max_uscg'
        ]:
            cmd += ' --%s %s' % (param.replace('_', '-'), params[param])

    for param in ['sensitivity', 'splits', 'composition']:
        cmd += ' --%s %s' % (param, params[param])

    for boolean in [
        'cleanup', 'restart', 'disable_find_proviruses',
        'disable_nn_classification', 'enable_score_calibration',
        'conservative_taxonomy', 'skip_integrase_identification',
        'skip_trna_identification', 'force_auto'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' %s' % fasta
    cmd += ' %s' % out_dir
    cmd += ' %s\n' % self.databases.paths['genomad']
    for suffix in [
        'aggregated_classification', 'annotate', 'find_proviruses',
        'marker_classification', 'nn_classification', 'summary'
    ]:
        cmd += 'tar cpfz %s/%s_%s.tar.gz -C %s %s_%s\n' % (
            out_dir, prefix, suffix, out_dir, prefix, suffix)
        cmd += 'rm -rf %s/%s_%s\n' % (out_dir, prefix, suffix)
    return cmd


def annotate_cmd(self):
    cmd = 'genomad annotate'  # metagenome.fna genomad_output genomad_db
    return cmd


def annotate(self):
    gene_prediction = ['prodigal']
    if self.soft.prev not in gene_prediction:
        sys.exit('[genomad_findproviruses] Only after MMseqs')

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_genomad(self, tech, fastas, group)


def find_proviruses_cmd(self):
    cmd = 'genomad find-proviruses'  # metagenome.fna genomad_output genomad_db
    return cmd


def marker_classification_cmd(self):
    cmd = 'genomad marker-classification'  # metagenome.fna genomad_output genomad_db
    return cmd


def nn_classification_cmd(self):
    cmd = 'genomad nn-classification'  # metagenome.fna genomad_output
    return cmd


def aggregated_classification_cmd(self):
    cmd = 'genomad aggregated-classification'  # metagenome.fna genomad_output
    return cmd


def score_calibration_cmd(self):
    cmd = 'genomad score-calibration'  # metagenome.fna genomad_output
    return cmd


def summary_cmd(self):
    cmd = 'genomad summary'  # metagenome.fna genomad_output
    return cmd


def get_genomad_end_to_end(
        self,
        tech: str,
        folders: dict,
        group: str
):
    """

    Parameters
    ----------
    self
    tech
    folders
    group
    """
    for genome, fastas in folders.items():
        key = genome_key(tech, group, genome)
        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        fasta = fastas[0]
        to_dos = status_update(self, tech, [fasta], group=group, genome=genome)

        prefix = splitext(basename(fasta.replace('.gz', '')))[0]
        fpo = '%s/%s_summary.log' % (out_dir, prefix)
        if self.config.force or to_do(fpo):
            cmd = genomad_end_to_end_cmd(self, tech, fasta, out_dir, prefix)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def genomad_(self) -> list:
    """Collect the name of the functions to run end-to-end geNomad.

    Parameters
    ----------
    self
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_genomad_end_to_end(self, tech, folders, group)


def genomad(self):
    """geNomad's primary goal is to identify viruses and plasmids in
    sequencing data (isolates, metagenomes, and metatranscriptomes). It also
    provides a couple of additional features that can help you in your analysis:
        - Taxonomic assignment of viral genomes.
        - Identification of viruses integrated in host genomes (proviruses).
        - Functional annotation of proteins.

    References
    ----------
    Camargo, A.P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain,
    P.S., Nayfach, S. and Kyrpides, N.C., 2023. Identification of mobile genetic
    elements with geNomad. Nature Biotechnology, pp.1-10.

    Notes
    -----
    GitHub  : https://github.com/apcamargo/genomad
    Paper   : https://doi.org/10.1038/s41587-023-01953-y
    Docs    : https://portal.nersc.gov/genomad/index.html

    Parameters
    ----------
    self
    """
    if 'genomad' not in self.databases.paths:
        sys.exit('[genomad] Path to "genomad" database folder needed')
    module_call = caller(self, __name__)
    module_call(self)


def rfplasmid(self):
    """Identifies the number of chromosomal marker genes, plasmid replication
    genes and plasmid typing genes using CheckM and DIAMOND Blast, and
    determines pentamer frequencies and contig sizes per contig. A prediction
    model was trained using Random Forest on an extensive set of plasmids and
    chromosomes from 19 different bacterial species and validated on separate
    test sets of known chromosomal and plasmid contigs of the different
    bacteria.

    References
    ----------
    van der Graaf-Van Bloois, L., Wagenaar, J.A. and Zomer, A.L.,
    2021. RFPlasmid: predicting plasmid sequences from short-read assembly
    data using machine learning. Microbial genomics, 7(11).

    Notes
    -----
    GitHub  : https://github.com/aldertzomer/RFPlasmid
    Paper   : https://doi.org/10.1099/mgen.0.000683
    Docs    : http://klif.uu.nl/rfplasmid/

    Parameters
    ----------
    self
    """
    pass


def deeplasmid(self):
    """Deeplasmid is a tool based on machine learning that separates plasmids
    from chromosomal sequences. It can identify plasmids in microbial isolate
    or metagenome assemblies. The input sequences are in the form of contigs
    and could have been produced from any sequencing technology or assembly
    algorithm.

    References
    ----------
    Andreopoulos, W.B., Geller, A.M., Lucke, M., Balewski, J., Clum, A.,
    Ivanova, N.N. and Levy, A., 2022. Deeplasmid: deep learning accurately
    separates plasmids from bacterial chromosomes. Nucleic acids research,
    50(3), pp.e17-e17.

    Notes
    -----
    GitHub  : https://github.com/wandreopoulos/deeplasmid
    Paper   : https://doi.org/10.1093/nar/gkab1115

    Parameters
    ----------
    self
    """
    pass


def plasflow(self):
    """PlasFlow is a set of scripts used for prediction of plasmid sequences
    in metagenomic contigs. It relies on the neural network models trained on
    full genome and plasmid sequences and is able to differentiate between
    plasmids and chromosomes with accuracy reaching 96%. It outperforms other
    available solutions for plasmids recovery from metagenomes and incorporates
    the thresholding which allows for exclusion of incertain predictions.

    References
    ----------
    Krawczyk, P.S., Lipinski, L. and Dziembowski, A., 2018. PlasFlow:
    predicting plasmid sequences in metagenomic data using genome signatures.
    Nucleic acids research, 46(6), pp.e35-e35.

    Notes
    -----
    GitHub  : https://github.com/smaegol/PlasFlow
    Paper   : https://doi.org/10.1093/nar/gkx1321

    Parameters
    ----------
    self
    """
    pass


def pprmeta(self):
    """PPR-Meta is designed to identify metagenomic sequences as phages,
    chromosomes or plasmids. The program calculate three score reflecting the
    likelihood of each input fragment as phage, chromosome or plasmid.
    PPR-Meta can run either on the virtual machine or physical host. For
    non-computer professionals, we recommend running the virtual machine
    version of PPR-Meta on local PC. In this way, users do not need to
    install any dependency package. If GPU is available, you can also choose
    to run the physical host version. This version can automatically speed up
    with GPU and is more suitable to handle large scale data.

    References
    ----------
    Fang, Z., Tan, J., Wu, S., Li, M., Xu, C., Xie, Z. and Zhu, H.,
    2019. PPR-Meta: a tool for identifying phages and plasmids from
    metagenomic fragments using deep learning. GigaScience, 8(6), p.giz066.

    Notes
    -----
    GitHub  : https://github.com/zhenchengfang/PPR-Meta
    Paper   : https://doi.org/10.1093/gigascience/giz066
    Docs    : http://cqb.pku.edu.cn/ZhuLab/PPR_Meta/

    Parameters
    ----------
    self
    """
    pass


def plasclass(self):
    """This module allows for easy classification of sequences as either
    plasmid or chromosomal. For example, it can be used to classify the
    contigs in a (metagenomic) assembly.

    References
    ----------
    Pellow D, Mizrahi I, Shamir R (2020) PlasClass improves plasmid sequence
    classification. PLoS Comput Biol 16(4): e1007781.

    Notes
    -----
    GitHub  : https://github.com/Shamir-Lab/PlasClass
    Paper   : https://doi.org/10.1371/journal.pcbi.1007781

    Parameters
    ----------
    self
    """
    pass
