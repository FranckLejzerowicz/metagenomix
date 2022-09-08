# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from metagenomix._io_utils import caller, io_update, to_do
from metagenomix._inputs import (group_inputs, genome_key, genome_out_dir,
                                 get_extension, get_assembler, add_folder)


def get_bin_folders(
        self,
        bin_paths: list
) -> dict:
    """Get the paths to the bin folder(s) per binning software/step key.

    Parameters
    ----------
    self
    bin_paths : list
        Paths to the bin folder(s) of the current binning software/step

    Returns
    -------
    bins : dict
        Paths to the bin folder(s) per binning software/step key
    """
    bins = {}
    if self.soft.prev == 'metawrap_refine':
        bins['refined_bins'] = bin_paths[0]
    elif self.soft.prev == 'metawrap_reassemble':
        bins = {}
        for group_path in bin_paths:
            bins['%s_reassembly' % group_path.split('_')[-1]] = group_path
    else:
        sys.exit('[drep] Not possible after "%s"' % self.soft.prev)
    return bins


def get_drep_bins(self) -> dict:
    """Get the genomes in an iterable format accommodating
    both the output of `metawrap_refine()` and of `metawrap_reassemble()`.

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software in the pipeline
        .pools : dict
            Samples per group as dispatched according to the pooling design
        .inputs : dict
            Path to the input files per pool

    Returns
    -------
    genomes : nested dict
        Path to folders with bins to dereplicate per technology and per binning:
        "stringent" or "permissive" after reassembly, empty for others.
    """
    genomes = {}
    for pool in self.pools:
        for (tech, group), bin_paths in self.inputs[pool].items():
            if (tech, pool) not in genomes:
                genomes[(tech, pool)] = {}
            bins_folders = get_bin_folders(self, bin_paths)
            for binning, bin_folder in bins_folders.items():
                genomes[(tech, pool)].setdefault(binning, []).append(bin_folder)
    return genomes


def get_drep_inputs(
        self,
        drep_dir: str,
        paths: list
) -> tuple:
    """Write the file containing the inputs bins to drep
    and the list of paths to these bins.

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    drep_dir :
        Path to the drep output folder
    paths : list
        Paths to the folders with the genome files to dereplicate

    Returns
    -------
    cmd : str
        Command to create the input
    drep_in : str
        File containing the paths corresponding to each bin
    bin_paths : int
        Genome files to use for dereplication
    """
    cmd = ''
    drep_in = '%s/input_genomes.txt' % drep_dir
    bin_paths = []
    for path in paths:
        bin_paths = glob.glob('%s/*fa' % path.replace('${SCRATCH_FOLDER}', ''))
        if self.config.dev:
            bin_paths = ['%s/a.fa' % path.replace('${SCRATCH_FOLDER}', ''),
                         '%s/b.fa' % path.replace('${SCRATCH_FOLDER}', '')]
        for bin_path in bin_paths:
            if cmd:
                cmd += 'echo "%s" >> %s\n' % (bin_path, drep_in)
            else:
                cmd += 'echo "%s" > %s\n' % (bin_path, drep_in)
    if cmd:
        cmd += 'envsubst < %s > %s.tmp\n' % (drep_in, drep_in)
        cmd += 'mv %s.tmp %s\n' % (drep_in, drep_in)
    return cmd, drep_in, bin_paths


def drep_cmd(
        self,
        algorithm: str,
        drep_in: str,
        drep_out: str,
        bin_paths: list,
        input_cmd: str
) -> str:
    """Collect dRep dereplicate command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    algorithm : str
        Alignment algorithm: 'fastANI', 'ANIn', 'ANImf', 'gANI', or 'goANI'
    drep_in : str
        File containing the paths corresponding to each bin
    drep_out : str
        Path to the output folder
    bin_paths : int
        Genome files to use for dereplication
    input_cmd : str
        Command to create the input

    Returns
    -------
    cmd : str
        dRep dereplicate command
    """
    cmd = '%s\ndRep dereplicate' % input_cmd
    cmd += ' %s' % drep_out
    cmd += ' --S_algorithm %s' % algorithm
    cmd += ' --ignoreGenomeQuality'
    cmd += ' --processors %s' % self.soft.params['cpus']
    if len(bin_paths) > self.soft.params['primary_chunksize']:
        cmd += ' --multiround_primary_clustering'
        cmd += ' --run_tertiary_clustering'
        if algorithm == 'fastANI':
            cmd += ' --greedy_secondary_clustering'
    else:
        for boolean in ['multiround_primary_clustering',
                        'greedy_secondary_clustering',
                        'run_tertiary_clustering']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean
    for boolean in ['SkipMash', 'SkipSecondary']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    cmd += ' --coverage_method %s' % self.soft.params['coverage_method']
    for param in [
        'S_ani', 'P_ani', 'warn_sim', 'warn_aln', 'n_PRESET',
        'warn_dist', 'cov_thresh', 'clusterAlg', 'MASH_sketch'
    ]:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    cmd += ' --genomes %s' % drep_in
    return cmd


def drep(self):
    """Dereplicate the binned genomes to obtain MAGs across samples using dRep.

    References
    ----------
    Olm, Matthew R., et al. "dRep: a tool for fast and accurate genomic
    comparisons that enables improved genome recovery from metagenomes
    through de-replication." The ISME journal 11.12 (2017): 2864-2868.

    Notes
    -----
    GitHub  : https://github.com/MrOlm/drep
    Docs    : https://drep.readthedocs.io/en/latest
    Paper   : https://doi.org/10.1038/ismej.2017.126

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for dRep
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    genomes = get_drep_bins(self)
    for (tech, pool), pool_paths in genomes.items():
        self.outputs['outs'][''] = {}
        for binning, paths in pool_paths.items():
            for algo in self.soft.params['S_algorithm']:

                bin_algo = '_'.join([binning, algo])
                pool_binning_algo = pool + '/' + bin_algo
                drep_out = '/'.join([self.dir, tech, pool_binning_algo])
                self.outputs['dirs'].append(drep_out)

                dereps = '%s/dereplicated_genomes' % drep_out
                self.outputs['outs'][''][(tech, pool_binning_algo)] = [dereps]

                cmd, drep_in, bin_paths = get_drep_inputs(self, drep_out, paths)
                if not bin_paths:
                    self.soft.add_status(tech, pool, paths, group=bin_algo,
                                         message='run previous')
                    if self.config.dev:
                        bin_paths = ['x'] * 5001

                out_dereps = '%s/*.fa' % dereps.replace('${SCRATCH_FOLDER}', '')
                if not self.config.force and glob.glob(out_dereps):
                    self.soft.add_status(tech, pool, 0, group=bin_algo)
                    continue
                key = (tech, pool_binning_algo)
                cmd = drep_cmd(self, algo, drep_in, drep_out, bin_paths, cmd)
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_d=paths, o_d=drep_out, key=key)
                self.soft.add_status(tech, pool, 1, group=bin_algo)


def tree_cmd(
        self,
        genomes_dir: str,
        tree_dir: str
) -> str:
    """Collect the command line for checkm tree.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    genomes_dir : str
        Path to the input folder
    tree_dir : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm tree command
    """
    cmd = '\ncheckm tree'
    cmd += ' --extension %s' % get_extension(self)
    cmd += ' --pplacer_threads %s' % self.soft.params['cpus']
    cmd += ' %s' % genomes_dir
    cmd += ' %s\n' % tree_dir
    return cmd


def tree(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm tree command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm tree
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    for genome, dirs in folders.items():
        genome_dir = dirs[0]
        out_dir = genome_out_dir(self, tech, genome_dir, group, genome)
        tree_dir = add_folder(self, 'checkm', out_dir, 'lineage_tree')

        self.outputs['dirs'].append(tree_dir)
        # outs = dirs + [tree_dir]
        outs = {genome: dirs + [tree_dir]}
        # self.outputs['outs'].setdefault((tech, group), []).extend(outs)
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        if self.soft.name != 'checkm' and to_do(folder=genome_dir):
            self.soft.add_status(
                tech, self.sam_pool, [genome_dir], group=group, genome=genome)

        if self.config.force or glob.glob('%s/*' % tree_dir):
            cmd = tree_cmd(self, genome_dir, tree_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_d=genome_dir, o_d=tree_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        elif self.soft.name == 'checkm':
            io_update(self, i_d=tree_dir, key=key)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def treeqa_cmd(
        self,
        tree_dir: str,
        out_dir: str
) -> str:
    """Collect the command line for checkm tree_qa.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    tree_dir : str
        Path to the input folder
    out_dir : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm tree_qa command
    """
    cmd = ''
    for (filename, out_format) in [
        ('tree_placement_summary.txt', '1'),
        ('tree_placement_lineage_stats.txt', '2'),
        ('genome_tree_IMG_genome_ids.nwk', '3'),
        ('genome_tree_taxonomy_strings.nwk', '4'),
        ('reference_genomes_and_bins.ali', '5')
    ]:
        tree_fpo = '%s/%s_%s' % (out_dir, out_format, filename)
        if self.config.force or to_do(tree_fpo):
            cmd += '\ncheckm tree_qa'
            cmd += ' --out_format %s' % out_format
            cmd += ' --file %s' % tree_fpo
            if self.soft.params['ali']:
                cmd += ' --ali'
            if self.soft.params['nt']:
                cmd += ' --nt'
            if self.soft.params['tab_table']:
                cmd += ' --tab_table'
            cmd += ' %s\n' % tree_dir
    return cmd


def treeqa(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm tree_qa command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm tree_qa
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    for genome, dirs in folders.items():
        tree_dir = dirs[-1]
        out_dir = genome_out_dir(self, tech, tree_dir, group, genome)
        qa_dir = out_dir
        if self.soft.name == 'checkm':
            qa_dir = add_folder(self, 'checkm', out_dir, 'lineage_tree_qa')
            tree_dir = add_folder(self, 'checkm', out_dir, 'lineage_tree')
        elif self.soft.prev != 'checkm_tree':
            sys.exit('[%s] Run "checkm_tree" first' % self.soft.name)

        self.outputs['dirs'].append(qa_dir)
        # outs = dirs + [qa_dir]
        outs = {genome: dirs + [qa_dir]}
        # self.outputs['outs'].setdefault((tech, group), []).extend(outs)
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        if self.soft.name != 'checkm' and to_do(folder=tree_dir):
            self.soft.add_status(
                tech, self.sam_pool, [tree_dir], group=group, genome=genome)

        out = '%s/1_tree_placement_summary.txt' % qa_dir
        if self.config.force or to_do(out):
            cmd = treeqa_cmd(self, tree_dir, qa_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_d=tree_dir, o_d=qa_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        elif self.soft.name == 'checkm':
            io_update(self, i_d=qa_dir, key=key)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def lineageset_cmd(
        self,
        tree_dir: str,
        lineage_ms: str
) -> str:
    """Collect the command line for checkm lineage_set.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    tree_dir : str
        Path to the input folder
    lineage_ms : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm lineage_set command
    """
    cmd = '\ncheckm lineage_set'
    cmd += ' --unique %s' % self.soft.params['unique']
    cmd += ' --multi %s' % self.soft.params['multi']
    if self.soft.params['force_domain']:
        cmd += ' --force_domain'
    if self.soft.params['no_refinement']:
        cmd += ' --no_refinement'
    cmd += ' %s' % tree_dir
    cmd += ' %s\n' % lineage_ms
    return cmd


def lineageset(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm lineage_set command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm lineage_set
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    for genome, dirs in folders.items():

        tree_dir = dirs[-1]
        out_dir = genome_out_dir(self, tech, tree_dir, group, genome)
        lineage_dir = out_dir
        if self.soft.name == 'checkm':
            tree_dir = add_folder(self, 'checkm', out_dir, 'lineage_tree')
            lineage_dir = add_folder(self, 'checkm', out_dir, 'lineage_set')
        elif self.soft.prev != 'checkm_tree':
            sys.exit('[%s] Run "checkm_tree" first' % self.soft.name)

        self.outputs['dirs'].append(lineage_dir)
        # outs = dirs + [lineage_dir]
        outs = {genome: dirs + [lineage_dir]}
        # self.outputs['outs'].setdefault((tech, group), []).extend(outs)
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        lineage_ms = '%s/lineage.ms' % lineage_dir
        if self.soft.name != 'checkm' and to_do(lineage_ms):
            self.soft.add_status(
                tech, self.sam_pool, [lineage_ms], group=group, genome=genome)

        if self.config.force or to_do(lineage_ms):
            cmd = lineageset_cmd(self, tree_dir, lineage_ms)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_d=tree_dir, o_d=lineage_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        elif self.soft.name == 'checkm':
            io_update(self, i_d=lineage_dir, key=key)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def analyze_cmd(
        self,
        genomes_dir: str,
        lineage_dir: str,
        analyze_dir: str
) -> str:
    """Collect the command line for checkm analyze.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    genomes_dir : str
        Path to the input folder
    lineage_dir : str
        Path to the folder resulting from lineage_set step
    analyze_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        checkm analyze command
    """
    cmd = '\ncheckm analyze'
    cmd += ' --extension %s' % get_extension(self)
    cmd += ' --threads %s' % self.soft.params['cpus']
    if self.soft.params['ali']:
        cmd += ' --ali'
    if self.soft.params['nt']:
        cmd += ' --nt'
    cmd += ' %s/lineage.ms' % lineage_dir
    cmd += ' %s' % genomes_dir
    cmd += ' %s\n' % analyze_dir
    return cmd


def analyze(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm analyze command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm analzye
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    for genome, dirs in folders.items():

        genomes_dir = dirs[0]
        analyze_dir = genome_out_dir(self, tech, dirs[-1], group, genome)
        if self.soft.name == 'checkm':
            analyze_dir = add_folder(
                self, 'checkm', analyze_dir, 'lineage_analyze')
            lineage_dir = add_folder(
                self, 'checkm', analyze_dir, 'lineage_set')
        else:
            if self.soft.prev != 'checkm_lineageset':
                sys.exit('[%s] Run "checkm_lineageset" first' % self.soft.name)
            lineage_dir = dirs[-1]

        self.outputs['dirs'].append(analyze_dir)
        # outs = dirs + [analyze_dir]
        outs = {genome: dirs + [analyze_dir]}
        # self.outputs['outs'].setdefault((tech, group), []).extend(outs)
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        out = '%s/lineage.ms' % lineage_dir
        if self.soft.name != 'checkm' and to_do(out):
            self.soft.add_status(
                tech, self.sam_pool, [out], group=group, genome=genome)

        if self.config.force or glob.glob('%s/*' % analyze_dir):
            cmd = analyze_cmd(self, genomes_dir, lineage_dir, analyze_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(
                self, i_d=[genomes_dir, lineage_dir], o_d=analyze_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        elif self.soft.name == 'checkm':
            io_update(self, i_d=analyze_dir, key=key)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def coverage_cmd(
        self,
        genome_dir: str,
        cov: str,
        bams: list
) -> str:
    """Collect the command line for checkm coverage.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    genome_dir : str
        Path to the input folder
    cov : str
        Path to the coverage output file
    bams : list
        Path to the input BAM file

    Returns
    -------
    cmd : str
        checkm coverage command
    """
    cmd = '\ncheckm coverage'
    cmd += ' --extension %s' % get_extension(self)
    if self.soft.params['all_reads']:
        cmd += ' --all_reads'
    for param in ['min_align', 'max_edit_dist', 'min_qc']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' %s' % genome_dir
    cmd += ' %s' % cov
    cmd += ' %s\n' % ' '.join(bams)
    return cmd


def get_bams(
        self,
        tech: str,
        group: str
) -> list:
    """Get the BAMs files for the assembler used before binning/drep.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm coverage
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    group : str
        Mame of a co-assembly pool's group

    Returns
    -------
    bams : list
        Path(s) to BAMs files
    """
    assembl = get_assembler(self)
    map_assembl = 'mapping_%s' % assembl
    if map_assembl in self.softs:
        bams = self.softs[map_assembl].outputs[self.sam_pool][(tech, group)]
        print(bams)
        print(bamsfds)
    elif self.config.dev:
        bams = ['/path/to/mapping1.bam', '/path/to/mapping2.bam']
    else:
        print('Run "%s %s" for checkm coverage' % (assembl, map_assembl))
        return []
    return bams


def coverage(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm coverage command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm coverage
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    if self.soft.params['coverage']:
        bams = get_bams(self, tech, group)
        if bams:
            for genome, dirs in folders.items():
                genome_dir = dirs[0]
                coverage_dir = genome_out_dir(self, tech, genome_dir,
                                              group, genome)
                if self.soft.name == 'checkm':
                    coverage_dir = add_folder(
                        self, 'checkm', coverage_dir, 'coverage')

                self.outputs['dirs'].append(coverage_dir)
                cov = '%s/coverage.tsv' % coverage_dir
                outs = {genome: cov}
                self.outputs['outs'].setdefault((tech, group), {}).update(outs)

                key = genome_key(tech, group, genome)
                if self.soft.name != 'checkm' and to_do(cov):
                    self.soft.add_status(
                        tech, self.sam_pool, [cov], group=group, genome=genome)

                if self.config.force or to_do(cov):
                    cmd = coverage_cmd(self, genome_dir, cov, bams)
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, i_f=bams, i_d=genome_dir, o_f=cov, key=key)
                    self.soft.add_status(
                        tech, self.sam_pool, 1, group=group, genome=genome)
                elif self.soft.name == 'checkm':
                    io_update(self, i_f=cov, key=key)
                else:
                    self.soft.add_status(
                        tech, self.sam_pool, 0, group=group, genome=genome)


def qa_cmd(
        self,
        lineage_ms: str,
        analyze_dir: str,
        qa_dir: str,
        cov: str
) -> str:
    """Collect the command line for checkm qa.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    lineage_ms : str
        Path to the lineage.ms file resulting from lineage_set step
    analyze_dir : str
        Path to the folder resulting from analyze step
    qa_dir : str
        Path to the output folder
    cov : str
        Path to the coverage output file

    Returns
    -------
    cmd : str
        checkm qa command
    """
    cmd = ''
    for idx, (filename, out_format) in enumerate([
        ('completeness_contamination.txt', 1),
        ('bin_statistics.txt', 2),
        ('bin_quality.txt', 3),
        ('marker_genes_and_counts.txt', 4),
        ('bin_id_marker_gene_id_gene_id.txt', 5),
        ('marker_genes_present_multiple_times_in_a_bin.txt', 6),
        ('marker_genes_present_multiple_times_on_same_scaffold.txt', 7),
        ('marker_gene_positions.txt', 8),
        ('marker_genes.fas', 9)
    ]):
        out = '%s/%s_qa_%s' % (qa_dir, out_format, filename)
        if self.config.force or to_do(out):
            cmd += '\ncheckm qa'
            cmd += ' --out_format %s' % out_format
            cmd += ' --file %s' % out
            cmd += ' --threads %s' % self.soft.params['cpus']
            if self.soft.params['tab_table']:
                cmd += ' --tab_table'
            if self.soft.params['individual_markers']:
                cmd += ' --individual_markers'
            if self.soft.params['skip_adj_correction']:
                cmd += ' --skip_adj_correction'
            if self.soft.params['skip_pseudogene_correction']:
                cmd += ' --skip_pseudogene_correction'
            if self.soft.params['ignore_thresholds']:
                cmd += ' --ignore_thresholds'
            if self.soft.params['coverage'] and cov:
                cmd += ' --coverage_file %s' % cov
            cmd += ' --aai_strain %s' % self.soft.params['aai_strain']
            cmd += ' --e_value %s' % self.soft.params['e_value']
            cmd += ' --length %s' % self.soft.params['length']
            if not idx:
                cmd += ' --alignment_file %s/multi_copy_genes_AAI.ali' % qa_dir
            cmd += ' %s' % lineage_ms
            cmd += ' %s\n' % analyze_dir
    return cmd


def get_cov(
        self,
        tech: str,
        group: str,
        genome: str
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm qa
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    group : str
        Mame of a co-assembly pool's group
    genome : str
        Name of the current genomes group

    Returns
    -------
    cov : str
        Path to the coverage output file
    """
    cov = None
    if self.soft.params['coverage'] and 'checkm_coverage' in self.softs:
        cov_soft = self.softs['checkm_coverage']
        run_on = (cov_soft.prev not in self.graph.paths[self.soft.name][0])
        if run_on or self.sam_pool not in cov_soft.outputs:
            sys.exit('[%s] Run on "%s" genomes to use coverage"' % (
                self.soft.name, cov_soft.prev))
        if (tech, group) in cov_soft.outputs[self.sam_pool]:
            cov = cov_soft.outputs[self.sam_pool][(tech, group)].get(genome)
    return cov


def qa(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm qa command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm qa
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    for genome, dirs in folders.items():

        analyze_dir = dirs[-1]
        qa_dir = genome_out_dir(self, tech, analyze_dir, group, genome)
        if self.soft.name == 'checkm':
            cov = '%s/coverage.tsv' % add_folder(
                self, 'checkm', analyze_dir, 'coverage')
            analyze_dir = add_folder(self, 'checkm', qa_dir, 'lineage_analyze')
            lineage_dir = add_folder(self, 'checkm', qa_dir, 'lineage_set')
            qa_dir = add_folder(self, 'checkm', qa_dir, 'lineage_qa')
        else:
            if self.soft.prev != 'checkm_analyze':
                sys.exit('[%s] Run "checkm_analyze" first' % self.soft.name)
            lineage_dir = dirs[1]
            cov = get_cov(self, tech, group, genome)

        self.outputs['dirs'].append(qa_dir)
        outs = {genome: [qa_dir]}
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        lineage_ms = '%s/lineage.ms' % lineage_dir
        if self.soft.name != 'checkm' and to_do(lineage_ms):
            self.soft.add_status(
                tech, self.sam_pool, [lineage_ms], group=group, genome=genome)

        out = '%s/1_completeness_contamination.txt' % qa_dir
        if self.config.force or to_do(out):
            cmd = qa_cmd(self, lineage_ms, analyze_dir, qa_dir, cov)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=[lineage_ms, cov], i_d=analyze_dir,
                      o_d=qa_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def unbinned_cmd(
        self,
        genomes_dir: str,
        contigs: str,
        unbinned_dir: str
) -> str:
    """Collect the command line for checkm unbinned.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    genomes_dir : str
        Path to the input folder
    contigs : str
        Path to the input file
    unbinned_dir : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm unbinned command
    """
    cmd = '\ncheckm unbinned'
    cmd += ' --extension %s' % get_extension(self)
    cmd += ' --min_seq_len %s' % self.soft.params['min_seq_len']
    cmd += ' %s' % genomes_dir
    cmd += ' %s' % contigs
    cmd += ' %s/unbinned.fa' % unbinned_dir
    cmd += ' %s/unbinned_stats.tsv\n' % unbinned_dir
    return cmd


def unbinned(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Get the checkm unbinned command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for checkm unbinned
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folder(s) containing the genome/MAGs
    group : str
        Mame of a co-assembly pool's group
    """
    if self.soft.name == 'checkm_unbinned' and self.sam_pool == '':
        sys.exit('[checkm_unbinned] Only per co-assembly (not after drep)')

    assembler = get_assembler(self)
    contigs = self.softs[assembler].outputs[self.sam_pool][(tech, group)][0]

    for genome, dirs in folders.items():
        genomes_dir = dirs[0]
        unbinned_dir = genome_out_dir(self, tech, genomes_dir, group, genome)
        if self.soft.name == 'checkm':
            unbinned_dir = add_folder(self, 'checkm', unbinned_dir, 'unbinned')

        self.outputs['dirs'].append(unbinned_dir)
        outs = {genome: unbinned_dir}
        self.outputs['outs'].setdefault((tech, group), {}).update(outs)

        key = genome_key(tech, group, genome)
        if self.soft.name != 'checkm' and not glob.glob('%s/*' % genomes_dir):
            self.soft.add_status(
                tech, self.sam_pool, [genomes_dir], group=group, genome=genome)

        if self.config.force or to_do('%s/unbinned.fa' % unbinned_dir):
            cmd = unbinned_cmd(self, genomes_dir, contigs, unbinned_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(
                self, i_f=contigs, i_d=genomes_dir, o_d=unbinned_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def tetra_cmd(
        self,
        fasta: str,
        out: str
) -> str:
    """Collect the command line for checkm tetra.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    fasta : str
        Path to the input fasta file
    out : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm tetra command
    """
    cmd = '\ncheckm tetra'
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' %s' % fasta
    cmd += ' %s\n' % out
    return cmd


def tetra(
        self,
        tech: str,
        fastas: dict,
        group: str
) -> None:
    """Collect the command to calculate tetranucleotide signature of
    the contig sequences using checkm tetra.

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
        Technology or technologies used for assembly
    fastas : dict
        Path(s) to the input fasta file(s) per genome/MAGs
    group : str
        Name of a co-assembly pool's group
    """
    for fasta in fastas.values():

        out_dir = genome_out_dir(self, tech, fasta[0], group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        key = genome_key(tech, group)
        if to_do(fasta[0]):
            self.soft.add_status(tech, self.sam_pool, [fasta[0]], group=group)

        out = '%s/tetra.txt' % out_dir
        if self.config.force or to_do(out):
            cmd = tetra_cmd(self, fasta[0], out)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta[0], o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def checkm_(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Perform all the key steps of the checkm analysis.

    Notes
    -----
    This is done if the pipeline specifies "drep checkm" and not specific
    modules of checkm, such as "drep checkm_qa" or "drep checkm_tree"

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
        Technology or technologies used for assembly
    folders : dict
        Paths to the input folders with the genomes/MAGs
    group : str
        Name of a co-assembly pool's group
    """
    for step in [tree, treeqa, lineageset, analyze, coverage, qa]:
        step(self, tech, folders, group)
    if self.sam_pool != '':
        unbinned(self, tech, folders, group)


def checkm(self) -> None:
    """CheckM provides a set of softwares for assessing the quality of genomes
    recovered from isolates, single cells, or metagenomes. It provides robust
    estimates of genome completeness and contamination by using collocated sets
    of genes that are ubiquitous and single-copy within a phylogenetic lineage.
    Assessment of genome quality can also be examined using plots depicting key
    genomic characteristics (e.g., GC, coding density) which highlight sequences
    outside the expected distributions of a typical genome. CheckM also provides
    softwares for identifying genome bins that are likely candidates for merging
    based on marker set compatibility, similarity in genomic characteristics,
    and proximity within a reference genome tree.

    References
    ----------
    Parks, Donovan H., et al. "CheckM: assessing the quality of microbial
    genomes recovered from isolates, single cells, and metagenomes." Genome
    research 25.7 (2015): 1043-1055.

    Notes
    -----
    GitHub  : https://github.com/Ecogenomics/CheckM
    Website : ecogenomics.github.io/checkm
    Docs    : https://github.com/Ecogenomics/CheckM/wiki
    Paper   : https://doi.org/10.1101/gr.186072.114

    Parameters
    ----------
    self : Commands class instance
        .soft.name : str
            Current software in the pipeline
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
    """
    __module_call__ = caller(self, __name__)
    assemblers = self.config.tools['assembling']
    if self.soft.name == 'checkm_tetra' and self.soft.prev not in assemblers:
        sys.exit('[checkm_tetra] Only after assembly: %s' % ''.join(assemblers))

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs, True)
            __module_call__(self, tech, folders, group)

    elif set(self.inputs) == {''}:
        for (tech, bin_algo), inputs in self.inputs[''].items():
            folders = group_inputs(self, inputs, True)
            __module_call__(self, tech, folders, bin_algo)


def checkm2_cmd(
        self,
        folder: str,
        out_dir: str
) -> str:
    """Collect the command line for checkm2 predict.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    folder : str
        Path to the input folder
    out_dir : str
        Path to the output file

    Returns
    -------
    cmd : str
        checkm2 predict command
    """
    cmd = '\ncheckm2 predict'
    cmd += ' --input %s' % folder
    cmd += ' --output-directory %s' % out_dir
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['lowmem', 'general', 'specific', 'allmodels',
                    'genes', 'force', 'dbg_cos', 'dbg_vectors']:
        if self.soft.params[boolean]:
            cmd += ' --%s %s' % (boolean, self.soft.params[boolean])
    return cmd


def get_checkm2(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """Collect the checkm2 predict command.

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for checkm tetra
        .pool : str
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
        Technology or technologies used for assembly
    folders : dict
        Path(s) to the input folders with the genome/MAGs
    group : str
        Name of a co-assembly pool's group
    """
    for genome, dirs in folders.items():

        folder = dirs[0]
        out_dir = genome_out_dir(self, tech, folder, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        key = genome_key(tech, group, genome)
        if to_do(folder=folder):
            self.soft.add_status(
                tech, self.sam_pool, [folder], group=group, genome=genome)

        out = '%s/quality_report.tsv' % out_dir
        if self.config.force or to_do(out):
            cmd = checkm2_cmd(self, folder, out_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=folder, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def checkm2(self) -> None:
    """Rapid assessment of genome bin quality using machine learning.

    Unlike CheckM1, CheckM2 has universally trained machine learning models
    it applies regardless of taxonomic lineage to predict the completeness
    and contamination of genomic bins. This allows it to incorporate many
    lineages in its training set that have few - or even just one -
    high-quality genomic representatives, by putting it in the context of all
    other organisms in the training set. As a result of this machine learning
    framework, CheckM2 is also highly accurate on organisms with reduced genomes
    or unusual biology, such as the Nanoarchaeota or Patescibacteria.

    CheckM2 uses two distinct machine learning models to predict genome
    completeness. The 'general' gradient boost model is able to generalize
    well and is intended to be used on organisms not well represented in
    GenBank or RefSeq (roughly, when an organism is novel at the level of
    order, class or phylum). The 'specific' neural network model is more
    accurate when predicting completeness of organisms more closely related
    to the reference training set (roughly, when an organism belongs to a
    known species, genus or family). CheckM2 uses a cosine similarity
    calculation to automatically determine the appropriate completeness model
    for each input genome, but you can also force the use of a particular
    completeness model, or get the prediction outputs for both. There is only
    one contamination model (based on gradient boost) which is applied
    regardless of taxonomic novelty and works well across all cases.

    References
    ----------
    Chklovski, Alex, et al. "CheckM2: a rapid, scalable and accurate tool for
    assessing microbial genome quality using machine learning." bioRxiv (2022).

    Notes
    -----
    GitHub  : https://github.com/chklovski/CheckM2
    Paper   : https://doi.org/10.1101/2022.07.11.499243

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Checkm2
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs, True)
            get_checkm2(self, tech, folders, group)

    elif set(self.inputs) == {''}:
        for (tech, bin_algo), inputs in self.inputs[''].items():
            folders = group_inputs(self, inputs, True)
            get_checkm2(self, tech, folders, bin_algo)
