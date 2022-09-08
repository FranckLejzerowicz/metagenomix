# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from metagenomix._io_utils import caller, io_update, to_do
from metagenomix._inputs import (group_inputs, genome_key, genome_out_dir,
                                 get_extension, get_reads, add_folder)


def lorikeet_cmd(
        self,
        genome: str,
        fasta_folder: str,
        group_reads: dict,
        out_dir: str,
        key: tuple,
        step: str
) -> str:
    """Collect lorikeet call, lorikeet consensus, or lorikeet genotype command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    genome : str
        MAGs/Genomes folder name or empty string (for assembly contigs)
    fasta_folder : str
        Path to aseembly fasta or to an input folder with the genomes/MAGs
    group_reads : dict
        Path(s) to the fastq files per sample
    out_dir : str
        Path to the output folder
    key : tuple
        Turrent tech, pool, group, genome, etc
    step : str
        name of the Lorikeet module ("call", "consensus", or "genotype")

    Returns
    -------
    cmd : str
        lorikeet command
    """
    tmp_dir = '$TMPDIR/lorikeet_%s_%s' % ('_'.join(key), step)
    cmd = 'mkdir -p %s\n' % tmp_dir

    cmd += 'lorikeet %s' % step
    cmd += ' --output-directory %s' % out_dir
    cmd += ' --output-prefix %s_output' % step
    cmd += ' --bam-file-cache-directory %s' % tmp_dir
    cmd += ' --threads %s' % self.soft.params['cpus']
    if genome:
        cmd += ' --reference %s' % fasta_folder
    else:
        cmd += ' --genome-fasta-directory %s' % fasta_folder
        cmd += ' --genome-fasta-extension %s' % get_extension(self)

    cmd += ' --method %s' % ' '.join(self.soft.params['method'])
    cmd_single, cmd_coupled, cmd_longreads = '', '', ''
    for sam, tech_sam_fastqs in group_reads.items():
        for tech_sam, fastqs in tech_sam_fastqs.items():
            if tech_sam[0] == 'illumina':
                if len(fastqs) == 1:
                    cmd_single += ' %s' % fastqs[0]
                elif len(fastqs) == 2:
                    cmd_coupled += ' %s' % ' '.join(fastqs)
            elif fastqs:
                cmd_longreads += ' %s' % fastqs[0]

    if cmd_single:
        cmd += ' --single%s' % cmd_single
    if cmd_coupled:
        cmd += ' --coupled%s' % cmd_coupled
    if cmd_longreads:
        cmd += ' --longreads%s' % cmd_longreads

    for param in [
        'min_read_aligned_length', 'min_read_percent_identity',
        'min_read_aligned_percent', 'min_read_aligned_length_pair',
        'min_read_percent_identity_pair', 'min_read_aligned_percent_pair',
        'min_covered_fraction', 'contig_end_exclusion', 'trim_min', 'trim_max',
        'mapper', 'longread_mapper', 'kmer_sizes', 'ploidy',
        'qual_by_depth_filter', 'qual_threshold', 'depth_per_sample_filter',
        'min_base_quality', 'min_mapq', 'base_quality_score_threshold',
        'max_input_depth', 'min_contig_size', 'min_sv_qual',
        'phred_scaled_global_read_mismapping_rate',
        'pair_hmm_gap_continuation_penalty',
        'heterozygosity', 'heterozygosity_stdev', 'indel_heterozygosity',
        'standard_min_confidence_threshold_for_calling',
        'active_probability_threshold', 'min_assembly_region_size',
        'max_assembly_region_size', 'assembly_region_padding',
        'min_dangling_branch_length', 'min_prune_factor', 'num_pruning_samples',
        'initial_error_rate_for_pruning', 'pruning_log_odds_threshold',
        'max_unpruned_variants', 'max_prob_propagation_distance',
        'max_mnp_distance', 'parallel_genomes'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])
    for boolean in [
        'discard_improper_pairs', 'discard_supplementary', 'include_secondary',
        'discard_unmapped', 'high_memory', 'splitbams',
        'minimap2_reference_is_index', 'calculate_fst', 'calculate_dnds',
        'features_vcf', 'do_not_call_svs', 'use_posteriors_to_calculate_qual',
        'annotate_with_num_discovered_alleles',
        'dont_increase_kmer_sizes_for_cycles', 'allow_non_unique_kmers_in_ref',
        'do_not_run_physical_phasing', 'recover_all_dangling_branches',
        'use_adaptive_pruning', 'graph_output', 'dont_use_soft_clipped_bases',
        'disable_optimizations', 'disable_avx'
    ]:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    for param in ['minimap2_params', 'bwa_params', 'ngmlr_params']:
        if self.soft.params[param]:
            cmd += ' --%s "%s"' % (
                param.replace('_', '-'), ' '.join(self.soft.params[param]))
    cmd += ' --force'
    return cmd


def call(
        self,
        tech: str,
        fastas_folders: dict,
        reads: dict,
        group: str
) -> None:
    """Call variants using local reassembly across multiple genomes and samples.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for lorikeet call
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
    fastas_folders : dict
        Path(s) to aseembly fasta or to the input folders with the genomes/MAGs
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    """
    get_lorikeet(self, tech, fastas_folders, reads, group, 'call')


def consensus(
        self,
        tech: str,
        fastas_folders: dict,
        reads: dict,
        group: str
) -> None:
    """Generate consensus genomes for each provided sample and genome.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for lorikeet consensus
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
    fastas_folders : dict
        Path(s) to aseembly fasta or to the input folders with the genomes/MAGs
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    """
    get_lorikeet(self, tech, fastas_folders, reads, group, 'consensus')


def genotype(
        self,
        tech: str,
        fastas_folders: dict,
        reads: dict,
        group: str
) -> None:
    """Resolves strain-level genotypes and abundance from metagenomes.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for lorikeet genotype
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
    fastas_folders : dict
        Path(s) to aseembly fasta or to the input folders with the genomes/MAGs
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    """
    get_lorikeet(self, tech, fastas_folders, reads, group, 'genotype')


def get_lorikeet(
        self,
        tech: str,
        fastas_folders: dict,
        reads: dict,
        group: str,
        step: str
) -> None:
    """Get the Lorikeet inputs and outputs for the current tech and
    MAGs/genomes and collect the commands and data structures.

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
    fastas_folders : dict
        Path(s) to aseembly fasta or to the input folders with the genomes/MAGs
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    step : str
        Name of the Lorikeet step
    """
    for genome, fasta_folder_ in fastas_folders.items():

        fasta_folder = fasta_folder_[0]
        out_dir = genome_out_dir(self, tech, fasta_folder, group, genome)
        if self.soft.name == 'lorikeet':
            out_dir = add_folder(self, 'lorikeet', out_dir, step)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        key = genome_key(tech, group, genome)
        if genome:
            condition = to_do(folder=fasta_folder)
        else:
            condition = to_do(fasta_folder)
        if condition:
            self.soft.add_status(
                tech, self.sam_pool, [fasta_folder], group=group, genome=genome)

        if self.config.force or not glob.glob('%s/*' % out_dir):

            if self.sam_pool:
                groups = self.pools[self.sam_pool][group]
                group_reads = dict(x for x in reads.items() if x[0] in groups)
            else:
                group_reads = reads

            cmd = lorikeet_cmd(
                self, genome, fasta_folder, group_reads, out_dir, key, step)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            if genome:
                io_update(self, i_d=fasta_folder, o_d=out_dir, key=key)
            else:
                io_update(self, i_f=fasta_folder, o_d=out_dir, key=key)


def lorikeet_(
        self,
        tech: str,
        fastas_folders: dict,
        reads: dict,
        group: str
) -> None:
    """Perform all the key steps of the lorikeet analysis.

    Notes
    -----
    This is done if the pipeline specifies "drep lorikeet" and not specific
    modules of checkm, such as "drep lorikeet_call" or "drep lorikeet_evolve"

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
    fastas_folders : dict
        Path(s) to aseembly fasta or to the input folders with the genomes/MAGs
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    """
    for step in ['call', 'consensus', 'genotype']:
        get_lorikeet(self, tech, fastas_folders, reads, group, step)


def lorikeet(self) -> None:
    """Lorikeet aims to call variants in metagenomes using local reassembly
    of haplotypes. Within-species variant analysis pipeline for metagenomic
    communities that utilizes both long and short reads. Lorikeet utilizes a
    re-implementaion of the GATK HaplotypeCaller algorithm, performing local
    re-assembly of potentially active regions within candidate genomes.
    Called variants can be clustered into likely strains using a combination
    of UMAP and HDBSCAN. Additional statistics, like consensus ANI,
    population ANI, and subpopulation ANI will also be calculated for each
    input geome providing values for each sample compared to the reference
    and also compared to all other samples.

    Lorikeet has a variety of subcommands with the main being call and
    genotype. The call pipeline will take any number of input genomes and
    samples and perform robust variant calling and ANI calculations. The
    genotype algorithm takes this a step further and attempts to reconstruct
    strain haplotypes from the called variants and return complete strain
    genomes.

    Input can either be reads and reference genome, or MAG. Or a BAM file
    and associated genome.

    References
    ----------
    In prep

    Notes
    -----
    GitHub  : https://github.com/rhysnewell/Lorikeet
    Docs    : https://rhysnewell.github.io/Lorikeet/
    Paper   : In prep

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
    reads = get_reads(self)

    __module_call__ = caller(self, __name__)
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas_folders = group_inputs(self, inputs, True)
            __module_call__(self, tech, fastas_folders, reads, group)

    elif set(self.inputs) == {''}:
        for (tech, bin_algo), inputs in self.inputs[''].items():
            folders = group_inputs(self, inputs, True)
            __module_call__(self, tech, folders, reads, bin_algo)
