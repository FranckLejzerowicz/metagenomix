# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import dirname
from metagenomix._io_utils import (
    caller, io_update, tech_specificity, status_update, get_assembly,
    to_do, get_assembly_contigs)
from metagenomix._inputs import (
    group_inputs, genome_key, genome_out_dir, get_extension, get_reads,
    get_group_reads, add_folder)


def lorikeet_cmd(
        self,
        is_folder: bool,
        contigs: list,
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
    is_folder : bool
        Whether fasta_folder is a folder or not
    contigs : list
        Path(s) to the input contigs fasta file(s)
    fasta_folder : str
        Path to assembly fasta or to an input folder with the genomes/MAGs
    group_reads : dict
        Path(s) to the fastq files per sample
    out_dir : str
        Path to the output folder
    key : tuple
        Current tech, pool, group, genome, etc
    step : str
        name of the Lorikeet module ("call", "consensus", or "genotype")

    Returns
    -------
    cmd : str
        lorikeet command
    """
    tmp_dir = '$TMPDIR/lorikeet_%s_%s' % ('_'.join(key), step)
    cmd_rm, cmd = '', 'mkdir -p %s\n' % tmp_dir
    for contig in contigs:
        if contig.endswith('.gz'):
            io_update(self, i_f=contig, key=key)
            cmd += 'gunzip -c %s > %s\n' % (contig, contig.rstrip('.gz'))
            cmd_rm += 'rm %s\n' % contig.rstrip('.gz')

    cmd += 'lorikeet %s' % step
    cmd += ' --output-directory %s' % out_dir
    # cmd += ' --graph-output %s/assembly_graph_debug.txt' % out_dir
    cmd += ' --bam-file-cache-directory %s' % tmp_dir
    if is_folder:
        cmd += ' --genome-fasta-directory %s' % fasta_folder
        cmd += ' --genome-fasta-extension %s' % get_extension(self)
    else:
        cmd += ' --reference %s' % ' '.join([c.rstrip('.gz') for c in contigs])

    if self.soft.params['features_vcf']:
        cmd += ' --features-vcf %s' % self.soft.params['features_vcf']

    cmd_1, cmd_2, cmd_single, cmd_coupled, cmd_longreads = '', '', '', '', ''
    for sam, tech_sam_fastqs in group_reads.items():
        for tech_sam, fastqs in tech_sam_fastqs.items():
            io_update(self, i_f=fastqs, key=key)
            if tech_sam[0] == 'illumina':
                if len(fastqs) == 1:
                    cmd_single += ' %s' % fastqs[0]
                elif len(fastqs) >= 2:
                    # cmd_1 += ' %s' % fastqs[0]
                    # cmd_2 += ' %s' % fastqs[1]
                    cmd_coupled += ' %s' % ' '.join(fastqs[:2])
                    if len(fastqs) == 3:
                        cmd_single += ' %s' % fastqs[-1]
            elif fastqs:
                cmd_longreads += ' %s' % fastqs[0]

    if cmd_1:
        cmd += ' -1%s' % cmd_1
    if cmd_2:
        cmd += ' -2%s' % cmd_2
    if cmd_single:
        cmd += ' --single%s' % cmd_single
    if cmd_coupled:
        cmd += ' --coupled%s' % cmd_coupled
    if cmd_longreads:
        cmd += ' --longreads%s' % cmd_longreads

    for pdx, param_ in enumerate(['min_read_aligned_length',
                                  'min_read_percent_identity',
                                  'min_read_aligned_percent']):
        param = param_.replace('_', '-')
        if (cmd_1 or cmd_coupled) and self.soft.params['proper_pairs_only']:
            if not pdx:
                cmd += ' --proper-pairs-only'
            param += '-pair'
        cmd += ' --%s %s' % (param, self.soft.params[param_])

    for param in [
        'mapper', 'longread_mapper', 'kmer_sizes', 'pcr_indel_model',
        'contig_end_exclusion', 'ploidy', 'qual_by_depth_filter',
        'qual_threshold', 'depth_per_sample_filter', 'min_long_read_size',
        'min_long_read_average_base_qual', 'min_base_quality', 'min_mapq',
        'base_quality_score_threshold', 'max_input_depth', 'min_contig_size',
        'min_sv_qual', 'phred_scaled_global_read_mismapping_rate',
        'pair_hmm_gap_continuation_penalty', 'min_assembly_region_size',
        'max_assembly_region_size', 'assembly_region_padding',
        'min_dangling_branch_length', 'min_prune_factor',
        'num_pruning_samples', 'max_unpruned_variants',
        'max_prob_propagation_distance', 'max_mnp_distance',
        # 'heterozygosity',
        'heterozygosity_stdev', 'indel_heterozygosity',
        'standard_min_confidence_threshold_for_calling',
        'initial_error_rate_for_pruning', 'trim_min', 'trim_max',
        'active_probability_threshold', 'pruning_log_odds_threshold',
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])

    for boolean in [
        'quiet', 'verbose', 'sharded', 'split_bams', 'calculate_fst',
        'calculate_dnds', 'do_not_call_svs', 'include_secondary',
        'exclude_supplementary', 'minimap2_reference_is_index',
        'use_posteriors_to_calculate_qual',
        'annotate_with_num_discovered_alleles',
        'dont_increase_kmer_sizes_for_cycles',
        'allow_non_unique_kmers_in_ref', 'recover_all_dangling_branches',
        'do_not_run_physical_phasing', 'disable_optimizations',
        'use_adaptive_pruning', 'discard_unmapped', 'disable_avx', 'force',
    ]:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    if self.soft.params['limiting_interval']:
        cmd += ' --limiting-interval %s' % self.soft.params['limiting_interval']

    if self.soft.params['dont_use_soft_clipped_bases']:
        if step in ['consensus', 'genotype']:
            cmd += ' --dont-use-soft-clipped-bases'

    for param in ['minimap2_params', 'ngmlr_params', 'bwa_params']:
        if self.soft.params[param]:
            cmd += ' --%s "%s"' % (
                param.replace('_', '-'), ' '.join(self.soft.params[param]))

    cmd += ' --parallel-genomes %s' % self.soft.params['cpus']
    cmd += ' --threads %s\n' % self.soft.params['cpus']
    cmd += cmd_rm
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
        group_reads: dict,
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
        Path(s) to assembly fasta or to the input folders with the genomes/MAGs
    group_reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    step : str
        Name of the Lorikeet step
    """
    for genome, fasta_folder_ in fastas_folders.items():

        fasta_folder = fasta_folder_[0]
        out_dir = genome_out_dir(self, tech, group, genome)
        if self.soft.name == 'lorikeet':
            out_dir = add_folder(self, 'lorikeet', out_dir, step)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        key = genome_key(tech, group, genome)
        if genome:
            is_folder = True
        else:
            is_folder = False
        to_dos = status_update(self, tech, [fasta_folder], folder=is_folder)

        assembly = get_assembly(self)
        contigs = get_assembly_contigs(self, tech, group, assembly)
        to_dos.extend(status_update(self, tech, contigs))

        if self.config.force or not glob.glob('%s/*' % out_dir):
            cmd = lorikeet_cmd(self, is_folder, contigs, fasta_folder,
                               group_reads, out_dir, key, step)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            if genome:
                io_update(self, i_d=fasta_folder, o_d=out_dir, key=key)
            else:
                io_update(self, i_f=fasta_folder, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


def lorikeet_(
        self,
        tech: str,
        fastas_folders: dict,
        group_reads: dict,
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
    group_reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    group : str
        Name of a co-assembly pool's group
    """
    for step in ['call', 'consensus', 'genotype']:
        get_lorikeet(self, tech, fastas_folders, group_reads, group, step)


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
    module_call = caller(self, __name__)
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas_folders = group_inputs(self, inputs, True)
            group_reads = get_group_reads(self, tech, group, reads)
            module_call(self, tech, fastas_folders, group_reads, group)


def get_sample_to_marker(
        self,
        tech,
        sam_dir: str,
        markers_dir: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
    sam_dir : str
    markers_dir : str
    """

    sam_dir_ = sam_dir.replace('${SCRATCH_FOLDER}', '')
    to_dos = status_update(self, tech, [sam_dir], folder=True)
    if self.config.force or len(glob.glob('%s/*.sam.bz2' % sam_dir_)):
        cmd = 'sample2markers.py'
        cmd += ' -i %s/*.sam.bz2' % sam_dir
        cmd += ' -o %s' % markers_dir
        cmd += ' -n %s' % self.soft.params['cpus']
        if to_dos:
            self.outputs['cmds'].setdefault((tech,), []).append(False)
        else:
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)


def extract_markers(
        self,
        tech,
        key,
        st: str,
        strain_dir: str,
        db_dir: str
) -> None:
    """Extract markers for the current species strains

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech
    key
    st : str
    strain_dir : str
    db_dir : str
    """
    to_dos = status_update(self, tech, [strain_dir], folder=True)
    if self.config.force or to_do('%s/%s.fna' % (strain_dir, st)):
        cmd = 'extract_markers.py'
        cmd += ' -c %s' % st
        cmd += ' -o %s' % db_dir
        if to_dos:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            self.outputs['cmds'].setdefault(key, []).append(cmd)


def get_strainphlan(
        self,
        tech: str,
        sam: str,
        inputs: dict
):
    bt2_fp, sam_fp = inputs[:2]
    sam_dir = dirname(sam_fp)
    db_dir = '%s/db_markers' % self.dir
    meta_dir = '%s/metadata' % self.dir
    markers_dir = '%s/consensus_markers' % self.dir
    wol_pd = self.databases.paths['wol']['taxonomy']

    self.outputs['outs'] = dict({})
    self.outputs['dirs'].extend([db_dir, meta_dir])
    get_sample_to_marker(self, tech, sam_dir, markers_dir)

    for strain_group, strains in self.config.strains.items():

        tmp = '$TMPDIR/strainphlan_%s' % strain_group
        key = (tech, strain_group)
        io_update(self, i_d=sam_dir, o_d=[markers_dir, db_dir, meta_dir],
                  key=key)

        strain_dir = '%s/%s' % (db_dir, strain_group)
        self.outputs['dirs'].append(strain_dir)
        for strain in [x.replace(' ', '_') for x in strains]:
            if strain[0] != 's':
                continue
            extract_markers(self, tech, key, strain, strain_dir, db_dir)
            odir = '%s/output/%s' % (self.dir, strain_group)
            self.outputs['dirs'].append(odir)
            self.outputs['outs'][(strain_group, strain)] = odir
            io_update(self, o_d=odir, key=key)
            tree = '%s/RAxML_bestTree.%s.StrainPhlAn3.tre' % (odir, strain)
            to_dos = status_update(self, tech, [markers_dir, db_dir],
                                   folder=True)
            if self.config.force or to_do(tree):
                # cmd = strainphlan_cmd(self, markers_dir, db_dir, strain,
                #                       wol_pd, key)
                cmd = 'strainphlan'
                cmd += ' -s %s/*.pkl' % markers_dir
                cmd += ' -m %s/%s.fna' % (db_dir, strain)
                if strain in wol_pd['species']:
                    fna = '%s/%s.fna.bz2' % (
                        self.databases.paths['wol']['fna'],
                        wol_pd.loc[wol_pd['species'] == strain, 0])
                    io_update(self, i_f=fna, key=key)
                    cmd += ' -r %s' % fna
                cmd += ' -o %s' % odir
                cmd += ' -n %s' % self.soft.params['cpus']
                cmd += ' -c %s' % strain
                cmd += ' --mutation_rates'
                cmd = 'mkdir -p %s\n' % tmp + cmd
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                self.soft.add_status(
                    tech, sam, 1, group=strain_group, genome=strain)
            else:
                self.soft.add_status(
                    tech, sam, 0, group=strain_group, genome=strain)


def strainphlan(self) -> None:
    """StrainPhlAn is a computational tool for tracking individual strains
    across a large set of samples. The input of StrainPhlAn is a set of
    metagenomic samples and for each species, the output is a multiple
    sequence alignment (MSA) file of all species strains reconstructed
    directly from the samples. From this MSA, StrainPhlAn calls (PhyloPhlAn
    3)[http://segatalab.cibio.unitn.it/tools/phylophlan3/index.html] to build
    the phylogenetic tree showing the strain evolution of the sample strains.

    References
    ----------
    Beghini, Francesco, et al. "Integrating taxonomic, functional,
    and strain-level profiling of diverse microbial communities with
    bioBakery 3." Elife 10 (2021): e65088.

    Notes
    -----
    Docs    : https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-3
    Website : http://segatalab.cibio.unitn.it/tools/phylophlan3/index.html
    Paper   : https://doi.org/10.7554/eLife.65088

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .databases
            All databases
        .config
            Configurations
    """
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, 'illumina', sam):
            continue
        get_strainphlan(self, tech, sam, inputs)


def instrain(self):
    """InStrain is a tool for analysis of co-occurring genome populations
    from metagenomes that allows highly accurate genome comparisons, analysis
    of coverage, microdiversity, and linkage, and sensitive SNP detection
    with gene localization and synonymous non-synonymous identification.

    References
    ----------
    Olm, M.R., Crits-Christoph, A., Bouma-Gregson, K., Firek, B.A., Morowitz,
    M.J. and Banfield, J.F., 2021. inStrain profiles population
    microdiversity from metagenomic data and sensitively detects shared
    microbial strains. Nature Biotechnology, 39(6), pp.727-736.

    Notes
    -----
    GitHub  : https://github.com/MrOlm/instrain
    Paper   : https://doi.org/10.1038/s41587-020-00797-0
    Docs    : https://instrain.readthedocs.io/en/latest/index.html

    Parameters
    ----------
    self
    """
    pass


def panphlan(self):
    """PanPhlAn is a strain-level metagenomic profiling tool for identifying
    the gene composition of individual strains in metagenomic samples.
    PanPhlAn’s ability for strain-tracking and functional analysis of unknown
    pathogens makes it an efficient tool for culture-free microbial
    population studies.

    References
    ----------
    Beghini, F., McIver, L.J., Blanco-Míguez, A., Dubois, L., Asnicar, F.,
    Maharjan, S., Mailyan, A., Manghi, P., Scholz, M., Thomas, A.M. and
    Valles-Colomer, M., 2021. Integrating taxonomic, functional,
    and strain-level profiling of diverse microbial communities with
    bioBakery 3. elife, 10, p.e65088.

    Notes
    -----
    GitHub  : https://github.com/segatalab/panphlan
    Paper   : https://elifesciences.org/articles/65088

    Parameters
    ----------
    self
    """
    pass


def strainsifter(self):
    """A straightforward bioinformatic pipeline for detecting the presence of
    a bacterial strain in one or more metagenome(s).
    StrainSifter is based on Snakemake. This pipeline allows you to output
    phylogenetic trees showing strain relatedness of input strains, as well
    as pairwise counts of single-nucleotide variants (SNVs) between input
    samples.

    References
    ----------
    Tamburini, F.B., Andermann, T.M., Tkachenko, E. et al. Precision
    identification of diverse bloodstream pathogens in the gut microbiome.
    Nat Med 24, 1809–1814 (2018). https://doi.org/10.1038/s41591-018-0202-8.

    Notes
    -----
    GitHub  : https://github.com/bhattlab/StrainSifter
    Paper   : https://doi.org/10.1038/s41591-018-0202-8

    Parameters
    ----------
    self
    """
    pass


def strainpro(self):
    """Characterizing the taxonomic diversity of a microbial community is
    very important to understand the roles of microorganisms. Next generation
    sequencing (NGS) provides great potential for investigation of a
    microbial community and leads to Metagenomic studies. NGS generates DNA
    sequences directly from microorganism samples, and it requires analysis
    tools to identify microbial species (or taxonomic composition) and
    estimate their relative abundance in the studied community. Here we
    developed a novel metagenomic analysis tool, called StrainPro, which is
    highly accurate both at characterizing microorganisms at strain-level and
    estimating their relative abundances. A unique feature of StrainPro is it
    identifies representative sequence segments from reference genomes. We
    generate three simulated datasets using known strain sequences and
    another three simulated datasets using unknown strain sequences.

    References
    ----------
    Lin, H.N., Lin, Y.L. and Hsu, W.L., 2019. StrainPro–a highly accurate
    Metagenomic strain-level profiling tool. bioRxiv, p.807149.

    Notes
    -----
    GitHub  : https://github.com/hsinnan75/StrainPro#download
    Paper   : https://doi.org/10.1101/807149

    Parameters
    ----------
    self
    """
    pass

