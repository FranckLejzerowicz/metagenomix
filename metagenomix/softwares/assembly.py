# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename, splitext
from metagenomix._inputs import genome_out_dir, genome_key
from metagenomix._io_utils import io_update, to_do, status_update
from metagenomix.core.parameters import tech_params


def quast_data(
        self,
        tech: str,
        pool: str,
        group_samples: dict
) -> tuple:
    """Collect the paths to the contigs for all the co-assembly groups of the
    current co-assembly pool name, and the respective metadata values
    qualifying these groups for the user-defined metadata columns.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pool : str
        Name of the co-assembly pool
    group_samples : dict
        Name of a group for the current co-assembly pool: list of pooled samples

    Returns
    -------
    contigs : list
        Paths to the contigs assembly files
    labels : list
        Labels for the co-assembly groups of the current pool
    """
    contigs, labels = [], []
    for group, assembly_outputs in sorted(self.inputs[pool].items()):
        if group[1] in self.soft.params['skip_samples']:
            continue
        if group[0] != tech:
            continue
        contigs.append(assembly_outputs[0])
        samples = group_samples[group[1]]
        label = '-'.join(group)
        if self.soft.params.get('label'):
            vals = self.config.meta.loc[samples, self.soft.params['label']]
            label += '-%s' % '__'.join(sorted(set(vals))).replace(' ', '_')
        labels.append(label)
    return contigs, labels


def quast_cmd(
        self,
        tech: str,
        pool: str,
        group_samples: dict,
        out_dir: str
) -> tuple:
    """Collect the QUAST command.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pool : str
        Name of the co-assembly pool
    group_samples : dict
        Name of a group for the current co-assembly pool: list of pooled samples
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        QUAST command
    contigs : list
        Paths to the contigs assembly files
    to_dos : list
    """
    contigs, labels = quast_data(self, tech, pool, group_samples)
    to_dos = status_update(self, tech, contigs, pool=pool)

    params = tech_params(self, tech)

    contigs_in = []
    cmd, cmd_rm = '', '\n'
    for contig in contigs:
        if contig.endswith('.gz'):
            contigs_in.append(contig.rstrip('.gz'))
            cmd += 'gunzip -c %s > %s\n' % (contig, contig.rstrip('.gz'))
            cmd_rm += 'rm %s\n' % contig.rstrip('.gz')
        else:
            contigs_in.append(contig)

    if params['circos_path']:
        cmd += 'export PATH=$PATH:%s\n' % params['circos_path']

    if params['path']:
        cmd += '%s/metaquast.py' % params['path']
    else:
        cmd += 'metaquast.py'

    # cmd += 'pe1'
    # cmd += 'pe2'
    # cmd += 'single'
    # cmd += 'pacbio'
    # cmd += 'nanopore'
    cmd += ' %s' % ' '.join(contigs_in)
    if labels:
        cmd += ' --labels "%s"' % ','.join(labels)
    cmd += ' --output-dir %s' % out_dir

    for param in [
        'sv_bedpe', 'sam', 'bam', 'ref_sam', 'ref_bam', 'plots_format',
        'est_insert_size', 'upper_bound_min_con', 'scaffold_gap_max_size',
        'unaligned_part_size', 'fragmented_max_indent', 'local_mis_size',
        'extensive_mis_size', 'ambiguity_score', 'min_identity', 'min_contig',
        'min_alignment', 'x_for_Nx', 'contig_thresholds', 'blast_db',
        'max_ref_number', 'operons', 'gene_thresholds', 'k_mer_size',
        'ambiguity_usage', 'r', 'references_list', 'g'
    ]:
        if params[param] or param in ['max_ref_number']:
            cmd += ' --%s %s' % (param.replace('_', '-'), params[param])

    if params['gene_finding']:
        cmd += ' --gene-finding'
    elif params['glimmer']:
        cmd += ' --glimmer'

    for boolean in [
        'circos', 'rna_finding', 'memory_efficient', 'conserved_genes_finding',
        'space_efficient', 'report_all_metrics', 'upper_bound_assembly',
        'skip_unaligned_mis_contigs', 'fragmented', 'strict_NA',
        'unique_mapping', 'use_all_alignments', 'reuse_combined_alignments',
        'use_input_ref_order', 'k_mer_stats', 'large', 'fungus', 'eukaryote',
        'split_scaffolds'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    if params['fast']:
        cmd += ' --fast'
    else:
        for boolean in ['no_check', 'no_plots', 'no_html', 'no_icarus',
                        'no_snps', 'no_gc', 'no_sv', 'no_read_stats']:
            if params[boolean]:
                cmd += ' --%s' % boolean.replace('_', '-')

    cmd += cmd_rm
    return cmd, contigs, to_dos


def quast(self) -> None:
    """QUAST stands for QUality ASsessment Tool. It evaluates
    genome/metagenome assemblies by computing various metrics. The current
    QUAST toolkit includes the general QUAST tool for genome assemblies,
    MetaQUAST, the extension for metagenomic datasets, QUAST-LG,
    the extension for large genomes (e.g., mammalians), and Icarus,
    the interactive visualizer for these softwares.

    References
    ----------
    Gurevich, Alexey, et al. "QUAST: quality assessment tool for genome
    assemblies." Bioinformatics 29.8 (2013): 1072-1075.

    Notes
    -----
    GitHub  : https://github.com/ablab/quast
    Docs    : http://quast.sourceforge.net/
    Paper   : https://doi.org/10.1093/bioinformatics/btt086

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for QUAST
        .pools : dict
            Pool name : { pool group : [pooled samples] }
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
    for pool, group_sams in self.pools.items():
        for tech in set([x[0] for x in self.inputs[pool]]):
            out_dir = '%s/%s/%s' % (self.dir, tech, pool)
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'][pool] = out_dir

            if self.config.force or to_do('%s/report.html' % out_dir):
                cmd, contigs, to_dos = quast_cmd(
                    self, tech, pool, group_sams, out_dir)
                key = (tech, pool)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=contigs, o_d=out_dir, key=key)
                self.soft.add_status(tech, pool, 1)
            else:
                self.soft.add_status(tech, pool, 0)


def spades_cmd(
        self,
        techs_inputs: dict,
        techs: list,
        tmp: str,
        out: str,
        log: str,
        hybrid: str,
        group: str
) -> tuple:
    """Create command lines for SPAdes.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    techs_inputs : dict
        Dict of the technology:paths to the input fasta/fastq.gz files
    techs : list
        Technologies to use for the (hybrid) assembly
    tmp : str
        Temporary folder for SPAdes
    out : str
        Path to the output folder
    log : str
        Path to the output log file
    hybrid : str
        Name of the hybrid assembly
    group : str
        Name of a group for the current co-assembly pool

    Returns
    -------
    cmd : str
        SPAdes.py command
    inputs : list
        Path to the input files
    to_dos : list
    """
    cmd = 'mkdir -p %s\n' % tmp
    cmd += 'spades.py'
    cmd += ' -o %s' % out
    cmd += ' -k %s' % ','.join(map(str, self.soft.params['k']))
    cmd += ' --memory %s' % self.soft.params['mem']
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --tmp-dir %s' % tmp
    for boolean in ['meta', 'careful', 'disable_gzip_output', 'disable_rr']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    cmd += ' --checkpoints %s' % self.soft.params['checkpoints']
    cmd += ' --cov-cutoff %s' % self.soft.params['cov_cutoff']

    if self.soft.params['only_assembler']:
        cmd += ' --only-assembler'
    elif self.soft.params['only_error_correction']:
        cmd += ' --only-error-correction'

    mode = self.soft.name.split('_')[-1]
    if mode == 'bio':
        cmd += ' --bio'
    elif mode == 'plasmid':
        cmd += ' --plasmid'
    elif mode == 'metaviral':
        cmd += ' --metaviral'

    inputs = []
    to_dos = []
    if self.soft.params['continue'] and not to_do(log):
        cmd += ' --continue'
    else:
        for tech, fastqs in techs_inputs.items():
            if tech not in techs:
                continue
            to_dos.extend(status_update(self, hybrid, fastqs, group=group))
            inputs.extend(fastqs)
            if tech != 'illumina':
                if len(fastqs):
                    cmd += ' --%s %s' % (tech, ' '.join(fastqs))
            else:
                for fastq in fastqs:
                    if len(fastq) == 1:
                        cmd += ' -s %s' % fastq
                    else:
                        if 'extendedFrags' in fastq:
                            cmd += ' --merge %s' % fastq
                        elif 'notCombined_1' in fastq:
                            cmd += ' -1 %s' % fastq
                        elif 'notCombined_2' in fastq:
                            cmd += ' -2 %s' % fastq
                        elif 'R1.fastq' in fastq:
                            cmd += ' -1 %s' % fastq
                        elif 'R2.fastq' in fastq:
                            cmd += ' -2 %s' % fastq

    cmd += '\nrm -rf %s\n' % tmp
    cmd += 'if [ -d %s/misc ]; then rm -rf %s/misc; fi\n' % (out, out)
    cmd += 'if [ -d %s/pipeline_state ];' % out
    cmd += ' then rm -rf %s/pipeline_state; fi\n' % out
    cmd += 'rm -rf %s/K*\n' % out
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out

    return cmd, inputs, to_dos


def hybridize_tech(
        self,
        techs_inputs: dict
) -> tuple:
    """Get the technologies to use for hybrid assembly and the folder name
    for this hybrid assembly.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    techs_inputs : dict
        Dict of the technology:paths to the input fasta/fastq.gz files

    Returns
    -------
    techs : list
        Technologies to use for the hybrid assembly
    hybrid : str
        Name of the hybrid assembly
    """
    techs = []
    for tech in self.soft.params['hybrid']:
        if tech in techs_inputs and techs_inputs[tech]:
            techs.append(tech)

    hybrid = '_'.join(techs)
    if len(techs) > 1:
        hybrid = 'hybrid__%s' % hybrid
    return techs, hybrid


def spades(self) -> None:
    """SPAdes - St. Petersburg genome assembler - is an assembly toolkit
    containing various assembly pipelines.

    References
    ----------
    Prjibelski, Andrey, et al. "Using SPAdes de novo assembler." Current
    protocols in bioinformatics 70.1 (2020): e102.
    Bankevich, Anton, et al. "SPAdes: a new genome assembly algorithm and its
    applications to single-cell sequencing." Journal of computational biology
    19.5 (2012): 455-477.

    Notes
    -----
    GitHub  : https://github.com/ablab/spades
    Docs    : https://cab.spbu.ru/software/spades/
    Paper   : https://doi.org/10.1002/cpbi.102

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SPAdes
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
    # iterate over inputs
    for group, techs_inputs in self.inputs[self.sam_pool].items():

        techs, hybrid = hybridize_tech(self, techs_inputs)

        # make output folder path
        tmp = '$TMPDIR/spades_%s_%s_%s' % (hybrid, group, self.sam_pool)
        out = '/'.join([self.dir, hybrid, self.sam_pool, group])
        # collect as folder to create under 'dirs'
        self.outputs['dirs'].append(out)

        # get expected output files
        contigs = '%s/contigs.fasta.gz' % out
        before_rr = '%s/before_rr.fasta.gz' % out
        first_pe = '%s/first_pe_contigs.fasta.gz' % out
        scaffolds = '%s/scaffolds.fasta.gz' % out
        log = '%s/spades.log.gz' % out
        outs = [contigs, before_rr, first_pe, scaffolds, log]
        self.outputs['outs'][(hybrid, group)] = outs

        # check if there is a need to run
        if self.config.force or to_do(contigs):
            cmd, inputs, to_dos = spades_cmd(self, techs_inputs, techs, tmp,
                                             out, log, hybrid, group)
            key = (hybrid, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=inputs, i_d=out, o_d=out, key=key)
            self.soft.add_status(hybrid, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(hybrid, self.sam_pool, 0, group=group)


def megahit_cmd(
        self,
        inputs: list,
        group: str,
        out: str,
        inter_dir: str,
        contigs: str,
        fastg: str,
        key: tuple
) -> str:
    """Create command lines for SPAdes.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    inputs : list
        Paths to the input fasta/fastq.gz files
    group : str
        Name of a group for the current co-assembly pool
    out : str
        Path to the output folder
    inter_dir : str
        Path to the intermediate_contigs folder
    contigs : str
        Path to the final contigs (gzipped)
    fastg : str
        Path to the longest k-mer contigs file (gzipped)
    key : tuple

    Returns
    -------
    cmd : str
        SPAdes.py command
    inputs : list
        Path to the input files
    """
    tmp = '$TMPDIR/megahit_%s_%s' % (self.sam_pool, group)
    cmd = 'mkdir -p %s\n' % tmp

    m_cmd = 'megahit'
    m_cmd += ' --tmp-dir %s' % tmp
    m_cmd += ' --out-prefix %s' % group
    m_cmd += ' --num-cpu-threads %s' % self.soft.params['cpus']
    for p in ['presets']:
        if self.soft.params[p]:
            m_cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
    for p in [
        'mem_flag', 'bubble_level', 'prune_level', 'k_min', 'k_max', 'k_step',
        'min_count', 'prune_depth', 'cleaning_rounds', 'min_contig_len',
        'merge_level', 'disconnect_ratio', 'low_local_ratio', 'memory'
    ]:
        if self.soft.params[p]:
            m_cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
    for boolean in [
        'no_mercy', 'no_local', 'kmin_1pass', 'no_hw_accel', 'keep_tmp_files'
    ]:
        if self.soft.params[boolean]:
            m_cmd += ' --%s' % boolean.replace('_', '-')
    m_cmd += ' --out-dir %s' % out
    if self.soft.params['continue'] and not to_do(inter_dir):
        m_cmd += ' --continue\n'
    else:
        if len(inputs) == 3:
            m_cmd += ' --read %s -1 %s -2 %s\n' % tuple(sorted(inputs))
        if len(inputs) == 2:
            m_cmd += ' -1 %s -2 %s\n' % tuple(inputs)
        if len(inputs) == 1:
            m_cmd += ' --read %s\n' % inputs[0]
    m_cmd += 'gzip %s\n' % contigs[:-3]

    if to_do(contigs) or self.config.force:
        cmd += 'rm -rf %s\n' % out
        cmd += m_cmd
        io_update(self, i_f=inputs, i_d=out, o_d=out, key=key)
    else:
        io_update(self, i_d=out, o_d=out, key=key)

    if to_do(fastg) or self.config.force:
        fasta = '%s/%s' % (inter_dir, basename(fastg).replace('.fg.gz', '.fa'))
        cmd += 'megahit_toolkit contig2fastg %s %s > %s\n' % (
            self.soft.params['k_max'], fasta, fastg[:-3])
        cmd += 'gzip %s\n' % (fastg[:-3])

    if to_do('%s.tar.gz' % inter_dir):
        cmd += '\ntar cpfz %s.tar.gz -C %s .\n' % (inter_dir, inter_dir)
        cmd += '\nrm -rf %s\n' % inter_dir

    cmd += '\nrm -rf %s\n' % tmp
    return cmd


def megahit(self) -> None:
    """MEGAHIT is an ultra-fast and memory-efficient NGS assembler. It is
    optimized for metagenomes, but also works well on generic single genome
    assembly (small or mammalian size) and single-cell assembly.

    References
    ----------
    Li, Dinghua, et al. "MEGAHIT: an ultra-fast single-node solution for
    large and complex metagenomics assembly via succinct de Bruijn graph."
    Bioinformatics 31.10 (2015): 1674-1676.

    Notes
    -----
    GitHub  : https://github.com/voutcn/megahit
    Docs    : https://github.com/voutcn/megahit/wiki
    Paper   : https://doi.org/10.1093/bioinformatics/btv033

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MEGAHIT
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
    tech = 'illumina'
    for group, techs_inputs in self.inputs[self.sam_pool].items():

        if tech not in techs_inputs:
            sys.exit('[megahit] Illumina reads are required (skipped)')
        inputs = techs_inputs[tech]
        to_dos = status_update(self, tech, inputs, group=group)

        out = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out)

        inter_dir = '%s/intermediate_contigs' % out
        contigs = '%s/%s.contigs.fa.gz' % (out, group)
        fastg = '%s/k%s.contigs.fg.gz' % (out, self.soft.params['k_max'])
        self.outputs['outs'][(tech, group)] = [contigs, out, fastg]

        if self.config.force or to_do(contigs) or to_do(fastg):
            key = (tech, group)
            cmd = megahit_cmd(
                self, inputs, group, out, inter_dir, contigs, fastg, key)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def plass_cmd(
        self,
        fastqs: list,
        contigs: str,
        tmp_dir: str
) -> str:
    """Collect PLASS command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastqs : list
        Paths to the fastqs input files
    contigs : str
        Path to the contigs output file
    tmp_dir : str
        Path to the temporary folder

    Returns
    -------
    cmd : str
        PLASS command
    """
    cmd = 'mkdir -p %s\n' % tmp_dir
    cmd += 'plass %s' % self.soft.params['type']
    cmd += ' %s' % ' '.join(fastqs)
    cmd += ' %s' % contigs
    cmd += ' %s' % tmp_dir
    cmd += ' --threads %s' % self.soft.params['cpus']
    for param in [
        'alph_size', 'k', 'split_memory_limit', 'e', 'c', 'a', 'cov_mode',
        'min_seq_id', 'min_aln_len', 'seq_id_mode', 'kmer_per_seq',
        'kmer_per_seq_scale', 'hash_shift', 'num_iterations', 'rescore_mode',
        'min_length', 'max_length', 'max_gaps', 'contig-start-mode',
        'contig_end_mode', 'orf_start_mode', 'forward_frames', 'reverse_frames',
        'translation_table', 'id_offset', 'protein_filter_threshold', 'dbtype',
        'db_load_mode', 'v', 'max_seq_len',
    ]:
        if len(param) == 1:
            cmd += ' -%s %s' % (param, self.soft.params[param])
        else:
            cmd += ' --%s %s' % (param, self.soft.params[param])

    for boolean in [
        'add_self_matches', 'spaced_kmer_mode', 'mask', 'mask_lower_case',
        'wrapped_scoring', 'adjust_kmer_len', 'include_only_extendable',
        'ignore_multi_kmer', 'translate', 'use_all_table_starts',
        'filter_proteins', 'shuffle', 'createdb_mode', 'compressed',
        'delete_tmp_inc', 'remove_tmp_files', 'filter_hits', 'sort_results',
        'create_lookup', 'write_lookup'
    ]:
        if self.soft.params[boolean]:
            cmd += ' --%s 1' % boolean
    cmd += '\n'
    cmd += 'rm -rf %s\n' % tmp_dir
    return cmd


def plass(self) -> None:
    """Plass (Protein-Level ASSembler) is a software to assemble short read
    sequencing data on a protein level. The main purpose of Plass is the
    assembly of complex metagenomic datasets. It assembles 10 times more
    protein residues in soil metagenomes than Megahit. Plass is GPL-licensed
    open source software that is implemented in C++ and available for Linux
    and macOS. The software is designed to run on multiple cores. Plass was
    used to create a Soil Reference Catalog (SRC) and a Marine Eukaryotic
    Reference Catalog (MERC).

    References
    ----------
    Steinegger, Martin, Milot Mirdita, and Johannes Söding. "Protein-level
    assembly increases protein sequence recovery from metagenomic samples
    manyfold." Nature methods 16.7 (2019): 603-606.

    Notes
    -----
    GitHub  : https://github.com/soedinglab/plass
    Docs    : https://ngs-docs.github.io/2018-cicese-metatranscriptomics/plass-paladin/
    Paper   : https://doi.org/10.1038/s41592-019-0437-4

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Plass
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech != 'illumina':
            continue
        out = '/'.join([self.dir, tech, self.sam_pool])
        self.outputs['dirs'].append(out)
        to_dos = status_update(self, tech, fastqs)

        if self.soft.params['type'] == 'nuclassemble':
            contigs = '%s/nucl_contigs.fasta' % out
        else:
            contigs = '%s/prot_contigs.fasta' % out
        self.outputs['outs'][(tech, sam)] = [contigs]

        if self.config.force or to_do(contigs):
            tmp_dir = '$TMPDIR/plass_%s' % self.sam_pool
            cmd = plass_cmd(self, fastqs, contigs, tmp_dir)
            if to_dos:
                self.outputs['cmds'].setdefault((tech,), []).append(False)
            else:
                self.outputs['cmds'].setdefault((tech,), []).append(cmd)
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def flye_cmd(
        self,
        tech: str,
        inputs: list,
        out: str
) -> str:
    """Create command lines for Flye.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'pacbio' or 'nanopore'
    inputs : list
        Path to the input files
    out : str
        Output folder for flye

    Returns
    -------
    cmd : str
        Flye command
    """
    params = tech_params(self, tech)
    cmd = 'flye'
    cmd += ' --%s %s' % (params['long_reads'], ' '.join(inputs))
    cmd += ' --out-dir %s' % out
    cmd += ' --threads %s' % params['cpus']
    if 'genome_size' in params:
        cmd += ' --genome-size %s' % params['genome_size']
    for boolean in ['meta', 'keep_haplotypes', 'scaffold', 'resume']:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    for param in ['iterations', 'min_overlap', 'asm_coverage', 'read_error',
                  'polish_target', 'stop_after', 'resume_from']:
        if params[param]:
            cmd += ' --%s %s\n' % (param.replace('_', '-'), params[param])
    return cmd


def flye(self) -> None:
    """Flye is a de novo assembler for single-molecule sequencing reads,
    such as those produced by PacBio and Oxford Nanopore Technologies. It is
    designed for a wide range of datasets, from small bacterial projects to
    large mammalian-scale assemblies. The package represents a complete
    pipeline: it takes raw PacBio / ONT reads as input and outputs polished
    contigs. Flye also has a special mode for metagenome assembly.
    Currently, Flye will produce collapsed assemblies of diploid genomes,
    represented by a single mosaic haplotype. To recover two phased
    haplotypes consider applying HapDup after the assembly.

    References
    ----------
    Kolmogorov, Mikhail, et al. "metaFlye: scalable long-read metagenome
    assembly using repeat graphs." Nature Methods 17.11 (2020): 1103-1110.

    Notes
    -----
    GitHub  : https://github.com/fenderglass/Flye
    Paper   : https://doi.org/10.1038/s41592-020-00971-x

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Flye
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
    for group, techs_inputs in self.inputs[self.sam_pool].items():

        if not sum([t in techs_inputs for t in ['pacbio', 'nanopore']]):
            sys.exit('[flye] Only for "pacbio" and/or "nanopore" (no input)')

        for tech, inputs in techs_inputs.items():
            if tech == 'illumina':
                continue

            out = '/'.join([self.dir, tech, self.sam_pool, group])
            self.outputs['dirs'].append(out)

            to_dos = status_update(self, tech, inputs, group=group)

            contigs = '%s/assembly.fasta' % out
            gfa = '%s/assembly_graph.gfa' % out
            gv = '%s/assembly_graph.gv' % out
            info = '%s/assembly_info.txt' % out
            self.outputs['outs'][(tech, group)] = [contigs, info, gfa, gv]

            if self.config.force or to_do(contigs):
                key = (tech, group)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    cmd = flye_cmd(self, tech, inputs, out)
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=inputs, i_d=out, o_d=out, key=key)
                self.soft.add_status(tech, self.sam_pool, 1, group=group)
            else:
                self.soft.add_status(tech, self.sam_pool, 0, group=group)


def canu_cmd(
        self,
        tech: str,
        inputs: list,
        group: str,
        out: str
) -> str:
    """Collect CANU command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'pacbio' or 'nanopore'
    inputs : list
        Input read fastx files
    group : str
        Name of the sample group within the pool
    out : str
        Output folder for CANU
    """
    params = tech_params(self, tech)
    cmd = 'if [ -d %s]; then rm -rf %s; fi\n' % (out, out)
    cmd += 'canu'
    cmd += ' -%s %s' % (tech, ' '.join(inputs))
    cmd += ' -d %s' % out
    cmd += ' -p %s' % group
    if 'specifications' in params:
        cmd += ' -s %s' % params['specifications']
    if 'genome_size' in params:
        cmd += ' genomeSize=%s' % params['genome_size']
    if params['processing']:
        cmd += ' -%s' % params['processing']
    if params['stage']:
        cmd += ' -%s' % params['stage']
    if params['rawErrorRate']:
        cmd += ' rawErrorRate=%s' % params['rawErrorRate']
    else:
        if tech == 'pacbio':
            cmd += ' rawErrorRate=0.300'
        else:
            cmd += ' rawErrorRate=0.500'
    if params['correctedErrorRate']:
        cmd += ' correctedErrorRate=%s' % params['correctedErrorRate']
    else:
        if tech == 'pacbio':
            cmd += ' correctedErrorRate=0.045'
        else:
            cmd += ' correctedErrorRate=0.144'

    for param in ['minReadLength', 'minOverlapLength', 'maxInputCoverage',
                  'corOutCoverage', 'corMhapSensitivity', 'corMinCoverage',
                  'redMemory', 'oeaMemory', 'batMemory']:
        cmd += ' %s=%s' % (param, params[param])

    for threads in ['bat', 'cns', 'cormhap', 'cormmap', 'corovl', 'cor',
                    'executive', 'hap', 'max', 'meryl', 'min', 'obtmhap',
                    'obtmmap', 'obtovl', 'oea', 'ovb', 'ovs', 'red',
                    'utgmhap', 'utgmmap', 'utgovl']:
        cmd += ' %sThreads=%s' % (threads, params['cpus'])
    return cmd


def canu(self) -> None:
    """Canu is a fork of the Celera Assembler designed for high-noise
    single-molecule sequencing (such as the PacBio RSII or Oxford Nanopore
    MinION).

    References
    ----------
    Koren, Sergey, et al. "Canu: scalable and accurate long-read assembly via
    adaptive k-mer weighting and repeat separation." Genome research 27.5
    (2017): 722-736.

    Notes
    -----
    GitHub  : https://github.com/marbl/canu
    Docs    : https://canu.readthedocs.io/en/latest/
    Paper   : https://doi.org/10.1101/gr.215087.116

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for CANU
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
    for group, techs_inputs in self.inputs[self.sam_pool].items():

        if not [t for t in ['pacbio', 'nanopore'] if techs_inputs.get(t)]:
            print('[canu] Only for "pacbio" and/or "nanopore" (no input)')
            continue

        for tech, inputs in techs_inputs.items():
            if tech == 'illumina':
                continue

            to_dos = status_update(self, tech, inputs, group=group)

            out = '/'.join([self.dir, tech, self.sam_pool, group])
            self.outputs['dirs'].append(out)

            report = '%s/%s.report' % (out, group)
            contigs = '%s/%s.contigs.fasta' % (out, group)
            self.outputs['outs'][(tech, group)] = [contigs, report]

            if self.config.force or to_do(contigs):
                cmd = canu_cmd(self, tech, inputs, group, out)
                key = (tech, group)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=inputs, i_d=out, o_d=out, key=key)
                self.soft.add_status(tech, self.sam_pool, 1, group=group)
            else:
                self.soft.add_status(tech, self.sam_pool, 0, group=group)


def unicycler_cmd(
        self,
        techs_inputs: dict,
        key: tuple,
        out: str,
        techs: list
) -> tuple:
    """Collect unicycler command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    techs_inputs : dict
        Path to the input fastqs files per technology
    key : tuple
        Technology and co-assembly group
    out : str
        Output folder for unicycler
    techs : list
        Pair or technologies

    Returns
    -------
    cmd : str
        Unicycler commands
    inputs : list
        Path to all input fastqs files
    """
    tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, self.sam_pool, '_'.join(key))
    cmd = 'unicycler'
    inputs = []
    for tech in techs:
        if tech == 'illumina':
            illumina = tuple(techs_inputs[tech])
            inputs.extend(techs_inputs[tech])
            if len(illumina) == 1:
                cmd += ' --unpaired %s' % illumina[0]
            elif len(illumina) == 2:
                cmd += ' --short1 %s --short2 %s' % illumina
            elif len(illumina) == 3:
                cmd += ' --short1 %s --short2 %s --unpaired %s' % illumina
        else:
            if len(techs_inputs[tech]):
                cmd += ' --long %s' % techs_inputs[tech][0]
                inputs.extend(techs_inputs[tech])
    if not inputs:
        return '', inputs
    cmd += ' --out %s' % out
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --spades_tmp_dir %s' % tmp
    for param in [
        'contamination', 'bcftools_path', 'java_path', 'bowtie2_path',
        'bowtie2_build_path', 'samtools_path', 'pilon_path', 'scores',
        'start_genes', 'existing_long_read_assembly', 'tblastn_path',
        'makeblastdb_path', 'spades_path', 'racon_pth', 'min_bridge_qual',
        'start_gene_cov', 'start_gene_id', 'max_kmer_frac', 'min_kmer_frac',
        'depth_filter', 'mode', 'min_component_size', 'min_dead_end_size',
        'min_polish_size', 'min_anchor_seg_len', 'min_fasta_length', 'keep',
        'kmer_count', 'linear_seqs', 'low_score', 'kmers', 'verbosity'
    ]:
        if param in self.soft.params:
            cmd += ' --%s %s' % (param, self.soft.params[param])
    for boolean in ['vcf', 'no_correct', 'largest_component',
                    'no_miniasm', 'no_rotate', 'no_pilon']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % self.soft.params[param]
    return cmd, inputs


def unicycler(self) -> None:
    """Unicycler is an assembly pipeline for bacterial genomes. It can
    assemble Illumina-only read sets where it functions as a
    SPAdes-optimiser. It can also assembly long-read-only sets (PacBio or
    Nanopore) where it runs a miniasm+Racon pipeline. For the best possible
    assemblies, give it both Illumina reads and long reads, and it will
    conduct a short-read-first hybrid assembly.

    References
    ----------
    Wick, Ryan R., et al. "Unicycler: resolving bacterial genome assemblies
    from short and long sequencing reads." PLoS computational biology 13.6
    (2017): e1005595.

    Notes
    -----
    GitHub  : https://github.com/rrwick/Unicycler
    Docs    : https://github.com/rrwick/Unicycler/wiki
    Paper   : https://doi.org/10.1371/journal.pcbi.1005595

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for unicycler
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
        fastqs = self.inputs[self.sam_pool]
    else:
        fastqs = {self.sam_pool: self.inputs[self.sam_pool]}

    for sam_group, techs_inputs in fastqs.items():
        if 'illumina' not in techs_inputs:
            sys.exit('[unicycler] Illumina reads are required (skipped)')

        techs = ['illumina']
        if self.soft.params['hybrid']:
            if not techs_inputs.get(self.soft.params['hybrid']):
                print('[unicycler] No long reads/hybrid assembly:', sam_group)
            techs += [self.soft.params['hybrid']]

        hybrid = '_'.join(techs)
        key = (hybrid, sam_group)
        if len(techs) > 1:
            hybrid = 'hybrid__%s' % hybrid

        out = '/'.join([self.dir, hybrid, self.sam_pool, sam_group])
        gfa = '%s/assembly.gfa' % out
        fasta = '%s/assembly.fasta' % out
        log = '%s/unicycler.log' % out
        self.outputs['dirs'].append(out)
        self.outputs['outs'][(hybrid, sam_group)] = [fasta, gfa, log]

        if self.config.force or to_do(fasta):
            cmd, inputs = unicycler_cmd(self, techs_inputs, key, out, techs)
            to_dos = status_update(self, hybrid, inputs, group=sam_group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=inputs, o_d=out, key=key)
            self.soft.add_status(hybrid, self.sam_pool, 1, group=sam_group)
        else:
            self.soft.add_status(hybrid, self.sam_pool, 0, group=sam_group)


# def miniasm_cmd(self):
#     cmd = '\nminimap2'
#     cmd += ' -x ava-pb'
#     cmd += ' -t8 pb-reads.fq pb-reads.fq'
#     cmd += ' | gzip -1 > reads.paf.gz'
#     cmd += ' > reads.paf.gz'
#     cmd += '\nminiasm'
#     cmd += ' -f reads.fq'
#     cmd += ' reads.paf.gz'
#     cmd += ' > reads.gfa'
#
#
def miniasm(self):
    """

    Parameters
    ----------
    self

    Returns
    -------

    """
    # cmd = miniasm_cmd(self)
    pass


def necat_cmd(
        self,
        group: str,
        fastxs: list,
        out: str
) -> str:
    """Collect NECAT command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    group : str
        Name of the sample group within the pool
    fastxs : list
        Paths to the input files
    out : str
        Output folder for NECAT

    Returns
    -------
    cmd : str
        NECAT command
    """
    reads_fp = '%s/reads.txt' % out
    cmd = necat_reads(fastxs, reads_fp)
    config_fp = '%s/config.txt' % out
    cmd += necat_config(self, reads_fp, config_fp, group)
    cmd += 'cd %s/..\n' % out
    cmd += 'export PATH="$PATH:%s"\n' % self.soft.params['path']
    cmd += 'necat.pl correct %s\n' % config_fp
    cmd += 'necat.pl assemble %s\n' % config_fp
    cmd += 'necat.pl bridge %s\n' % config_fp
    return cmd


def necat_reads(
        fastxs: list,
        reads_fp: str
) -> str:
    """Make the `reads.txt` input file.

    Parameters
    ----------
    fastxs : list
        Paths to the reads input files
    reads_fp : str
        Path to the `reads.txt` input file

    Returns
    -------
    cmd : str
        echo commands to make the `reads.txt` input file
    """
    cmd = ''
    for fdx, fastx in enumerate(fastxs):
        if fdx:
            cmd += 'echo -e "%s\\n" >> %s\n' % (fastx, reads_fp)
        else:
            cmd += 'echo -e "%s\\n" > %s\n' % (fastx, reads_fp)
    return cmd


def necat_config(
        self,
        reads_fp: str,
        config_fp: str,
        group: str
) -> str:
    """Create command lines for NECAT.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    reads_fp : str
        Path to the `reads.txt` input file
    config_fp : str
        Path to the `config.txt` input file
    group : str
        Name of the sample group within the pool

    Returns
    -------
    cmd : str
        echo commands to make the `config.txt` input file
    """

    opts = {x: [y for y in self.soft.params if y.startswith(y.rsplit('_')[0])]
            for x in self.soft.params if 'options' in x}

    config = list()
    config.append('PROJECT=%s' % group)
    config.append('ONT_READ_LIST=%s' % reads_fp)
    if self.soft.params['genome_size']:
        config.append('GENOME_SIZE=%s' % self.soft.params['genome_size'])
    config.append('THREADS=%s' % self.soft.params['cpus'])
    config.append('MIN_READ_LENGTH=%s' % self.soft.params['min_read_length'])
    config.append('PREP_OUTPUT_COVERAGE=%s' % self.soft.params[
        'prep_output_coverage'])
    for opt, opt_params in opts.items():
        params = ['-%s %s' % (x[-1], self.soft.params[x]) for x in opt_params]
        config.append('%s=%s' % (opt[:-2].upper(), ' '.join(params)))
    config.append('NUM_ITER=%s' % self.soft.params['num_iter'])
    config.append('CNS_OUTPUT_COVERAGE=%s' % self.soft.params[
        'cns_output_coverage'])
    config.append('CLEANUP=%s' % self.soft.params['cleanup'])
    config.append('POLISH_CONTIGS=%s' % self.soft.params['polish_contig'])

    cmd = ''
    for cdx, conf in enumerate(config):
        if cdx:
            cmd += 'echo -e "%s\\n" >> %s\n' % (conf, config_fp)
        else:
            cmd += 'echo -e "%s\\n" > %s\n' % (conf, config_fp)
    return cmd


def necat(self) -> None:
    """NECAT is an error correction and de-novo assembly tool for Nanopore
    long noisy reads.

    References
    ----------
    Chen, Ying, et al. "Efficient assembly of nanopore reads via highly
    accurate and intact error correction." Nature Communications 12.1
    (2021): 1-10.

    Notes
    -----
    GitHub  : https://github.com/xiaochuanle/necat
    Docs    : http://www.tgsbioinformatics.com/necat
    Paper   : https://doi.org/10.1038/s41467-020-20236-7

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for NECAT
        .pool : str
            Pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    tech = 'nanopore'
    for group, fastxs in self.inputs[self.sam_pool].get(tech, {}).items():
        # define the output folder
        out = '/'.join([self.dir, tech, self.sam_pool, group])

        to_dos = status_update(self, tech, fastxs, group=group)

        # get some of the expected outputs files
        final = '%s/1-consensus/cns_final.fasta' % out
        contigs = '%s/4-fsa/contigs.fasta' % out
        bridged = '%s/6-bridge_contigs/bridged_contigs.fasta' % out
        polish = ''
        if self.soft.params['polish_contigs']:
            polish = '%s/6-bridge_contigs/polished_contigs.fasta' % out

        # [machinery] Collect the folder(s) to be created
        self.outputs['dirs'].append(out)
        # [machinery] Collect the outputs files in a standard / step-wise way
        self.outputs['outs'][(tech, group)] = [contigs, final, bridged, polish]

        # if user in command line want to force re-creation of tool's command
        # or if the file "contigs" is to do (i.e. not yet written by the tool)
        if self.config.force or to_do(contigs):
            cmd = necat_cmd(self, group, fastxs, out)
            # [machinery] Collect the command line
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fastxs, o_d=out, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def metamic_cmd(
        self,
        tech: str,
        bam: str,
        contigs: str,
        assembler: str,
        out_dir: str
) -> str:
    """Collect the command line for metaMIC.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    bam : str
        Path to the input mapping alignment
    contigs : str
        Path to the input contigs fasta.(gz) file
    assembler : str
        Name of the assembler
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        metaMIC commands
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    if contigs.endswith('.fa.gz') or contigs.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % contigs.rstrip('.gz')
        contigs = contigs.rstrip('.gz')

    cmd += 'samtools faidx %s\n' % contigs

    pileup = '%s.pileup' % splitext(bam)[0]
    cmd += '\nsamtools mpileup'
    cmd += ' -C 50 -A -f %s %s' % (contigs, bam)
    cmd += " | awk '$3 != \"N\"' > %s\n" % pileup

    cmd += '\nmetaMIC extract_feature'
    cmd += ' --bam %s' % bam
    cmd += ' --contig %s' % contigs
    cmd += ' --pileup %s' % pileup
    cmd += ' --output %s' % out_dir
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --mlen %s' % params['mlen']
    cmd += ' --mode meta\n'

    cmd += '\nmetaMIC predict'
    cmd += ' --contig %s' % contigs
    cmd += ' --output %s' % out_dir
    if assembler == 'megahit':
        cmd += ' --assembler MEGAHIT'
    elif assembler == 'idba_ud':
        cmd += ' --assembler IDBA_UD'
    elif assembler == 'spades':
        cmd += ' --assembler metaSPAdes'
    for param in ['mlen', 'slen', 'nb', 'rb', 'at', 'st']:
        if params[param] is not None:
            cmd += ' --%s %s' % (param, params[param])
    cmd += ' --mode meta\n'
    cmd += cmd_rm

    return cmd


def get_metamic(
        self,
        bam: str,
        bam_infos: list,
        assembler: str
) -> None:
    """Get the info and collect the command for metaMIC.

    Parameters
    ----------
    self
    bam : str
        Path to the input BAM file
    bam_infos : list
        [tech, sample, aligner, coassembly, coassembly group, contig path]
    assembler : str
        Name of the assembler
    """
    tech, sample, aligner, _, group, contigs = bam_infos
    out_dir = genome_out_dir(self, tech, group) + '/' + aligner
    self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
    self.outputs['dirs'].append(out_dir)

    contigs_gz = contigs + '.gz'
    to_dos = status_update(
        self, tech, [bam, contigs_gz], self.sam_pool, group=group)

    key = genome_key(tech, group, aligner)
    out_fp = '%s/output' % out_dir
    if self.config.force or to_do(out_fp):
        if to_dos:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            cmd = metamic_cmd(self, tech, bam, contigs_gz, assembler, out_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_f=[bam, contigs_gz], o_d=out_dir, key=key)
        self.soft.add_status(tech, self.sam_pool, 1, group=group)
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group=group)


def metamic(self) -> None:
    """metaMIC is a fully automated tool for identifying and correcting
    misassemblies of (meta)genomic assemblies with the following three steps.
    Firstly, metaMIC extracts various types of features from the alignment
    between paired-end sequencing reads and the assembled contigs. Secondly,
    the features extracted in the first step will be used as input of a
    random forest classifier for identifying misassembled metagenomic
    assemblies. Thirdly, metaMIC will localize misassembly breakpoints for
    each misassembled contig and then correct misassemblies by splitting
    into parts at the breakpoints.

    References
    ----------
    Lai, S., Pan, S., Sun, C., Coelho, L.P., Chen, W.H. and Zhao, X.M.,
    2022. metaMIC: reference-free Misassembly Identification and Correction
    of de novo metagenomic assemblies. Genome Biology, 23(1), p.242.

    Notes
    -----
    GitHub  : https://github.com/ZhaoXM-Lab/metaMIC
    Paper   : https://doi.org/10.1186/s13059-022-02810-y

    Parameters
    ----------
    self : Commands class instance
    """
    if not self.soft.prev.startswith('mapping'):
        sys.exit("[%s] Only run after a mapping_* command" % self.soft.name)

    assembler = self.softs[self.soft.prev][
        self.hashes[tuple(self.soft.path[:-1])]].prev
    assemblers = ['megahit', 'idba_ud', 'spades']
    if assembler not in assemblers:
        sys.exit("[%s] mapping_* not done after %s" % (
            self.soft.name, ', or '.join(assemblers)))

    if self.sam_pool in self.pools:
        for bam, bam_infos in self.inputs[self.sam_pool].items():
            get_metamic(self, bam, bam_infos, assembler)


def trycycler(self) -> None:
    """Trycycler is a tool for generating consensus long-read assemblies for
    bacterial genomes. I.e. if you have multiple long-read assemblies for the
    same isolate, Trycycler can combine them into a single assembly that is
    better than any of your inputs.

    References
    ----------
    Wick, R.R., Judd, L.M., Cerdeira, L.T., Hawkey, J., Méric, G., Vezina,
    B., Wyres, K.L. and Holt, K.E., 2021. Trycycler: consensus long-read
    assemblies for bacterial genomes. Genome biology, 22, pp.1-17.

    Notes
    -----
    GitHub  : https://github.com/rrwick/Trycycler
    Paper   : https://doi.org/10.1186/s13059-021-02483-z

    Parameters
    ----------
    self : Commands class instance
    """
    pass
