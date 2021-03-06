# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
from os.path import basename

from metagenomix._io_utils import io_update, to_do
from metagenomix.parameters import tech_params


def quast_data(
        self,
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
        contigs.append(assembly_outputs[1])
        if self.soft.params['label']:
            samples = group_samples[group]
            vals = self.config.meta.loc[samples, self.soft.params['label']]
            labels.append('__'.join(sorted(set(vals))))
    return contigs, labels


def quast_cmd(
        self,
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
    """
    contigs, labels = quast_data(self, pool, group_samples)

    cmd = 'metaquast.py'
    cmd += ' %s' % ' '.join(contigs)
    if labels:
        cmd += ' --labels %s' % ','.join(labels)
    cmd += ' --output-dir %s' % out_dir
    cmd += ' --min-contig %s' % self.soft.params['min_contig']
    cmd += ' --ambiguity-usage %s' % self.soft.params['ambiguity_usage']
    for boolean in ['circos', 'glimmer', 'rna_finding',
                    'conserved_genes_finding', 'space_efficient']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    if self.soft.params['fast']:
        cmd += ' --fast'
    else:
        for boolean in ['no_check', 'no_plots', 'no_html', 'no_icarus',
                        'no_snps', 'no_gc', 'no_sv', 'no_read_stats']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean.replace('_', '-')
    return cmd, contigs


def quast(self) -> None:
    """Create command lines for QUAST.

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
    for pool, group_samples in self.pools.items():

        out_dir = '%s/%s' % (self.dir, pool)
        if self.soft.params['label']:
            out_dir += '_%s' % self.soft.params['label']
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][pool] = out_dir

        if self.config.force or not to_do('%s/report.html' % out_dir):
            cmd, contigs = quast_cmd(self, pool, group_samples, out_dir)
            self.outputs['cmds'].setdefault(pool, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out_dir, key=pool)


def spades_cmd(
        self,
        techs_inputs: dict,
        techs: list,
        tmp: str,
        out: str,
        log: str,
        mode: str
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
        Technologies to use for the hybrid assembly
    tmp : str
        Temporary folder for SPAdes
    out : str
        Path to the output folder
    log : str
        Path to the output log file
    mode : str
        Spades mode indicated after the underscore of the tool

    Returns
    -------
    cmd : str
        SPAdes.py command
    inputs : list
        Path to the input files
    """
    cmd = 'mkdir -p %s\n' % tmp
    cmd += 'spades.py'
    cmd += ' -k %s' % ','.join(map(str, self.soft.params['k']))
    cmd += ' --memory %s' % self.soft.params['mem']
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --tmp-dir %s' % tmp
    for boolean in ['meta', 'careful', 'checkpoints', 'disable_gzip_output',
                    'disable_rr', 'cov_cutoff']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    if self.soft.params['only_assembler']:
        cmd += ' --only-assembler'
    elif self.soft.params['only_error_correction']:
        cmd += ' --only-error-correction'

    if mode == 'bio':
        cmd += ' --bio'
    elif mode == 'plasmid':
        cmd += ' --plasmid'
    elif mode == 'metaviral':
        cmd += ' --metaviral'

    inputs = []
    cmd += ' -o %s' % out
    if self.soft.params['continue'] and not to_do(log):
        cmd += ' --continue'
    else:
        for tech, fastqs in techs_inputs.items():
            if tech not in techs:
                continue
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
                            cmd += ' -2 %s' % fastq
                        elif 'R2.fastq' in fastq:
                            cmd += ' -2 %s' % fastq
    cmd += '\nrm -rf %s\n' % tmp
    return cmd, inputs


def hybridize(
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
        if tech in techs_inputs:
            techs.append(tech)

    hybrid = '_'.join(techs)
    if len(techs) > 1:
        hybrid = 'hybrid__%s' % hybrid
    return techs, hybrid


def spades(self) -> None:
    """Create command lines for SPAdes.

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
    mode = self.soft.name.split('_')[-1]
    for group, techs_inputs in self.inputs[self.pool].items():
        techs, hybrid = hybridize(self, techs_inputs)
        mode_group = '%s/%s/%s' % (mode, hybrid, group)

        tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, mode_group, self.pool)
        out = '%s/%s/%s/%s/%s' % (self.dir, mode, hybrid, self.pool, group)
        self.outputs['dirs'].append(out)

        before_rr = '%s/before_rr.fasta' % out
        contigs = '%s/contigs.fasta' % out
        first_pe = '%s/first_pe_contigs.fasta' % out
        scaffolds = '%s/scaffolds.fasta' % out
        log = '%s/spades.log' % out
        outs = [before_rr, contigs, first_pe, scaffolds, log]
        self.outputs['outs'][group] = outs

        if self.config.force or to_do(contigs):
            cmd, inputs = spades_cmd(
                self, techs_inputs, techs, tmp, out, log, mode)
            self.outputs['cmds'].setdefault(mode_group, []).append(cmd)
            io_update(self, i_f=inputs, i_d=out, o_d=out, key=mode_group)


def megahit_cmd(
        self,
        inputs: list,
        group: str,
        tmp: str,
        out: str,
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
    tmp : str
        Temporary folder for SPAdes
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        SPAdes.py command
    inputs : list
        Path to the input files
    """
    cmd = 'mkdir -p %s\n' % tmp
    cmd += 'megahit'
    cmd += ' --tmp-dir %s' % tmp
    cmd += ' --out-prefix %s' % group
    cmd += ' --num-cpu-threads %s' % self.soft.params['cpus']
    for p in ['presets']:
        if self.soft.params[p]:
            cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
    for p in [
        'mem_flag', 'bubble_level', 'prune_level', 'k_min', 'k_max', 'k_step',
        'min_count', 'prune_depth', 'cleaning_rounds', 'min_contig_len',
        'merge_level', 'disconnect_ratio', 'low_local_ratio', 'memory'
    ]:
        if self.soft.params[p]:
            cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
    for boolean in [
        'no_mercy', 'no_local', 'kmin_1pass', 'no_hw_accel', 'keep_tmp_files'
    ]:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --out-dir %s' % out
    out_folder = '%s/intermediate_contigs' % out
    if self.soft.params['continue'] and not to_do(folder=out_folder):
        cmd += ' --continue'
    else:
        if len(inputs) == 3:
            cmd += ' --read %s -1 %s -2 %s' % tuple(inputs)
        if len(inputs) == 2:
            cmd += ' -1 %s -2 %s' % tuple(inputs)
        if len(inputs) == 1:
            cmd += ' --read %s' % inputs[0]
    cmd += '\nrm -rf %s\n' % tmp
    return cmd


def megahit(self) -> None:
    """Create command lines for MEGAHIT.

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
    for group, techs_inputs in self.inputs[self.pool].items():

        if 'illumina' not in techs_inputs:
            sys.exit('[megahit] Illumina reads are required (skipped)')
        inputs = techs_inputs['illumina']

        tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, self.pool, group)
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        self.outputs['dirs'].append(out)

        contigs = '%s/%s.contigs.fa' % (out, group)
        self.outputs['outs'][group] = [out, contigs]

        if self.config.force or to_do(contigs):
            cmd = megahit_cmd(self, inputs, group, tmp, out)
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=inputs, i_d=out, o_d=out, key=group)


def plass(self) -> None:
    """Create command lines for Plass.

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
    out = '%s/%s' % (self.dir, self.sam)
    self.outputs['dirs'].append(out)
    inputs = self.inputs[self.sam]
    if self.soft.prev == 'flash':
        inputs = [x.replace('.fasta', '.fastq') for x in inputs if
                  'notCombined_' in basename(x)]
    out_fp = '%s/plass_%sassembly.fasta' % (out, self.soft.params.type)
    self.outputs['outs'].append(out_fp)
    tmp_dir = '$TMPDIR/plass_%s' % self.sam
    cmd = 'plass %sassemble' % self.soft.params.type
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --num-iterations 6'
    cmd += ' --min-seq-id 0.8'
    cmd += ' -c 0.8'
    cmd += ' %s' % ' '.join(sorted(inputs))
    cmd += ' %s' % out_fp
    cmd += ' %s' % tmp_dir
    if self.config.force or to_do(out_fp):
        cmd = 'mkdir -p %s\n%s\nrm -rf %s' % (tmp_dir, cmd, tmp_dir)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, i_d=tmp_dir, o_f=out_fp)


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
    """Create command lines for Flye.

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
    for group, techs_inputs in self.inputs[self.pool].items():

        if not sum([t in techs_inputs for t in ['pacbio', 'nanopore']]):
            sys.exit('[flye] Only for "pacbio" and/or "nanopore" (no input)')

        for tech, inputs in techs_inputs.items():
            if tech == 'illumina':
                continue

            out = '%s/%s/%s/%s' % (self.dir, tech, self.pool, group)
            self.outputs['dirs'].append(out)

            contigs = '%s/assembly.fasta' % out
            gfa = '%s/assembly_graph.gfa' % out
            gv = '%s/assembly_graph.gv' % out
            info = '%s/assembly_info.txt' % out
            self.outputs['outs'][(tech, group)] = [info, contigs, gfa, gv]

            if self.config.force or to_do(contigs):
                cmd = flye_cmd(self, tech, inputs, out)
                tech_group = '%s_%s' % (tech, group)
                self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
                io_update(self, i_f=inputs, i_d=out, o_d=out, key=tech_group)


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
    if 'path' in params:
        cmd += 'export PATH=$PATH:%s\n' % params['path']
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
    """Create command lines for CANU.

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
    for group, techs_inputs in self.inputs[self.pool].items():

        if not [t for t in ['pacbio', 'nanopore'] if techs_inputs.get(t)]:
            sys.exit('[flye] Only for "pacbio" and/or "nanopore" (no input)')

        for tech, inputs in techs_inputs.items():
            if tech == 'illumina':
                continue

            out = '%s/%s/%s' % (self.dir, self.pool, group)
            self.outputs['dirs'].append(out)

            report = '%s/%s.report' % (out, group)
            contigs = '%s/%s.contigs.fasta' % (out, group)
            self.outputs['outs'][(tech, group)] = [report, contigs]

            if self.config.force or to_do(contigs):
                cmd = canu_cmd(self, tech, inputs, group, out)
                tech_group = '%s_%s' % (tech, group)
                self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
                io_update(self, i_f=inputs, i_d=out, o_d=out, key=tech_group)


def unicycler_cmd(
        self,
        techs_inputs: dict,
        tmp: str,
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
    tmp : str
        Path to a temporary directory
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
                cmd += ' --unpaired %s --short1 %s --short2 %s' % illumina
        else:
            cmd += ' --long %s' % techs_inputs[tech][0]
            inputs.extend(techs_inputs[tech])
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
    """Create command lines for unicycler.

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
    if self.pool in self.pools:
        fastqs = self.inputs[self.pool]
    else:
        fastqs = self.inputs

    for sam_group, techs_inputs in fastqs.items():
        if 'illumina' not in techs_inputs:
            sys.exit('[unicycler] Illumina reads are required (skipped)')

        long_techs = [t for t in ['pacbio', 'nanopore'] if techs_inputs.get(t)]
        if not long_techs and self.soft.params['hybrid']:
            print('[unicycler] No long reads/hybrid assembly w/ %s' % sam_group)

        techs_pairs = [['illumina']]
        if self.soft.params['hybrid']:
            techs_pairs += [['illumina', tech] for tech in long_techs]

        for techs in techs_pairs:
            hybrid = '_'.join(techs)
            key = '%s_%s' % (hybrid, sam_group)
            print(sam_group)
            if len(techs) > 1:
                hybrid = 'hybrid__%s' % hybrid

            tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, self.pool, key)
            out = '%s/%s/%s/%s' % (self.dir, hybrid, self.pool, sam_group)
            gfa = '%s/assembly.gfa' % out
            fasta = '%s/assembly.fasta' % out
            log = '%s/unicycler.log' % out
            if self.pool in self.pools:
                self.outputs['outs'][sam_group] = [gfa, fasta, log]
            else:
                self.outputs['outs'] = [gfa, fasta, log]
            if self.config.force or to_do(fasta):
                cmd, inputs = unicycler_cmd(self, techs_inputs, tmp, out, techs)
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                if self.pool in self.pools:
                    io_update(self, i_f=inputs, i_d=tmp, o_d=out, key=key)
                else:
                    io_update(self, i_f=inputs, i_d=tmp, o_d=out, key=key)


def miniasm_cmd(self):
    cmd = '\nminimap2'
    cmd += ' -x ava-pb'
    cmd += ' -t8 pb-reads.fq pb-reads.fq'
    cmd += ' | gzip -1 > reads.paf.gz'
    cmd += ' > reads.paf.gz'
    cmd += '\nminiasm'
    cmd += ' -f reads.fq'
    cmd += ' reads.paf.gz'
    cmd += ' > reads.gfa'


def miniasm(self):
    cmd = miniasm_cmd(self)
    pass


def necat_cmd(self, group: str, fastxs: list, out: str):
    """Create command lines for NECAT.

    Parameters
    ----------
    config_fp : str
        Config file for NECAT
    out : str
        Output folder for NECAT
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


def necat_reads(fastxs: list, reads_fp: str):
    cmd = ''
    for fdx, fastx in enumerate(fastxs):
        if fdx:
            cmd += 'echo -e "%s\\n" >> %s\n' % (fastx, reads_fp)
        else:
            cmd += 'echo -e "%s\\n" > %s\n' % (fastx, reads_fp)
    return cmd


def necat_config(self, reads_fp: str, config_fp: str, group: str):

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
    """Create command lines for NECAT.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for NECAT
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
    tech = 'nanopore'
    for group, fastxs in self.inputs[self.pool].get(tech, {}).items():
        # define the output folder
        out = '%s/%s/%s' % (self.dir, self.pool, group)

        # get some of the expected outputs files
        cns_final = '%s/1-consensus/cns_final.fasta' % out
        contigs = '%s/4-fsa/contigs.fasta' % out
        bridged = '%s/6-bridge_contigs/bridged_contigs.fasta' % out
        polished = ''
        if self.soft.params['polish_contigs']:
            polished = '%s/6-bridge_contigs/polished_contigs.fasta' % out

        # [machinery] Collect the folder(s) to be created
        self.outputs['dirs'].append(out)
        # [machinery] Collect the outputs files in a standard / step-wise way
        self.outputs['outs'][group] = [cns_final, contigs, bridged, polished]

        # if user in command line want to force re-creation of tool's command
        # or if the file "contigs" is to do (i.e. not yet written by the tool)
        if self.config.force or to_do(contigs):
            cmd = necat_cmd(self, group, fastxs, out)
            # [machinery] Collect the command line
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=fastxs, o_d=out, key=group)
