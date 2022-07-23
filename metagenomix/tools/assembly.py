# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename, splitext
from metagenomix._io_utils import io_update, todo


def get_quast_contigs_labels(self, pool, group_samples):
    contigs = []
    labels = []
    for group, fps in sorted(self.inputs[pool].items()):
        contigs.append(fps[1])
        if self.soft.params['label']:
            samples = group_samples[group]
            vals = self.config.meta.loc[samples, self.soft.params['label']]
            labels.append('__'.join(sorted(set(vals))))
    return contigs, labels


def quast(self):
    for pool, group_samples in self.pools.items():
        contigs, labels = get_quast_contigs_labels(self, pool, group_samples)

        out_dir = '%s/%s' % (self.dir, pool)
        if self.soft.params['label']:
            out_dir += '_%s' % self.soft.params['label']

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

        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][pool] = out_dir
        self.outputs['cmds'].setdefault(pool, []).append(cmd)
        io_update(self, i_f=contigs, o_d=out_dir, key=pool)


def spades_cmd(fastas: list, params: dict, tmp: str, out: str):
    """Create command lines for spades.

    Parameters
    ----------
    fastas : list
        Input fasta/fastq.gz files
    params : dict
        Parameters for Spades
    tmp : str
        Temporary folder for Spades
    out : str
        Output folder for Spades
    """
    cmd = 'spades.py'
    cmd += ' -m %s' % params['mem']
    cmd += ' -k %s' % ','.join(map(str, params['k']))
    cmd += ' -t %s' % params['cpus']
    if params['bio']:
        cmd += ' --bio'
    if params['meta']:
        cmd += ' --meta'
    if params['plasmid']:
        cmd += ' --plasmid'
    if params['only_assembler']:
        cmd += ' --only-assembler'
    cmd += ' --tmp-dir %s -o %s' % (tmp, out)
    for fasta in fastas[:3]:
        if 'extendedFrags' in fasta:
            cmd += ' --merge %s' % fasta
        elif 'notCombined_1' in fasta:
            cmd += ' -1 %s' % fasta
        elif 'notCombined_2' in fasta:
            cmd += ' -2 %s' % fasta
    return cmd


def spades(self) -> None:
    """Create command lines for spades.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
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
    for group, fastas in self.inputs[self.pool].items():
        tmp = '$TMPDIR/spades_%s_%s' % (self.pool, group)
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = 'mkdir -p %s\n' % tmp
        cmd += '%s\n' % spades_cmd(fastas, self.soft.params, tmp, out)
        cmd += 'rm -rf %s\n' % tmp
        before_rr = '%s/before_rr.fasta' % out
        contigs = '%s/contigs.fasta' % out
        first_pe = '%s/first_pe_contigs.fasta' % out
        scaffolds = '%s/scaffolds.fasta' % out
        log = '%s/spades.log' % out
        outs = [before_rr, contigs, first_pe, scaffolds, log]
        if self.config.force or todo(contigs):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(
                self, i_f=fastas[:3], i_d=tmp, o_f=outs, o_d=out, key=group)
        self.outputs['outs'][group] = outs


def viralverify(self):
    if self.soft.prev != 'spades':
        sys.exit('[viralVerify] can only be run on assembly output')
    pfam = '%s/Pfam-A.hmm' % self.databases.paths.get('pfam')
    if not self.config.dev and todo(pfam):
        sys.exit('[viralVerify] Needs the Pfam .hmm database in database yaml')
    for group, spades_outs in self.inputs[self.pool].items():
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = 'export PATH=$PATH:%s' % self.soft.params['path']
        cmd += '\nviralverify'
        cmd += ' -f %s' % spades_outs[1]
        cmd += ' -o %s' % out
        if self.soft.params['db']:
            cmd += ' --db'
        if self.soft.params['p']:
            cmd += ' --p'
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -thr %s' % self.soft.params['thr']
        cmd += ' --hmm %s' % pfam
        out_fp = '%s_result_table.csv' % splitext(spades_outs[1])[0]
        if self.config.force or todo(out_fp):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out, key=group)
        self.outputs['outs'][group] = out_fp


def plass(self):
    """Create command lines for Plass.

    self : Commands class instance
        .dir : str
            Path to pipeline output folder for plass
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
    if self.config.force or todo(out_fp):
        cmd = 'mkdir -p %s\n%s\nrm -rf %s' % (tmp_dir, cmd, tmp_dir)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, i_d=tmp_dir, o_f=out_fp)


def get_merges(self) -> list:
    """Get the extension of the file to merge to make the pools.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev
            Software run just before pooling

    Returns
    -------
    merges : list
        Extension of the outputs to merge
    """
    merges = ['']
    if self.soft.prev == 'flash':
        merges = ['extendedFrags', 'notCombined_1', 'notCombined_2']
    return merges


def pool_cmd(self, pool: str, paths: list,
             fasta: str, group: str) -> None:
    """Write the pooling command and collect the output and io for FLASh.

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    pool : str
        Name of the pool
    paths : list
        Path the input fasta files
    fasta : str
        Path to an output fasta file
    group : str
        Name of the sample group within the pool
    """
    if self.config.force or todo(fasta):
        for pdx, path in enumerate(paths):
            if pdx:
                cmd = 'cat %s >> %s' % (path, fasta)
            else:
                cmd = 'cat %s > %s' % (path, fasta)
            if pool not in self.cmds:
                self.cmds[pool] = {}
            self.cmds[pool].setdefault(group, []).append(cmd)


def get_pools(self, pool: str, group: str, sams: list) -> list:
    """Get the output fasta files returned by FLASh and write the commands.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for pools
        .inputs : dict
            Input files
        .cmds
            Command lines
        .soft.dirs : set
            All directories to create
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    pool : str
        Name of the pool
    group : str
        Name of the sample group within the pool
    sams : list
        Samples to pool

    Returns
    -------
    fasta_fps : list
        Paths to the output fasta files
    """
    out_dir = self.dir + '/' + pool
    self.soft.dirs.add(out_dir.replace('${SCRATCH_FOLDER}', ''))
    fasta_fps = []
    print()
    print(self.inputs)
    print()
    print(self.longs)
    for merge in get_merges(self):
        paths = [fp for sam in sams for fp in self.inputs[sam] if merge in fp]
        print()
        print(paths)
        print(selfinputs)
        # only pool if there is min 2 samples being merged
        if len(sams) > 1:
            fasta = '%s/%s.%s.fasta' % (out_dir, group, merge)
            fasta = fasta.replace(' ', '_').replace('..', '.')
            self.soft.io[(pool, group)][('I', 'f')].update(paths)
            pool_cmd(self, pool, paths, fasta, group)
        else:
            if len(paths) > 1:
                sys.exit('[pooling] Error in co-assembly %s:%s' % (pool, group))
            fasta = paths[0]
        fasta_fps.append(fasta)
    self.soft.io[(pool, group)][('O', 'd')].add(out_dir)
    self.soft.io[(pool, group)][('O', 'f')].update(fasta_fps)
    return fasta_fps


def pooling(self, pool):
    """Create command lines for pooling samples' reads using FLASh.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    pool : str
        Name of the pool
    """
    for group, group_pd in self.config.meta.groupby(pool):
        # get the full list of samples and the samples per pooling group
        sams = group_pd.index.tolist()
        self.pools[pool][group] = sams
        self.soft.io[(pool, group)] = {('I', 'd'): set(), ('I', 'f'): set(),
                                       ('O', 'd'): set(), ('O', 'f'): set()}
        # get the outputs for the current group and collect pooling commands
        self.soft.outputs[pool][group] = get_pools(self, pool, group, sams)


def flye_cmd(fastx: list, params: dict, out: str):
    """Create command lines for flye.

    Parameters
    ----------
    fastx : list
        Input read fasta/fastq.gz files
    params : dict
        Parameters for flye
    out : str
        Output folder for flye
    """
    cmd = 'flye'
    cmd += ' --%s %s' % (params['long_reads'], ' '.join(fastx))
    cmd += ' --out-dir %s' % out
    cmd += ' --threads %s' % params['cpus']
    if 'genome_size' in params:
        cmd += ' --genome-size %s' % params['genome_size']
    for boolean in ['meta', 'keep_haplotypes', 'scaffold']:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    for param in ['iterations', 'min_overlap', 'asm_coverage',
                  'read_error', 'polish_target', 'resume']:
        if params[param]:
            cmd += ' --%s %s\n' % (param.replace('_', '-'), params[param])
    return cmd


def flye(self) -> None:
    """Create command lines for flye.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for flye
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
    for group, fastx in self.inputs[self.pool].items():
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = flye_cmd(fastx, self.soft.params, out)
        contigs = '%s/assembly.fasta' % out
        gfa = '%s/assembly_graph.gfa' % out
        gv = '%s/assembly_graph.gv' % out
        info = '%s/assembly_info.txt' % out
        self.outputs['outs'][group] = [info, contigs, gfa, gv]
        if self.config.force or todo(contigs):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=fastx, i_d=out, o_d=out, key=group)


def canu_cmd(fastx: list, params: dict, group: str, out: str):
    """Create command lines for flye.

    Parameters
    ----------
    fastx : list
        Input read fastx files
    params : dict
        Parameters for flye
    group : str
        Name of the sample group within the pool
    out : str
        Output folder for flye
    """
    cmd = 'if [ -d %s]; then rm -rf %s; fi\n' % (out, out)
    if 'path' in params:
        cmd += 'export PATH=$PATH:%s\n' % params['path']
    cmd += 'canu'
    cmd += ' -%s %s' % (params['technology'], ' '.join(fastx))
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
        if 'pacbio' in params['technology']:
            cmd += ' rawErrorRate=0.300'
        else:
            cmd += ' rawErrorRate=0.500'
    if params['correctedErrorRate']:
        cmd += ' correctedErrorRate=%s' % params['correctedErrorRate']
    else:
        if 'pacbio' in params['technology']:
            cmd += ' correctedErrorRate=0.045'
        else:
            cmd += ' correctedErrorRate=0.144'
    cmd += ' minReadLength=%s' % params['minReadLength']
    cmd += ' minOverlapLength=%s' % params['minOverlapLength']
    for threads in ['bat', 'cns', 'cormhap', 'cormmap', 'corovl', 'cor',
                    'executive', 'hap', 'max', 'meryl', 'min', 'obtmhap',
                    'obtmmap', 'obtovl', 'oea', 'ovb', 'ovs', 'red',
                    'utgmhap', 'utgmmap', 'utgovl']:
        cmd += ' %sThreads=%s' % (threads, params['cpus'])
    return cmd


def canu(self) -> None:
    """Create command lines for flye.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for flye
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
    for group, fastx in self.inputs[self.pool].items():
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        report = '%s/%s.report' % (out, group)
        contigs = '%s/%s.contigs.fasta' % (out, group)
        self.outputs['outs'][group] = [report, contigs]
        if self.config.force or todo(contigs):
            cmd = canu_cmd(fastx, self.soft.params, group, out)
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=fastx, o_d=out, key=group)


def unicycler_cmd(fastqs: list, params: dict, tmp: str, out: str, long=[]):
    """Create command lines for unicycler.

    Parameters
    ----------
    fastqs : list
        Input read fastqs files
    params : dict
        Parameters for flye
    tmp : str
        Path to a temporary directory
    out : str
        Output folder for unicycler
    """
    cmd = 'unicycler'
    if len(fastqs) == 2:
        cmd += ' --short1 %s --short2 %s' % tuple(fastqs)
    # --unpaired
    # --long
    cmd += ' --out %s' % out
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
        if param in params:
            cmd += ' --%s %s' % (param, params[param])
    for boolean in ['vcf', 'no_correct', 'largest_component',
                    'no_miniasm', 'no_rotate', 'no_pilon']:
        if params[boolean]:
            cmd += ' --%s' % params[param]
    return cmd


def get_fastqs(self):
    if self.pool in self.pools:
        fastqs = self.inputs[self.pool]
    else:
        fastqs = {self.sam: self.inputs[self.sam]}
    return fastqs


def get_longs(self):
    if self.pool in self.pools:
        groups = self.inputs[self.pool]
        longs = {group: [self.longs[sam] for sams
                         in self.pools[group] for sam in sams]
                 for group in groups}
    else:
        longs = {self.sam: self.long[self.sam]}
    return longs


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
    fastqs = get_fastqs(self)
    longs = get_longs(self)
    for sam_group, fqs in fastqs.items():
        tmp = '$TMPDIR/unicycler_%s_%s' % (self.pool, sam_group)
        out = '%s/%s/%s' % (self.dir, self.pool, sam_group)
        report = '%s/%s.' % (out, sam_group)
        contigs = '%s/%s.' % (out, sam_group)
        if self.pool in self.pools:
            self.outputs['outs'][sam_group] = [report, contigs]
        else:
            self.outputs['outs'] = [report, contigs]
        if self.config.force or todo(contigs):
            cmd = unicycler_cmd(fqs, self.soft.params, tmp, out, long)
            if self.pool in self.pools:
                self.outputs['cmds'].setdefault(sam_group, []).append(cmd)
                io_update(self, i_f=fqs, i_d=tmp, o_d=out, key=sam_group)
            else:
                self.outputs['cmds'].setdefault(sam_group, []).append(cmd)
                io_update(self, i_f=fqs, i_d=tmp, o_d=out)


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

# USE CALLER LIKE METAWRAP
# def circlator_cmd(fastqs: list, params: dict, tmp: str, out: str):
#     """Create command lines for circlator.
#
#     Parameters
#     ----------
#     fastqs : list
#         Input read fastqs files
#     params : dict
#         Parameters for flye
#     tmp : str
#         Path to a temporary directory
#     out : str
#         Output folder for flye
#     """
#     cmd = 'circlater %s' % analysis
#     return cmd
#
#
# def circlator(self) -> None:
#     """Create command lines for circlator.
#
#     Parameters
#     ----------
#     self : Commands class instance
#         .dir : str
#             Path to pipeline output folder for circlator
#         .pool : str
#             Pool name.
#         .inputs : dict
#             Input files
#         .outputs : dict
#             All outputs
#         .soft.params
#             Parameters
#         .config
#             Configurations
#     """
#     if self.pool in self.pools:
#         fastqs_d = self.inputs[self.pool]
#     else:
#         fastqs_d = {self.sam: self.inputs[self.sam]}
#     for sam_group, fastqs in fastqs_d.items():
#         tmp = '$TMPDIR/spades_%s_%s' % (self.pool, sam_group)
#         out = '%s/%s/%s' % (self.dir, self.pool, sam_group)
#         report = '%s/%s.' % (out, sam_group)
#         contigs = '%s/%s.' % (out, sam_group)
#         if self.pool in self.pools:
#             self.outputs['outs'][sam_group] = [report, contigs]
#         else:
#             self.outputs['outs'] = [report, contigs]
#         if self.config.force or todo(contigs):
#             cmd = circlator_cmd(fastqs, self.soft.params, tmp, out)
#             if self.pool in self.pools:
#                 self.outputs['cmds'].setdefault(sam_group, []).append(cmd)
#                 io_update(self, i_f=fastqs, i_d=tmp, o_d=out, key=sam_group)
#             else:
#                 self.outputs['cmds'].setdefault(sam_group, []).append(cmd)
#                 io_update(self, i_f=fastqs, i_d=tmp, o_d=out)
