# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
import pkg_resources
from os.path import basename, splitext
from metagenomix._io_utils import io_update, status_update
from metagenomix.core.parameters import tech_params
from metagenomix._io_utils import to_do
from metagenomix._inputs import (
    group_inputs, genome_key, genome_out_dir, get_reads, get_group_reads)
from metagenomix.softwares.alignment import (
    bowtie2_cmd,
    minimap2_cmd,
    # bbmap_cmd,
    # bwa_cmd
)

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def get_mapping_target(self) -> tuple:
    source = self.soft.name.split('_', 1)[-1]
    if '_' not in self.soft.name or source == 'mapping':
        sys.exit('[mapping] Add "_<soft>" to specify what to map to "%s" data' %
                 self.soft.prev)
    step = self.config.tools.get(source)
    if step == 'preprocessing':
        func = raw
    elif step in ['paired-read merging', 'coassembly setup']:
        func = merged
    elif step in ['MAG', 'binning']:
        func = genomes
    elif step == 'assembling':
        func = assembly
    else:
        sys.exit('[mapping] Not possible to map "%s" on "%s"' % (
            source, self.soft.prev))
    return func, source, step


# def bwa():
#     cmd = 'bwa index -p %s/index %s\n' % (out_dir, contigs)
#     cmd += 'bwa mem -p %s/index %s' % (out_dir, fastqs)
#     cmd += ' | samtools view -bh'
#     cmd += ' | samtools sort -o %s/%s\n' % (out_dir, bam)
#     cmd += 'samtools index %s %s\n' % (bam, bai)


def assembly(self, func, fastas, tech, group):
    pass


def merged(self, func, fastas, tech, group):
    pass


def genomes(self, func, fastas, tech, group):
    pass


def get_minimap2_db_cmd(
        fasta: str,
        sam_tech_dir: str,
        params: dict
) -> tuple:
    """Commands to build a minimap2 database indices for the current fasta
    reference and paths to these indices.

    Parameters
    ----------
    fasta : str
    sam_tech_dir : str
    params : dict

    Returns
    -------
    cmd : str
        minimap2-build command
    cmd_rm : str
        minimap2-build command removers
    dbs : dict
        Databases indices
    """
    cmd, cmd_rm = '', ''
    dbs = {}
    return cmd, cmd_rm, dbs


def get_bowtie2_db_cmd(
        out: str,
        fastas: list,
        params: dict,
) -> tuple:
    """Commands to build a Bowtie2 database indices for the
    current fasta reference and paths to these indices.

    Parameters
    ----------
    out : str
    fastas : list
    params : dict

    Returns
    -------
    cmd : str
        bowtie2-build command
    cmd_rm : str
        bowtie2-build command removers
    dbs : dict
        Databases indices
    """
    dbs = {}
    cmd, cmd_gz, cmd_rm = '', '', ''
    for fasta_ in fastas:
        fasta = fasta_
        if fasta.endswith('.gz'):
            cmd_gz += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
            fasta = fasta.rstrip('.gz')
            cmd_rm += 'rm %s\n' % fasta

        base = splitext(basename(fasta))[0]
        bam_dir = out + '/' + base
        bam = '%s/alignment.bowtie2.sorted.bam' % bam_dir
        dbs[base] = (fasta, bam, bam_dir)
        if to_do(bam):
            cmd += 'mkdir -p %s/%s\n' % (out, base)
            cmd += 'mkdir -p %s/dbs/%s\n' % (out, base)
            cmd += 'bowtie2-build'
            cmd += ' --threads %s' % params['cpus']
            cmd += ' %s %s/dbs/%s\n' % (fasta, out, base)
    if cmd:
        cmd = cmd_gz + cmd
        cmd_rm += '\nrm -rf %s/dbs\n' % out
    return cmd, cmd_rm, dbs


def get_cmds(
        self,
        ref_group: str,
        reads_tech: str,
        sam: str,
        out: str,
        fastqs: list,
        fastas: list,
        ali: str,
        params: dict
) -> tuple:
    """

    Parameters
    ----------
    self
    ref_group : str
    reads_tech : str
    sam : str
    out : str
    fastqs : list
    fastas : list
    ali : str
    params : dict

    Returns
    -------
    cmd : str
    bams : list
    bam_dirs : list
    fastas_bams : dict
    """
    bams, bam_dirs, fastas_bams = [], [], {}
    cmd, cmd_rm, dbs = globals()['get_%s_db_cmd' % ali](out, fastas, params)
    for db, (fasta, bam, bam_dir) in dbs.items():
        db_path = '%s/dbs/%s' % (out, db)
        if to_do(bam):
            cmd += globals()['%s_cmd' % ali](
                sam, fastqs, db, db_path, bam_dir, params)
            cmd += ' | samtools view -b -'
            cmd += ' | samtools sort -o %s\n' % bam
        else:
            bams.append(bam)
        if to_do('%s.bai' % bam):
            cmd += 'samtools index %s\n' % bam
        fastas_bams[bam] = [
            reads_tech, sam, ali, self.sam_pool, ref_group, fasta]
        bam_dirs.append(bam_dir)
    if cmd:
        cmd += cmd_rm
    return cmd, bams, bam_dirs, fastas_bams


def raw(
        self,
        source: str,
        ref_tech: str,
        ref_group: str,
        reads: dict,
        fastas: list,
        key: list,
        out_dir: str,
        to_dos: list
) -> None:
    """

    Parameters
    ----------
    self
    source : str
    ref_tech : str
    ref_group : str
    reads : dict
    fastas : list
    key : list
    out_dir : str
    to_dos: list
    """
    for sam, reads_tech_fastqs in reads.items():
        for (reads_tech, _), fastqs in reads_tech_fastqs.items():
            reads_to_dos = status_update(self, reads_tech, fastqs)
            for ali in self.soft.params['aligners']:
                params = tech_params(self, reads_tech, ali)
                out = '/'.join([out_dir, sam, reads_tech, ali])
                cmd, bams, bam_dirs, fastas_bams = get_cmds(
                    self, ref_group, reads_tech, sam, out, fastqs, fastas,
                    ali, params)
                self.outputs['dirs'].extend(bam_dirs)
                self.outputs['outs'].update(fastas_bams)
                to_do_list = [to_do(x) for x in fastas_bams.keys()]
                if self.config.force or sum(to_do_list):
                    if to_dos or reads_to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                    if bams:
                        io_update(self, i_f=bams, o_d=out, key=key)
                    else:
                        io_update(self, i_f=(fastqs + fastas), o_d=out, key=key)
                    self.soft.add_status(
                        ref_tech, self.sam_pool, 1, group=ref_group)
                else:
                    self.soft.add_status(
                        ref_tech, self.sam_pool, 0, group=ref_group)


def get_mapping(
        self,
        func,
        source: str,
        ref_tech: str,
        ref_group: str,
        reads: dict,
        references: dict,
) -> None:
    """

    Parameters
    ----------
    self
    func
    source : str
    ref_tech : str
    ref_group : str
    reads : dict
    references : dict
    """
    for genome, fastas in references.items():
        key = genome_key(ref_tech, ref_group, genome)
        out_dir = genome_out_dir(self, ref_tech, ref_group, genome)
        to_dos = status_update(
            self, ref_tech, fastas, group=ref_group, genome=genome)
        func(self, source, ref_tech, ref_group,
             reads, fastas, key, out_dir, to_dos)


def mapping(self):
    """Mapping would be a rather specific process that generically consists
    of aligning reads that can be raw, filtered, or for a specific scope
    (e.g., that were used to re-assemble a MAGs or specific proteins).

    It can be done using different aligners, such as BWA, bowtie2, minimap2,
    bbmap, or others, but all should yield a .sam file that can be piped in
    samtools in order to obtain one .bam file and one .bai file per mapping.
    These files are then possibly used by several other softwares.

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    func, source, step = get_mapping_target(self)
    all_reads = get_reads(self, soft=source)
    if self.sam_pool in self.pools:
        for (ref_tech, ref_group), inputs in self.inputs[self.sam_pool].items():
            refs = group_inputs(self, inputs)
            reads = get_group_reads(self, ref_tech, ref_group, all_reads)
            get_mapping(self, func, source, ref_tech, ref_group, reads, refs)


def get_pysam_target(self) -> tuple:

    if self.soft.name.startswith('pysam_'):
        target = self.soft.name.split('_', 1)[1]
        if target not in self.softs:
            sys.exit("[%s] %s has not been run" % (self.soft.name, target))
    else:
        target = self.softs[self.soft.prev].prev

    classif = self.config.tools[target].split()[0]
    if classif == 'assembling':
        func = assembling
    elif classif == 'binning':
        func = binning
    elif classif == 'MAG':
        func = mag
    elif classif == 'annotation':
        func = annotation
    else:
        sys.exit("[%s] Counting '%s' not possible" % (self.soft.name, target))

    return target, func


def binning(
        target: str,
        prev: str,
        maps: dict,
        fas: list,
        out_dir: str,
) -> tuple:
    """

    Parameters
    ----------
    target
    prev
    maps
    fas
    out_dir

    Returns
    -------
    cmd : str
        pysam commands for bins sequence quantification
    """
    pass


def mag(
        target: str,
        prev: str,
        maps: dict,
        fas: list,
        out_dir: str,
) -> str:
    """

    Parameters
    ----------
    target
    prev
    maps
    fas
    out_dir

    Returns
    -------
    cmd : str
        pysam commands for MAGs quantification
    """
    pass


def annotation(
        target: str,
        prev: str,
        maps: dict,
        fas: list,
        out_dir: str,
) -> tuple:
    """

    Parameters
    ----------
    target
    prev : str
    maps : dict
    fas : list
    out_dir : str

    Returns
    -------
    cmd : str
        pysam commands for annotated sequence quantification
    sam_bams : list
        per sample bam files for the current contigs counting
    aas : list
        Input files
    """
    fastas = {}
    aas, sam_bams = [], []
    cmd, cmd_rm = '', ''
    for fa in fas:
        aa = '%s/protein.translations.fasta.gz' % fa
        aas.append(aa)
        if aa.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (aa, aa.rstrip('.gz'))
            aa = aa.rstrip('.gz')
            coassembly, group = fa.rstrip('/').rsplit('/', 2)[-2:]
            bams = [[x] + y[:-1] for x, y in maps.items()
                    if y[-3] == coassembly and y[-2] == group]
            sam_bams.extend([x[0] for x in bams])
            sam_bams.extend(['%s.bai' % x[0] for x in bams])
            fastas[aa] = bams

    cmd += get_pysam_cmd(prev, target, fastas, out_dir, 'annotation')
    cmd += cmd_rm
    return cmd, sam_bams, aas


def assembling(
        target: str,
        prev: str,
        maps: dict,
        fas: list,
        out_dir: str,
) -> tuple:
    """

    Parameters
    ----------
    target : str
    prev : str
    maps : str
    fas : list
    out_dir : str

    Returns
    -------
    cmd : str
        pysam commands for the assembly sequence quantification
    sam_bams : list
        per sample bam files for the current contigs counting
    fas : list
        Input files
    """
    fastas = {}
    sam_bams = []
    cmd, cmd_rm = '', ''
    for fa in fas:
        if fa.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (fa, fa.rstrip('.gz'))
            fa = fa.rstrip('.gz')
            cmd_rm += 'rm %s\n' % fa
            bams = [[x] + y[:-1] for x, y in maps.items() if y[-1] == fa]
            sam_bams.extend([x[0] for x in bams])
            sam_bams.extend(['%s.bai' % x[0] for x in bams])
            fastas[fa] = bams

    cmd += get_pysam_cmd(prev, target, fastas, out_dir, 'assembling')
    cmd += cmd_rm
    return cmd, sam_bams, fas


def get_pysam_cmd(
        prev: str,
        target: str,
        fastas: dict,
        out_dir: str,
        classif='assembling'
) -> str:
    """

    Parameters
    ----------
    prev : str
        Name of the previous software
    target : str
    fastas : dict
        Pairs of fastas files to counts per reference bam file
    out_dir : str
        Output folder
    classif : str
        target to count the reads for ("assembling", "annotation")

    Returns
    -------
    cmd : str
        pysam counting command lines
    """
    cmd = get_pysam_inputs(prev, target, fastas, out_dir, classif)
    cmd += 'python3 %s/pysam_count.py' % RESOURCES
    cmd += ' -i %s/inputs.txt' % out_dir
    cmd += ' -o %s/reads.txt\n' % out_dir
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    return cmd


def get_pysam_inputs(
        prev: str,
        target: str,
        fastas: dict,
        out_dir: str,
        classif: str
) -> str:
    """

    Parameters
    ----------
    prev : str
    target : str
    fastas : dict
    out_dir : str
    classif : str

    Returns
    -------
    cmd : str
        Commands to make the inputs file
    """
    echo = 'fasta\\tbam\\ttech\\tsample\\tali'
    echo += '\\tcoassembly\\tgroup\\tprev\\ttarget\\tclassif'
    cmd = 'echo -e "%s" > %s/inputs.txt\n' % (echo, out_dir)
    for fa, bams in fastas.items():
        for bs in bams:
            echo = '%s' % '\\t'.join(([fa] + bs + [prev, target, classif]))
            cmd += 'echo -e "%s" >> %s/inputs.txt\n' % (echo, out_dir)
    cmd += 'envsubst < %s/inputs.txt > %s/inputs.tmp\n' % (out_dir, out_dir)
    cmd += 'mv %s/inputs.tmp %s/inputs.txt\n' % (out_dir, out_dir)
    return cmd


def pysam_cmd(
        self,
        tech: str,
        group: str,
        fas: list,
        sam_bams: list,
        key: tuple,
        to_dos: list,
        out_dir: str,
        out: str,
        cmd: str
) -> None:
    if self.config.force or to_do(out):
        if to_dos and not self.config.dev:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_f=(fas + sam_bams), o_d=out_dir, key=key)
        self.soft.add_status(tech, self.sam_pool, 1, group=group)
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group=group)


def get_pysam(
        self,
        func,
        mappings,
        tech: str,
        group: str,
        target: str,
        references: dict
) -> None:
    """

    Parameters
    ----------
    self
    func
        Function to use to get the mapping couting adapted to:
            assembly data: func "assembling"
            binning data: func "binning"
            MAGs data: func "mag"
            annotation data: func "annotation"
    mappings
        Name of the mapping softwares
    tech : str
        Name of the technology
    group : str
        Name of the co-assembly group
    target : str
        Name of the reads to map
    references : dict
        Name of the references to map against
    """
    prev = mappings.prev
    maps = mappings.outputs[self.sam_pool]
    for genome, fas in references.items():
        key = genome_key(tech, group, genome)
        out_dir = '/'.join([genome_out_dir(self, tech, group, genome),
                            'map_%s' % prev, 'count_%s' % target])
        out = '%s/reads.txt.gz' % out_dir
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault(key, []).append(out)

        cmd, bams, fs = func(target, prev, maps, fas, out_dir)
        to_dos = status_update(self, tech, bams, group=group, genome=genome)
        to_dos.extend(status_update(self, tech, fs, group=group, genome=genome))
        pysam_cmd(self, tech, group, fs, bams, key, to_dos, out_dir, out, cmd)


def pysam(self):
    """Uses a pysam python script to count the reads aligned using a mapping_
    command for each of the mapped-onto sequences. Whether these mapped-onto
    sequences are from a single file (typically, the contigs of an assembly),
    or from multiple files (typically, contigs binned in different MAGs), is
    detected automatically and all counts are consigned into a table including
    the file information.

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    if not self.soft.prev.startswith('mapping'):
        sys.exit("[%s] Only run after a mapping_* command" % self.soft.name)
    mappings = self.softs[self.soft.prev]
    target, func = get_pysam_target(self)
    if self.sam_pool in self.pools:
        pool_inputs = self.softs[target].outputs[self.sam_pool]
        for (tech, group), inputs in pool_inputs.items():
            references = group_inputs(self, inputs, target=target)
            get_pysam(self, func, mappings, tech, group, target, references)


def get_salmon_reads_cmd(
        self,
        tech: str,
        fastx: str,
        fastqs: list,
        out_dir: str,
        key: tuple
) -> str:
    """Get the Salmon quant command (using reads).

    Parameters
    ----------
    self
    tech : str
    fastx : str
    fastqs : list
    out_dir : str
    key: tuple

    Returns
    -------
    cmd : str
        Salmon quant command (using reads)
    """
    params = tech_params(self, tech)
    cmd = 'salmon quant'
    cmd += ' --libType %s'
    cmd += ' --index %s' % fastx
    if len(fastqs) == 3:
        cmd += ' --unmatedReads %s' % fastqs[0]
        cmd += ' --mates1 %s' % fastqs[1]
        cmd += ' --mates2 %s' % fastqs[2]
    elif len(fastqs) == 2:
        cmd += ' --mates1 %s' % fastqs[0]
        cmd += ' --mates2 %s' % fastqs[1]
    elif len(fastqs) == 1:
        cmd += ' --unmatedReads %s' % fastqs[0]
    cmd += ' --output %s' % out_dir

    for param in [
        'sigDigits', 'thinningFactor', 'numBootstraps', 'numGibbsSamples',
        'rangeFactorizationBins', 'numBiasSamples', 'numAuxModelSamples',
        'numPreAuxModelSamples', 'maxOccsPerHit', 'maxReadOcc',
        'maxRecoverReadOcc', 'minAssignedFrags', 'reduceGCMemory',
        'biasSpeedSamp', 'fldMax', 'fldMean', 'fldSD', 'minAssignedFrags',
        'kmerLen', 'filterSize', 'scoreExp', 'numErrorBins',
        'mappingCacheMemoryLimit', 'mismatchSeedSkip', 'ma', 'mp', 'go',
        'ge', 'bandwidth', 'vbPrior', 'minAlnProb', 'decoyThreshold',
        'consensusSlack', 'preMergeChainSubThresh', 'postMergeChainSubThresh',
        'orphanChainSubThresh', 'minScoreFraction', 'incompatPrior',
        'forgettingFactor'
    ]:
        cmd += ' --%s %s' % (param, params[param])
    return cmd


def get_salmon_sam_cmd(
        self,
        tech: str,
        fastx: str,
        sam: str,
        key: tuple
) -> str:
    """Get the Salmon command.

    Parameters
    ----------
    self
    tech : str
    fastx : str
    sam : str
    key : tuple

    Returns
    -------
    cmd : str
        Salmon quant command
    """
    params = tech_params(self, tech)
    cmd = 'salmon quant'
    for param in [
        'sigDigits', 'thinningFactor', 'numBootstraps', 'numGibbsSamples',
        'rangeFactorizationBins', 'numBiasSamples', 'numAuxModelSamples',
        'numPreAuxModelSamples', 'maxOccsPerHit', 'maxReadOcc',
        'maxRecoverReadOcc', 'minAssignedFrags', 'reduceGCMemory',
        'biasSpeedSamp', 'fldMax', 'fldMean', 'fldSD', 'minAssignedFrags',
        'kmerLen', 'filterSize', 'scoreExp', 'numErrorBins',
        'mappingCacheMemoryLimit', 'mismatchSeedSkip', 'ma', 'mp', 'go',
        'ge', 'bandwidth', 'vbPrior', 'minAlnProb', 'decoyThreshold',
        'consensusSlack', 'preMergeChainSubThresh', 'postMergeChainSubThresh',
        'orphanChainSubThresh', 'minScoreFraction', 'incompatPrior',
        'forgettingFactor'
    ]:
        cmd += ' --%s %s' % (param, params[param])
    return cmd


def get_salmon_indexing_cmd(
        self,
        tech: str,
        fasta: str,
        fastx: str,
        key: tuple
) -> str:
    """Get the Salmon indexing command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for humann
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fasta : str
        Path to the reference fasta file
    fastx : str
        Path to the indexed reference fasta file
    key : tuple
        Current input unit keys

    Returns
    -------
    cmd : str
        Salmon indexing command
    """
    tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, self.sam_pool, '_'.join(key))
    params = tech_params(self, tech)
    cmd = 'salmon index'
    cmd += ' --transcripts %s' % fasta
    cmd += ' --index %s' % fastx
    cmd += ' --threads %s' % params['cpus']
    for boolean in [
        'gencode', 'features', 'keepDuplicates', 'keepFixedFasta', 'sparse'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean
    for param in [
        'kmerLen', 'filterSize', 'decoys', 'type'
    ]:
        cmd += ' --%s %s' % (param, params[param])
    cmd += ' --tmpdir %s' % tmp
    return cmd


def mapdamage2_cmd(
        self,
        tech: str,
        bam: str,
        contigs: str,
        out_dir: str
) -> str:
    """Collect the command line for filtering.

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
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        filtering commands
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    if contigs.endswith('.fa.gz') or contigs.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % contigs.rstrip('.gz')
        contigs = contigs.rstrip('.gz')

    cmd += '\nmapDamage'
    cmd += ' --input=%s' % bam
    cmd += ' --reference=%s' % contigs
    cmd += ' --folder=%s' % out_dir
    for param in [
        'downsample_seed', 'length', 'around', 'min_basequal', 'ymax',
        'readplot', 'refplot', 'rand', 'burn', 'adjust', 'iter', 'seq_length',
        'rescale_length_5p', 'rescale_length_3p'
    ]:
        cmd += ' --%s=%s' % (param.replace('_', '-'), params[param])

    for boolean in [
        'merge_reference_sequences', 'fasta', 'plot_only', 'quiet', 'verbose',
        'mapdamage_modules', 'forward', 'reverse', 'var_disp', 'jukes_cantor',
        'diff_hangs', 'fix_nicks', 'use_raw_nick_freq', 'single_stranded',
        'theme_bw', 'stats_only', 'no_stats', 'check_R_packages', 'rescale',
        'rescale_only'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    if params['downsample']:
        cmd += ' --downsample %s' % params['downsample']
    cmd += '\n%s\n' % cmd_rm
    return cmd


def get_mapdamage2(
        self,
        bam: str,
        bam_infos: list
) -> None:
    """Get the info and collect the commasnd for mapDamage2.

    Parameters
    ----------
    self
    bam : str
        Path to the input BAM file
    bam_infos : list
        [tech, sample, aligner, coassembly, coassembly group, contig path]
    """
    tech, sample, aligner, _, group, contigs = bam_infos
    out_dir = genome_out_dir(self, tech, group) + '/' + aligner
    self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
    self.outputs['dirs'].append(out_dir)

    contigs_gz = contigs + '.gz'
    to_dos = status_update(
        self, tech, [bam, contigs_gz], self.sam_pool, group=group)

    key = genome_key(tech, group, aligner)
    pdf = '%s/Fragmisincorporation_plot.pdf' % out_dir
    if self.config.force or to_do(pdf):
        if to_dos:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            cmd = mapdamage2_cmd(self, tech, bam, contigs_gz, out_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_f=[bam, contigs_gz], o_d=out_dir, key=key)
        self.soft.add_status(tech, self.sam_pool, 1, group=group)
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group=group)


def mapdamage2(self) -> None:
    """.

    Notes
    -----
    GitHub  : https://ginolhac.github.io/mapDamage
    paper   : https://doi.org/10.1093/bioinformatics/btt193

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MapDamage2
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
    if not self.soft.prev.startswith('mapping'):
        sys.exit("[%s] Only run after a mapping_* command" % self.soft.name)
    if self.sam_pool in self.pools:
        for bam, bam_infos in self.inputs[self.sam_pool].items():
            get_mapdamage2(self, bam, bam_infos)


# def salmon(self) -> None:
#     """Salmon is a tool for wicked-fast transcript quantification from
#     RNA-seq data. It requires a set of target transcripts (either from a
#     reference or de-novo assembly) to quantify. All you need to run Salmon
#     is a FASTA file containing your reference transcripts and a (set of)
#     FASTA/FASTQ file(s) containing your reads. Optionally, Salmon can make
#     use of pre-computed alignments (in the form of a SAM/BAM file) to the
#     transcripts rather than the raw reads.
#
#     References
#     ----------
#     Patro, Rob, et al. "Salmon provides fast and bias-aware quantification of
#     transcript expression." Nature methods 14.4 (2017): 417-419.
#
#     Notes
#     -----
#     GitHub  : https://github.com/COMBINE-lab/salmon
#     Docs    : https://salmon.readthedocs.io/en/latest/salmon.html
#     Paper   : https://doi.org/10.1038/nmeth.4197
#
#     Parameters
#     ----------
#     self
#     """
#     if self.soft.params['useAlignments']:
#         print(self.config.__dict__.keys())
#         print(self.config.softs)
#         print(selfconfigsofts)
#
#     if self.sam_pool in self.pools:
#         for (tech, group), inputs in self.inputs[self.sam_pool].items():
#             reads = get_reads(self, tech, group)
#             references = group_inputs(self, inputs)
#             get_mapping(self, func, tech, group, reads, references)
#             # get_salmon_indexing_cmd(self, tech, fasta, fastx, key)
#
#     idx = '%s.si' % (out_dir, basename(splitext(fasta)[0]))
#
#     for (tech, sample), fastxs in self.inputs[self.sam_pool].items():
#         if tech_specificity(self, fastxs, tech, sample):
#             continue
#
#         out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
#         self.outputs['outs'][(tech, self.sam_pool)] = dict()
#
#         for db, db_path in self.soft.params['databases'].items():
#
#             db_out = '%s/%s' % (out, db)
#             params = tech_params(self, tech)
#             cmd, sam = get_minimap2_cmd(fastxs, db_path, db_out, params)
#             self.outputs['outs'][(tech, self.sam_pool)][(db, 'minimap2')] = sam
#
#             if self.config.force or to_do(sam):
#                 if status_update(self, tech, fastxs):
#                     self.outputs['cmds'].setdefault((tech,), []).append(False)
#                 else:
#                     self.outputs['cmds'].setdefault((tech,), []).append(cmd)
#                 io_update(self, i_f=fastxs, i_d=db_out, o_d=db_out, key=tech)
#                 self.soft.add_status(tech, self.sam_pool, 1)
#             else:
#                 self.soft.add_status(tech, self.sam_pool, 0)
#             self.outputs['dirs'].append(db_out)
#
#
# def kallisto(self) -> None:
#     """kallisto is a program for quantifying abundances of transcripts from
#     RNA-Seq data, or more generally of target sequences using high-throughput
#     sequencing reads. It is based on the novel idea of pseudoalignment for
#     rapidly determining the compatibility of reads with targets, without the
#     need for alignment. On benchmarks with standard RNA-Seq data, kallisto
#     can quantify 30 million human bulk RNA-seq reads in less than 3 minutes
#     on a Mac desktop computer using only the read sequences and a
#     transcriptome index that itself takes than 10 minutes to build.
#     Pseudoalignment of reads preserves the key information needed for
#     quantification, and kallisto is therefore not only fast, but also
#     comparably accurate to other existing quantification tools. In fact,
#     because the pseudoalignment procedure is robust to errors in the reads,
#     in many benchmarks kallisto significantly outperforms existing tools.
#
#     References
#     ----------
#     Bray, N.L., Pimentel, H., Melsted, P. and Pachter, L., 2016. Near-optimal
#     probabilistic RNA-seq quantification. Nature biotechnology, 34(5),
#     pp.525-527.
#
#     Notes
#     -----
#     GitHub  : https://github.com/pachterlab/kallisto
#     Docs    : http://pachterlab.github.io/kallisto/manual.html
#     Paper   : https://doi.org/10.1038/nbt.3519
#
#     Parameters
#     ----------
#     self
#     """
#     pass