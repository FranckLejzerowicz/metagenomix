# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
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
    # bwa_cmd,
)


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


def bwa():
    cmd = 'bwa index -p %s/index %s\n' % (out_dir, contigs)
    cmd += 'bwa mem -p %s/index %s' % (out_dir, fastqs)
    cmd += ' | samtools view -bh'
    cmd += ' | samtools sort -o %s/%s\n' % (out_dir, bam)
    cmd += 'samtools index %s %s\n' % (bam, bai)


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
    db_cmds : tuple
        (Command to create the fasta db, Command to delete the fasta db)
    dbs : dict
        Databases indices
    """
    db_cmds = ''
    dbs = {}

    return db_cmds, dbs


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
        db = splitext(basename(fasta))[0]
        bam_dir = out + '/' + db
        bam = '%s/alignment.bowtie2.bam' % bam_dir
        bam_sorted = '%s.sorted.bam' % splitext(bam)[0]
        dbs[db] = (bam, bam_sorted)
        if to_do(bam_sorted):
            if to_do(bam):
                cmd += 'mkdir -p %s\n' % bam_dir
                cmd += 'mkdir -p %s/dbs/%s\n' % (out, db)
                cmd += 'bowtie2-build'
                cmd += ' --threads %s' % params['cpus']
                cmd += ' %s %s/dbs/%s\n' % (fasta, out, db)
    if cmd:
        cmd = cmd_gz + cmd
        cmd_rm += '\nrm -rf %s/dbs\n' % out
    return cmd, cmd_rm, dbs


def get_cmds(
        sam: str,
        out: str,
        fastqs: list,
        fastas: list,
        aligner: str,
        params: dict
) -> tuple:
    """

    Parameters
    ----------
    sam : str
    out : str
    fastqs : list
    fastas : list
    aligner : str
    params : dict

    Returns
    -------
    cmd : str
    bams : list
    bams_sorted : list
    """
    bams, bams_sorted = [], []
    cmd, cmd_rm, dbs = globals()['get_%s_db_cmd' % aligner](out, fastas, params)
    for db_, (bam, bam_sorted) in dbs.items():
        db = '%s/dbs/%s' % (out, db_)
        if to_do(bam):
            cmd += globals()['%s_cmd' % aligner](sam, fastqs, db,
                                                 out, bam, params)
        else:
            cmd_rm = ''
        if to_do(bam_sorted):
            cmd += 'samtools sort %s > %s\n' % (bam, bam_sorted)
            cmd += 'samtools index %s\n' % bam_sorted
        bams.append(bam)
        bams_sorted.append(bam_sorted)
    if cmd:
        cmd += cmd_rm
    return cmd, bams, bams_sorted


def raw(
        self,
        tech: str,
        group: str,
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
    tech : str
    group : str
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
                cur_key = tuple(list(key) + [sam, reads_tech])
                out = '/'.join([out_dir, sam, reads_tech, ali])
                self.outputs['dirs'].append(out)
                cmd, bams, bams_sorted = get_cmds(
                    sam, out, fastqs, fastas, ali, params)
                self.outputs['outs'].setdefault(cur_key, []).append(out)
                if self.config.force or to_do(bams_sorted[-1]):
                    if to_dos or reads_to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                    if bams:
                        io_update(self, i_f=bams, o_d=out, key=key)
                    else:
                        io_update(self, i_f=(fastqs + fastas), o_d=out, key=key)
                    self.soft.add_status(tech, self.sam_pool, 1, group=group)
                else:
                    self.soft.add_status(tech, self.sam_pool, 0, group=group)


def get_mapping(
        self,
        func,
        tech: str,
        group: str,
        reads: dict,
        references: dict,
) -> None:
    """

    Parameters
    ----------
    self
    func
    tech : str
    group : str
    reads : dict
    references : dict
    """
    for genome, fastas in references.items():
        key = genome_key(tech, group, genome)
        out_dir = genome_out_dir(self, tech, group, genome)
        to_dos = status_update(self, tech, fastas, group=group, genome=genome)
        func(self, tech, group, reads, fastas, key, out_dir, to_dos)


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
        for (ref_tech, group), inputs in self.inputs[self.sam_pool].items():
            references = group_inputs(self, inputs)
            reads = get_group_reads(self, ref_tech, group, all_reads)
            get_mapping(self, func, ref_tech, group, reads, references)


# def prep_map__spades_prodigal(self):
#     if 'prodigal' not in self.softs or 'mapping' not in self.softs:
#         return None
#     if self.softs['prodigal'].prev != 'spades':
#         return None
#     if self.softs['mapping'].prev != 'spades':
#         return None
#     prodigals_fps = self.softs['prodigal'].outputs
#     sams_fps = self.softs['mapping'].outputs
#     group_fps = self.inputs[self.pool]
#     self.outputs['outs'] = {}
#     for group, fps in group_fps.items():
#         self.outputs['outs'][group] = {}
#         for sam in self.pools[self.pool][group]:
#             bam = sams_fps[self.pool][group][sam]
#             prot = prodigals_fps[self.pool][group][1]
#             out_dir = '%s/%s/%s' % (self.dir, self.pool, sam)
#             out = '%s/reads.txt' % out_dir
#             if not isfile(out):
#                 cmd = 'pysam_reads_to_prodigal.py \\\n'
#                 cmd += '-prodigal %s \\\n' % prot
#                 cmd += '-bam %s \\\n' % bam
#                 cmd += '-out %s\n' % out
#                 self.outputs['cmds'].setdefault(
#                     (self.pool, group), []).append(cmd)
#             self.outputs['outs'][group][sam] = out
#             io_update(self, i_f=[prot, bam, out])


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