# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import sys
from os.path import basename, splitext
from metagenomix._io_utils import io_update, status_update, to_do
from metagenomix.core.parameters import tech_params
from metagenomix._inputs import (
    group_inputs, genome_key, genome_out_dir, get_reads, get_group_reads)
from metagenomix.softwares.alignment import (
    bowtie2_cmd,
    minimap2_cmd,
    # bbmap_cmd,
    # bwa_cmd,
)

# Keep line because read mapping alignments will be python scripts
# scripts = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


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
        fastas: list,
        out_dir: str,
        params: dict,
        cmd: str
) -> dict:
    """Commands to build a Bowtie2 database indices for the current fasta
    reference and paths to these indices.

    Parameters
    ----------
    fastas : list
    out_dir : str
    params : dict
    cmd: str

    Returns
    -------
    dbs : dict
        Databases indices
    """
    cmd_rm, dbs = '', {}
    for fasta_ in fastas:
        fasta = fasta_
        if fasta.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
            fasta = fasta.rstrip('.gz')
            cmd_rm += 'rm %s\n' % fasta
        db = splitext(basename(fasta))[0]
        dbs[db] = '%s/dbs/%s' % (out_dir, db)
        cmd += '\nbowtie2-build'
        cmd += ' --threads %s' % params['cpus']
        cmd += ' %s %s\n' % (fasta, dbs[db])
        cmd += 'rm %s\n' % fasta
    cmd += cmd_rm
    return dbs


def get_cmds(
        sam: str,
        fastqs: list,
        fastas: list,
        out_dir: str,
        aligner: str,
        params: dict
) -> tuple:
    ali_db_cmd = globals()['get_%s_db_cmd' % aligner]
    ali_cmd = globals()['%s_cmd' % aligner]
    cmd, bams = '', []
    dbs = ali_db_cmd(fastas, out_dir, params, cmd)
    for db, db_index in dbs.items():
        bam = '%s/%s/alignment.bowtie2.bam' % (out_dir, db)
        cmd += ali_cmd(sam, fastqs, db_index, out_dir, bam, params)
        bams.append(bam)
    return cmd, bams


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
        # print("sample:", sam)
        for (reads_tech, _), fastqs in reads_tech_fastqs.items():
            reads_to_dos = status_update(self, reads_tech, fastqs)
            for ali in self.soft.params['aligners']:
                params = tech_params(self, reads_tech, ali)
                cur_key = tuple(list(key) + [sam, reads_tech])
                out = '/'.join([out_dir, sam, reads_tech, ali])
                self.outputs['dirs'].append(out)
                self.outputs['outs'].setdefault(cur_key, []).append(out)
                # print("reads_tech:", reads_tech)
                # print("cur_key:", cur_key)
                # print("sam_tech_dir:", out)
                # print("to_dos:", to_dos)
                # print("reads_to_dos:", reads_to_dos)
                # print("fastqs:", fastqs)
                if self.config.force or glob.glob(
                        '%s/*.bam' % out.replace('${SCRATCH_FOLDER}', '')):
                    cmd, bams = get_cmds(sam, fastqs, fastas, out, ali, params)
                    if to_dos or reads_to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
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
        # print('genome:', genome)
        # print('fastas:', fastas)
        # print('key:', key)
        # print('out_dir:', out_dir)
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
    # print("func:", func)
    # print("source:", source)
    # print("step:", step)
    # print("all_reads:", all_reads)
    if self.sam_pool in self.pools:
        for (ref_tech, group), inputs in self.inputs[self.sam_pool].items():
            references = group_inputs(self, inputs)
            reads = get_group_reads(self, ref_tech, group, all_reads)
            # print("ref_tech:", ref_tech)
            # print("group:", group)
            # print("inputs:", inputs)
            # print("references:", references)
            # print("reads:", reads)
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