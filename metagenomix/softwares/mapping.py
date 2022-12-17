# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import splitext
from metagenomix._io_utils import io_update, status_update, to_do
from metagenomix.core.parameters import tech_params
from metagenomix._inputs import (
    sample_inputs, group_inputs, genome_key, genome_out_dir)
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
        sys.exit('[mapping] Add "_<soft>" to specify what to map to "%s"' %
                 self.soft.prev)
    category = self.config.tools.get(source)
    if category == 'preprocessing':
        func = raw
    elif category in ['paired read merging', 'coassembly setup']:
        func = merged
    elif category in ['MAG', 'binning']:
        func = genomes
    elif category == 'assembling':
        func = assembly
    else:
        sys.exit('[mapping] Not possible to map "%s" on "%s"' % (
            source, self.soft.prev))
    return func, source


def bwa():
    cmd = 'bwa index -p %s/index %s\n' % (out_dir, contigs)
    cmd += 'bwa mem -p %s/index %s' % (out_dir, fastqs)
    cmd += ' | samtools view -bh'
    cmd += ' | samtools sort -o %s/%s\n' % (out_dir, bam)
    cmd += 'samtools index %s %s\n' % (bam, bai)


def get_reads(self, source, tech, group) -> dict:
    sams = self.pools[self.sam_pool][group]
    same_tech = self.soft.params.get('per_tech', False)
    reads = {}
    for sam in sams:
        tech_sam = self.softs[source].outputs[sam]
        if same_tech:
            reads[sam] = {tech: tech_sam[(tech, sam)]}
        else:
            reads[sam] = {t: paths for (t, sam), paths in tech_sam.items()}
    return reads


def assembly(self, func, fastas, tech, group):
    pass


def merged(self, func, fastas, tech, group):
    pass


def genomes(self, func, fastas, tech, group):
    pass


def get_bowtie2_db_cmd(
        fasta: str,
        sam_tech_dir: str,
        params: dict
) -> tuple:
    """Commands to build a Bowtie2 database indices for the current fasta
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


def raw(
        self,
        tech: str,
        group: str,
        reads: dict,
        fasta: list,
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
    fasta : list
    key : list
    out_dir : str
    to_dos: list
    """
    for sam, reads_tech_fastqs in reads.items():
        for reads_tech, fastqs in reads_tech_fastqs.items():
            reads_to_dos = status_update(self, reads_tech, fastqs)
            cur_key = tuple(key + [sam, reads_tech])
            sam_tech_dir = '/'.join([
                self.dir, tech, self.sam_pool, group, sam, reads_tech])
            self.outputs['dirs'].append(sam_tech_dir)
            self.outputs['outs'].setdefault(cur_key, []).append(sam_tech_dir)
            out = '%s/ali' % sam_tech_dir
            print()
            print(to_dos)
            print(reads_to_dos)
            print(sam_tech_dir)
            print(fasta)
            print(fastqs)
            print(out)
            # print(tech_fastqsfds)
            if self.config.force or to_do(out):
                params = tech_params(self, reads_tech)
                for aligner in params['aligners']:
                    get_aligner_db_cmd = 'get_%s_db_cmd' % aligner
                    get_aligner_cmd = '%s_cmd' % aligner
                    db_cmds, dbs = globals()[get_aligner_db_cmd](
                        fasta, sam_tech_dir, params)
                    for db, db_index in dbs.items():
                        cmd, sam = globals()[get_aligner_cmd](
                            sam, fastqs, db, sam_tech_dir, params)
                        bam = '%s.bam' % splitext(sam)[0]
                        bai = '%s.bai' % splitext(sam)[0]
                        cmd += ' | samtools view -bh'
                        # cmd += ' | samtools sort -o %s/%s\n' %
                        cmd += 'samtools index %s %s\n' % (bam, bai)


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
    for genome, fasta in references.items():
        key = list(genome_key(tech, group, genome))
        out_dir = genome_out_dir(self, tech, group, genome)
        to_dos = status_update(self, tech, fasta, group=group, genome=genome)
        func(self, tech, group, reads, fasta, key, out_dir, to_dos)


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
    func, source = get_mapping_target(self)
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            reads = get_reads(self, source, tech, group)
            references = group_inputs(self, inputs)
            get_mapping(self, func, tech, group, reads, references)

    elif set(self.inputs) == {''}:
        print("self.sam_pool")
        print(self.sam_pool)
        for (tech, mags), inputs in self.inputs[''].items():
            reads = get_reads(self, source, tech, mags)
            references = group_inputs(self, inputs)
            get_mapping(self, func, tech, mags, reads, references)


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


def salmon(self) -> None:
    """Salmon is a tool for wicked-fast transcript quantification from
    RNA-seq data. It requires a set of target transcripts (either from a
    reference or de-novo assembly) to quantify. All you need to run Salmon
    is a FASTA file containing your reference transcripts and a (set of)
    FASTA/FASTQ file(s) containing your reads. Optionally, Salmon can make
    use of pre-computed alignments (in the form of a SAM/BAM file) to the
    transcripts rather than the raw reads.

    References
    ----------
    Patro, Rob, et al. "Salmon provides fast and bias-aware quantification of
    transcript expression." Nature methods 14.4 (2017): 417-419.

    Notes
    -----
    GitHub  : https://github.com/COMBINE-lab/salmon
    Docs    : https://salmon.readthedocs.io/en/latest/salmon.html
    Paper   : https://doi.org/10.1038/nmeth.4197

    Parameters
    ----------
    self
    """
    if self.soft.params['useAlignments']:
        print(self.config.__dict__.keys())
        print(self.config.softs)
        print(selfconfigsofts)

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            reads = get_reads(self, tech, group)
            references = group_inputs(self, inputs)
            get_mapping(self, func, tech, group, reads, references)
            # get_salmon_indexing_cmd(self, tech, fasta, fastx, key)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            reads = get_reads(self, tech, mags)
            references = group_inputs(self, inputs)
            get_mapping(self, func, tech, mags, reads, references)
            # get_salmon_indexing_cmd(self, tech, fasta, fastx, key)

    idx = '%s.si' % (out_dir, basename(splitext(fasta)[0]))

    for (tech, sample), fastxs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastxs, tech, sample):
            continue

        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['outs'][(tech, self.sam_pool)] = dict()

        for db, db_path in self.soft.params['databases'].items():

            db_out = '%s/%s' % (out, db)
            params = tech_params(self, tech)
            cmd, sam = get_minimap2_cmd(fastxs, db_path, db_out, params)
            self.outputs['outs'][(tech, self.sam_pool)][(db, 'minimap2')] = sam

            if self.config.force or to_do(sam):
                if status_update(self, tech, fastxs):
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                io_update(self, i_f=fastxs, i_d=db_out, o_d=db_out, key=tech)
                self.soft.add_status(tech, self.sam_pool, 1)
            else:
                self.soft.add_status(tech, self.sam_pool, 0)
            self.outputs['dirs'].append(db_out)


def kallisto(self) -> None:
    """kallisto is a program for quantifying abundances of transcripts from
    RNA-Seq data, or more generally of target sequences using high-throughput
    sequencing reads. It is based on the novel idea of pseudoalignment for
    rapidly determining the compatibility of reads with targets, without the
    need for alignment. On benchmarks with standard RNA-Seq data, kallisto
    can quantify 30 million human bulk RNA-seq reads in less than 3 minutes
    on a Mac desktop computer using only the read sequences and a
    transcriptome index that itself takes than 10 minutes to build.
    Pseudoalignment of reads preserves the key information needed for
    quantification, and kallisto is therefore not only fast, but also
    comparably accurate to other existing quantification tools. In fact,
    because the pseudoalignment procedure is robust to errors in the reads,
    in many benchmarks kallisto significantly outperforms existing tools.

    References
    ----------
    Bray, N.L., Pimentel, H., Melsted, P. and Pachter, L., 2016. Near-optimal
    probabilistic RNA-seq quantification. Nature biotechnology, 34(5),
    pp.525-527.

    Notes
    -----
    GitHub  : https://github.com/pachterlab/kallisto
    Docs    : http://pachterlab.github.io/kallisto/manual.html
    Paper   : https://doi.org/10.1038/nbt.3519

    Parameters
    ----------
    self
    """
    pass