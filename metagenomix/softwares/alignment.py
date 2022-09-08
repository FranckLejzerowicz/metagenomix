# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext
from metagenomix.core.parameters import tech_params
from metagenomix._io_utils import (io_update, to_do, tech_specificity,
                                   not_paired, status_update)

# Keep line because read mapping alignments will be python scripts
# scripts = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


def grep_sam_tail_cmd(
        cmd: str,
        reads: str,
        sam: str
) -> str:
    """Grep the end of the alignment and compare is to end of input to make
    sure that the alignment was not aborted and hence that the output file is
    not incomplete.

    Notes
    -----
    It is possible that the last sequences of the input file do not align,
    but it is safer to recompute the alignment than to assume that and
    potentially leave many sequences out.

    Parameters
    ----------
    cmd : str
        Alignment commands
    reads : str
        last reads of the input sequences file
    sam : str
        Path to the output .sam alignment file

    Returns
    -------
    grep_cmd : str
        Command to grep the end of the alignment and compare is to end of input
    """
    grep_cmd = '\nlast_reads=`%s`\n' % reads
    grep_cmd += 'GREPOK=0\n'
    grep_cmd += 'for i in $listVar\n'
    grep_cmd += 'do\n'
    grep_cmd += '    if grep -q "$i" %s\n' % sam
    grep_cmd += '        then\n'
    grep_cmd += '            GREPOK=1\n'
    grep_cmd += '            break\n'
    grep_cmd += '    fi\n'
    grep_cmd += 'done\n'
    grep_cmd += "if [ $GREPOK = 0 ]\n"
    grep_cmd += 'then\n'
    grep_cmd += '%s\n' % cmd
    grep_cmd += 'fi\n'
    return grep_cmd


def condition_ali_command(
        fastxs: list,
        cmd: str,
        sam: str
) -> str:
    """Add the conditional statements checking that the alignment
    already available was not aborted, in which case the alignment
    will be re-run.

    Parameters
    ----------
    fastxs : list
        Path to the input fasta/fastq/fastq.gz files
    cmd : str
        Alignment command
    sam : str
        Path to the alignment file

    Returns
    -------
    cmd : str
        Alignment command potentially decorated with unix check to rerun
    """
    if fastxs[-1].endswith('fastq.gz'):
        reads = 'zcat %s | tail -n 80 | grep "^@"' % fastxs[-1]
    elif fastxs[-1].endswith('fastq'):
        reads = 'tail -n 80 %s | grep "^@"' % fastxs[-1]
    elif fastxs[-1].endswith('fasta'):
        reads = 'tail -n 40 %s | grep "^>"' % fastxs[-1]
    else:
        return cmd

    cmd = grep_sam_tail_cmd(cmd, reads, sam)
    return cmd


# def get_burst_cmd(self, tech: str, fastx: str, db_path: str, out: str) -> str:
#     """Get the burst alignment command.
#
#     Parameters
#     ----------
#     tech : str
#         Technology: 'illumina', 'pacbio', or 'nanopore'
#     fastx : str
#         Path to the input fasta file
#     db_path : str
#         Path to the burst database
#     out : str
#         Path to the output folder
#
#     Returns
#     -------
#     cmd : str
#         Alignment command
#     """
#     params = tech_params(self, tech)
#     cmd = 'burst'
#     cmd += ' -q %s' % fastx
#     cmd += ' -r %s.edx' % db_path
#     cmd += ' -a %s.acx' % db_path
#     cmd += ' -sa'
#     cmd += ' -fr'
#     cmd += ' -i 0.98'
#     cmd += ' -t %s' % params['cpus']
#     cmd += ' -o %s' % out
#     return cmd


def get_bowtie2_cmd(
        self,
        tech: str,
        fastx: list,
        db_path: str,
        db_out: str
) -> tuple:
    """Get the bowtie2 alignment command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for humann
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastx : list
        Path to the input fasta/fastq/fastq.gz files
    db_path : str
        Path to the bowtie2 database
    db_out: str
        Path to the output folder

    Returns
    -------
    cmd : str
        Alignment command
    sam : str
        Alignment .sam file
    """
    params = tech_params(self, tech)
    cmd = 'bowtie2'
    cmd += ' -p %s' % params['cpus']
    cmd += ' -x %s' % db_path
    sam = '%s/alignment.bowtie2' % db_out
    if len(fastx) == 2 and params['pairing'] == 'paired':
        cmd += ' -1 %s' % fastx[0]
        cmd += ' -2 %s' % fastx[1]
        if not params['discordant']:
            cmd += ' --no-discordant'
            sam += '.nodiscordant'
    else:
        cmd += ' -U %s' % ','.join(fastx)
    sam += '.sam'
    cmd += ' -S %s' % sam
    cmd += ' --seed 12345'
    cmd += ' --very-sensitive'
    cmd += ' -k %s' % params['k']
    cmd += ' --np %s' % params['np']
    cmd += ' --mp "%s"' % params['mp']
    cmd += ' --rdg "%s"' % params['rdg']
    cmd += ' --rfg "%s"' % params['rfg']
    cmd += ' --score-min "%s"' % params['score-min']
    cmd += ' --met-file %s.met' % splitext(sam)[0]
    cmd += ' --met 240'
    cmd += ' --no-unal'
    return cmd, sam


def get_alignment_cmd(
        fastxs: list,
        cmd: str,
        ali: str
) -> str:
    """Get the alignment command with ot without unix checks to decide
    whether the alignment (if exists) was aligned to completion, and not
    aborted (e.g., for lack of memory reasons).

    Parameters
    ----------
    fastxs : list
        Path to the input fasta/fastq/fastq.gz files
    cmd : str
        Alignment command
    ali : str
        Path to the alignment file

    Returns
    -------
    cmd : str
        Alignment command potentially decorated with unix check to rerun
    """
    if to_do(ali):
        return cmd
    else:
        return condition_ali_command(fastxs, cmd, ali)


def bowtie2(self) -> None:
    """Get the full command lines to perform sequence alignment using
     bowtie2 and potentially, re-run if the output was existing but
     possibly aborted.

    References
    ----------
    Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with
    Bowtie 2. Nature methods, 9(4), 357-359.

    Notes
    -----
    GitHub  : https://github.com/BenLangmead/bowtie2
    Docs    : http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    Paper   : https://doi.org/10.1038/nmeth.1923

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam_pool : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .config
            Configurations
    """
    for (tech, sample), fastxs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastxs, tech, sample):
            continue
        status_update(self, tech, fastxs)
        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['outs'][(tech, self.sam_pool)] = dict()
        for db, db_path in self.soft.params['databases'].items():
            db_out = '%s/%s/%s' % (out, db, self.soft.params['pairing'])
            cmd, sam = get_bowtie2_cmd(self, tech, fastxs, db_path, db_out)
            self.outputs['outs'][(tech, self.sam_pool)][(db, 'bowtie2')] = sam
            if self.config.force or to_do(sam):
                cmd = get_alignment_cmd(fastxs, cmd, sam)
                self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                io_update(self, i_f=fastxs, i_d=db_out, o_d=db_out, key=tech)
                self.soft.add_status(tech, self.sam_pool, 1)
            else:
                self.soft.add_status(tech, self.sam_pool, 0)
            self.outputs['dirs'].append(db_out)


def flash_cmd(
        self,
        tech: str,
        fastqs: list,
        out: str
) -> str:
    """Get the FLASh merging command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for FLASh
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Paths to the input files
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        FLASh command
    """
    params = tech_params(self, tech)
    cmd = 'flash %s' % ' '.join(fastqs)
    cmd += ' --min-overlap %s' % params['min_overlap']
    cmd += ' --max-overlap %s' % params['max_overlap']
    cmd += ' --max-mismatch-density %s' % params['mismatch']
    cmd += ' --output-directory %s' % out
    cmd += ' --output-prefix %s' % self.sam_pool
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --compress'
    return cmd


def flash(self) -> None:
    """FLASH (Fast Length Adjustment of SHort reads) is a very fast and
    accurate software tool to merge paired-end reads from next-generation
    sequencing experiments. FLASH is designed to merge pairs of reads when
    the original DNA fragments are shorter than twice the length of reads.
    The resulting longer reads can significantly improve genome assemblies.
    They can also improve transcriptome assembly when FLASH is used to merge
    RNA-seq data.

    References
    ----------
    Magoč, Tanja, and Steven L. Salzberg. "FLASH: fast length adjustment of
    short reads to improve genome assemblies." Bioinformatics 27.21 (2011):
    2957-2963.

    Notes
    -----
    SourceForge : https://sourceforge.net/projects/flashpage/files/
    GitHub      : https://github.com/genome-vendor/FLASH
    Docs        : https://ccb.jhu.edu/software/FLASH/
    Paper       : https://doi.org/10.1093/bioinformatics/btr507

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for FLASh
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for FLASh
        .config
            Configurations
    """
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        status_update(self, tech, fastqs)

        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        rad = out + '/' + self.sam_pool
        ext = '%s.extendedFrags.fastq.gz' % rad
        nc1 = '%s.notCombined_1.fastq.gz' % rad
        nc2 = '%s.notCombined_2.fastq.gz' % rad
        outs = [ext, nc1, nc2]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            cmd = flash_cmd(self, tech, fastqs, out)
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def ngmerge_cmd(
        self,
        tech: str,
        fastqs: list,
        ext: str,
        log: str,
        fail: str) -> str:
    """Get the NGmerge merging command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for NGmerge
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Paths to the input files
    ext : str
        Path to the merged reads output file
    log : str
        Path to the reads merging log file
    fail : str
        Path to the unmerged reads output file

    Returns
    -------
    cmd : str
        NGmerge command
    """
    params = tech_params(self, tech)
    cmd = 'NGmerge'
    cmd += ' -1 %s' % fastqs[0]
    cmd += ' -2 %s' % fastqs[1]
    cmd += ' -o %s' % ext.replace('.gz', '')
    if params['a']:
        cmd += ' -c %s.log' % log
    else:
        cmd += ' -j %s.log' % log
        cmd += ' -f %s' % fail

    for boolean in ['a', 'd', 's', 'i', 'g']:
        if params[boolean]:
            cmd += ' -%s' % boolean
    for param in ['p', 'm', 'e', 'q', 'u']:
        cmd += ' -%s %s' % (param, params[param])
    cmd += ' -n %s' % params['cpus']
    cmd += ' -z'
    return cmd


def ngmerge(self) -> None:
    """NGmerge operates on paired-end high-throughput sequence reads in two
    distinct modes.
    - In the default stitch mode, NGmerge combines paired-end reads that
      overlap into a single read that spans the full length of the original
      DNA fragment (Fig. 1A). The ends of the merged read are defined by the
      5' ends of the original reads. Reads that fail the stitching process
      (due to a lack of sufficient overlap, or excessive sequencing errors)
      are placed into secondary output files, if the user requires them.
    - The alternative adapter-removal mode returns the original reads as
      pairs, removing the 3' overhangs of those reads whose valid stitched
      alignment has this characteristic (Fig. 1B). Reads whose alignments do
      not have such overhangs (or do not align at all) will also be printed
      to the output files, unmodified.

    References
    ----------
    Gaspar, John M. "NGmerge: merging paired-end reads via novel
    empirically-derived models of sequencing errors." BMC bioinformatics 19.1
    (2018): 1-9.

    Notes
    -----
    GitHub  : https://github.com/harvardinformatics/NGmerge
    Paper   : https://doi.org/10.1186/s12859-018-2579-2

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for NGmerge
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for NGmerge
        .config
            Configurations
    """
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        status_update(self, tech, fastqs)

        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        rad = out + '/' + self.sam_pool
        if self.soft.params['a']:
            ext = '%s.trim' % rad
            log = '%s_trim' % rad
        else:
            ext = '%s.extendedFrags.fastq.gz' % rad
            log = '%s_merging' % rad
        fail = '%s.notCombined' % rad
        nc1 = '%s.notCombined_1.fastq.gz' % rad
        nc2 = '%s.notCombined_2.fastq.gz' % rad
        outs = [ext, nc1, nc2]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            cmd = ngmerge_cmd(self, tech, fastqs, ext, log, fail)
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def pear_cmd(
        self,
        tech: str,
        fastqs: list,
        rad: str,
        outs: list,
        outs_: list,
        na: str) -> str:
    """Get the PEAR merging command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for PEAR
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Paths to the input files
    rad : str
        Radical of the path to the output files
    outs : list
        Paths to the reads merging output files (after renaming)
    outs_ : list
        Paths to the reads merging output files (before renaming)
    na : str
        Path to the unmerged reads output file

    Returns
    -------
    cmd : str
        PEAR command
    """
    params = tech_params(self, tech)
    cmd = 'pear'
    cmd += ' --forward-fastq %s' % fastqs[0]
    cmd += ' --reverse-fastq %s' % fastqs[1]
    cmd += ' --output %s' % rad
    for param in [
        'min_assembly_length', 'max_assembly_length', 'quality_threshold',
        'min_trim_length', 'score_method', 'min_overlap', 'test_method',
        'phred_base', 'cap', 'p_value', 'max_uncalled_base'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    for boolean in ['empirical_freqs', 'keep_original', 'stitch', 'nbase']:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --threads %s\n' % params['cpus']
    for (src, snk) in zip(outs_, outs):
        cmd += 'gzip %s\nmv %s %s\n' % (src.replace('.gz', ''), src, snk)
    cmd += 'gzip %s\n' % na.replace('.gz', '')
    return cmd


def pear(self) -> None:
    """PEAR assembles Illumina paired-end reads if the DNA fragment sizes are
    smaller than twice the length of reads. PEAR can assemble 95% of reads
    with 35-bp mean overlap with a false-positive rate of 0.004. PEAR also
    works with multiplexed data sets where the true underlying DNA fragment
    size varies. PEAR has an extremely low false-positive rate of 0.0003 on
    data sets where no overlap exists between the two reads (i.e. when DNA
    fragment sizes are larger than twice the read length).

    References
    ----------
    Zhang, Jiajie, et al. "PEAR: a fast and accurate Illumina Paired-End
    reAd mergeR." Bioinformatics 30.5 (2014): 614-620.

    Notes
    -----
    GitHub  : https://github.com/tseemann/PEAR
    Docs    : http://www.exelixis-lab.org/web/software/pear
    Paper   : https://doi.org/10.1093/bioinformatics/btt593

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for PEAR
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for PEAR
        .config
            Configurations
    """
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        status_update(self, tech, fastqs)

        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        rad = out + '/' + self.sam_pool
        ext_ = '%s.assembled.fastq.gz' % rad
        nc1_ = '%s.unassembled.forward.fastq.gz' % rad
        nc2_ = '%s.unassembled.reverse.fastq.gz' % rad
        outs_ = [ext_, nc1_, nc2_]
        ext = '%s.extendedFrags.fastq.gz' % rad
        nc1 = '%s.notCombined_1.fastq.gz' % rad
        nc2 = '%s.notCombined_2.fastq.gz' % rad
        outs = [ext, nc1, nc2]
        na = '%s.discarded.fastq.gz' % rad
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            cmd = pear_cmd(self, tech, fastqs, rad, outs, outs_, na)
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def bbmerge_cmd(
        self,
        tech: str,
        fastqs: list,
        os: list
) -> str:
    """Get the bbmerge.sh merging command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for bbmerge.sh
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Paths to the input files
    os : list
        Paths to the output files

    Returns
    -------
    cmd : str
        bbmerge.sh command
    """
    params = tech_params(self, tech)
    cmd = 'bbmerge.sh'
    cmd += ' in1=%s in2=%s' % tuple(fastqs)
    cmd += ' out=%s outu1=%s outu2=%s outinsert=%s outc=%s ihist=%s' % tuple(os)
    if params['strictness']:
        cmd += ' %s=t' % params['strictness']
    else:
        for param in [
            'reads', 'ziplevel', 'trimq','minlength', 'maxlength',
            'minavgquality', 'maxexpectederrors', 'forcetrimleft',
            'forcetrimright', 'forcetrimright2', 'forcetrimmod',
            'mininsert', 'mininsert0', 'minoverlap', 'minoverlap0', 'minq',
            'maxq', 'efilter', 'pfilter', 'kfilter', 'maxratio',
            'ratiomargin', 'ratiooffset', 'maxmismatches',
            'ratiominoverlapreduction', 'minsecondratio', 'margin',
            'mismatches', 'k', 'extend', 'extend2', 'iterations',
            'mindepthseed', 'mindepthextend', 'branchmult1', 'branchmult2',
            'branchlower', 'prealloc', 'prefilter', 'filtermem', 'minprob',
            'minapproxoverlap'
        ]:
            cmd += ' %s=%s' % (param, params[param])
        for boolean in [
            'interleaved', 'nzo', 'showhiststats', 'ordered', 'mix',
            'qtrim', 'qtrim2', 'tbo', 'ooi', 'trimpolya', 'usejni', 'merge',
            'ecco', 'trimnonoverlapping', 'useoverlap', 'entropy', 'ouq',
            'owq', 'usequality', 'iupacton', 'ratiomode', 'forcemerge',
            'flatmode', 'requireratiomatch', 'trimonfailure', 'rem', 'rsem',
            'ecctadpole', 'reassemble', 'removedeadends', 'removebubbles',
            'ibb', 'eccbloom', 'testmerge', 'eoom', 'da',
        ]:
            if params[boolean]:
                cmd += ' %s=t' % boolean
            else:
                cmd += ' %s=f' % boolean
    return cmd


def bbmerge(self) -> None:
    """BBMerge is designed to merge two overlapping paired reads into a
    single read. For example, a 2x150bp read pair with an insert size of
    270bp would result in a single 270bp read. This is useful in amplicon
    studies, as clustering and consensus are far easier with single reads
    than paired reads, and also in assembly, where longer reads allow the
    use of longer kmers (for kmer-based assemblers) or fewer comparisons
    (for overlap-based assemblers). And in either case, the quality of the
    overlapping bases is improved. BBMerge is also capable of error-correcting
    the overlapping portion of reads without merging them, as well as merging
    nonoverlapping reads, if enough coverage is available. BBMerge is the
    fastest and by far the most accurate overlap-based read merger currently
    in existence.

    References
    ----------
    Bushnell, Brian, Jonathan Rood, and Esther Singer. "BBMerge–accurate
    paired shotgun read merging via overlap." PloS one 12.10 (2017): e0185056.

    Notes
    -----
    SourceForge : https://sourceforge.net/projects/bbmap/
    Docs        : https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/
    Paper       : https://doi.org/10.1371/journal.pone.0185056

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for BBmerge
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for BBmerge
        .config
            Configurations
    """
    # iterate over the inputs
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        status_update(self, tech, fastqs)

        # make the output directory
        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        # get the expected names of some of the ouptuts:
        # - those you want to collect in 'outs' because they will be future inpt
        # - at least one that will help knowing whether the software already run
        rad = out + '/' + self.sam_pool
        ext = '%s.extendedFrags.fastq.gz' % rad
        nc1 = '%s.notCombined_1.fastq.gz' % rad
        nc2 = '%s.notCombined_2.fastq.gz' % rad
        ins = '%s.inserts.txt' % rad
        kmer = '%s.kmer_cardinality.txt' % rad
        ihist = '%s.insert_length.hist' % rad
        outs = [ext, nc1, nc2]
        outs_cmd = outs + [ins, kmer, ihist]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        # check if the tool already run (or if --force) to allow getting command
        if self.config.force or sum([to_do(x) for x in outs]):
            # collect the command line
            cmd = bbmerge_cmd(self, tech, fastqs, outs_cmd)
            # add is to the 'cmds'
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


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
    print()