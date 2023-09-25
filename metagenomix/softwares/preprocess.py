# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import basename
from metagenomix._io_utils import (io_update, to_do, tech_specificity,
                                   status_update)
from metagenomix.core.parameters import tech_params
from metagenomix.softwares.alignment import (
    bowtie2_cmd,
    minimap2_cmd,
    # bbmap_cmd,
    # bwa_cmd
)


def count_cmd(
        self,
        idx: int,
        fastx: str,
        out: str
) -> str:
    """Get the command to count the reads in a fasta, fastq or fastq.gz file.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .soft.params
            Parameters
    idx : int
        Unique, incremental numeric identifier.
    fastx : str
        Path to the input file name.
    out : str
        Path to the output file name.

    Returns
    -------
    cmd : str
        Command to count reads and output the counts.
    """
    if fastx.endswith('fastq.gz'):
        cmd = "n%s=`zcat %s | wc -l | " % (idx, fastx)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
    elif fastx.endswith('fasta'):
        cmd = "n%s=`wc -l %s | " % (idx, fastx)
        cmd += "sed 's/ //g' | awk '{x=$1/2; print x}'`\n"
    elif fastx.endswith('fastq'):
        cmd = "n%s=`wc -l %s | " % (idx, fastx)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
    else:
        raise IOError("Input sequence file invalid %s" % fastx)
    if idx:
        cmd += 'echo "%s,2,${n%s}" >> %s\n' % (self.sam_pool, idx, out)
    else:
        cmd += 'echo "%s,1,${n%s}" > %s\n' % (self.sam_pool, idx, out)
    return cmd


def count(self) -> None:
    """Create command lines for counting the reads.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MIDAS
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
    outs = {}
    for (tech, sam), fastxs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastxs, tech, sam):
            continue
        to_dos = status_update(self, tech, fastxs)

        out_dir = '%s/%s/%s' % (self.dir, tech, sam)
        self.outputs['dirs'].append(out_dir)

        out = '%s/read_count.tsv' % out_dir
        outs[(tech, self.sam_pool)] = out

        if self.config.force or to_do(out):
            for idx, fastx in enumerate(fastxs):
                cmd = count_cmd(self, idx, fastx, out)
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
            io_update(self, i_f=fastxs, i_d=out_dir, o_f=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)
    self.outputs['outs'] = outs


def fastqc(self) -> None:
    """FastQC is a program designed to spot potential problems in high
    througput sequencing datasets. It runs a set of analyses on one or more
    raw sequence files in fastq or bam format and produces a report which
    summarises the results.

    References
    ----------
    Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput
    Sequence Data [Online].

    Notes
    -----
    GitHub  : https://github.com/s-andrews/FastQC
    Docs    : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for FastQC
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
    outs = {}
    for (tech, sam), fastxs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastxs, tech, sam):
            continue
        to_dos = status_update(self, tech, fastxs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        out = ['%s/%s_fastqc.html' % (
            out_dir, basename(x).rsplit('.fastq', 1)[0]) for x in fastxs]
        outs[(tech, self.sam_pool)] = out

        if self.config.force or sum([to_do(x) for x in out]):
            cmd = 'fastqc %s -o %s' % (' '.join(fastxs), out_dir)
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastxs, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)

    self.outputs['outs'] = outs


def get_fastp_cmd(
        self,
        tech: str,
        fastqs: list,
        out_dir: str,
        outs: list,
        paired: bool
) -> str:
    """Collect the Fastp command.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Path to the input files
    out_dir : str
        Path to the output folder
    outs : list
        Path to the output files
    paired : bool
        Whether the inputs are paired reads'

    Returns
    -------
    cmd : str
        Fastp command
    """
    params = tech_params(self, tech)
    cmd = 'fastp'
    for fdx, fastq in enumerate(fastqs):
        if paired:
            if fdx:
                i, o = 'I', 'O'
                out = '%s/%s_2.fastq.gz' % (out_dir, self.sam_pool)
            else:
                i, o = 'i', 'o'
                out = '%s/%s_1.fastq.gz' % (out_dir, self.sam_pool)
        else:
            i, o = 'i', 'o'
            out = '%s/%s.fastq.gz' % (out_dir, self.sam_pool)
        outs.append(out)
        cmd += '  -%s %s -%s %s' % (i, fastq, o, out)
    cmd += ' --thread %s' % params['cpus']
    cmd += ' --report_title="%s"' % self.sam_pool
    cmd += ' --json=%s/%s.json' % (out_dir, self.sam_pool)
    cmd += ' --html=%s/%s.html' % (out_dir, self.sam_pool)

    if params['split_by_lines']:
        cmd += ' --split_by_lines'
        cmd += ' --split_prefix_digits %s' % params[
            'split_prefix_digits']
    elif params['split']:
        cmd += ' --split'
        cmd += ' --split_prefix_digits %s' % params[
            'split_prefix_digits']

    if params['overrepresentation_analysis']:
        cmd += ' --overrepresentation_analysis'
        cmd += ' --overrepresentation_sampling %s' % params[
            'overrepresentation_sampling']

    if params['correction']:
        cmd += ' --correction'
        cmd += ' --overlap_len_require %s' % params[
            'overlap_len_require']
        cmd += ' --overlap_diff_limit %s' % params[
            'overlap_diff_limit']

    if params['disable_length_filtering']:
        cmd += ' --disable_length_filtering'
    else:
        cmd += ' --length_required %s' % params['length_required']
        cmd += ' --length_limit %s' % params['length_limit']

    if params['disable_adapter_trimming']:
        cmd += ' --disable_adapter_trimming'
    else:
        if params['detect_adapter_for_pe']:
            cmd += ' --detect_adapter_for_pe'
        if params['adapter_fasta']:
            cmd += ' --adapter_fasta %s' % params['adapter_fasta']
        else:
            for adapter in ['adapter_sequence', 'adapter_sequence_r2']:
                if params[adapter]:
                    cmd += ' --%s %s' % (adapter, params[adapter])

    if params['disable_quality_filtering']:
        cmd += ' --disable_quality_filtering'
    else:
        for param in ['qualified_quality_phred', 'unqualified_percent_limit',
                      'n_base_limit', 'average_qual']:
            cmd += ' --%s %s' % (param, params[param])

    for param in ['compression', 'trim_front1', 'trim_tail1', 'max_len1',
                  'trim_front2', 'trim_tail2', 'max_len2',
                  'dup_calc_accuracy', 'cut_window_size', 'cut_mean_quality']:
        if params[param]:
            cmd += ' --%s %s' % (param, params[param])

    for boolean in ['dont_overwrite', 'phred64', 'verbose', 'dedup',
                    'dont_eval_duplication', 'trim_poly_g', 'trim_poly_x',
                    'cut_front', 'cut_tail', 'cut_right']:
        if params[boolean]:
            if boolean == 'trim_poly_g' and not params['disable_trim_poly_g']:
                cmd += ' --%s' % boolean
                cmd += ' --poly_g_min_len %s' % params['poly_g_min_len']
            else:
                cmd += ' --%s' % boolean
                if boolean == 'trim_poly_x':
                    cmd += ' --poly_x_min_len %s' % params['poly_x_min_len']
                elif boolean.startswith('cut_'):
                    window_size = '%s_window_size' % boolean
                    mean_quality = '%s_mean_quality' % boolean
                    cmd += ' --%s %s' % (window_size, params[window_size])
                    cmd += ' --%s %s' % (mean_quality, params[mean_quality])
            if boolean == 'cut_front':
                cmd += ' --cut_front_mean_quality %s' % (
                    params['cut_front_mean_quality'])
                cmd += ' --cut_front_window_size %s' % (
                    params['cut_front_window_size'])
            if boolean == 'cut_tail':
                cmd += ' --cut_tail_mean_quality %s' % (
                    params['cut_tail_mean_quality'])
                cmd += ' --cut_tail_window_size %s' % (
                    params['cut_tail_window_size'])
            if boolean == 'cut_right':
                cmd += ' --cut_right_mean_quality %s' % (
                    params['cut_right_mean_quality'])
                cmd += ' --cut_right_window_size %s' % (
                    params['cut_right_window_size'])

    return cmd


def fastp_cmd(
        self,
        tech: str,
        fastqs: list,
        out_dir: str
) -> tuple:
    """Collect the Fastp command.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Path to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Fastp command
    outs : list
        Path to the output files
    """
    cmd = ''
    outs = []

    paired = [x for x in fastqs if '_1.fastq' in x or '_2.fastq' in x
              or '_R1.fastq' in x or '_R2.fastq' in x]
    if paired:
        cmd += get_fastp_cmd(self, tech, paired, out_dir, outs, True)

    unpair = [x for x in fastqs if x not in paired]
    if unpair:
        cmd += get_fastp_cmd(self, tech, unpair, out_dir, outs, False)

    return cmd, outs


def fastp(self) -> None:
    """A tool designed to provide fast all-in-one preprocessing for FastQ
    files. This tool is developed in C++ with multithreading supported to
    afford high performance.

    References
    ----------
    Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor."
    Bioinformatics 34.17 (2018): i884-i890.

    Notes
    -----
    GitHub  : https://github.com/OpenGene/fastp
    Paper   : https://doi.org/10.1093/bioinformatics/bty560

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Fastp
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
        if tech_specificity(self, fastqs, tech, sam):
            continue
        to_dos = status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outs = fastp_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def cutadapt_cmd(
        self,
        tech: str,
        fastqs: list,
        out_dir: str
) -> tuple:
    """Collect the Cutadapt command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Path to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Cutadapt command
    outs : list
        Path to the output files
    """
    params = tech_params(self, tech)
    cmd = 'cutadapt'
    cmd += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    cmd += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    for boolean in [
        'trim_n', 'revcomp', 'zero_cap', 'no_indels', 'pair_adapters',
        'discard_trimmed', 'discard_untrimmed', 'match_read_wildcards',
        'no_match_adapter_wildcards'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    for param in [
        'times', 'overlap', 'error_rate', 'quality_base', 'report',
        'pair-filter', 'action'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    cmd += ' --minimum-length 10'
    cmd += ' --cores %s' % params['cpus']

    r1_o = '%s/%s.R1.fastq.gz' % (out_dir, self.sam_pool)
    cmd += ' --output %s' % r1_o
    outs = [r1_o]
    if len(fastqs) == 2:
        r2_o = '%s/%s.R2.fastq.gz' % (out_dir, self.sam_pool)
        cmd += ' --paired-output %s ' % r2_o
        outs.append(r2_o)
    cmd += ' '.join(fastqs)
    return cmd, outs


def cutadapt(self) -> None:
    """Cutadapt finds and removes adapter sequences, primers, poly-A tails
    and other types of unwanted sequence from your high-throughput sequencing
    reads.

    Cleaning your data in this way is often required: Reads from small-RNA
    sequencing contain the 3’ sequencing adapter because the read is longer
    than the molecule that is sequenced. Amplicon reads start with a primer
    sequence. Poly-A tails are useful for pulling out RNA from your sample,
    but often you don’t want them to be in your reads.

    Cutadapt helps with these trimming tasks by finding the adapter or primer
    sequences in an error-tolerant way. It can also modify and filter
    single-end and paired-end reads in various ways. Adapter sequences can
    contain IUPAC wildcard characters. Cutadapt can also demultiplex your reads.

    References
    ----------
    Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput
    sequencing reads." EMBnet. journal 17.1 (2011): 10-12.

    Notes
    -----
    GitHub  : https://github.com/marcelm/cutadapt
    Docs    : https://cutadapt.readthedocs.io/en/stable/develop.html
    Paper   : https://doi.org/10.14806/ej.17.1.200

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Cutadapt
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
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        to_dos = status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outs = cutadapt_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_f=outs, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def atropos_outs(
        self,
        fastqs,
        out_dir: str
) -> list:
    """Collect the output file paths for Atropos.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
    fastqs : list
        Path to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    outs : list
        Path to the output files
    """
    if len(fastqs) == 1:
        outs = ['%s/%s.fastq' % (out_dir, self.sam_pool)]
    else:
        outs = ['%s/%s_R1.fastq' % (out_dir, self.sam_pool),
                '%s/%s_R2.fastq' % (out_dir, self.sam_pool)]
    outs = ['%s.gz' % outs[idx] if x.endswith('.gz') else x
            for idx, x in enumerate(fastqs)]
    return outs


def atropos_cmd(
        self,
        tech: str,
        inputs: list,
        out_dir: str,
        outputs: list
) -> str:
    """Collect the command line for Atropos filtering and trimming.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Path to the input files
    out_dir : str
        Path to the output folder
    outputs : list
        Path to the output files

    Returns
    -------
    cmd : str
        Atropos commands
    """
    params = tech_params(self, tech)
    cmd = 'atropos trim'
    for a in params['a']:
        cmd += ' -a %s' % a
    for A in params['A']:
        cmd += ' -A %s' % A
    cmd += ' -q %s' % params['q']
    cmd += ' --minimum-length %s' % params['minimum_length']
    cmd += ' --nextseq-trim %s' % params['nextseq_trim']
    cmd += ' --pair-filter %s' % params['pair_filter']
    cmd += ' --overlap %s' % params['overlap']
    cmd += ' --max-reads %s' % params['max_reads']
    cmd += ' --indel-cost %s' % params['indel_cost']
    cmd += ' --error-rate %s' % params['error_rate']
    cmd += ' --quality-cutoff %s' % params['quality_cutoff']
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --report-file %s/%s.log' % (out_dir, self.sam_pool)
    cmd += ' --report-formats %s' % params['report_formats']
    if len(inputs) == 1:
        cmd += ' -se %s' % inputs[0]
        cmd += ' -o %s' % outputs[0]
    else:
        cmd += ' -pe1 %s' % inputs[0]
        cmd += ' -pe2 %s' % inputs[1]
        cmd += ' -o %s' % outputs[0]
        cmd += ' -p %s' % outputs[1]
    return cmd


def atropos(self) -> None:
    """Atropos is tool for specific, sensitive, and speedy trimming of NGS
    reads. It is a fork of the venerable Cutadapt read trimmer.

    References
    ----------
    Didion, John P., Marcel Martin, and Francis S. Collins. "Atropos:
    specific, sensitive, and speedy trimming of sequencing reads." PeerJ 5 (
    2017): e3720.

    Notes
    -----
    GitHub  : https://github.com/jdidion/atropos
    Docs    : https://atropos.readthedocs.io/en/1.1/
    Paper   : https://doi.org/10.7717/peerj.3720

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Atropos
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
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        to_dos = status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        outs = atropos_outs(self, fastqs, out_dir)

        cmd = atropos_cmd(self, tech, fastqs, out_dir, outs)
        if self.config.force or sum([to_do(x) for x in outs]):
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_f=outs, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)

        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)


def hifiadapterfilt_cmd(
        self,
        fastqs: list,
        out_dir: str,
        sam: str
) -> str:
    """Collect the command line for HiFiAdapterFilt.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .soft.params
            Parameters
    fastqs : list
        Path to the input files
    out_dir : str
        Path to the output folder
    sam : str
        Name of the sample (prefix)

    Returns
    -------
    cmd : str
        HiFiAdapterFilt commands
    """
    cmd = ''
    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s\n' % self.soft.params['path']
        cmd += 'export PATH=$PATH:%s/DB\n' % self.soft.params['path']

    cmd += 'cd %s\n' % out_dir
    cmd += 'cp %s %s/.\n' % (' '.join(fastqs), out_dir)

    if self.soft.params['path']:
        cmd += 'bash %s/pbadapterfilt.sh' % self.soft.params['path']
    else:
        cmd += 'pbadapterfilt.sh'
    cmd += ' -p %s' % sam
    cmd += ' -l %s' % self.soft.params['l']
    cmd += ' -m %s' % self.soft.params['m']
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -o %s\n' % out_dir

    for fastq in fastqs:
        cmd += 'rm %s/%s\n' % (out_dir, basename(fastq))
    return cmd


def hifiadapterfilt(self) -> None:
    """Convert .bam to .fastq and remove reads with remnant PacBio adapter
    sequences.

    References
    ----------
    Sim, Sheina B., et al. "HiFiAdapterFilt, a memory efficient read
    processing pipeline, prevents occurrence of adapter sequence in PacBio
    HiFi reads and their negative impacts on genome assembly." BMC genomics
    23.1 (2022): 1-7.

    Notes
    -----
    GitHub  : https://github.com/sheinasim/HiFiAdapterFilt
    Paper   : https://doi.org/10.1186/s12864-022-08375-1

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for HiFiAdapterFilt
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
        if tech_specificity(self, fastqs, tech, sam, ['pacbio']):
            continue
        to_dos = status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        out = '%s/%s.filt.fastq.gz' % (out_dir, sam)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(out)

        if self.config.force or to_do(out):
            cmd = hifiadapterfilt_cmd(self, fastqs, out_dir, sam)
            key = (tech, sam)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fastqs, o_d=out_dir, key=key)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def kneaddata_cmd(
        self,
        tech: str,
        inputs: list,
        out_dir: str
) -> tuple:
    """Collect the Kneaddata commands.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .outputs : dict
            All outputs
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Path to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Kneaddata commands
    outputs : list
        Path to the output files
    """
    params = tech_params(self, tech)
    cmd = 'kneaddata'
    for inp in inputs:
        cmd += ' --input %s' % inp
    for database in params['databases']:
        cmd += ' --reference-db %s' % database
    cmd += ' --output %s' % out_dir
    cmd += ' --processes %s' % params['cpus']
    cmd += ' --max-memory %s%s' % (params['mem'], params['mem_dim'][0])
    if params['trimmomatic']:
        cmd += ' --trimmomatic %s' % params['trimmomatic']
    if params['bowtie2']:
        cmd += ' --bowtie2 %s' % params['bowtie2']
    cmd += ' --bowtie2-options="--very-sensitive-local -p %s"' % params['cpus']
    cmd += ' --remove-intermediate-output'
    cmd += ' --log %s/%s.log\n' % (out_dir, self.sam_pool)
    # write these additional output files manipulations commands
    rad = basename(inputs[0]).replace('.fastq.gz', '')
    outputs = []
    for typ in ['paired', 'unmatched']:
        for orient in ['1', '2']:
            out = '%s/%s_kneaddata_%s_%s.fastq' % (out_dir, rad, typ, orient)
            r = out.replace('.R1', '')
            self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(r)
            if out == r:
                continue
            outputs.append(r)
            cmd += 'mv %s %s\n' % (out, r)
    if params['purge']:
        cmd += 'rm %s/%s*\n' % (out_dir, rad)
    else:
        cmd += 'tar cpfz %s/%s_kneadfiles.tar.gz %s/%s* --remove-files\n' % (
            out_dir, self.sam_pool, out_dir, rad)
    return cmd, outputs


def kneaddata(self) -> None:
    """KneadData is a tool designed to perform quality control on metagenomic
    and metatranscriptomic sequencing data, especially data from microbiome
    experiments. In these experiments, samples are typically taken from a
    host in hopes of learning something about the microbial community on the
    host. However, sequencing data from such experiments will often contain a
    high ratio of host to bacterial reads. This tool aims to perform
    principled in silico separation of bacterial reads from these
    "contaminant" reads, be they from the host, from bacterial 16S sequences,
    or other user-defined sources. Additionally, KneadData can be used for
    other filtering tasks. For example, if one is trying to clean data
    derived from a human sequencing experiment, KneadData can be used to
    separate the human and the non-human reads.

    Notes
    -----
    GitHub  : https://github.com/biobakery/kneaddata
    Docs    : https://huttenhower.sph.harvard.edu/kneaddata

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Kneaddata
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
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        to_dos = status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outputs = kneaddata_cmd(self, tech, fastqs, out_dir)

        if self.config.force or sum([to_do(x) for x in outputs]):
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def filtering_cmd(
        self,
        sam: str,
        tech: str,
        fastqs: list,
        out_dir: str,
        outs: list,
) -> str:
    """Collect the command line for filtering.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    sam : str
        Name of the sample
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Path to the input fastq files
    out_dir : str
        Path to the output folder
    outs : list
        Path to the output fastq.gz files

    Returns
    -------
    cmd : str
        filtering commands
    """
    dbs = self.soft.params['databases']
    aligner = self.soft.params['aligner']
    params = tech_params(self, tech, aligner)
    inputs = sorted(fastqs)
    cmd = ''
    for dx, db in enumerate(dbs):
        db_path = self.databases.builds[db][aligner]
        bam = '%s/%s.bam' % (out_dir, db)
        cmd += globals()['%s_cmd' % aligner](
            sam, inputs, db, db_path, out_dir, params)
        cmd += ' | samtools view -b -f 12 -F 256 -'
        cmd += ' | samtools sort -@ %s -n - > %s\n' % (params['cpus'], bam)
        cmd += 'bedtools bamtofastq -i %s' % bam
        out1 = '%s/%s_R1.fastq' % (out_dir, db)
        cmd += ' -fq %s' % out1
        fastqs = [out1]
        if len(inputs) > 1:
            out2 = '%s/%s_R2.fastq' % (out_dir, db)
            cmd += ' -fq2 %s\n' % out2
            fastqs.append(out2)
        cmd += 'rm %s\n' % bam
        if dx:
            cmd += 'rm %s\n' % ' '.join(inputs)
        inputs = list(fastqs)

    if fastqs:
        cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
        cmd += 'mv %s.gz %s\n' % (fastqs[0], outs[0])
        if len(fastqs) > 1:
            cmd += 'mv %s.gz %s\n' % (fastqs[1], outs[1])
    return cmd


def filtering(self):
    """Perform the filtering of fastq files so that every sequences that
    aligns to any of the passed databases is discarded.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for filtering
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
        if tech_specificity(self, fastqs, tech, sam):
            continue
        to_dos = status_update(self, tech, fastqs)
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        if len(fastqs) == 1:
            outs = ['%s/%s.fastq.gz' % (out_dir, self.sam_pool)]
        else:
            outs = ['%s/%s_R1.fastq.gz' % (out_dir, self.sam_pool),
                    '%s/%s_R2.fastq.gz' % (out_dir, self.sam_pool)]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)
        if self.config.force or sum([to_do(x) for x in outs]):
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                cmd = filtering_cmd(self, sam, tech, fastqs, out_dir, outs)
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, i_d=out_dir, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def berokka(self):
    """Trim, circularise, orient & filter long read bacterial genome assemblies.

    Useful if you only have the contig files and do not have the corrected
    reads anymore, if your contigs are simple cases with clear overhang and
    could be done manually with BLAST, and if circlator fails on your data even
    after troubleshooting.

    Notes
    -----
    GitHub  : https://github.com/tseemann/berokka

    Parameters
    ----------
    self : Commands class instance
    """
