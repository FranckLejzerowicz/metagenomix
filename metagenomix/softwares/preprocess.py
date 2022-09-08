# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from os.path import basename
from metagenomix._io_utils import (io_update, to_do, tech_specificity,
                                   status_update)
from metagenomix.core.parameters import tech_params


def get_cat_zcat(
        fastq_fp: str
) -> str:
    """Get the command `cat` or `gzcat` to run for fastq editing.

    Parameters
    ----------
    fastq_fp : str
        Fastq file path.

    Returns
    -------
    cat : str
        The command to use (`cat` or `zcat`).
    """
    cat = "cat"
    if fastq_fp.endswith('.fastq.gz'):
        cat = "zcat"
    return cat


def edit_fastq_cmd(
        fastq_fp: str,
        num: int
) -> str:
    """Get the unix command to run on each fastq file that
    needs editing, depending on the source of fastq file.

    Parameters
    ----------
    fastq_fp : str
        Fastq file path.
    num : int
        Can be 1 or 2 for the forward and reverse read, respectively.

    Returns
    -------
    cmd : str
        Unix command to run to edit a fastq file.
    """
    cat = get_cat_zcat(fastq_fp)

    cmd = 'bioawk_exe=$(which bioawk)\n'
    cmd += 'if [ -x "$bioawk_exe" ] ; then\n'
    cmd += '%s %s | ' % (cat, fastq_fp)
    cmd += "bioawk -c fastx "
    cmd += "'{print \"@\"$1\"/%s\\n\"$2\"\\n+\\n\"$3}' " % str(num)
    if fastq_fp.endswith('.fastq'):
        cmd += " > %s_renamed\n" % fastq_fp
    else:
        cmd += " | gzip > %s_renamed\n" % fastq_fp
    cmd += "mv %s_renamed %s\n" % (fastq_fp, fastq_fp)
    cmd += 'else\n'
    cmd += 'echo "Need to load bioawk module or environment with bioawk"\n'
    cmd += 'exit 1\n'
    cmd += 'fi\n'
    return cmd


def get_fastq_header(
        fastq_path: str
) -> str:
    """Get the first line of the fastq file.

    Parameters
    ----------
    fastq_path : str
        Path to a fastq.gz or fastq file.

    Returns
    -------
    fastq_line : str
        First line of the fastq file.
    """
    fastq_line = ''
    if fastq_path.endswith('fastq.gz'):
        with gzip.open(fastq_path) as f:
            for line in f:
                fastq_line = line.decode().strip()
                break
    elif fastq_path.endswith('fastq'):
        with open(fastq_path) as f:
            for line in f:
                fastq_line = line.strip()
                break
    return fastq_line


def get_edit_cmd(
        self,
        num: int,
        tech: str
) -> None:
    """Get the unix command to run on each fastq file that needs editing.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .config
            Configurations
    num : int
        Can be 1 or 2 for the forward and reverse read, respectively.
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    """
    fastqs_fps = self.inputs[self.sam_pool][(tech, self.sam_pool)]
    status_update(self, tech, fastqs_fps)

    fastq_fp = fastqs_fps[num - 1]
    line = get_fastq_header(
        self.config.fastq[self.sam_pool][(tech, self.sam_pool)][num - 1])
    line_split = line.split()
    if not line_split[0].endswith('/%s' % num):
        if len(line_split) > 1:
            cmd = edit_fastq_cmd(fastq_fp, num)
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)
            self.soft.add_status(tech, self.sam_pool, 1)
        else:
            self.soft.add_status(tech, self.sam_pool, 0)
    io_update(self, i_f=fastqs_fps, o_f=fastqs_fps, key=tech)


def edit(self) -> None:
    """Create command lines for editing the fasta if need be so that
    the fastq header all contaim the trailing "/1" and "/2" that allow
    identifying paired reads.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .config
            Configurations
    """
    tech = 'illumina'
    for r in [1, 2]:
        if len(self.config.fastq[self.sam_pool][(tech, self.sam_pool)]) < r:
            continue
        get_edit_cmd(self, r, tech)


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
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
    elif fastx.endswith('fasta'):
        cmd = "n%s=`wc -l %s | " % (idx, fastx)
        cmd += "sed 's/ //g' | awk '{x=$1/2; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/2; print x}'`"
    elif fastx.endswith('fastq'):
        cmd = "n%s=`wc -l %s | " % (idx, fastx)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
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
        status_update(self, tech, fastxs)

        out_dir = '%s/%s' % (self.dir, tech)
        self.outputs['dirs'].append(out_dir)

        out = '%s/%s_read_count.tsv' % (out_dir, self.sam_pool)
        outs[(tech, self.sam_pool)] = out

        if self.config.force or to_do(out):
            for idx, fastx in enumerate(fastxs):
                cmd = count_cmd(self, idx, fastx, out)
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
        status_update(self, tech, fastxs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        out = ['%s_fastqc.html' % x.rsplit('.fastq', 1)[0] for x in fastxs]
        outs[(tech, self.sam_pool)] = out

        if self.config.force or not sum([to_do(x) for x in out]):
            cmd = 'fastqc %s -o %s' % (' '.join(fastxs), out_dir)
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastxs, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)

    self.outputs['outs'] = outs


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
    params = tech_params(self, tech)
    outs = []
    cmd = 'fastp'
    for fdx, fastq in enumerate(fastqs):
        if fdx:
            i, o = 'I', 'O'
            out = '%s/%s_2.fastq' % (out_dir, self.sam_pool)
        else:
            i, o = 'i', 'o'
            out = '%s/%s_1.fastq' % (out_dir, self.sam_pool)
        if fastq.endswith('.gz'):
            out += '.gz'
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
                    'average_qual', 'cut_front', 'cut_tail', 'cut_right']:
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
        status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outs = fastp_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
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
        status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outs = cutadapt_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
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
        status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        outs = atropos_outs(self, fastqs, out_dir)

        cmd = atropos_cmd(self, tech, fastqs, out_dir, outs)
        if self.config.force or sum([to_do(x) for x in outs]):
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
    cmd = 'export PATH=$PATH:%s\n' % self.soft.params['path']
    cmd += 'export PATH=$PATH:%s/DB\n' % self.soft.params['path']
    cmd += 'cd %s\n' % out_dir
    cmd += 'cp %s %s/.\n' % (' '.join(fastqs), out_dir)
    cmd += 'bash %s/pbadapterfilt.sh' % self.soft.params['path']
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
        status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        out = '%s/%s.filt.fastq.gz' % (out_dir, sam)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(out)

        if self.config.force or to_do(out):
            cmd = hifiadapterfilt_cmd(self, fastqs, out_dir, sam)
            key = (tech, sam)
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
        status_update(self, tech, fastqs)

        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outputs = kneaddata_cmd(self, tech, fastqs, out_dir)

        if self.config.force or sum([to_do(x) for x in outputs]):
            self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def filtering_cmd(
        self,
        tech: str,
        db_path: str,
        inputs: list,
        out1: str,
        out2: str
) -> tuple:
    """Collect the command line for filtering.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    db_path : str
        Path to the reference bowtie2 database
    inputs : list
        Path to the input fastq files
    out1 : str
        Path to the output fastq.gz files
    out2 : str
        Path to the output fastq.gz files

    Returns
    -------
    cmd : str
        filtering commands
    fastqs : list
        Path to the output fastq files
    """
    params = tech_params(self, tech)
    cmd = '\nbowtie2'
    cmd += ' -p %s' % params['cpus']
    cmd += ' -x %s' % db_path
    cmd += ' --very-sensitive'
    if len(inputs) == 1:
        cmd += ' -U %s' % inputs[0]
    else:
        cmd += ' -1 %s' % inputs[0]
        cmd += ' -2 %s' % inputs[1]
    cmd += ' | samtools view -f 12 -F 256'
    cmd += ' | samtools sort -@ %s -n' % params['cpus']
    cmd += ' | samtools view -bS'
    cmd += ' | bedtools bamtofastq -i -'
    cmd += ' -fq %s' % out1.replace('fastq.gz', 'fastq')
    fastqs = [out1.replace('fastq.gz', 'fastq')]
    if len(inputs) > 1:
        cmd += ' -fq2 %s' % out2.replace('fastq.gz', 'fastq')
        fastqs.append(out2.replace('fastq.gz', 'fastq'))
    cmd += '\n'
    return cmd, fastqs


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
    for (tech, sam), fastqs_ in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs_, tech, sam):
            continue
        status_update(self, tech, fastqs_)

        fastqs = fastqs_
        cmds = ''
        out_dirs = []
        databases = self.soft.params['databases']
        for dx, (db, db_path) in enumerate(databases.items()):
            inputs = list(fastqs)
            out_dir = '%s/%s/%s_%s/%s' % (self.dir, tech, dx, db, self.sam_pool)
            out_dirs.append(out_dir)
            self.outputs['dirs'].append(out_dir)
            if len(inputs) == 1:
                out1 = '%s/%s.fastq' % (out_dir, self.sam_pool)
            else:
                out1 = '%s/%s_R1.fastq' % (out_dir, self.sam_pool)
            out2 = '%s/%s_R2.fastq' % (out_dir, self.sam_pool)
            if inputs[0].endswith('.gz'):
                out1 += '.gz'
                out2 += '.gz'
            cmd, fastqs = filtering_cmd(self, tech, db_path, inputs, out1, out2)
            if self.config.force or to_do(out1):
                cmds += cmd
                if dx:
                    cmds += '\nrm %s' % ' '.join(inputs)
        if fastqs[0].endswith('.fastq'):
            for fastq in fastqs:
                cmds += '\ngzip -f %s' % fastq
        fastqs_gz = ['%s.gz' % x for x in fastqs]
        self.outputs['outs'].setdefault(
            (tech, self.sam_pool), []).extend(fastqs_gz)
        if (self.config.force or sum([to_do(x) for x in fastqs_gz])) and cmds:
            self.outputs['cmds'][(tech,)] = [cmds]
            io_update(self, i_f=fastqs_, i_d=out_dirs, o_f=fastqs_gz, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


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
    pass
