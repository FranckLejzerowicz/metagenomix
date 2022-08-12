# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from os.path import basename
from metagenomix._io_utils import io_update, to_do, tech_specificity
from metagenomix.parameters import tech_params


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
    cmd = '%s %s | ' % (cat, fastq_fp)
    cmd += "bioawk -c fastx "
    cmd += "'{print \"@\"$1\"/%s\\n\"$2\"\\n+\\n\"$3}' " % str(num)
    if fastq_fp.endswith('.fastq'):
        cmd += " > %s_renamed\n" % fastq_fp
    else:
        cmd += " | gzip > %s_renamed\n" % fastq_fp
    cmd += "mv %s_renamed %s" % (fastq_fp, fastq_fp)
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
    fastq_fp = fastqs_fps[num - 1]
    line = get_fastq_header(
        self.config.fastq[self.sam_pool][(tech, self.sam_pool)][num - 1])
    line_split = line.split()
    if not line_split[0].endswith('/%s' % num):
        if len(line_split) > 1:
            cmd = edit_fastq_cmd(fastq_fp, num)
            self.outputs['cmds'].setdefault(tech, []).append(cmd)
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
        out_dir = '%s/%s' % (self.dir, tech)
        self.outputs['dirs'].append(out_dir)
        out = '%s/%s_read_count.tsv' % (out_dir, self.sam_pool)
        outs[(tech, self.sam_pool)] = out
        if self.config.force or to_do(out):
            for idx, fastx in enumerate(fastxs):
                cmd = count_cmd(self, idx, fastx, out)
                self.outputs['cmds'].setdefault(tech, []).append(cmd)
            io_update(self, i_f=fastxs, o_f=out, key=tech)
    self.outputs['outs'] = outs


def fastqc(self) -> None:
    """Check quality based on the FastQC software.

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
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        out = ['%s_fastqc.html' % x.rsplit('.fastq', 1)[0] for x in fastxs]
        outs[(tech, self.sam_pool)] = out
        if self.config.force or not sum([to_do(x) for x in out]):
            cmd = 'fastqc %s -o %s' % (' '.join(fastxs), out_dir)
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastxs, o_d=out_dir, key=tech)
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
    """Filter the fastq files based on the Fastp software.

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
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        cmd, outs = fastp_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)
        if self.config.force or sum([to_do(x) for x in outs]):
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)


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
    """Filter the fastq files based on the Cutadapt software.

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
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)

        cmd, outs = cutadapt_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or sum([to_do(x) for x in outs]):
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastqs, o_f=outs, key=tech)


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
    """Filter the fastq files based on the Atropos software.

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
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        outs = atropos_outs(self, fastqs, out_dir)
        cmd = atropos_cmd(self, tech, fastqs, out_dir, outs)
        if self.config.force or sum([to_do(x) for x in outs]):
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastqs, o_f=outs, key=tech)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)


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
    """Filter the fastq files based on the BioBakery Kneaddata software.

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
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out_dir)
        cmd, outputs = kneaddata_cmd(self, tech, fastqs, out_dir)
        if self.config.force or sum([to_do(x) for x in outputs]):
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)


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
                cmds += '\ngzip %s' % fastq
        fastqs_gz = ['%s.gz' % x for x in fastqs]
        self.outputs['outs'].setdefault(
            (tech, self.sam_pool), []).extend(fastqs_gz)
        if cmds:
            self.outputs['cmds'][tech] = [cmds]
            io_update(self, i_f=fastqs_, i_d=out_dirs, o_f=fastqs_gz, key=tech)
