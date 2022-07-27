# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from os.path import basename
from metagenomix._io_utils import io_update, to_do
from metagenomix.parameters import tech_params


def get_cat_zcat(fastq_fp: str) -> str:
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


def edit_fastq_cmd(fastq_fp: str, num: int) -> str:
    """Get the unix command to run on each fastq file that
    needs editing, depending on the source of fastq file.

    Parameters
    ----------
    fastq_fp : str
        Fastq file path.
    num : int
        Can be 1 or 2 for the forward and reverse read, respectively.
    source : str
        Source of the fastq that defines format and editing type.

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


def get_fastq_header(fastq_path: str) -> str:
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


def get_edit_cmd(self, num: int, tech: str) -> None:
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
        Technology
    """
    fastqs_fps = self.inputs[self.sam][tech]
    fastq_fp = fastqs_fps[num - 1]
    line = get_fastq_header(self.config.fastq[self.sam][tech][num - 1])
    line_split = line.split()
    if not line_split[0].endswith('/%s' % num):
        if len(line_split) > 1:
            cmd = edit_fastq_cmd(fastq_fp, num)
            # self.outputs['cmds'].append(cmd)
            self.outputs['cmds'].setdefault(tech, []).append(cmd)
    # io_update(self, i_f=fastqs_fps, o_f=fastqs_fps)
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
        if len(self.config.fastq[self.sam][tech]) < r:
            continue
        get_edit_cmd(self, r, tech)


def count_cmd(self, idx: int, input_fp: str, out: str) -> str:
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
    input_fp : str
        Path to the input file name.
    out : str
        Path to the output file name.

    Returns
    -------
    cmd : str
        Command to count reads and output the counts.
    """
    if input_fp.endswith('fastq.gz'):
        cmd = "n%s=`zcat %s | wc -l | " % (idx, input_fp)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
    elif input_fp.endswith('fasta'):
        cmd = "n%s=`wc -l %s | " % (idx, input_fp)
        cmd += "sed 's/ //g' | awk '{x=$1/2; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/2; print x}'`"
    elif input_fp.endswith('fastq'):
        cmd = "n%s=`wc -l %s | " % (idx, input_fp)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
    else:
        raise IOError("Input sequence file invalid %s" % input_fp)
    if idx:
        cmd += 'echo "%s,2,${n%s}" >> %s\n' % (self.sam, idx, out)
    else:
        cmd += 'echo "%s,1,${n%s}" > %s\n' % (self.sam, idx, out)
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
    for tech, inputs in self.inputs[self.sam].items():
        out_dir = '%s/%s' % (self.dir, tech)
        self.outputs['dirs'].append(out_dir)
        out = '%s/%s_read_count.tsv' % (out_dir, self.sam)
        outs[tech] = out
        if self.config.force or to_do(out):
            for idx, input_path in enumerate(inputs):
                cmd = count_cmd(self, idx, input_path, out)
                # self.outputs['cmds'].append(cmd)
                self.outputs['cmds'].setdefault(tech, []).append(cmd)
            # io_update(self, i_f=inputs, o_f=out)
            io_update(self, i_f=inputs, o_f=out, key=tech)
    self.outputs['outs'] = outs


def fastqc(self) -> None:
    outs = {}
    for tech, inputs in self.inputs[self.sam].items():
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
        self.outputs['dirs'].append(out_dir)
        out = ['%s_fastqc.html' % x.rsplit('.fastq', 1)[0] for x in inputs]
        outs[tech] = out
        if self.config.force or not sum([to_do(x) for x in out]):
            cmd = 'fastqc %s -o %s' % (' '.join(inputs), out_dir)
            # self.outputs['cmds'].append(cmd)
            if tech in self.outputs['cmds']:
                print(self.outputs['cmds'][tech])
                print(selfoutputscmdstech)
            self.outputs['cmds'][tech] = [cmd]
            # io_update(self, i_f=inputs, o_d=out_dir)
            io_update(self, i_f=inputs, o_d=out_dir, key=tech)
    self.outputs['outs'] = outs


def fastp_cmd(self, tech, fastqs, out_dir):
    params = tech_params(self, tech)
    outs = []
    cmd = 'fastp'
    for fdx, fastq in enumerate(fastqs):
        if fdx:
            i, o = 'I', 'O'
            out = '%s/%s_2.fastq' % (out_dir, self.sam)
        else:
            i, o = 'i', 'o'
            out = '%s/%s_1.fastq' % (out_dir, self.sam)
        if fastq.endswith('.gz'):
            out += '.gz'
        outs.append(out)
        cmd += '  -%s %s -%s %s' % (i, fastq, o, out)
    cmd += ' --thread %s' % params['cpus']
    cmd += ' --report_title="%s"' % self.sam
    cmd += ' --json=%s/%s.json' % (out_dir, self.sam)
    cmd += ' --html=%s/%s.html' % (out_dir, self.sam)

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
    for tech, fastqs in self.inputs[self.sam].items():
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
        self.outputs['dirs'].append(out_dir)
        cmd, outs = fastp_cmd(self, tech, fastqs, out_dir)
        self.outputs['outs'].setdefault(tech, []).extend(outs)
        if self.config.force or sum([to_do(x) for x in outs]):
            # self.outputs['cmds'].append(cmd)
            if tech in self.outputs['cmds']:
                print(self.outputs['cmds'][tech])
                print(selfoutputscmdstech)
            self.outputs['cmds'][tech] = [cmd]
            # io_update(self, i_f=fastqs, o_d=out_dir)
            io_update(self, i_f=fastqs, o_d=out_dir, key=tech)


def cutadapt(self) -> None:
    tech = 'illumina'
    out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
    self.outputs['dirs'].append(out_dir)
    r1_o = '%s/%s.R1.fastq.gz' % (out_dir, self.sam)
    r2_o = '%s/%s.R2.fastq.gz' % (out_dir, self.sam)
    cmd = 'cutadapt'
    cmd += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    cmd += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    cmd += ' -m 10'
    cmd += ' -o %s' % r1_o
    cmd += ' -p %s ' % r2_o
    cmd += ' '.join(self.inputs[self.sam][tech])
    outs = [r1_o, r2_o]
    self.outputs['outs'].setdefault(tech, []).extend(outs)
    if self.config.force or to_do(r1_o) or to_do(r2_o):
        # self.outputs['cmds']] = list([cmd])
        if tech in self.outputs['cmds']:
            print(self.outputs['cmds'][tech])
            print(selfoutputscmdstech)
        self.outputs['cmds'][tech] = [cmd]
        # io_update(self, i_f=self.inputs[self.sam][tech], o_f=outs)
        io_update(self, i_f=self.inputs[self.sam][tech], o_f=outs, key=tech)


def atropos_outs(self, inputs, out_dir) -> list:
    if len(inputs) == 1:
        out1 = '%s/%s.fastq' % (out_dir, self.sam)
    else:
        out1 = '%s/%s_R1.fastq' % (out_dir, self.sam)
    out2 = '%s/%s_R2.fastq' % (out_dir, self.sam)
    if inputs[0].endswith('.gz'):
        out1 += '%s.gz'
        out2 += '%s.gz'
    return [out1, out2]


def atropos_cmd(self, tech, inputs, out_dir, outputs):
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
    cmd += ' --report-file %s/%s.log' % (out_dir, self.sam)
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


def atropos(self):
    tech = 'illumina'
    out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
    self.outputs['dirs'].append(out_dir)
    inputs = self.inputs[self.sam][tech]
    outputs = atropos_outs(self, inputs, out_dir)
    cmd = atropos_cmd(self, tech, inputs, out_dir, outputs)
    if self.config.force or sum([to_do(x) for x in outputs]):
        # self.outputs['cmds'].append(cmd)
        if tech in self.outputs['cmds']:
            print(self.outputs['cmds'][tech])
            print(selfoutputscmdstech)
        self.outputs['cmds'][tech] = [cmd]
        # io_update(self, i_f=inputs, o_f=outputs)
        io_update(self, i_f=inputs, o_f=outputs, key=tech)
    self.outputs['outs'].setdefault(tech, []).extend(outputs)


def kneaddata_cmd(self, tech, inputs, out_dir) -> tuple:
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
    cmd += ' --log %s/%s.log\n' % (out_dir, self.sam)
    # write these additional output files manipulations commands
    rad = basename(inputs[0]).replace('.fastq.gz', '')
    outputs = []
    for typ in ['paired', 'unmatched']:
        for r in ['1', '2']:
            out = '%s/%s_kneaddata_%s_%s.fastq' % (out_dir, rad, typ, r)
            renamed = out.replace('.R1', '')
            self.outputs['outs'].setdefault(tech, []).append(renamed)
            if out == renamed:
                continue
            outputs.append(renamed)
            cmd += 'mv %s %s\n' % (out, renamed)
    if params['purge']:
        cmd += 'rm %s/%s*\n' % (out_dir, rad)
    else:
        cmd += 'tar cpfz %s/%s_kneadfiles.tar.gz %s/%s* --remove-files\n' % (
            out_dir, self.sam, out_dir, rad)
    return cmd, outputs


def kneaddata(self):
    tech = 'illumina'
    out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
    self.outputs['dirs'].append(out_dir)
    inputs = self.inputs[self.sam][tech]
    cmd, outputs = kneaddata_cmd(self, tech, inputs, out_dir)
    if self.config.force or sum([to_do(x) for x in outputs]):
        # self.outputs['cmds'].append(cmd)
        if tech in self.outputs['cmds']:
            print(self.outputs['cmds'][tech])
            print(selfoutputscmdstech)
        self.outputs['cmds'][tech] = [cmd]
        # io_update(self, i_f=inputs, o_d=out_dir)
        io_update(self, i_f=inputs, o_d=out_dir, key=tech)


def filtering_cmd(self, tech, db_path, inputs, out1, out2) -> tuple:
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
    outs = [out1.replace('fastq.gz', 'fastq')]
    if len(inputs) > 1:
        cmd += ' -fq2 %s' % out2.replace('fastq.gz', 'fastq')
        outs.append(out2.replace('fastq.gz', 'fastq'))
    cmd += '\n'
    return cmd, outs


def filtering(self):
    for tech in self.inputs[self.sam]:
        outs = self.inputs[self.sam][tech]
        if not outs:
            continue
        cmds = ''
        databases = self.soft.params['databases']
        for ddx, (db, db_path) in enumerate(databases.items()):
            inputs = list(outs)
            out_dir = '%s/%s/%s_%s/%s' % (self.dir, tech, ddx, db, self.sam)
            self.outputs['dirs'].append(out_dir)
            if len(inputs) == 1:
                out1 = '%s/%s.fastq' % (out_dir, self.sam)
            else:
                out1 = '%s/%s_R1.fastq' % (out_dir, self.sam)
            out2 = '%s/%s_R2.fastq' % (out_dir, self.sam)
            if inputs[0].endswith('.gz'):
                out1 += '.gz'
                out2 += '.gz'
            cmd, outs = filtering_cmd(self, tech, db_path, inputs, out1, out2)
            if self.config.force or to_do(out1):
                cmds += cmd
                if ddx:
                    cmds += '\nrm %s' % ' '.join(inputs)
        if outs[0].endswith('.fastq'):
            for output in outs:
                cmds += '\ngzip %s' % output
        outputs_gz = ['%s.gz' % x for x in outs]
        self.outputs['outs'].setdefault(tech, []).extend(outputs_gz)
        if cmds:
            # self.outputs['cmds'].append(cmds)
            if tech in self.outputs['cmds']:
                print(self.outputs['cmds'][tech])
                print(selfoutputscmdstech)
            self.outputs['cmds'][tech] = [cmds]
            # io_update(self, o_f=outputs_gz)
            io_update(
                self, i_f=self.inputs[self.sam][tech], o_f=outputs_gz, key=tech)
