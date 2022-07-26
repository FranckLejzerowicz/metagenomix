# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
from os.path import basename
from metagenomix._io_utils import io_update, todo


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
            self.outputs['cmds'].append(cmd)
    io_update(self, i_f=fastqs_fps, o_f=fastqs_fps)


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
        cmd = "n%s=`%s %s | wc -l | " % (idx, self.soft.params['cat'], input_fp)
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
        if self.config.force or todo(out):
            for idx, input_path in enumerate(inputs):
                cmd = count_cmd(self, idx, input_path, out)
                self.outputs['cmds'].append(cmd)
            io_update(self, i_f=inputs, o_f=out)
    self.outputs['outs'] = outs


def fastqc(self) -> None:
    outs = {}
    for tech, inputs in self.inputs[self.sam].items():
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
        self.outputs['dirs'].append(out_dir)
        out = ['%s_fastqc.html' % x.rsplit('.fastq', 1)[0] for x in inputs]
        outs[tech] = out
        if self.config.force or not sum([todo(x) for x in out]):
            cmd = 'fastqc %s -o %s' % (' '.join(inputs), out_dir)
            self.outputs['cmds'].append(cmd)
            io_update(self, i_f=inputs, o_d=out_dir)
    self.outputs['outs'] = outs


def fastp(self) -> None:
    for tech, fastqs in self.inputs[self.sam].items():
        out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
        self.outputs['dirs'].append(out_dir)
        out_json = '%s/%s.json' % (out_dir, self.sam)
        out_html = '%s/%s.html' % (out_dir, self.sam)
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
        cmd += ' --verbose'
        if self.soft.params['trim_poly_g']:
            cmd += ' --trim_poly_g'
        cmd += ' --json=%s' % out_json
        cmd += ' --html=%s' % out_html
        cmd += ' --thread %s' % self.soft.params['cpus']
        cmd += ' --report_title="%s"' % self.sam
        cmd += ' --average_qual=%s' % self.soft.params['average_qual']

        if self.soft.params['split_by_lines']:
            cmd += ' --split_by_lines'
            cmd += ' --split_prefix_digits %s' % self.soft.params[
                'split_prefix_digits']
        elif self.soft.params['split']:
            cmd += ' --split'
            cmd += ' --split_prefix_digits %s' % self.soft.params[
                'split_prefix_digits']

        if self.soft.params['overrepresentation_analysis']:
            cmd += ' --overrepresentation_analysis'
            cmd += ' --overrepresentation_sampling %s' % self.soft.params[
                'overrepresentation_sampling']

        if self.soft.params['correction']:
            cmd += ' --correction'
            cmd += ' --overlap_len_require %s' % self.soft.params[
                'overlap_len_require']
            cmd += ' --overlap_diff_limit %s' % self.soft.params[
                'overlap_diff_limit']

        if self.soft.params['disable_length_filtering']:
            cmd += ' --disable_length_filtering'
        else:
            cmd += ' --length_required %s' % self.soft.params['length_required']
            cmd += ' --length_limit %s' % self.soft.params['length_limit']


        self.outputs['outs'].setdefault(tech, []).extend(outs)
        if self.config.force or sum([todo(x) for x in outs]):
            self.outputs['cmds'].append(cmd)
            io_update(self, i_f=fastqs, o_d=out_dir)


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
    if self.config.force or todo(r1_o) or todo(r2_o):
        self.outputs['cmds'] = list([cmd])
        io_update(self, i_f=self.inputs[self.sam][tech], o_f=outs)


def atropos(self):
    tech = 'illumina'
    inputs = self.inputs[self.sam][tech]
    out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
    self.outputs['dirs'].append(out_dir)
    if len(inputs) == 1:
        out_1 = '%s/%s.fastq' % (out_dir, self.sam)
    else:
        out_1 = '%s/%s_R1.fastq' % (out_dir, self.sam)
    out_2 = '%s/%s_R2.fastq' % (out_dir, self.sam)
    log = '%s/%s.log' % (out_dir, self.sam)
    if inputs[0].endswith('.gz'):
        out_1 += '%s.gz'
        out_2 += '%s.gz'
    cmd = 'atropos trim'
    for a in self.soft.params['a']:
        cmd += ' -a %s' % a
    for A in self.soft.params['A']:
        cmd += ' -A %s' % A
    cmd += ' -q %s' % self.soft.params['q']
    cmd += ' --minimum-length %s' % self.soft.params['minimum_length']
    cmd += ' --nextseq-trim %s' % self.soft.params['nextseq_trim']
    cmd += ' --pair-filter %s' % self.soft.params['pair_filter']
    cmd += ' --overlap %s' % self.soft.params['overlap']
    cmd += ' --max-reads %s' % self.soft.params['max_reads']
    cmd += ' --indel-cost %s' % self.soft.params['indel_cost']
    cmd += ' --error-rate %s' % self.soft.params['error_rate']
    cmd += ' --quality-cutoff %s' % self.soft.params['quality_cutoff']
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --report-file %s' % log
    cmd += ' --report-formats %s' % self.soft.params['report_formats']
    if len(inputs) == 1:
        cmd += ' -se %s' % inputs[0]
        cmd += ' -o %s' % out_1
    else:
        cmd += ' -pe1 %s' % inputs[0]
        cmd += ' -pe2 %s' % inputs[1]
        cmd += ' -o %s' % out_1
        cmd += ' -p %s' % out_2
    if self.config.force or todo(out_1) or todo(out_2):
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, o_f=[out_1, out_2])
    self.outputs['outs'].setdefault(tech, []).extend([out_1, out_2])


def kneaddata(self):
    tech = 'illumina'
    inputs = self.inputs[self.sam][tech]
    out_dir = '%s/%s/%s' % (self.dir, tech, self.sam)
    self.outputs['dirs'].append(out_dir)

    cmd = 'kneaddata'
    for inp in inputs:
        cmd += ' --input %s' % inp
    for database in self.soft.params['databases']:
        cmd += ' --reference-db %s' % database
    cmd += ' --output %s' % out_dir
    # cmd += ' --output-prefix %s' % self.sam
    cpus = self.soft.params['cpus']
    cmd += ' --processes %s' % cpus
    cmd += ' --max-memory %s%s' % (self.soft.params['mem'],
                                   self.soft.params['mem_dim'][0])
    if self.soft.params['trimmomatic']:
        cmd += ' --trimmomatic %s' % self.soft.params['trimmomatic']
    if self.soft.params['bowtie2']:
        cmd += ' --bowtie2 %s' % self.soft.params['bowtie2']
    cmd += ' --bowtie2-options="--very-sensitive-local -p %s"' % cpus
    cmd += ' --remove-intermediate-output'
    cmd += ' --log %s/%s.log' % (out_dir, self.sam)
    self.outputs['cmds'].append(cmd)
    # write these additional output files manipulations commands
    rad = basename(inputs[0]).replace('.fastq.gz', '')
    for typ in ['paired', 'unmatched']:
        for r in ['1', '2']:
            out = '%s/%s_kneaddata_%s_%s.fastq' % (out_dir, rad, typ, r)
            renamed = out.replace('.R1', '')
            self.outputs['outs'].setdefault(tech, []).append(renamed)
            if out == renamed:
                continue
            self.outputs['cmds'].append('mv %s %s' % (out, renamed))

    if self.soft.params['purge']:
        tar_cmd = 'rm %s/%s*' % (out_dir, rad)
    else:
        tar_cmd = 'tar cpfz %s/%s_kneadfiles.tar.gz %s/%s* --remove-files' % (
            out_dir, self.sam, out_dir, rad)

    self.outputs['cmds'].append(tar_cmd)
    io_update(self, i_f=inputs, o_d=out_dir)


def filtering(self):
    for tech in self.inputs[self.sam]:
        outputs = self.inputs[self.sam][tech]
        if not outputs:
            continue
        io_update(self, i_f=self.inputs[self.sam][tech])
        cmds = ''
        databases = self.soft.params['databases']
        for ddx, (db, db_path) in enumerate(databases.items()):
            inputs = list(outputs)
            outputs = []
            out_dir = '%s/%s/%s_%s/%s' % (self.dir, tech, ddx, db, self.sam)
            self.outputs['dirs'].append(out_dir)
            if len(inputs) == 1:
                out_1 = '%s/%s.fastq' % (out_dir, self.sam)
            else:
                out_1 = '%s/%s_R1.fastq' % (out_dir, self.sam)
            out_2 = '%s/%s_R2.fastq' % (out_dir, self.sam)
            if inputs[0].endswith('.gz'):
                out_1 += '.gz'
                out_2 += '.gz'

            cmd = '\nbowtie2'
            cmd += ' -p %s' % self.soft.params['cpus']
            cmd += ' -x %s' % db_path
            cmd += ' --very-sensitive'
            if len(inputs) == 1:
                cmd += ' -U %s' % inputs[0]
            else:
                cmd += ' -1 %s' % inputs[0]
                cmd += ' -2 %s' % inputs[1]
            cmd += ' | samtools view -f 12 -F 256'
            cmd += ' | samtools sort -@ %s -n' % self.soft.params['cpus']
            cmd += ' | samtools view -bS'
            cmd += ' | bedtools bamtofastq -i -'
            cmd += ' -fq %s' % out_1.replace('fastq.gz', 'fastq')
            outputs.append(out_1.replace('fastq.gz', 'fastq'))
            if len(inputs) > 1:
                cmd += ' -fq2 %s' % out_2.replace('fastq.gz', 'fastq')
                outputs.append(out_2.replace('fastq.gz', 'fastq'))
            if ddx:
                cmd += '\nrm %s\n' % ' '.join(inputs)
            if self.config.force or todo(out_1):
                cmds += cmd
        if outputs[0].endswith('.fastq'):
            for output in outputs:
                cmds += '\ngzip %s' % output
        if cmds:
            self.outputs['cmds'].append(cmds)
        outputs_gz = ['%s.gz' % x for x in outputs]
        self.outputs['outs'].setdefault(tech, []).extend(outputs_gz)
        io_update(self, o_f=outputs_gz)
