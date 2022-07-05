# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import gzip
from os.path import basename, isfile
from metagenomix._io_utils import io_update


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


def edit_fastq_cmd(fastq_fp: str, num: int, source: str) -> str:
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
    cmd += "awk '{ if (NR%s4==1) " % '%'
    if source == 'illumina':
        cmd += "{ print $1\"/%s\" } " % str(num)
    elif source == 'ebi':
        cmd += "{ gsub(/ .*/,\"/%s\",$0); print } " % num
    cmd += "else if (NR%s2 == 1) { print \"+\" } " % '%'
    cmd += "else { print } }'"
    if fastq_fp.endswith('.gz'):
        cmd += " | gzip"
    cmd += " > %s_renamed\n" % fastq_fp
    cmd += "mv %s_renamed %s" % (fastq_fp, fastq_fp)
    if fastq_fp.endswith('.fastq'):
        cmd += "\n"
    else:
        cmd += ".gz\n"
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


def get_edit_cmd(self, num: int) -> None:
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
    """
    fastqs_fps = self.inputs[self.sam]
    fastq_fp = fastqs_fps[num - 1]
    line = get_fastq_header(self.config.fastq[self.sam][num - 1])
    line_split = line.split()
    if not line_split[0].endswith('/%s' % num):
        if len(line_split) > 1:
            if line_split[1].startswith('%s:N:' % num):
                cmd = edit_fastq_cmd(fastq_fp, num, 'illumina')
            elif line_split[0].endswith('.%s' % num):
                cmd = edit_fastq_cmd(fastq_fp, num, 'ebi')
            else:
                sys.exit('Fastq file headers not editable:\n%s' % line)
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
    for r in [1, 2]:
        if len(self.config.fastq[self.sam]) < r:
            continue
        get_edit_cmd(self, r)


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
        cmd = "n%s=`%s %s | wc -l | " % (self.soft.params['cat'], idx, input_fp)
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
    out = '%s/%s_read_count.tsv' % (self.dir, self.sam)
    self.outputs['outs'].append(out)
    if self.config.force or not isfile(out):
        inputs = self.inputs[self.sam]
        for idx, input_path in enumerate(inputs):
            cmd = count_cmd(self, idx, input_path, out)
            self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, o_f=out)


def fastqc(self) -> None:
    out_dir = '%s/%s' % (self.dir, self.sam)
    ins = self.inputs[self.sam]
    outs = ['%s_fastqc.html' % x.rsplit('.fastq', 1)[0] for x in ins]
    self.outputs['outs'].append(out_dir)
    if self.config.force or [not isfile(x) for x in outs]:
        cmd = 'fastqc %s -o %s' % (' '.join(ins), out_dir)
        self.outputs['dirs'].append(out_dir)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=ins, o_d=out_dir)


def fastp(self) -> None:
    fastqs = self.inputs[self.sam]
    out_dir = '%s/%s' % (self.dir, self.sam)
    out_json = '%s/%s.json' % (out_dir, self.sam)
    out_html = '%s/%s.html' % (out_dir, self.sam)
    outs = []
    self.outputs['outs'].append(out_dir)
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
        outs.append(isfile(out))
        cmd += '  -%s %s -%s %s' % (i, out, o, out)
    cmd += ' --verbose'
    if self.soft.params['trim_poly_g']:
        cmd += ' --trim_poly_g'
    cmd += ' --json=%s' % out_json
    cmd += ' --html=%s' % out_html
    cmd += ' -w %s' % self.soft.params['cpus']
    cmd += ' --report_title="fastp report %s"' % self.sam
    cmd += ' --average_qual=%s' % self.soft.params['average_qual']
    if self.config.force or not sum(outs):
        self.outputs['dirs'].append(out_dir)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=fastqs, o_d=out_dir)


def cutadapt(self) -> None:
    r1_o = '%s/%s.R1.fastq.gz' % (self.dir, self.sam)
    r2_o = '%s/%s.R2.fastq.gz' % (self.dir, self.sam)
    cmd = 'cutadapt'
    cmd += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    cmd += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    cmd += ' -m 10'
    cmd += ' -o %s' % r1_o
    cmd += ' -p %s ' % r2_o
    cmd += ' '.join(self.inputs[self.sam])
    outs = [r1_o, r2_o]
    self.outputs['outs'].extend(outs)
    if self.config.force or not isfile(r1_o) or not isfile(r2_o):
        self.outputs['cmds'] = list([cmd])
        io_update(self, i_f=self.inputs[self.sam], o_f=outs)


def atropos(self):
    inputs = self.inputs[self.sam]
    if len(inputs) == 1:
        out_1 = '%s/%s/%s.fastq' % (self.dir, self.sam, self.sam)
    else:
        out_1 = '%s/%s/%s_R1.fastq' % (self.dir, self.sam, self.sam)
    out_2 = '%s/%s/%s_R2.fastq' % (self.dir, self.sam, self.sam)
    log = '%s/%s/%s.log' % (self.dir, self.sam, self.sam)
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
    if self.config.force or not isfile(out_1) or not isfile(out_2):
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, o_f=[out_1, out_2])
    self.outputs['outs'].extend([out_1, out_2])


def kneaddata(self):
    inputs = self.inputs[self.sam]
    out_dir = '%s/%s' % (self.dir, self.sam)
    log = '%s/%s.log' % (out_dir, self.sam)
    cmd = 'kneaddata'
    for inp in inputs:
        cmd += ' --input %s' % inp
    for database in self.soft.params['databases']:
        cmd += ' --reference-db %s' % database
    cmd += ' --output %s' % out_dir
    # cmd += ' --output-prefix %s' % self.sam
    cpus = self.soft.params['cpus']
    cmd += ' --processes %s' % cpus
    cmd += ' --max-memory %s%s' % (self.soft.params['mem_num'],
                                   self.soft.params['mem_dim'][0])
    if self.soft.params['trimmomatic']:
        cmd += ' --trimmomatic %s' % self.soft.params['trimmomatic']
    if self.soft.params['bowtie2']:
        cmd += ' --bowtie2 %s' % self.soft.params['bowtie2']
    cmd += ' --bowtie2-options="--very-sensitive-local -p %s"' % cpus
    cmd += ' --remove-intermediate-output'
    cmd += ' --log %s' % log
    self.outputs['cmds'].append(cmd)
    # write these additional output files manipulations commands
    rad = basename(inputs[0]).replace('.fastq.gz', '')
    for typ in ['paired', 'unmatched']:
        for r in ['1', '2']:
            out = '%s/%s_kneaddata_%s_%s.fastq' % (out_dir, rad, typ, r)
            renamed = out.replace('.R1', '')
            self.outputs['outs'].append(renamed)
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
    outputs = self.inputs[self.sam]
    io_update(self, i_f=self.inputs[self.sam])
    cmds = ''
    databases = self.soft.params['databases']
    for ddx, db_path in enumerate(databases):
        inputs = list(outputs)
        outputs = []
        db = db_path.split('/')[-1]
        if len(inputs) == 1:
            out_1 = '%s/%s_%s/%s.fastq' % (self.dir, ddx, db, self.sam)
        else:
            out_1 = '%s/%s_%s/%s_R1.fastq' % (self.dir, ddx, db, self.sam)
        out_2 = '%s/%s_%s/%s_R2.fastq' % (self.dir, ddx, db, self.sam)
        if inputs[0].endswith('.gz'):
            out_1 += '%s.gz'
            out_2 += '%s.gz'

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
        if self.config.force or not isfile(out_1):
            cmds += cmd
    if outputs[0].endswith('.fastq'):
        for output in outputs:
            cmds += '\ngzip %s' % output
    if cmds:
        self.outputs['cmds'].append(cmds)
    self.outputs['outs'].extend(outputs)
    io_update(self, o_f=outputs)
