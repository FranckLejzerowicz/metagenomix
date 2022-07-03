# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import pkg_resources
from os.path import isfile, splitext
from metagenomix._io_utils import io_update

# Keep line because read mapping alignments will be python scripts
# scripts = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


def condition_ali_command(fastx: list, cmd: str, sam: str) -> str:
    """Add the conditional statements checking that the alignment
    already available was not aborted, in which case the alignment
    will be re-run.

    Parameters
    ----------
    fastx : list
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
    if fastx[-1].endswith('fastq.gz'):
        reads = 'zcat %s | tail -n 80 | grep "^@"' % fastx[-1]
    elif fastx[-1].endswith('fastq'):
        reads = 'tail -n 80 %s | grep "^@"' % fastx[-1]
    elif fastx[-1].endswith('fasta'):
        reads = 'tail -n 40 %s | grep "^>"' % fastx[-1]
    else:
        return cmd

    cmd = grep_sam_tail_cmd(cmd, reads, sam)
    return cmd


def grep_sam_tail_cmd(cmd, reads, sam) -> str:
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


def get_burst_cmd(fastx: str, db_path: str, params: dict, out: str) -> str:
    """Get the burst alignment command.

    Parameters
    ----------
    fastx : str
        Path to the input fasta file
    db_path : str
        Path to the burst database
    params : dict
        Run parameters.
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Alignment command
    """
    cmd = 'burst'
    cmd += ' -q %s' % fastx
    cmd += ' -r %s.edx' % db_path
    cmd += ' -a %s.acx' % db_path
    cmd += ' -sa'
    cmd += ' -fr'
    cmd += ' -i 0.98'
    cmd += ' -t %s' % params['cpus']
    cmd += ' -o %s' % out
    return cmd


def get_bowtie2_cmd(self, fastx: list, db_path: str, db_out: str) -> tuple:
    """Get the bowtie2 alignment command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters for humann
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
    cmd = 'bowtie2'
    cmd += ' -p %s' % self.soft.params['cpus']
    cmd += ' -x %s' % db_path
    sam = '%s/alignment.bowtie2' % db_out
    if len(fastx) == 2 and self.soft.params['pairing'] == 'paired':
        cmd += ' -1 %s' % fastx[0]
        cmd += ' -2 %s' % fastx[1]
        if not self.soft.params['discordant']:
            cmd += ' --no-discordant'
            sam += '.nodiscordant'
    else:
        cmd += ' -U %s' % ','.join(fastx)
    sam += '.sam'
    cmd += ' -S %s' % sam
    cmd += ' --seed 12345'
    cmd += ' --very-sensitive'
    cmd += ' -k %s' % self.soft.params['k']
    cmd += ' --np %s' % self.soft.params['np']
    cmd += ' --mp "%s"' % self.soft.params['mp']
    cmd += ' --rdg "%s"' % self.soft.params['rdg']
    cmd += ' --rfg "%s"' % self.soft.params['rfg']
    cmd += ' --score-min "%s"' % self.soft.params['score-min']
    cmd += ' --met-file %s.met' % splitext(sam)[0]
    cmd += ' --met 240'
    cmd += ' --no-unal'
    return cmd, sam


def get_alignment_cmd(fastx: list, cmd: str, ali: str) -> str:
    """Get the alignment command with ot without unix checks to decide
    whether the alignment (if exists) was aligned to completion, and not
    aborted (e.g., for lack of memory reasons).

    Parameters
    ----------
    fastx : list
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
    if not isfile(ali):
        return cmd
    else:
        return condition_ali_command(fastx, cmd, ali)


def bowtie2(self) -> None:
    """Get the full command lines to perform sequence alignment using
     bowtie2 and potentially, re-run if the output was existing but
     possibly aborted.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
    """
    out = self.dir + '/' + self.sam
    self.outputs['outs'] = dict()
    fastx = self.inputs[self.sam]
    for db, db_path in self.soft.params['databases'].items():
        db_out = '%s/%s/%s' % (out, db, self.soft.params['pairing'])
        cmd, sam = get_bowtie2_cmd(self, fastx, db_path, db_out)
        cmd = get_alignment_cmd(fastx, cmd, sam)
        key = (db, self.soft.params['pairing'])
        self.outputs['outs'][key] = sam
        self.outputs['cmds'].append(cmd)
        self.outputs['dirs'].append(db_out)
        io_update(self, i_f=fastx, o_d=db_out)


def check_bowtie_k_np(soft, params, defaults, let_go):
    for opt in ['k', 'np']:
        if opt in params and not str(params[opt]).isdigit():
            sys.exit('[bowtie2] "%s" option invalid' % opt)
        else:
            soft.params[opt] = defaults[opt]
        let_go.append(opt)


def check_bowtie_mp_rdg_rfg(soft, params, defaults, let_go):
    for opt in ['mp', 'rdg', 'rfg']:
        if opt in params:
            if len([x.isdigit() for x in str(params[opt]).split(',')]) != 2:
                sys.exit('[bowtie2] "%s" option invalid' % opt)
        else:
            soft.params[opt] = defaults[opt]
        let_go.append(opt)


def check_bowtie_score_min(soft, params, defaults, let_go):
    if 'score-min' in params:
        s = params['score-min'].split(',')
        if len(s) != 3 or s[0] not in 'CLSG':
            sys.exit('[bowtie2] "score-min" option invalid')
        else:
            for r in [1, 2]:
                try:
                    float(s[r])
                except ValueError:
                    sys.exit('[bowtie2] "score-min" option invalid')
    else:
        soft.params['score-min'] = defaults['score-min']
    let_go.append('score-min')


def flash(self):
    """Create command lines for FLASh

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
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
    min_overlap = self.soft.params['min_overlap']
    max_overlap = self.soft.params['max_overlap']
    mismatch = self.soft.params['mismatch']
    out = '%s/m%s_M%s_e%s/%s' % (
        self.dir, min_overlap, max_overlap, mismatch, self.sam)
    rad = out + '/' + self.sam
    ext_fq = '%s.extendedFrags.fastq' % rad
    ext = ext_fq.replace('.fastq', '.fasta')
    nc1_fq = '%s.notCombined_1.fastq' % rad
    nc1 = nc1_fq.replace('.fastq', '.fasta')
    nc2_fq = '%s.notCombined_2.fastq' % rad
    nc2 = nc2_fq.replace('.fastq', '.fasta')
    out_fps = [ext, nc1, nc2]
    self.outputs['dirs'].append(out)
    self.outputs['outs'].extend(out_fps)
    if self.config.force or len([not isfile(x) for x in out_fps]):
        cmd = 'flash %s' % ' '.join(self.inputs[self.sam])
        cmd += ' -m %s' % min_overlap
        cmd += ' -M %s' % max_overlap
        cmd += ' -x %s' % mismatch
        cmd += ' -d %s' % out
        cmd += ' -o %s' % self.sam
        cmd += ' -t %s\n' % self.soft.params['cpus']
        cmd += 'seqtk seq -A %s > %s\n' % (ext_fq, ext)
        cmd += 'seqtk seq -A %s > %s\n' % (nc1_fq, nc1)
        cmd += 'seqtk seq -A %s > %s\n' % (nc2_fq, nc2)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=self.inputs[self.sam], o_f=out_fps)
