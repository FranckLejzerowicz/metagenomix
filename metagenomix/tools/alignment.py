# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from os.path import isfile, splitext

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


def get_bowtie2_cmd(params: dict, fastx: list, db_path: str,
                    db_out: str) -> tuple:
    """Get the bowtie2 alignment command.

    Parameters
    ----------
    params : dict
        Run parameters
    fastx : list
        Path to the input fasta/fastq/fastq.gz files
    db : str
        Database name
    db_path : str
        Path to the bowtie2 database
    db_out: str
        Path to the output folder

    Returns
    -------
    cmd : str
        Alignment command
    """
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


def bowtie(out_dir, sample, inputs, params) -> dict:
    """Get the full command lines to perform sequence alignment using
     bowtie2 and potentially, re-run if the output was existing but
     possibly aborted.

    Parameters
    ----------
    out_dir : str
        Path to output folder
    sample : str
        Sample name
    inputs : dict
        Input files
    params : dict
        Software parameters

    Returns
    -------
    outputs : dict
        All outputs
    """
    outputs = {'io': {'I': {'d': set(), 'f': set()}, 'O': {'d': set()}},
               'cmds': [], 'dirs': [], 'outs': dict({})}
    out = out_dir + '/' + sample
    fastx = inputs[sample]
    outputs['io']['I']['f'].update(fastx)
    for db, db_path in params['databases'].items():
        db_out = '%s/%s/%s' % (out, db, params['pairing'])
        outputs['io']['O']['d'].add(db_out)
        outputs['dirs'].append(db_out)
        cmd, sam = get_bowtie2_cmd(params, fastx, db_path, db_out)
        outputs['outs'][(db, params['pairing'])] = sam
        outputs['cmds'].append(get_alignment_cmd(fastx, cmd, sam))
    return outputs