# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import basename, isfile, splitext

from metagenomix._io_utils import check_min_lines_count


def append_cmd(config, cmd: list, ret: int, tab: str, cmds: list):
    """Add or not the current SHOGUN command to the list of commands to run
    depending on whether the file already exists or exists but is only a header.

    Parameters
    ----------
    config
        Configuration.
    cmd : list
        List of current SHOGUN command lines.
    ret : int
        Return value of `check_combined` function telling if input is valid.
    tab : str
        Path to the output table of the current command.
    cmds : list
        List of command to possibly expand with the current command.
    """
    if ret == 1 or config.force:
        cmds.extend(cmd)
    elif not isfile(tab) or not check_min_lines_count(tab):
        cmds.append('file="%s"' % tab)
        cmds.append('if [ -f "$file" ]')
        cmds.append('then')
        cmds.append("FILESIZE=`ls -l $file | cut -d ' ' -f5`")
        cmds.append("if [[ $FILESIZE -lt 1000 ]]")
        cmds.append('then')
        cmds.extend(cmd)
        cmds.append('fi')


def get_functional_cmd(ret: int, fun_tab: str, fun: str,
                       out_dir: str, cmds: list, config) -> None:
    """

    Parameters
    ----------
    ret : int
        Return value of `check_combined` function telling if input is valid.
    fun_tab : str
        Path to the input table file.
    fun : str
        Path to the functional database.
    out_dir : str
        Path to the output functional table.
    cmds : list
        List of command to possibly expand with the current command.
    config
        Configuration.
    """
    cmd = 'shogun functional'
    cmd += ' -i %s' % fun_tab
    cmd += ' -d %s' % fun
    cmd += ' -o %s' % out_dir
    cmd += ' -l genus'
    out_file = '%s/%s.genus.normalized.txt' % (
        out_dir, splitext(basename(fun_tab))[0].replace('_taxt', '_funt'))
    append_cmd(config, list([cmd]), ret, out_file, cmds)


def get_coverage_cmd(ret: int, ali: str, tax: str,
                     level: str, cov_tab: str, cmds: list, config) -> None:
    """Get the SHOGUN command that calculates the coverage per taxon.

    Parameters
    ----------
    ret : int
        Return value of `check_combined` function telling if input is valid.
    ali : str
        Path to the alignment.
    tax : str
        Path to the taxonomic database.
    level : str
        Taxonomic level for coverage calculation.
    cov_tab : str
        Path to the output coverage table file.
    cmds : list
        List of command to possibly expand with the current command.
    config
        Configuration.
    """
    cmd = 'shogun coverage'
    cmd += ' -i %s' % ali
    cmd += ' -d %s' % tax
    cmd += ' -l %s' % level
    cmd += ' -o %s' % cov_tab
    append_cmd(config, list([cmd]), ret, cov_tab, cmds)


def get_redistribution_command(
        ret: int, tax: str, tax_norm: str, cmds: list, config) -> list:
    """Get the SHOGUN command to redistribute a taxonomic profile
    at various taxonomic levels.

    Parameters
    ----------
    args : dict
        Arguments incl. 'softs', 'soft', 'config', 'inputs' and 'databases'.
    ret : int
        Return value of `check_combined` function telling if input is valid.
    tax : str
        Path to the taxonomic database.
    tax_norm : str
        Path to the input taxonomic table file.
    cmds : list
        List of command to possibly expand with the current command.
    config
        Configuration.

    Returns
    -------
    redists : list
        Paths to the redistributed taxonomic table files.
    """
    redists, redist_cmds = [], []
    for level in ['phylum', 'genus', 'strain']:
        out_redist = '%s.redist.%s.tsv' % (splitext(tax_norm)[0], level)
        redists.append(out_redist)
        cmd = 'shogun redistribute'
        cmd += ' -d %s' % tax
        cmd += ' -l %s' % level
        cmd += ' -i %s' % tax_norm
        cmd += ' -o %s' % out_redist
        redist_cmds.append(cmd)
    redist_out = '%s.redist.strain.tsv' % splitext(tax_norm)[0]
    append_cmd(config, redist_cmds, ret, redist_out, cmds)
    return redists


def assign_taxonomy_cmd(aligner: str, ali: str, tax: str,
                        out_tab: str, sub_db: str) -> str:
    """Get the assign taxonomy command.

    Parameters
    ----------
    aligner : str
        Name of the aligner.
    ali : str
        Path to the alignment.
    tax : str
        Path to the taxonomic database.
    out_tab : str
        Path to the output table.
    sub_db : str
        Suffix to the database.

    Returns
    -------
    cmd : str
        assign taxonomy command.
    """
    cmd = 'shogun assign_taxonomy'
    cmd += ' -a %s' % aligner
    cmd += ' -i %s' % ali
    cmd += ' -d %s%s' % (tax, sub_db)
    cmd += ' -o %s\n' % out_tab
    return cmd


def get_normalize_cmd(input_fp: str, output_fp: str) -> str:
    """Get the SHOGUN command to normalize table.

    Parameters
    ----------
    input_fp : str
        Path to the input table to normalize.
    output_fp : str
        Path to the output normalized table.

    Returns
    -------
    cmd : str
        Command line for table normalization.
    """
    cmd = 'shogun normalize -i %s -o %s\n' % (input_fp, output_fp)
    return cmd


def get_assign_cmd(ret: int, aligner: str, ali: str, tax: str,
                   out_tab: str, cmds: list, sub_db: str, config) -> None:
    """Get the taxonomic assignment and the taxonomic table
    normalization SHOGUN commands.

    Parameters
    ----------
    ret : int
        Return value of `check_combined` function telling if input is valid.
    aligner : str
        Name of the aligner.
    ali : str
        Path to the alignment.
    tax : str
        Path to the taxonomic database.
    out_tab : str
        Path to the output table.
    cmds : list
        List of command to possibly expand with the current command.
    sub_db : str
        Suffix to the database.
    config
        Configuration.
    """
    out_norm = '%s_norm.tsv' % splitext(out_tab)[0]
    cmd = assign_taxonomy_cmd(aligner, ali, tax, out_tab, sub_db)
    cmd += get_normalize_cmd(out_tab, out_norm)
    append_cmd(config, list([cmd]), ret, out_norm, cmds)


def condition_ali_command(fasta: str, ali_out: str, cmd: str) -> list:
    """Add the conditional statements checking that the alignment
    already available was not aborted, in which case the alignment
    will be re-run.

    Parameters
    ----------
    fasta : str
        Path to the input fasta file.
    ali_out
        Path to the aligment file.
    cmd : str
        Alingment command using SHOGUN.

    Returns
    -------
    cmds : list
        Full list of commands.
    """
    tail = [x[1:] for x in os.popen('tail -n 40 %s' % fasta).read().split('\n')
            if '>' in x]
    cmds = list()
    cmds.append('\nlistVar="%s"' % ' '.join(tail))
    cmds.append('GREPOK=0')
    cmds.append('for i in $listVar')
    cmds.append('do')
    cmds.append('    if grep -q "$i" %s' % ali_out)
    cmds.append('        then')
    cmds.append('            GREPOK=1')
    cmds.append('            break')
    cmds.append('    fi')
    cmds.append('\ndone')
    cmds.append("if [ $GREPOK = 0 ]")
    cmds.append('then')
    cmds.append(cmd)
    cmds.append('fi')
    return cmds


def get_ali_cmd(params: dict, aligner: str, fasta: str, tax: str,
                out: str) -> tuple:
    """Get the SHOGUN alignment command.

    Parameters
    ----------
    params : dict
        Run parameters.
    aligner : str
        Name of the aligner.
    fasta : str
        Path to the input fasta file.
    tax : str
        Path to the taxonomic database.
    db : str
        Name of the taxonomic database.
    out : str
        Path to the output folder.

    Returns
    -------
    cmd : str
        Alignment command using SHOGUN aligner.
    ali_out : str
        Path to the output alignment file.
    """
    if aligner == 'burst':
        ali_out = '%s/alignment.burst.b6' % out
        # edx_acx = args['databases'].shogun
        # cmd = 'burst'
        # cmd += ' -q %s' % fasta
        # cmd += ' -r %s' % edx_acx[db]['edx']
        # cmd += ' -a %s' % edx_acx[db]['acx']
        # cmd += ' -sa'
        # cmd += ' -fr'
        # cmd += ' -i 0.98'
        # cmd += ' -t %s' % args['soft'].params['cpus']
        # cmd += ' -o %s' % ali_out
    else:
        ali_out = '%s/alignment.bowtie2.sam' % out

    cmd = 'shogun align -a %s' % aligner
    cmd += ' -i %s' % fasta
    cmd += ' -d %s' % tax
    cmd += ' -t %s' % params['cpus']
    cmd += ' -o %s' % out
    return cmd, ali_out


def get_full_ali_cmd(
        params: dict, ret: int, aligner: str, fasta: str,
        tax: str, out: str) -> list:
    """Get the full command lines to perform sequence alignment using
    SHOGUN aligner and potentially, re-run if the output was existing but
    possibly aborted.

    Parameters
    ----------
    params : dict
        Run parameters.
    ret : int
        Return value of `check_combined` function telling if input is valid.
    aligner : str
        Name of the aligner.
    fasta : str
        Path to the input fasta file.
    tax : str
        Path to the taxonomic database.
    out : str
        Path to the output folder.

    Returns
    -------
    ali_cmds : list
        Alignment commands.
    """
    cmd, ali_out = get_ali_cmd(params, aligner, fasta, tax, out)
    if ret or not isfile(ali_out):
        ali_cmds = [cmd]
    else:
        ali_cmds = condition_ali_command(fasta, ali_out, cmd)
    return ali_cmds


def get_out_paths(out: str, ali_base: str, tax_or_fun: str) -> tuple:
    """Get output path for alignment, table and normalized table.

    Parameters
    ----------
    out : str
        Path to the output folder.
    ali_base : str
        Basename of the alignment as outputted by burst or bowtie2.
    tax_or_fun : str
        "tax" of "fun"

    Returns
    -------
    ali : str
        Path to alignment file.
    tab : str
        Path to raw table file.
    norm : str
        Path to normalized table file.
    """
    ali = '%s/%s' % (out, ali_base)
    tab = '%s/%stable.tsv' % (out, tax_or_fun)
    norm = '%s/%stable_norm.tsv' % (out, tax_or_fun)
    return ali, tab, norm


def get_alignment_basename(aligner):
    if aligner == 'bowtie2':
        ali_base = 'alignment.bowtie2.sam'
    elif aligner == 'burst':
        ali_base = 'alignment.burst.b6'
    return ali_base


def get_combine_cmd(sam: str, inputs: dict, io: dict) -> tuple:
    """

    Parameters
    ----------
    sam : str
        Sample name.
    inputs : dict
        Input files.
    io : dict
        Inputs and outputs to potentially move to scratch and back.

    Returns
    -------
    fastas : list
        Paths to fasta files to combine for SHOGUN.
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files.
    """
    fastas, combine_cmds = [], []
    # get the fastq versions of input file
    for path in inputs[sam]:
        if path.endswith('fastq') or path.endswith('fastq.gz'):
            # replace non-fasta by fasta extensions
            fas_path = path.replace('.fastq', '.fasta').replace('.gz', '')
            # Prepare the extraction command and collect it
            if not isfile(fas_path):
                to_fasta_cmd = 'seqtk seq -A %s > %s' % (path, fas_path)
                combine_cmds.append(to_fasta_cmd)
            fastas.append(fas_path)
        elif path.endswith('fasta'):
            fastas.append(path)
        io['I']['f'].append(path)
    return fastas, combine_cmds


def check_combined(fasta: str, fastas: list) -> int:
    """Check whether the last sequences of the combined fasta file
    contains the same sequences as the last fasta file to be combined in.

    Parameters
    ----------
    fasta : str
        Path to the combined sequences fasta file.
    fastas : list
        Paths to fasta files to combine for SHOGUN.

    Returns
    -------
    ret : int
        0: combining ended up properly, with the same last sequences.
        1: combining did not happen.
        2: combining did not end up properly.
    """
    ret = 1
    if isfile(fasta):
        combined_size = os.path.getsize(fasta)
        not_combined_tail = os.popen('tail -n 4 %s' % list(fastas)[-1]).read()
        combined_tail = os.popen('tail -n 4 %s' % fasta).read()
        # Checks that the combining was done in order
        if not_combined_tail == combined_tail:
            ret = 0
        elif combined_size:
            ret = 2
    return ret


def combine_inputs(sam: str, inputs: dict, out_dir: str, io: dict) -> tuple:
    """Combine the fastq/a of the sample to run shogun.

    Parameters
    ----------
    sam : str
        Sample name.
    inputs : dict
        Input files.
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    io : dict
        Inputs and outputs to potentially move to scratch and back.

    Returns
    -------
    combine_cmds : list
        Command lines to perform the combining.
    comb_path : str
        Path to the combined sequences fasta file.
    ret : int
        0: combining ended up properly, with the same last sequences.
        0: combining did not happen.
        2: combining did not end up properly.
    """
    fastas, combine_cmds = get_combine_cmd(sam, inputs, io)
    fasta = '%s/combined.fasta' % out_dir
    ret = check_combined(fasta, list(fastas))
    if ret == 1:
        combine_cmds.append('cat %s > %s' % (' '.join(fastas), fasta))
    return combine_cmds, fasta, ret


def add_to_dirs_io(out, dirs, io):
    dirs.append(out)
    io['O']['d'].append(out)


def shogun(out_dir: str, sam: str, inputs: dict, params: dict, shogun: dict,
           config):
    """

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for SHOGUN.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    params : dict
        Run parameters.
    shogun : dict
        SHOGUN databases.
    config
        Configuration.

    Returns
    -------
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All SHOGUN command lines.
    outputs : list
        All outputs paths.
    dirs : list
        List of output folders to create.
    """
    io = {'I': {'f': [], 'd': []}, 'O': {'f': [], 'd': []}}
    outputs, dirs, cmds = {}, [], []
    out_dir = out_dir + '/' + sam
    add_to_dirs_io(out_dir, dirs, io)

    combine_cmds, fasta, ret = combine_inputs(sam, inputs, out_dir, io)
    io['O']['f'].append(fasta)

    ali_cmds = []
    for aligner in ['bowtie2', 'burst']:
        ali_base = get_alignment_basename(aligner)
        for db, tax in shogun.items():
            key = (aligner, 'tax', db)
            if aligner == 'burst' and db != 'wol':
                continue
            out = out_dir + '/' + aligner + '/tax/' + db
            add_to_dirs_io(out, dirs, io)

            ali, tax_tab, tax_norm = get_out_paths(out, ali_base, 'tax')
            outputs.setdefault(key, []).extend([ali, tax_tab, tax_norm])
            io['O']['f'].extend([ali, tax_tab, tax_norm])

            ali_cmds = get_full_ali_cmd(params, ret, aligner, fasta, tax, out)
            get_assign_cmd(ret, aligner, ali, tax, tax_tab, cmds, '', config)
            redists = get_redistribution_command(
                ret, tax, tax_norm, cmds, config)
            io['O']['f'].extend(redists)

            if aligner == 'burst':
                cov_tab = '%s/%s_coverage.tsv' % (out, sam)
                outputs[key].append(cov_tab)
                io['O']['f'].append(cov_tab)
                get_coverage_cmd(ret, ali, tax, 'strain', cov_tab, cmds, config)

        for (fun, db) in shogun.items():
            key = (aligner, 'fun', db)
            out = out_dir + '/' + aligner + '/fun/' + db
            ali, fun_tab, fun_norm = get_out_paths(out, ali_base, 'fun')

            if db == 'shogun':
                if aligner == 'burst':
                    continue
                for sub_db in ['kegg', 'refseq', 'uniprot']:
                    out_db = '%s-%s' % (out, sub_db)
                    sub_fun_tab = fun_tab.replace(out, out_db)
                    sub_fun_norm = fun_norm.replace(out, out_db)
                    outputs.setdefault(key, []).append(sub_fun_tab)
                    outputs[key].append(sub_fun_norm)
                    io['O']['f'].extend([sub_fun_tab, sub_fun_norm])
                    io['O']['d'].append(out_db)
                    ali = ali.replace('/fun/', '/tax/')
                    get_assign_cmd(ret, aligner, ali, fun, sub_fun_tab,
                                   cmds, '/functions-%s' % sub_db, config)

            key = (aligner, 'fun_algo', db)
            out_db = out_dir + '/' + aligner + '/fun_algo/' + db
            fun_tab = '%s/%s_taxtable_norm.tsv' % (
                out_db.replace('/fun_algo/', '/tax/'), sam)
            out_dir = '%s/%s_fun_table_norm_results' % (out_db, sam)
            get_functional_cmd(ret, fun_tab, fun, out_dir, cmds, config)
            outputs.setdefault(key, []).append(out_dir)
            io['O']['d'].append(out_dir)
            io['O']['f'].append(fun_tab)

    if cmds:
        cmds = list(combine_cmds) + list(ali_cmds) + list(cmds)

    return io, cmds, outputs, dirs
