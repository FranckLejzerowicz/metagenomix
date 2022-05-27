# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
from os.path import basename, isfile, splitext

from metagenomix._io_utils import check_min_lines_count

scripts = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


def append_cmd(config, cmd: list, tab: str, cmds: list):
    """Add or not the current SHOGUN command to the list of commands to run
    depending on whether the file already exists or exists but is only a header.

    Parameters
    ----------
    config
        Configuration.
    cmd : list
        List of current SHOGUN command lines.
    tab : str
        Path to the output table of the current command.
    cmds : list
        List of command to possibly expand with the current command.
    """
    if config.force:
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


def get_functional_cmd(
        fun_tab: str, fun: str, out_dir: str, cmds: list, config) -> None:
    """

    Parameters
    ----------
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
    append_cmd(config, list([cmd]), out_file, cmds)


def get_coverage_cmd(ali: str, tax: str, level: str,
                     cov_tab: str, cmds: list, config) -> None:
    """Get the SHOGUN command that calculates the coverage per taxon.

    Parameters
    ----------
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
    append_cmd(config, list([cmd]), cov_tab, cmds)


def get_redistribution_command(
        tax: str, tax_norm: str, cmds: list, config) -> list:
    """Get the SHOGUN command to redistribute a taxonomic profile
    at various taxonomic levels.

    Parameters
    ----------
    args : dict
        Arguments incl. 'softs', 'soft', 'config', 'inputs' and 'databases'.
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
    append_cmd(config, redist_cmds, redist_out, cmds)
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


def get_assign_cmd(aligner: str, ali: str, tax: str,
                   out_tab: str, cmds: list, sub_db: str, config) -> None:
    """Get the taxonomic assignment and the taxonomic table
    normalization SHOGUN commands.

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
    append_cmd(config, list([cmd]), out_norm, cmds)


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
        if params.get('b2') == 'paired':
            cmd = 'bowtie2'
            cmd += ' -p %s' % params['cpus']
            cmd += ' -x %s' % tax
            cmd += ' -S %s' % out
            cmd += ' -1 02_fastp/SRR17458627_R1.fastq.gz'
            cmd += ' -2 02_fastp/SRR17458627_R2.fastq.gz'
            cmd += ' --seed 12345'
            cmd += ' --very-sensitive'
            cmd += ' -k 16'
            cmd += ' --np 1'
            cmd += ' --mp "1,1"'
            cmd += ' --rdg "0,1"'
            cmd += ' --rfg "0,1"'
            cmd += ' --score-min "L,0,-0.05"'
            cmd += ' --no-head'
            cmd += ' --no-unal'
    cmd = 'shogun align -a %s' % aligner
    cmd += ' -i %s' % fasta
    cmd += ' -d %s' % tax
    cmd += ' -t %s' % params['cpus']
    cmd += ' -o %s' % out
    return cmd, ali_out


def get_full_ali_cmd(
        params: dict, aligner: str, fasta: str,
        tax: str, out: str) -> list:
    """Get the full command lines to perform sequence alignment using
    SHOGUN aligner and potentially, re-run if the output was existing but
    possibly aborted.

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
    out : str
        Path to the output folder.

    Returns
    -------
    ali_cmds : list
        Alignment commands.
    """
    cmd, ali_out = get_ali_cmd(params, aligner, fasta, tax, out)
    if not isfile(ali_out):
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


def get_orients(inputs: list):
    if len(inputs) == 1:
        orients = ['']
    elif len(inputs) == 2:
        orients = ['1', '2']
    elif len(inputs) == 3:
        orients = ['', '1', '2']
    else:
        raise IOError('Input to shogun must be 1, 2 or 3 fasta files')
    return orients


def get_combine_cmd(sam: str, inputs: dict, out_dir: str, io: dict) -> tuple:
    """

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
    fastas : list
        Paths to fasta files to combine for SHOGUN.
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files.
    """
    fastas, combine_cmds = [], []
    # get the fastq versions of input file
    orients = get_orients(inputs[sam])
    for pdx, path_ in enumerate(sorted(inputs[sam])):
        io['I']['f'].append(path_)
        path = path_
        if path_.endswith('fastq') or path_.endswith('fastq.gz'):
            # replace non-fasta by fasta extensions
            path = path_.replace('.fastq', '.fasta').replace('.gz', '')
            # Prepare the extraction command and collect it
            if isfile(path_) and not isfile(path):
                to_fasta_cmd = 'seqtk seq -A %s > %s' % (path_, path)
                combine_cmds.append(to_fasta_cmd)

        path_out = '%s/%s' % (out_dir, basename(path))
        edit_fasta = 'python3 %s/fasta4shogun.py -i %s -o %s -s %s' % (
            scripts, path, path_out, sam)
        orient = orients[pdx]
        if orient:
            edit_fasta += ' -r %s' % orient
        combine_cmds.append(edit_fasta)
        fastas.append(path_out)
    return fastas, combine_cmds


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
    """
    fastas, combine_cmds = get_combine_cmd(sam, inputs, out_dir, io)
    fasta = '%s/combined.fasta' % out_dir
    combine_cmds.append('cat %s > %s' % (' '.join(fastas), fasta))
    combine_cmds.append('rm %s' % ' '.join(fastas))
    return combine_cmds, fasta


def add_to_dirs_io(out, dirs, io):
    dirs.append(out)
    io['O']['d'].append(out)


def shogun(
        out_dir: str,
        sam: str,
        inputs: dict,
        params: dict,
        shogun: dict,
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
    prev : str
        Previous software.
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
    combine_cmds, fasta = combine_inputs(sam, inputs, out_dir, io)
    io['O']['f'].append(fasta)

    ali_cmds = []
    # for aligner in ['bowtie2', 'burst']:
    for aligner in ['bowtie2']:
        ali_base = get_alignment_basename(aligner)
        for db, tax in shogun.items():
            key = (aligner, 'tax', db)
            if aligner == 'burst' and db != 'wol':
                continue
            out = out_dir + '/' + aligner + '/tax/' + db #+ '/%s' % params['b2']
            add_to_dirs_io(out, dirs, io)

            ali, tax_tab, tax_norm = get_out_paths(out, ali_base, 'tax')
            outputs.setdefault(key, []).extend([ali, tax_tab, tax_norm])
            io['O']['f'].extend([ali, tax_tab, tax_norm])

            ali_cmds = get_full_ali_cmd(params, aligner, fasta, tax, out)
            get_assign_cmd(aligner, ali, tax, tax_tab, cmds, '', config)
            redists = get_redistribution_command(
                tax, tax_norm, cmds, config)
            io['O']['f'].extend(redists)

            if aligner == 'burst':
                cov_tab = '%s/%s_coverage.tsv' % (out, sam)
                outputs[key].append(cov_tab)
                io['O']['f'].append(cov_tab)
                get_coverage_cmd(ali, tax, 'strain', cov_tab, cmds, config)

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
                    get_assign_cmd(aligner, ali, fun, sub_fun_tab,
                                   cmds, '/functions-%s' % sub_db, config)

            key = (aligner, 'fun_algo', db)
            out_db = out_dir + '/' + aligner + '/fun_algo/' + db
            fun_tab = '%s/%s_taxtable_norm.tsv' % (
                out_db.replace('/fun_algo/', '/tax/'), sam)
            out_dir = '%s/%s_fun_table_norm_results' % (out_db, sam)
            get_functional_cmd(fun_tab, fun, out_dir, cmds, config)
            outputs.setdefault(key, []).append(out_dir)
            io['O']['d'].append(out_dir)
            io['O']['f'].append(fun_tab)

    if cmds:
        cmds = list(combine_cmds) + list(ali_cmds) + list(cmds)

    return io, cmds, outputs, dirs
