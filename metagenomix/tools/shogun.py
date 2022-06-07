# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from os.path import basename, isfile, splitext

from metagenomix._io_utils import check_min_lines_count
from metagenomix.tools.alignment import get_alignment_cmd

scripts = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


def append_cmd(cmds: list, tab: str, outputs: dict, config):
    """Add or not the current SHOGUN command to the list of commands to run
    depending on whether the file already exists or exists but is only a header.

    Parameters
    ----------
    cmds : list
        List of current SHOGUN command lines
    tab : str
        Path to the output table of the current command
    outputs : dict
        All outputs
    config
        Configuration
    """
    if config.force:
        outputs['cmds'].extend(cmds)
    elif not isfile(tab) or not check_min_lines_count(tab):
        cmd = 'file="%s"\n' % tab
        cmd += 'if [ -f "$file" ]\n'
        cmd += 'then\n'
        cmd += "FILESIZE=`ls -l $file | cut -d ' ' -f5`\n"
        cmd += "if [[ $FILESIZE -lt 1000 ]]\n"
        cmd += 'then\n'
        cmd += '\n'.join(cmds)
        cmd += 'fi\n'
        outputs['cmds'].append(cmd)


def redistribute(db: str, aligner: str, tax_norm: str, db_path: str,
                 outputs: dict, config) -> None:
    """Get the SHOGUN command to redistribute a taxonomic profile
    at various taxonomic levels.

    Parameters
    ----------
    db : str
        Database name
    aligner : str
        Name of the aligner
    tax_norm : str
        Path to the input taxonomic table file
    db_path : str
        Path to the taxonomic database
    outputs : dict
        All outputs
    config
        Configuration
    """
    redist_cmds = []
    for level in ['phylum', 'genus', 'strain']:
        redist = '%s.redist.%s.tsv' % (splitext(tax_norm)[0], level)
        outputs['io']['O']['f'].add(redist)
        outputs['outs'].setdefault((db, aligner), []).append(redist)
        cmd = 'shogun redistribute'
        cmd += ' -i %s' % tax_norm
        cmd += ' -d %s' % db_path
        cmd += ' -l %s' % level
        cmd += ' -o %s' % redist
        redist_cmds.append(cmd)
    redist_out = '%s.redist.strain.tsv' % splitext(tax_norm)[0]
    append_cmd(redist_cmds, redist_out, outputs, config)


def assign_taxonomy(aligner: str, ali: str, tax: str, db_path: str,
                    sub_db: str) -> str:
    """Get the assign taxonomy command.

    Parameters
    ----------
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    tax : str
        Path to the output table
    db_path : str
        Path to the database
    sub_db : str
        Suffix to the database

    Returns
    -------
    cmd : str
        assign taxonomy command
    """
    cmd = 'shogun assign_taxonomy'
    cmd += ' -a %s' % aligner
    cmd += ' -i %s' % ali
    cmd += ' -d %s%s' % (db_path, sub_db)
    cmd += ' -o %s\n' % tax
    return cmd


def normalize(input_fp: str, output_fp: str) -> str:
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


def get_ali_cmd(fasta: str, db_path: str, out: str, aligner: str,
                params: dict) -> str:
    """Get the SHOGUN alignment command.

    Parameters
    ----------
    fasta : str
        Path to the input fasta file
    db_path : str
        Path to the shogun database
    out : str
        Path to the output folder
    aligner : str
        Name of the aligner
    params : dict
        Run parameters

    Returns
    -------
    cmd : str
        Alignment command using SHOGUN aligner
    """
    cmd = 'shogun align -a %s' % aligner
    cmd += ' -i %s' % fasta
    cmd += ' -d %s' % db_path
    cmd += ' -t %s' % params['cpus']
    cmd += ' -o %s' % out
    return cmd


def get_alignment_basename(aligner):
    if aligner == 'bowtie2':
        ali_base = 'alignment.bowtie2.sam'
    elif aligner == 'burst':
        ali_base = 'alignment.burst.b6'
    else:
        raise ValueError('No output basename for aligner "%s"' % aligner)
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


def get_combine_cmd(sample: str, inputs: dict, out_dir: str,
                    combine_cmds: list, outputs : dict) -> list:
    """

    Parameters
    ----------
    sample : str
        Sample name
    inputs : dict
        Input files
    out_dir : str
        Path to pipeline output folder for SHOGUN
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files
    outputs : dict
        All outputs

    Returns
    -------
    fastas : list
        Paths to fasta files to combine for SHOGUN
    """
    fastas = []
    # get the fastq versions of input file
    orients = get_orients(inputs[sample])
    for pdx, path_ in enumerate(sorted(inputs[sample])):
        outputs['io']['I']['f'].add(path_)
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
            scripts, path, path_out, sample)
        orient = orients[pdx]
        if orient:
            edit_fasta += ' -r %s' % orient
        combine_cmds.append(edit_fasta)
        fastas.append(path_out)
    return fastas


def combine_inputs(sample: str, inputs: dict, out_dir: str,
                   combine_cmds: list, outputs: dict) -> str:
    """Combine the fastq/a of the sample to run shogun.

    Parameters
    ----------
    sample : str
        Sample name
    inputs : dict
        Input files
    out_dir : str
        Path to pipeline output folder for SHOGUN
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files
    outputs : dict
        All outputs

    Returns
    -------
    fasta : str
        Path to the combined sequences fasta file
    """
    fastas = get_combine_cmd(sample, inputs, out_dir, combine_cmds, outputs)
    fasta = '%s/combined.fasta' % out_dir
    combine_cmds.extend([
        'cat %s > %s' % (' '.join(fastas), fasta), 'rm %s' % ' '.join(fastas)])
    return fasta


def align_cmds(fasta: str, out: str, db: str, db_path: str, aligner: str,
               outputs: dict, ali_cmds: list, params: dict) -> str:
    """

    Parameters
    ----------
    fasta : str
        Path to the combined sequences fasta file
    out : str
        Path to pipeline output folder for SHOGUN
    db : str
        Database name
    db_path : str
    aligner : str
    outputs : dict
    ali_cmds : list
    params : dict

    Returns
    -------
    ali : str
        Alignment output file
    """
    out_dir = '%s/%s/%s' % (out, db, aligner)
    ali = '%s/%s' % (out_dir, get_alignment_basename(aligner))

    outputs['dirs'].append(out_dir)
    outputs['outs'].setdefault((db, aligner), []).append(ali)
    outputs['io']['O']['d'].add(out_dir)

    cmd = get_ali_cmd(fasta, db_path, out_dir, aligner, params)
    cmd = get_alignment_cmd([fasta], cmd, ali)
    ali_cmds.append(cmd)
    if cmd.startswith('shogun'):
        outputs['io']['I']['f'].add(ali)
    return ali


def format_sam(sam_: str, ali_cmds: list, sample: str) -> str:
    sam = '%s_formatted.sam' % splitext(sam_)[0]
    cmd = '%s/sam4shogun.py -i %s -o %s -s %s' % (scripts, sam_, sam, sample)
    ali_cmds.append(cmd)
    return sam


def get_dir(out_dir: str, db: str, aligner: str) -> str:
    """

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for SHOGUN
    db : str
        Database name
    aligner : str
        Name of the aligner

    Returns
    -------
    out : str
        Path to pipeline output folder for SHOGUN
    """
    out = '%s/%s' % (out_dir, db)
    if aligner:
        # i.e., if the previous software is not bowtie2
        out += '/%s' % aligner
    return out


def assign_normalize(tab: str, aligner: str, ali: str, db_path: str,
                     outputs: dict, sub_db: str, config):
    """Get the assignment and the normalization SHOGUN commands for
    the taxonomy or functional classifications.

    Parameters
    ----------
    tab : str
        Path to the output table
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    db_path : str
        Path to the database
    outputs : dict
        All outputs
    sub_db : str
        Suffix to the database
    config
        Configuration
    """
    cmd = ''
    tab_norm = '%s_norm.tsv' % splitext(tab)[0]
    if config.force or not isfile(tab):
        cmd += assign_taxonomy(aligner, ali, tab, db_path, sub_db)
    if config.force or not isfile(tab_norm):
        cmd += normalize(tab, tab_norm)
    if cmd:
        append_cmd(list([cmd]), tab_norm, outputs, config)


def get_paths(out: str, aligner: str, db: str,
              outputs: dict, tax_fun: str) -> tuple:
    """

    Parameters
    ----------
    out : str
        Path to pipeline output folder for SHOGUN
    aligner : str
        Name of the aligner
    db : str
        Database name
    outputs : dict
        All outputs
    tax_fun : str
        "taxonomy" or "function"

    Returns
    -------
    tax : str
        Path to the output table
    tax_norm : str
        Path to the output table (relative abundances)
    """
    out_dir = get_dir(out, db, aligner)
    outputs['io']['I']['d'].add(out_dir)
    outputs['dirs'].append(out_dir)

    tab = '%s/%s.tsv' % (out_dir, tax_fun)
    tab_norm = '%s/%s_norm.tsv' % (out_dir, tax_fun)
    outputs['outs'].setdefault((db, aligner), []).extend([tab, tab_norm])
    outputs['io']['O']['f'].update([tab, tab_norm])
    outputs['io']['O']['f'].update([tab, tab_norm])
    return tab, tab_norm


def get_coverage_cmd(ali: str, db_path: str, level: str, cov_tab: str,
                     outputs: dict, config) -> None:
    """Get the SHOGUN command that calculates the coverage per taxon.

    Parameters
    ----------
    ali : str
        Path to the alignment
    db_path : str
        Path to the taxonomic database
    level : str
        Taxonomic level for coverage calculation
    cov_tab : str
        Path to the output coverage table file
    outputs : dict
        All outputs, incl. commands
    config
        Configuration
    """
    cmd = 'shogun coverage'
    cmd += ' -i %s' % ali
    cmd += ' -d %s' % db_path
    cmd += ' -l %s' % level
    cmd += ' -o %s' % cov_tab
    append_cmd(list([cmd]), cov_tab, outputs, config)


def coverage(out: str, aligner: str, ali: str, db: str, db_path: str,
             outputs: dict, config) -> None:
    if aligner == 'burst':
        cov_tab = '%s/coverage.tsv' % out
        outputs['outs'].setdefault((db, aligner), []).append(cov_tab)
        if config.force or not isfile(cov_tab):
            outputs['io']['O']['f'].append(cov_tab)
            get_coverage_cmd(ali, db_path, 'strain', cov_tab, outputs, config)


def taxonomy(out: str, aligner: str, ali: str, db: str, db_path: str,
             outputs: dict, config) -> str:
    """

    Parameters
    ----------
    out : str
        Path to pipeline output folder for SHOGUN
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    db : str
        Database name
    db_path : str
        Path to the taxonomic database
    outputs : dict
        All outputs
    config
        Configuration

    Returns
    -------
    norm : str
        Path to the normalized output table
    """
    tab, norm = get_paths(out, aligner, db, outputs, 'taxonomy')
    assign_normalize(tab, aligner, ali, db_path, outputs, '', config)
    redistribute(db, aligner, norm, db_path, outputs, config)
    return norm


def function(out: str, aligner: str, ali: str, db: str, db_path: str,
             outputs: dict, norm: str, config):
    for sub in ['kegg', 'refseq', 'uniprot']:
        fun, f_norm = get_paths(out, aligner, db, outputs, 'function-%s' % sub)
        assign_normalize(fun, aligner, ali, db_path, outputs,
                         '/functions-%s' % sub, config)
        redistribute(db, aligner, f_norm, db_path, outputs, config)

    out_dir = out + '/functional'
    for level in ['genus', 'species']:
        functional(norm, fun, out_dir, outputs, level, config)


def functional(norm: str, db_path: str, out_dir: str,
               outputs: dict, level: str, config) -> None:
    """

    Parameters
    ----------
    norm : str
        Path to the input table file
    db_path : str
        Path to the functional database
    out_dir : str
        Path to the output functional table
    outputs : dict
        All outputs
    level : str
        Taxonomic level
    config
        Configuration
    """
    cmd = 'shogun functional'
    cmd += ' -i %s' % norm
    cmd += ' -o %s' % out_dir
    cmd += ' -d %s' % db_path
    cmd += ' -l %s' % level
    base = splitext(basename(norm))[0]
    out = '%s/%s.%s.normalized.txt' % (out_dir, base, level)
    append_cmd(list([cmd]), out, outputs, config)


def shogun(prev: str, out_dir: str, sample: str, inputs: dict, params: dict,
           databases, config):
    """Get the SHOGUN commands that consist of the classifications using the
    databases built for shogun and the alignment in shogun, or not as this
    may be run after bowtie2, in which case the classifications are performed
    on the alignment obtained before, provided that it was done on a database
    formatted for SHOGUN.

    Parameters
    ----------
    prev : str
        Previous software
    out_dir : str
        Path to pipeline output folder for SHOGUN
    sample : str
        Sample name
    inputs : dict
        Input files
    params : dict
        Run parameters
    databases
        All databases class instance
    config
        Configuration class instance

    Returns
    -------
    outputs : list
        All outputs
    """
    outputs = {
        'io': {'I': {'d': set(), 'f': set()}, 'O': {'d': set(), 'f': set()}},
        'cmds': [], 'dirs': [], 'outs': dict({})}
    combine_cmds, ali_cmds = [], []

    out = out_dir + '/' + sample
    outputs['dirs'].append(out)
    outputs['io']['O']['d'].add(out)

    if prev == 'bowtie2':
        for (db, pairing), sam_ in inputs[sample].items():
            path = '%s/shogun' % databases.paths[db]
            ali = format_sam(sam_, ali_cmds, sample)
            taxonomy(out, '', ali, db, path, outputs, config)
    elif params['databases']:
        fasta = combine_inputs(sample, inputs, out, combine_cmds, outputs)
        for db, aligners in params['databases'].items():
            path = '%s/shogun' % databases.paths[db]
            for aligner in aligners:
                ali = align_cmds(
                    fasta, out, db, path, aligner, outputs, ali_cmds, params)
                norm = taxonomy(out, aligner, ali, db, path, outputs, config)
                coverage(out, aligner, ali, db, path, outputs, config)
                function(out, aligner, ali, db, path, outputs, norm, config)

    if outputs['cmds']:
        outputs['cmds'] = combine_cmds + ali_cmds + outputs['cmds']

    return outputs
