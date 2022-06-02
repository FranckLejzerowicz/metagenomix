# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import isdir, isfile


def get_kraken2_db(db, databases, config):
    if db == 'default':
        db_path = databases.paths['kraken2']
    elif db in databases.paths:
        db_path = '%s/kraken2' % databases.paths[db]
        if not config.dev and not isdir(db_path):
            sys.exit('[kraken2] Not a database: %s' % db_path)
    else:
        sys.exit('[kraken2] Database not found: %s' % db)
    return db_path


def kraken2(out_dir: str, sam: str, inputs: dict,
            params: dict, databases, config) -> tuple:
    """Create command lines for kraken2.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    params : dict
        Kraken2 parameters.
    databases
        Path to the database.

    Returns
    -------
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.
    outputs : list
        All outputs paths.
    """
    io, cmds, outputs = {'I': [], 'O': []}, [], []
    for db in params['databases']:
        o = '%s/%s/%s' % (out_dir, sam, db)
        report = '%s/report.tsv' % o
        result = '%s/result.tsv' % o
        if config.force or not isfile(result):
            io['I'].extend(inputs)
            io['O'].extend([report, result])
            outputs.append(result)
            db_path = get_kraken2_db(db, databases, config)
            cmd = 'kraken2 '
            cmd += ' -db %s' % db_path
            cmd += ' --threads %s' % params['cpus']
            cmd += ' --report %s/report.tsv' % report
            cmd += ' --gzip-compressed'
            cmd += ' --confidence 0.5'
            if len(inputs[sam]) > 1:
                unclass = ['%s/unclassified_%s.fastq' % (o, r) for r in [1, 2]]
                classif = ['%s/classified_%s.fastq' % (o, r) for r in [1, 2]]
                cmd += ' --unclassified-out %s/unclassified#.fastq' % o
                cmd += ' --classified-out %s/classified#.fastq' % o
                cmd += ' --paired'
            else:
                unclass = ['%s/unclassified.fastq' % o]
                classif = ['%s/classified.fastq' % o]
                cmd += ' --unclassified-out %s/unclassified.fastq' % o
                cmd += ' --classified-out %s/classified.fastq' % o
            if inputs[sam][0].endswith('.gz'):
                cmd += ' --gzip-compressed'
            cmd += ' %s > %s' % (' '.join(inputs[sam]), result)
            cmds.append(cmd)
    return io, cmds, outputs


def get_species_select(db: str, species_list: str) -> set:
    """Get the species names for which there is a reference in the database.

    Parameters
    ----------
    db : str
        Database name.
    species_list : str
        Path to file containing list of species to focus on.

    Returns
    -------
    select : set
        Species names for which there is a reference in the database.
    """
    select = set()
    if species_list:
        midas_species = [x.strip() for x in open(species_list).readlines()]
        with open('%s/species_info.txt' % db) as f:
            for line in f:
                for genus_species in midas_species:
                    if genus_species.replace(' ', '_') in line:
                        select.add(line.split('\t')[0])
    return select


def get_cmd(
        out: str,
        inputs: list,
        db: str,
        cpus: str,
        analysis: str,
        select: set = None) -> str:
    """Build the command line for MIDAS analysis.

    Parameters
    ----------
    out : str
        Path to the output folder.
    inputs : list
        Input files.
    db : str
        Path to MIDAS database.
    cpus : str
        Number of CPUs.
    analysis : str
        MIDAS analysis (any of "species", "genes" or "snps").
    select : set
        Species names for which there is a reference in the database.

    Returns
    -------
    cmd : str
        Midas command line for the species level.
    """
    cmd = 'run_midas.py %s' % analysis
    cmd += ' %s' % out
    cmd += ' -1 %s' % inputs[0]
    if len(inputs) > 1:
        cmd += ' -2 %s' % inputs[1]
    cmd += ' -d %s' % db
    cmd += ' -t %s' % cpus
    cmd += ' --remove_temp'
    if analysis != 'species':
        cmd += ' --species_cov 1'
    if select:
        cmd += ' --species_id %s' % ','.join(list(select))
    return cmd
