# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import isfile


def midas(out_dir: str, sam: str, inputs: dict, path: str, cpus: str,
          focus: str, db_species: tuple) -> tuple:
    """Create command lines for MIDAS for the current database's species focus.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    path : str
        Path to the midas database.
    cpus : str
        Number of cpus.
    focus : str

    db_species : tuple
        (Database name, Path to file containing list of species to focus on.)

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
    db, species_list = db_species

    focus_dir = '%s/%s/%s' % (out_dir, focus, sam)
    if db == path:
        species_out = '%s/species' % focus_dir
        if isfile('%s/species_profile.txt' % species_out):
            io['I'].append(species_out)
        else:
            cmds.append(get_cmd(focus_dir, inputs[sam], db, cpus, 'species'))
            io['O'].append(species_out)
        outputs.append(species_out)

    select = set(get_species_select(db, species_list))

    genes_out = '%s/genes' % focus_dir
    if not isfile('%s/readme.txt' % genes_out):
        cmds.append(get_cmd(focus_dir, inputs[sam], db, cpus, 'genes', select))
        io['O'].append(genes_out)
    outputs.append(genes_out)

    snps_out = '%s/snps' % focus_dir
    if not isfile('%s/readme.txt' % snps_out):
        cmds.append(get_cmd(focus_dir, inputs[sam], db, cpus, 'snps', select))
        io['O'].append(snps_out)
    outputs.append(snps_out)

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
