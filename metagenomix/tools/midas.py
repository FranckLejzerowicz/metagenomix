# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import isfile
from metagenomix._io_utils import io_update


def get_cmd(
        self,
        focus_dir: str,
        db: str,
        analysis: str,
        select: set = None) -> str:
    """Build the command line for MIDAS analysis.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .soft.params
            Parameters
    focus_dir : str
        Path to the output folder.
    db : str
        Path to MIDAS database.
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
    cmd += ' %s' % focus_dir
    cmd += ' -1 %s' % self.inputs[self.sam][0]
    if len(self.inputs[self.sam]) > 1:
        cmd += ' -2 %s' % self.inputs[self.sam][1]
    cmd += ' -d %s' % db
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' --remove_temp'
    if analysis != 'species':
        cmd += ' --species_cov 1'
    if select:
        cmd += ' --species_id %s' % ','.join(list(select))
    return cmd


def midas_species(self, focus_dir, db) -> None:
    if db == self.databases.paths['midas']:
        species_out = '%s/species' % focus_dir
        species_profile = '%s/species_profile.txt' % species_out
        if not self.config.force and isfile(species_profile):
            io_update(self, i_d=species_out)
        else:
            self.outputs['cmds'].append(get_cmd(self, focus_dir, db, 'species'))
            io_update(self, i_f=self.inputs[self.sam], o_d=species_out)
        self.outputs['outs'].append(species_out)
        self.outputs['dirs'].append(species_out)


def midas_genus(self, focus_dir, genes_out, db, select):
    if self.config.force or not isfile('%s/readme.txt' % genes_out):
        self.outputs['cmds'].append(
            get_cmd(self, focus_dir, db, 'genes', select))
        io_update(self, o_d=genes_out)
    self.outputs['outs'].append(genes_out)
    self.outputs['dirs'].append(genes_out)


def midas_snps(self, focus_dir, snps_out, db, select):
    if self.config.force or not isfile('%s/readme.txt' % snps_out):
        self.outputs['cmds'].append(
            get_cmd(self, focus_dir, db, 'snps', select))
        io_update(self, o_d=snps_out)
    self.outputs['outs'].append(snps_out)
    self.outputs['dirs'].append(snps_out)


def get_species_select(self, db: str, species_list: str) -> set:
    """Get the species names for which there is a reference in the database.

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
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
        if self.config.dev and not isfile(species_list):
            select.add('Escherichia')
            return select
        else:
            species = [x.strip() for x in open(species_list).readlines()]
        with open('%s/species_info.txt' % db) as f:
            for line in f:
                for genus_species in species:
                    if genus_species.replace(' ', '_') in line:
                        select.add(line.split('\t')[0])
    return select


def midas(self) -> None:
    """Create command lines for MIDAS for the current database's species focus.

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
        .databases
            All databases class instance
        .config
            Configurations
    """
    for focus, (db, species_list) in self.soft.params['focus'].items():

        focus_dir = '%s/%s/%s' % (self.dir, focus, self.sam)
        midas_species(self, focus_dir, db)

        select = set(get_species_select(self, db, species_list))

        genes_out = '%s/genes' % focus_dir
        midas_genus(self, focus_dir, genes_out, db, select)

        snps_out = '%s/snps' % focus_dir
        midas_snps(self, focus_dir, snps_out, db, select)
