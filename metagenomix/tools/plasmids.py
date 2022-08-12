# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import dirname

from metagenomix._inputs import (sample_inputs, group_inputs,
                                 genome_key, genome_out_dir)
from metagenomix._io_utils import io_update, to_do


def plasforest_cmd(
        self,
        in_fp: str,
        out_fp: str,
        out_dir: str
) -> str:
    """Collect plasforest command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    in_fp : str
        Path to the input file
    out_fp : str
        Path to the output file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        plasforest command
    """
    binary = self.soft.params['binary']
    cmd = 'cd %s\n' % out_dir
    cmd += 'cp %s/plasforest.sav %s/.\n' % (dirname(binary), out_dir)
    cmd += 'cp %s/*.fasta* %s/.\n' % (dirname(binary), out_dir)
    cmd += 'python3 %s' % binary
    cmd += ' -i %s' % in_fp
    cmd += ' -o %s' % out_fp
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    cmd += '\nrm %s/plasforest.sav\n' % out_dir
    cmd += 'rm %s/*.fasta*\n' % out_dir
    return cmd


def get_plasforest(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample or co-assembly name
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, fasta[0], sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        out_fp = '%s/plasmids.csv' % out_dir
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_fp)

        if self.config.force or to_do(out_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = plasforest_cmd(self, fasta[0], out_fp, out_dir)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta[0], i_d=out_dir, o_f=out_fp, key=key)


def plasmidfinder_cmd(
        self,
        fasta: str,
        out_dir: str,
        key: str
) -> str:
    """Collect PlasmidFinder command.

    Parameters
    ----------
    self : Commands class instance
        .pool : str
            Pool name.
        .soft.params
            PlasmidFinder parameters
    fasta : str
        Path to the input file
    out_dir : str
        Path to the output folder for the current sample/MAG
    key : str
        Concatenation of the technology and/or co-assembly pool group name

    Returns
    -------
    cmd : str
        PlasmidFinder command
    """
    tmp_dir = '$TMPDIR/plasmidfinder_%s' % key
    cmd = 'mkdir -p %s\n' % tmp_dir
    cmd += 'plasmidfinder.py'
    if len(fasta) == 2:
        cmd += ' --infile %s' % ' '.join(fasta)
    else:
        cmd += ' --infile %s' % fasta[0]
    cmd += ' --outputPath %s' % out_dir
    cmd += ' --tmp_dir %s' % tmp_dir
    cmd += ' --methodPath %s' % self.soft.params['methodPath']
    cmd += ' --databasePath %s' % self.soft.params['databasePath']
    if 'databases' in self.soft.params:
        cmd += ' --databases %s' % self.soft.params['databases']
    cmd += ' --mincov %s' % self.soft.params['mincov']
    cmd += ' --threshold %s' % self.soft.params['threshold']
    if self.soft.params['extented_output']:
        cmd += ' --extented_output'
    cmd += '\nrm -rf %s\n' % tmp_dir
    return cmd


def get_plasmidfinder(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample or co-assembly name
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, fasta[0], sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)

        json_fp = '%s/data.json' % out_dir
        if self.config.force or to_do(json_fp):
            key = genome_key(tech, sam_group, genome)
            cmd = plasmidfinder_cmd(self, fasta, out_dir, key)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, i_d=out_dir, o_d=out_dir, key=key)


def dispatch(self) -> None:
    """Classify assembly contigs or genomes/MAGs as plasmids or not.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
    """
    __plasmid_tool__ = getattr(sys.modules[__name__], 'get_%s' % self.soft.name)

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            __plasmid_tool__(self, tech, fastas, group)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            __plasmid_tool__(self, tech, fastas, mags)

    else:
        if self.soft.name == 'plasmidfinder':
            tech_fastas = sample_inputs(self, raw=True)
        else:
            tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            __plasmid_tool__(self, tech, fastas, self.sam_pool)


def plasmidfinder(self) -> None:
    """Detect plasmids in fasta files of contigs or binned/dereplicated MAGs
    using PlasmidFinder.

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for PlasmidFinder
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
        .sam
            Sample name
    """
    dispatch(self)


def plasforest(self) -> None:
    """Classify assembly contigs as plasmids or not using PlasForest.

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for PlasForest
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
        .sam
            Sample name
    """
    dispatch(self)
