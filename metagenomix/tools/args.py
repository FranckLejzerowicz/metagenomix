# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import basename
from metagenomix._io_utils import (caller, io_update, to_do, tech_specificity,
                                   not_paired, get_genomes_fastas)


def predict_cmd(
        self,
        typ: str,
        fasta: str,
        out: str,
) -> str:
    """Collect deeparg predict commands.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    typ : str
        'nucl' or 'prot' for fasta from assemblers or from prodigal/plass
    fasta : str
        Path to the input fasta file
    out : str
        Paths to the output file's prefix

    Returns
    -------
    cmd : str
        deeparg predict commands
    """
    cmd = 'deeparg predict'
    cmd += ' --input-file %s' % fasta
    cmd += ' --output-file %s' % out
    cmd += ' --data-path %s' % self.soft.params['database']

    for param in [
        'min_prob', 'arg_alignment_overlap', 'arg_alignment_evalue',
        'arg_alignment_identity', 'arg_num_alignments_per_entry',
        'model_version'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])

    if typ in ['nucl', 'prot']:
        cmd += ' --type %s' % typ
    else:
        cmd += ' --type nucl'

    if self.soft.prev == 'plass':
        cmd += ' --model SS'
    else:
        cmd += ' --model LS'
    return cmd


def get_predict(
        self,
        fastas: dict,
        tech: str,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    fastas : dict
        Fasta files per 'nucl' or 'prot' (type of data)
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
    """
    out = '%s/%s/%s' % (self.dir, tech, sam_group)
    self.outputs['dirs'].append(out)

    for typ, fasta in fastas.items():
        prefix = '%s/%s_%s' % (out, typ, sam_group)
        arg = '%s.ARG' % prefix
        pot_arg = '%s.potential.ARG' % prefix
        outs = [arg, pot_arg]
        self.outputs['outs'].setdefault((tech, sam_group), []).extend(outs)

        # check if the tool already run (or if --force) to allow getting command
        if self.config.force or to_do(arg):
            # collect the command line
            cmd = predict_cmd(self, typ, fasta, prefix)
            # add is to the 'cmds'
            key = '_'.join([tech, sam_group])
            self.outputs['cmds'][key] = [cmd]
            io_update(self, i_f=fasta, o_d=out, key=key)


def predict(self) -> None:
    """Run deeparg predict module.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name
        .pools : dict
            Pools
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group) in self.inputs[self.pool]:
            fastas = get_genomes_fastas(self, tech, group)
            get_predict(self, fastas, tech, group)
    else:
        for (tech, sam), contigs in self.inputs[self.sam].items():
            if self.soft.prev == 'plass' and tech != 'illumina':
                continue
            fastas = {basename(contigs).split('_')[0]: contigs}
            get_predict(self, fastas, tech, sam)


def short_cmd(
        self,
        fastqs: list,
        prefix: str
) -> str:
    """Collect deeparg short_reads_pipeline commands.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastqs : list
        Paths to the input files
    prefix : str
        Paths to the output file's prefix

    Returns
    -------
    cmd : str
        deeparg short_reads_pipeline commands
    """
    cmd = 'deeparg short_reads_pipeline'
    cmd += ' --forward_pe_file %s --reverse_pe_file %s' % tuple(fastqs)
    cmd += ' --output_file %s' % prefix
    for param in [
        'deeparg_data_path', 'deeparg_identity', 'deeparg_probability',
        'deeparg_evalue', 'gene_coverage', 'bowtie_16s_identity'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])
    return cmd


def short(self) -> None:
    """Run deeparg short_reads_pipeline module.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    # iterate over the inputs
    for (tech, sam), fastqs in self.inputs[self.sam].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, fastqs):
            continue

        # make the output directory
        out = '%s/%s/%s' % (self.dir, tech, self.sam)
        self.outputs['dirs'].append(out)

        # get the expected names of some of the ouptuts:
        # - those you want to collect in 'outs': will be used as future inputs
        # - at least one: will help to know whether the software already run
        prefix = out + '/' + self.sam
        arg = '%s.ARG' % prefix
        pot_arg = '%s.potential.ARG' % prefix
        outs = [arg, pot_arg]
        self.outputs['outs'].setdefault((tech, self.sam), []).extend(outs)

        # check if the tool already run (or if --force) to allow getting command
        if self.config.force or to_do(arg):
            # collect the commmand line
            cmd = short_cmd(self, fastqs, prefix)
            # add is to the 'cmds'
            self.outputs['cmds'][tech] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)


def deeparg(self) -> None:
    """Create command lines for metaWRAP

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    # This function splits the name of the software and calls as function
    # the last underscore-separated field (which is in this module)
    module_call = caller(self, __name__)
    module_call(self)
