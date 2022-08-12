# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
from os.path import basename, dirname, splitext
from metagenomix._inputs import group_inputs, genome_key, genome_out_dir
from metagenomix._io_utils import (caller, io_update, to_do, tech_specificity,
                                   not_paired)


def predict_cmd(
        self,
        fasta: str,
        out: str,
        typ: str
) -> str:
    """Collect deeparg predict commands.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fasta : str
        Path to the input fasta file
    out : str
        Paths to the output file's prefix
    typ : str
        Type of input data: 'nucl' (DNA) or 'prot' (protein)

    Returns
    -------
    cmd : str
        deeparg predict commands
    """
    cmd = 'deeparg predict'
    cmd += ' --input-file %s' % fasta
    cmd += ' --output-file %s' % out
    cmd += ' --data-path %s' % self.soft.params['db_dir']
    cmd += ' --type %s' % typ
    for param in [
        'min_prob', 'arg_alignment_overlap', 'arg_alignment_evalue',
        'arg_alignment_identity', 'arg_num_alignments_per_entry',
        'model_version'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])
    if self.soft.prev == 'plass':
        cmd += ' --model SS'
    else:
        cmd += ' --model LS'
    return cmd


def predict_inputs(
        self,
        fasta: str
) -> dict:
    """Get input files per sequence type.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
    fasta : str
        Path to the input file or folder

    Returns
    -------
    typ_seqs : dict
        Paths to the input fasta file per sequence type
    """
    if self.soft.prev == 'prodigal':
        typ_seqs = {'prot': '%s/protein.translations.fasta' % fasta,
                    'nucl': '%s/nucleotide.sequences.fasta' % fasta}
    elif self.soft.prev == 'plass':
        typ_seqs = {'prot': '%s/prot_contigs.fasta' % dirname(fasta),
                    'nucl': '%s/nucl_contigs.fasta' % dirname(fasta)}
    else:
        sys.exit('[%s] Only avail after prodigal or plass' % self.soft.name)
    return typ_seqs


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
        Paths to the input fasta files per genome/MAG
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str

    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, fasta[0], sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        typ_seqs = predict_inputs(self, fasta[0])
        for typ, seq in typ_seqs.items():

            key = genome_key(tech, sam_group, genome)
            if to_do(seq):
                self.soft.status.add('Run %s (%s)' % (self.soft.prev, key))

            base = splitext(basename(seq))[0]
            prefix = '%s/%s' % (out_dir, base)
            arg = '%s.mapping.ARG' % prefix
            pot_arg = '%s.potential.ARG' % prefix
            outs = [arg, pot_arg]
            self.outputs['outs'].setdefault((tech, sam_group), []).append(outs)

            # check if tool already run (or if --force) to allow getting command
            if self.config.force or to_do(arg):
                # collect the command line
                cmd = predict_cmd(self, seq, prefix, typ)
                # add is to the 'cmds'
                key += '_' + typ
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta[0], o_d=out_dir, key=key)


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
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_predict(self, fastas, tech, group)

    elif self.soft.prev == 'drep':
        for (tech, bin_algo), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_predict(self, fastas, tech, bin_algo)
    else:
        for (tech, sam), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
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
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, fastqs):
            continue

        # make the output directory
        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        # get the expected names of some of the ouptuts:
        # - those you want to collect in 'outs': will be used as future inputs
        # - at least one: will help to know whether the software already run
        prefix = out + '/' + self.sam_pool
        arg = '%s.mapping.ARG' % prefix
        pot_arg = '%s.potential.ARG' % prefix
        outs = [arg, pot_arg]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

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
