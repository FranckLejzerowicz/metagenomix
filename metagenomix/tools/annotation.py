# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys

from os.path import basename, dirname, splitext
from metagenomix._io_utils import (caller, get_out_dir, write_hmms,
                                   io_update, to_do, get_genomes_fastas)
from metagenomix.parameters import tech_params


def get_input(
        self,
        index: int = 0
) -> dict:
    """Get input files for Prodigal or an annotation tool taking Prodigal
    output as input.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to software's output folder
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    index : int
        Index in the list of input files (i.e. outputs from the previous step)

    Returns
    -------
    inputs : dict
        Path to the assembly contig file per technology
    """
    tool = self.soft.name
    if tool == 'prodigal':
        if self.soft.prev == 'plass' and 'illumina' in self.inputs[self.sam]:
            if self.softs['plass'].params['type'] == 'assembly':
                sys.exit('[%s] Can only predict for plass nuclassembly' % tool)
            inputs = {'illumina': self.inputs[self.sam]['illumina']}
        else:
            inputs = {t: self.inputs[self.sam][t][index] for t in [
                'pacbio', 'nanopore'] if self.inputs[self.sam].get(t)}
    else:
        inputs = {t: self.inputs[self.sam][t][index] for t in [
            'illumina', 'pacbio', 'nanopore'] if self.inputs[self.sam].get(t)}
    return inputs


def prodigal_cmd(
        self,
        contig_fp: str,
        outputs: list,
        tech: str
) -> str:
    """Create command lines for Prodigal.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prodigal
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    contig_fp : str
        Path to the input contigs file
    outputs : list
        Paths to the output files
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    --------
    cmd : str
        Prodigal command
    """
    params = tech_params(self, tech)
    cmd = ''
    contig_fa = contig_fp
    if contig_fp.endswith('fastq.gz'):
        contig_fa = contig_fp.replace('fastq.gz', 'fasta')
        cmd += 'seqtk seq -a %s > %s\n' % (contig_fp, contig_fa)
    cmd = 'prodigal'
    cmd += ' -i %s' % contig_fa
    cmd += ' -a %s' % outputs[0]
    cmd += ' -d %s' % outputs[1]
    cmd += ' -s %s' % outputs[2]
    cmd += ' -o %s' % outputs[3]
    for param in ['f', 'p']:
        cmd += ' -%s %s' % (param, params[param])
    for boolean in ['c', 'm', 'n', 'q']:
        if params[boolean]:
            cmd += ' -%s' % params[boolean]
    return cmd


def prodigal_outputs(
        self,
        out: str
) -> list:
    """Get the output files of Prodigal.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    out : str
        Path to the output folder

    Returns
    -------
    outputs : list
        Paths to the Prodigal output files
    """
    proteins = '%s/protein.translations.fasta' % out
    nucleotides = '%s/nucleotide.sequences.fasta' % out
    genes = '%s/potential.starts.fasta' % out
    gbk = '%s/gene.coords.%s' % (out, self.soft.params['f'])
    outputs = [proteins, nucleotides, genes, gbk]
    return outputs


def get_prodigal(
        self,
        contig_fp: str,
        out: str,
        tech: str,
        sam_group: str,
):
    """Get the prodigal command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    contig_fp : str
        Path to the input contigs file
    out : str
        Path to the output folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    """
    self.outputs['dirs'].append(out)
    outputs = prodigal_outputs(self, out)
    self.outputs['outs'][(tech, sam_group)] = outputs
    tech_sam_group = '_'.join([tech, sam_group])
    if self.config.force or to_do(outputs[0]):
        cmd = prodigal_cmd(self, contig_fp, outputs, tech)
        self.outputs['cmds'].setdefault(tech_sam_group, []).append(cmd)
        io_update(self, i_f=contig_fp, o_f=outputs, o_d=out, key=tech_sam_group)


def prodigal(self) -> None:
    """Create command lines for Prodigal.

    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prodigal
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group), inputs in self.inputs[self.pool].items():
            out = '/'.join([self.dir, tech, self.pool, group])
            get_prodigal(self, inputs[1], out, tech, group)
    else:
        tech_contigs = get_input(self)
        for tech, contigs in tech_contigs.items():
            out = '/'.join([self.dir, tech, self.sam])
            get_prodigal(self, contigs, out, tech, self.sam)


def macsyfinder_cmd(
        self,
        fp: str,
        out_dir: str,
        models: str,
        tech: str,
) -> tuple:
    """Collect command for MacSyFinder.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    fp : str
        Path to the input fasta file
    out_dir : str
        Path to the output folder
    models : str
        Path to the models folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    -------
    cmd : str
        MacSyFinder command
    outs : list
        Path to the output folders
    """
    cmd = ''
    input_fp = fp
    params = tech_params(self, tech)
    if tech in ['illumina', 'pacbio', 'nanopore']:
        input_fp = '%s_edit.fasta' % fp.replace('.fasta', '')
        cmd += 'header_space_replace.py -i %s -o %s --n\n' % (fp, input_fp)

    outs = []
    for model in params['models']:
        model_dir = '%s/%s' % (out_dir, model)
        outs.append(model_dir)
        self.outputs['dirs'].append(model_dir)
        res = '%s/macsyfinder.log' % model_dir
        if self.config.force or to_do(res):
            cmd += 'macsyfinder'
            for param in [
                'db_type', 'replicon_topology', 'e_value_search',
                'i_evalue_sel', 'coverage_profile', 'mandatory_weight',
                'accessory_weight', 'exchangeable_weight', 'redundancy_penalty',
                'out_of_cluster'
            ]:
                cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
            cmd += ' --out-dir %s' % model_dir
            cmd += ' --res-search-suffix _hmm.tsv'
            cmd += ' --res-extract-suffix _out.tsv'
            cmd += ' --worker %s' % params['cpus']
            cmd += ' --sequence-db %s' % input_fp
            cmd += ' --models-dir %s/data/models' % models
            cmd += ' --models %s all' % model
            cmd += ' --verbosity\n'
    return cmd, outs


def get_macsyfinder(
        self,
        fp: str,
        out_dir: str,
        tech: str,
        sam_group: str,
) -> None:
    """Get the MacSyFinder command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Databases
        .outputs : dict
            All outputs
    fp : str
        Path to the input fasta file
    out_dir : str
        Path to the output folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    """
    self.outputs['dirs'].append(out_dir)
    models = self.databases.paths['macsyfinder']
    cmd, outs = macsyfinder_cmd(self, fp, out_dir, models, tech)
    self.outputs['outs'][(tech, sam_group)] = outs
    if cmd:
        tech_sam_group = '_'.join([tech, sam_group])
        self.outputs['cmds'].setdefault(tech_sam_group, []).append(cmd)
        io_update(self, i_f=fp, i_d=models, o_d=out_dir, key=tech_sam_group)


def macsyfinder(self) -> None:
    """Detect different types of secretion systems using MacSyFinder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MacSyFinder
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group) in self.inputs[self.pool]:
            out_dir, fp = get_out_dir(self, self.pool, tech, group)
            get_macsyfinder(self, fp, out_dir, tech, group)
    else:
        tech_proteins = get_input(self)
        for tech, _ in tech_proteins.items():
            out_dir, fp = get_out_dir(self, self.sam, tech)
            get_macsyfinder(self, fp, out_dir, tech, self.sam)


def integronfinder_cmd(
        self,
        fp: str,
        o_dir: str,
        tech: str,
        group: str
) -> str:
    """Collect command for integron_finder.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Path to the reference databases
    fp : str
        Path to the input fasta file
    o_dir : str
        Path to the output folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    group : str
        Group for the current co-assembly

    Returns
    -------
    cmd : str
        integron_finder command
    """
    params = tech_params(self, tech)
    fp_out = fp.replace('.fasta', '_%s.fasta' % params['min_length'])
    cmd, hmms = write_hmms(self)
    cmd += 'filter_on_length.py'
    cmd += ' -i %s' % fp
    cmd += ' -o %s' % fp_out
    cmd += ' -t %s\n' % params['min_length']
    cmd += 'integron_finder'
    for boolean in ['local_max', 'promoter_attI', 'mute', 'pdf',
                    'gbk', 'union_integrases']:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --verbose'
    cmd += ' --outdir %s' % o_dir
    cmd += ' --cpu %s' % params['cpus']
    if self.databases.hmms_dias:
        cmd += ' --func-annot'
        cmd += ' --path-func-annot %s' % hmms
        io_update(self, i_f=hmms, key='_'.join([tech, group]))
    if params['prot_file']:
        cmd += ' --prot-file %s' % params['prot_file']
    if params['attc_model']:
        cmd += ' --attc-model %s' % params['attc_model']
        cmd += ' --evalue-attc %s' % params['evalue_attc']
        cmd += ' --max-attc-size %s' % params['max_attc_size']
        cmd += ' --min-attc-size %s' % params['min_attc_size']
    if params['topology_file']:
        cmd += ' --topology-file %s' % params['topology_file']
    else:
        cmd += ' --%s' % params['topology']
    cmd += ' %s' % fp_out
    return cmd


def get_integronfinder(
        self,
        fp: str,
        out_dir: str,
        tech: str,
        sam_group: str,
) -> None:
    """Get the MacSyFinder command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Databases
        .outputs : dict
            All outputs
    fp : str
        Path to the input fasta file
    out_dir : str
        Path to the output folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    """
    self.outputs['dirs'].append(out_dir)
    self.outputs['outs'].setdefault((tech, self.pool), []).append(out_dir)
    out = '%s/Results_Integron_Finder_mysequences/mysequences.summary' % out_dir
    if self.config.force or to_do(out):
        tech_sam_group = '_'.join([tech, sam_group])
        cmd = integronfinder_cmd(self, fp, out_dir, tech, sam_group)
        self.outputs['cmds'].setdefault(tech_sam_group, []).append(cmd)
        io_update(self, i_f=fp, o_d=out_dir, key=tech_sam_group)


def integronfinder(self) -> None:
    """Finder integrons, integration sites and cassettes with integron_finder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for integron_finder
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group) in self.inputs[self.pool]:
            out_dir, fp = get_out_dir(self, self.pool, tech, group)
            get_integronfinder(self, fp, out_dir, tech, group)
    else:
        tech_proteins = get_input(self)
        for tech, _ in tech_proteins.items():
            out_dir, fp = get_out_dir(self, self.sam, tech)
            get_integronfinder(self, fp, out_dir, tech, self.sam)


# def custom_cmd(
#         self,
#         input_file,
#         out_dir,
#         dia_hmm
# ) -> list:
#     outs = {}
#     for target, gene_hmm_dia in self.databases.hmms_dias.items():
#         if target == 'cazy':
#             continue
#         for gene, (hmm, dia) in gene_hmm_dia.items():
#             sam_dir = '%s/%s' % (out_dir, gene)
#             self.outputs['dirs'].append(sam_dir)
#             out = '%s/%s.tsv' % (sam_dir, self.sam)
#             if to_do(out):
#                 if dia_hmm == 'hmm':
#                     stdout = '%s_hmmer.out' % splitext(out)[0]
#                     io_update(self, i_f=hmm, o_f=[out, stdout], key=target)
#                     cmd = 'hmmsearch'
#                     cmd += ' --cut_tc'
#                     cmd += ' --tblout %s' % out
#                     cmd += ' -o %s' % stdout
#                     cmd += ' --cpu %s' % self.soft.params['cpus']
#                     cmd += ' %s %s' % (hmm, input_file)
#                 else:
#                     io_update(self, i_f=dia, o_f=out, key=target)
#                     tmp = '%s/%s_tmp' % (sam_dir, self.sam)
#                     cmd = 'mkdir -p %s\n' % tmp
#                     cmd += 'diamond blastp'
#                     cmd += ' -d %s' % dia
#                     cmd += ' -q %s' % input_file
#                     cmd += ' -o %s' % out
#                     cmd += ' -k 1 --id 80'
#                     cmd += ' -p %s' % self.soft.params['cpus']
#                     cmd += ' -t %s' % tmp
#                 self.outputs['cmds'].setdefault(target, []).append(cmd)
#             else:
#                 io_update(self, i_f=out, key=target)
#             out2 = out.replace('.tsv', '_contigs.tsv')
#             if to_do(out2):
#                 io_update(self, o_f=out2, key=target)
#                 cmd = 'extract_custom_searched_contigs.py'
#                 cmd += ' -i %s' % out
#                 cmd += ' -o %s' % out2
#                 self.outputs['cmds'].setdefault(target, []).append(cmd)
#             outs.setdefault(target, []).extend([out, out2])
#     return outs


def search_cmd(
        self,
        out_dir: str,
        fp: str,
        tech: str,
        sam_group: str
) -> list:
    """Dispatch the command line creation for hmmer or diamond and
    for the different reference databases.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample or co-assembly pool name
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    out_dir : str
        Path to the main output folder
    fp : str
        Path to the input fasta file
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly

    Returns
    -------
    outs : list
        Paths to the output files
    """
    outs = []
    params = tech_params(self, tech)
    tmp = '$TMPDIR/%s_%s_%s_%s' % (self.soft.name, self.sam, tech, sam_group)
    module_call = caller(self, __name__)
    for db, db_paths in params['databases'].items():
        for db_path in db_paths:
            db_base = splitext(basename(db_path))[0]
            db_out_dir = '%s/%s' % (out_dir, db_base)
            self.outputs['dirs'].append(db_out_dir)
            out = '%s/%s.tsv' % (db_out_dir, db_base)
            outs.append(db_out_dir)
            if self.config.force or to_do(out):
                tech_sam_group = '_'.join([tech, sam_group])
                cmd = module_call(params, fp, db_path, out, tmp)
                io_update(self, i_f=fp, o_d=db_out_dir, key=tech_sam_group)
                self.outputs['cmds'].setdefault(tech_sam_group, []).append(cmd)
    return outs


def diamond(
        params: dict,
        fp: str,
        db_path: str,
        out: str,
        tmp_dir: str
) -> str:
    """Collect the command line for DIAMOND.

    Parameters
    ----------
    params : dict
        Parameters for diamond and the current technology
    fp : str
        Path to the input fasta file
    db_path : str
        Path to the reference hmm database
    out : str
        Path to the output file
    tmp_dir : str
        Path to the temporary folder

    Returns
    -------
    cmd : str
        diamond command
    """
    cmd = 'mkdir -p %s\n' % tmp_dir
    cmd += 'diamond blastp'
    cmd += ' --db %s' % db_path
    cmd += ' --query %s' % fp
    cmd += ' --out %s' % out
    if params['mode']:
        cmd += ' --%s' % params['mode']
    if params['top']:
        cmd += ' --top %s' % params['top']
    else:
        cmd += ' -k %s' % params['max_target_seqs']
    cmd += ' --threads %s' % params['cpus']
    for param in ['id', 'strand', 'masking', 'evalue',
                  'query_cover', 'subject_cover']:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    cmd += ' --tmpdir %s\n' % tmp_dir
    return cmd


def hmmer(
        params: dict,
        fp: str,
        db_path: str,
        out: str,
        _: str
) -> str:
    """Collect the command line for HMMER.

    Parameters
    ----------
    params : dict
        Parameters for hmmer and the current technology
    fp : str
        Path to the input fasta file
    db_path : str
        Path to the reference hmm database
    out : str
        Path to the output file
    _ : str
        Unused argument

    Returns
    -------
    cmd : str
        hmmsearch command
    """
    stdout = '%s_per-sequence.out' % splitext(out)[0]
    domout = '%s_per-domain.out' % splitext(out)[0]
    cmd = 'hmmsearch'
    cmd += ' -o %s' % stdout
    cmd += ' --tblout %s' % out
    cmd += ' --domtblout %s' % domout
    for boolean in ['nobias', 'noali', 'max']:
        if params[boolean]:
            cmd += ' --%s' % boolean
    for cut in ['cut_ga', 'cut_nc', 'cut_tc']:
        if params[cut]:
            cmd += ' --%s' % cut
            break
    for param in ['E', 'Z', 'domE', 'domZ', 'textw']:
        cmd += ' --%s %s' % (param, params[param])
    cmd += ' --cpu %s' % params['cpus']
    cmd += ' %s' % db_path
    cmd += ' %s' % fp
    return cmd


def get_search(
        self,
        fp: str,
        out_dir: str,
        tech: str,
        sam_group: str
) -> None:
    """Get the DIAMOND or HMMER search command
    and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Databases
        .outputs : dict
            All outputs
    out_dir : str
        Path to the output folder
    fp : str
        Path to the input fasta file
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    """
    outs = search_cmd(self, out_dir, fp, tech, sam_group)
    if outs:
        self.outputs['outs'].setdefault((tech, sam_group), []).extend(outs)


def search(self) -> None:
    """Create command lines for custom searches based on DIAMOND or HMMER.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for searches (DIAMOND or HMMER)
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .soft.prev : str
            Previous software
        .config
            Configurations
    """
    # This function splits the name of the software and calls as function (
    # from this module) the part of the software name that follows the first
    # underscore, e.g. software "search_diamond" would call `diamond()`
    if self.pool in self.pools:
        for (tech, group) in self.inputs[self.pool]:
            out_dir, fp = get_out_dir(self, self.pool, tech, group)
            get_search(self, fp, out_dir, tech, group)
    else:
        tech_proteins = get_input(self)
        for tech, _ in tech_proteins.items():
            out_dir, fp = get_out_dir(self, self.sam, tech)
            get_search(self, fp, out_dir, tech, self.sam)


def prokka_cmd(
        self,
        contigs: str,
        out_dir: str,
        pref: str,
        config: dict,
        cols: list
) -> str:
    """Collect the command line for Prokka.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    contigs : str
        Path to the contigs input file
    out_dir : str
        Path to the output folder
    pref : str
        Prefix to add to the output files
    config : dict
        One field#:taxon dict obtained by parsing the user-defined config
    cols : list
        Columns of the taxonomic config: 'genus', 'species', 'strain', 'plasmid'

    Returns
    -------
    cmd : str
        Prokka command
    """
    cmd = 'prokka'
    cmd += ' --mincontiglen %s' % self.soft.params['mincontiglen']
    cmd += ' --cpus %s' % self.soft.params['cpus']
    for boolean in ['metagenome', 'notrna', 'norrna']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --force'
    for col in cols:
        if col and config[col]:
            cmd += ' --%s %s' % (col, config[col])
    cmd += ' --prefix %s' % pref
    cmd += ' --outdir %s' % out_dir
    cmd += ' %s' % contigs
    return cmd


def prokka_config(
        self,
        cols: list
) -> list:
    """Parse the user-defined taxonomic config file to collect the taxonomic
    levels to annotate more precisely using Prokka.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    cols : list
        Columns of the taxonomic config: 'genus', 'species', 'strain', 'plasmid'

    Returns
    -------
    configs : list
        List of field#:taxon dicts obtained by parsing the user-defined config
    """
    configs = []
    cols += ['proteins']
    if 'taxa' in self.soft.params:
        with open(self.soft.params['taxa']) as f:
            for ldx, line in enumerate(f):
                ls = line.strip('\n').split('\t')
                if not ldx:
                    continue
                configs.append({cols[idx]: val for idx, val in enumerate(ls)})
    return configs


def prokka_configs(self) -> tuple:
    """Get the config files for the different run of Prokka at
    different taxonomic levels.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters

    Returns
    -------
    configs : list
        Dicts of config files per taxonomic level
    cols : list
        Taxonomic levels
    """
    if self.soft.params['config']:
        cols = ['genus', 'species', 'strain', 'plasmid']
        configs = prokka_config(self, cols)
    else:
        cols = ['']
        configs = [{'': 'no_config'}]
    return configs, cols


def get_prokka(
        self,
        inputs: list,
        out_dir: str
) -> str:
    """Get the Prokka commands for the different configs
    (at each taxonomic level).

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .config
            Configurations
    inputs : list
        Paths to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Prokka command
    """
    cmd = ''
    configs, cols = prokka_configs(self)
    for config in configs:
        pref = '_'.join([config[x] for x in cols if config[x]])
        if config.get('proteins'):
            pref += '_%s' % splitext(basename(config['proteins']))[0]
        file_out = '%s/%s.out' % (out_dir, pref)
        if self.config.force or to_do(file_out):
            cmd += prokka_cmd(self, inputs[1], out_dir, pref, config, cols)
    return cmd


def prokka(self) -> None:
    """Perform prokaryotic genome annotation using Prokka and possibly,
    for specific prokaryotic taxa as per the user-defined config file.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prokka
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group), inputs in self.inputs[self.pool].items():
            out_dir = '/'.join([self.dir, tech, self.pool, group])
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'].setdefault(group, []).append(out_dir)

            cmd = get_prokka(self, inputs, out_dir)
            if cmd:
                tech_group = '_'.join([tech, group])
                self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
                io_update(self, i_f=inputs[1], o_d=out_dir, key=tech_group)


def barrnap_cmd(
        self,
        out: str,
        fasta: str,
        tech: str,
) -> str:
    """Collect the command line for Barrnap.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    out : str
        Path to the output file
    fasta : str
        Path to the input fasta file
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    -------
    cmd : str
        Barrnap command
    """
    params = tech_params(self, tech)
    cmd = 'barrnap'
    cmd += ' --kingdom %s' % params['kingdom']
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --reject %s' % params['reject']
    cmd += ' --lencutoff %s' % params['lencutoff']
    cmd += ' --evalue %s' % params['evalue']
    if params['incseq']:
        cmd += ' --incseq'
    cmd += ' --outseq %s' % out
    cmd += ' %s' % fasta
    return cmd


def ccmap_cmd(
        self,
        out: str,
        fasta: str,
        tech: str,
) -> str:
    """Collect the command line for CCmap.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .config
            Configurations
    out : str
        Path to the output file
    fasta : str
        Path to the input fasta file
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    -------
    cmd : str
        CCMap command
    """
    params = tech_params(self, tech)
    out_dir = dirname(out)
    cmd = ''
    if self.config.force:
        cmd += 'if [ -d "%s" ]; then rm -rf %s; fi\n' % (out_dir, out_dir)
    cmd += 'ccfind %s %s' % (fasta, out_dir)
    cmd += ' --terminal-fragment-size %s' % params['terminal_fragment_size']
    cmd += ' --min-percent-identity %s' % params['min_percent_identity']
    cmd += ' --min-aligned-length %s' % params['min_aligned_length']
    cmd += ' --ncpus %s' % params['cpus']
    if params['preserve_tmpdir']:
        cmd += ' --preserve-tmpdir'
    return cmd


def get_generic_on_fasta(
        self,
        fasta: str,
        tech: str,
        sam_group: str
) -> None:
    """Get the Barrnap or CCmap command
    and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prokka
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    fasta : str
        Path to the input fasta file
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or name of a co-assembly pool's group
    """
    out_dir = '/'.join([self.dir, tech, self.pool])
    if self.pool != sam_group:
        out_dir += '/%s' % sam_group
    self.outputs['dirs'].append(out_dir)

    out = '%s/%s.fas' % (out_dir, splitext(basename(fasta))[0])
    self.outputs['outs'].setdefault((tech, sam_group), []).append(out)

    if self.soft.name == 'barrnap':
        condition = not to_do(out)
    elif self.soft.name == 'ccmap':
        condition = not to_do(folder=dirname(out))
    else:
        print(error)

    if self.config.force or not condition:
        soft_cmd = '%s_cmd' % self.soft.name
        cmd = globals()[soft_cmd](self, out, fasta, tech)
        tech_sam_group = '_'.join([tech, sam_group])
        self.outputs['cmds'].setdefault(tech_sam_group, []).append(cmd)
        io_update(self, i_f=fasta, o_d=out, key=tech_sam_group)


def generic_on_fasta(self) -> None:
    """Run either Barrnap or CCmap on the fasta files resulting from
    the assembly or genome dereplication step.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for filtering
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    if self.pool in self.pools:
        for (tech, group) in self.inputs[self.pool]:
            fastas_d = get_genomes_fastas(self, tech, group)
            for _, fastas in fastas_d.items():
                for fasta in fastas:
                    get_generic_on_fasta(self, fasta, tech, group)
    else:
        tech_inputs = get_input(self)
        for tech, fasta in tech_inputs.items():
            get_generic_on_fasta(self, fasta, tech, self.sam)


def barrnap(self) -> None:
    """Dispatch the collection of Barrnap commands using generic function.

    Parameters
    ----------
    self : Commands class instance
    """
    generic_on_fasta(self)


def ccmap(self) -> None:
    """Dispatch the collection of CCmap commands using generic function.

    Parameters
    ----------
    self : Commands class instance
    """
    generic_on_fasta(self)


def antismash_cmd(
        self,
        fa: str,
        base: str,
        out: str
) -> str:
    """Collect the command line for Antismash.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .soft.prev : str
            Previous software
    fa : str
        Path to the input genome fasta file
    base : str
        Basename for the output files
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Antismash command
    """
    cmd = '\nantismash'
    for boolean in ['rre', 'asf', 'cassis', 'pfam2go', 'tigrfam', 'fullhmmer',
                    'cc-mibig', 'cb_general', 'smcog_trees', 'clusterhmmer',
                    'cb_subclusters', 'cb_knownclusters']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --tta-threshold %s' % self.soft.params['--tta-threshold']
    cmd += ' --genefinding-tool %s' % self.soft.params['genefinding_tool']
    cmd += ' --html-description after_%s' % self.soft.prev
    cmd += ' --taxon %s' % self.soft.params['taxon']
    cmd += ' --cpus %s' % self.soft.params['cpus']
    cmd += ' --output-basename %s' % base
    cmd += ' --html-title "%s: %s"' % base
    cmd += ' --output-dir %s' % out
    cmd += ' %s\n' % fa
    return cmd


def antismash(self) -> None:
    """Perform metabolomic profiling of dereplicated genomes using Antismash.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Antismash
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .soft.prev : str
            Previous software
        .config
            Configurations
    """
    if self.soft.prev != 'drep':
        sys.exit('[antismash] Possible only after drep (!=%s)' % self.soft.prev)
    for pool in self.pools:
        for algo, folder in self.inputs[pool].items():
            out_dir = '/'.join([self.dir, pool, algo])
            io_update(self, o_d=out_dir, key=pool)
            for fa in glob.glob('%s/*.fa' % folder):
                io_update(self, i_f=fa, key=pool)
                base = fa.split('/')[-1].replace('.fa', '')
                out = '%s/%s' % (out_dir, base)
                self.outputs['dirs'].append(out)
                outs = glob.glob('%s/*' % out.replace('${SCRATCH_FOLDER}', ''))
                if self.config.force or not outs:
                    cmd = antismash_cmd(self, fa, base, out)
                    self.outputs['cmds'].setdefault(
                        (pool, algo), []).append(cmd)


def tiara(self) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for tiara
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
    if self.soft.prev != 'spades':
        sys.exit('[tiara] can only be run on assembly output')
    for (tech, group), spades_outs in self.inputs[self.pool].items():
        out_dir = '/'.join([self.dir, tech, self.pool, group])
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir
        out_fp = '%s/classifications.txt' % out_dir
        if self.config.force or to_do(out_fp):
            cmd = 'cd %s\n' % out_dir
            cmd += 'tiara'
            cmd += ' --input %s' % spades_outs[1]
            cmd += ' --output %s' % out_fp
            cmd += ' --threads %s' % self.soft.params['cpus']
            for boolean in ['probabilities', 'verbose', 'gzip']:
                if self.soft.params[boolean]:
                    cmd += ' --%s' % boolean
            for param in ['min_len', 'first_stage_kmer', 'second_stage_kmer']:
                cmd += ' --%s %s' % (param, self.soft.params[param])
            for param in ['prob_cutoff', 'to_fasta']:
                cmd += ' --%s %s' % (param, ' '.join(self.soft.params[param]))
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out_dir, key=group)


# def write_dbcan_subset(
#         self,
#         taxa: list,
#         folder: str,
#         name: str
# ) -> None:
#     """Write the fasta file suset to the target features and
#     the command to make it a diamond database.
#
#     Parameters
#     ----------
#     self : Commands class instance
#         .databases : dict
#             Path to the reference databases
#         .outputs : dict
#             All outputs
#     taxa : list
#         Taxa from the user file.
#     folder : str
#         dbCAN-Seq folder.
#     name : str
#         Current subset name.
#
#     Returns
#     -------
#     cmd : str
#         Command to make the diamond db from the subsets fasta file.
#     """
#     path = self.databases['dbcan']
#     cmd = ""
#     fas_fp = '%s/%s.fa' % (folder, name)
#     dia_fp = '%s.dmnd' % splitext(fas_fp)[0]
#     io_update(self, i_f=fas_fp, o_f=dia_fp)
#     if to_do(dia_fp):
#         with open(fas_fp, 'w') as o:
#             taxon_found = False
#             for taxon in taxa:
#                 meta_taxon_pd = self.databases.dbcan_meta.loc[
#                     self.databases.dbcan_meta.genome_name.str.contains(taxon),:]
#                 if not meta_taxon_pd.shape[0]:
#                     continue
#                 gcf_dir = '%s/dbCAN-seq/CAZyme_seq_list' % path
#                 if to_do(folder=gcf_dir):
#                     os.makedirs(gcf_dir)
#                 for gcf in set(meta_taxon_pd.index):
#                     gcf_fas = '%s/%s.fasta' % (gcf_dir, gcf)
#                     if not to_do(gcf_fas):
#                         for e in read(gcf_fas, 'fasta'):
#                             o.write('>%s\n%s\n' % (e.metadata['id'], e))
#                         taxon_found = True
#         if taxon_found:
#             cmd = "diamond makedb --in %s -d %s\n" % (fas_fp, dia_fp)
#     self.outputs['cmds'].append(cmd)
#
#
# def set_dbcan_taxa(self) -> None:
#     """
#
#     Parameters
#     ----------
#     self : Commands class instance
#         .soft.params : dict
#             Parameters
#     """
#     for name, fp in self.soft.params['taxa'].items():
#         taxa = list(reads_lines(fp))
#         if not taxa:
#             continue
#         folder = '%s/subsets' % self.dir
#         if to_do(folder=folder):
#             os.makedirs(folder)
#         self.databases.cazys[name] = folder
#         write_dbcan_subset(self, taxa, folder, name)
#
#
# def cazy(self) -> None:
#     """
#
#     Parameters
#     ----------
#     self : Commands class instance
#     """
#     set_dbcan_taxa(self)
#
#
# def prepare_ioncom_inputs(
#         self,
#         ion_com: str,
#         fp: str
# ) -> str:
#     """Create the preformatting command lines before running IonCom.
#
#     Parameters
#     ----------
#     self : Commands class instance
#         .dir : str
#             Path to pipeline output folder for IonCom
#         .sam : str
#             Sample name
#         .inputs : dict
#             Input files
#         .outputs : dict
#             All outputs
#         .soft.params : dict
#             Parameters
#         .config
#             Configurations
#     tech : str
#         Technology: 'illumina', 'pacbio', or 'nanopore'
#     ion_com : str
#         Path to the IonCom standalone
#     fp : str
#         Path to the output file
#     """
#     cmd = 'prepare_ioncom_inputs.py'
#     cmd += ' -d %s -p %s -n 1' % (ion_com, fp)
#     cmd += '%s/run_IonCom.pl' % ion_com
#     return cmd
#
#
# def ioncom(self) -> None:
#     """Create command lines for IonCom.
#
#     Parameters
#     ----------
#     self : Commands class instance
#         .dir : str
#             Path to pipeline output folder for IonCom
#         .sam : str
#             Sample name
#         .inputs : dict
#             Input files
#         .outputs : dict
#             All outputs
#         .soft.params : dict
#             Parameters
#         .config
#             Configurations
#     """
#     # -------------------------------------------------------
#     # --- WILL NEED TO BE REVISED AFTER I-TASSER DOWNLOAD ---
#     # -------------------------------------------------------
#     i_tasser_libs = self.databases.ioncom['itasser']
#     i_tasser = self.soft.params['itasser']
#     ion_com = '%s/IonCom_standalone' % self.databases.paths['ioncom']
#
#     if self.pool in self.pools:
#         for group in self.pools[self.pool]:
#             o_dir, fp = get_out_dir(self, self.pool, group)
#             cmd = prepare_ioncom_inputs(self, ion_com, fp)
#             self.outputs['cmds'].setdefault(group, []).append(cmd)
#             output = '%s/output' % ion_com
#             cmd = 'for i in %s/*rehtml\n' % output
#             cmd += 'do\n'
#             cmd += '    mv %s/*rehtml %s/.\n' % (output, o_dir)
#             cmd += 'done\n'
#             cmd += 'for i in %s/*/*/prob_*.txt\n' % output
#             cmd += 'do\n'
#             cmd += '    mkdir -p "$i"\n'
#             cmd += 'done'
#             self.outputs['cmds'].setdefault(group, []).append(cmd)
#             self.outputs['dir'].append(o_dir)
#             self.outputs['outs'].setdefault(self.pool, []).append(o_dir)
#             io_update(self, i_f=fp, i_d=[ion_com, i_tasser, i_tasser_libs],
#                       o_d=o_dir, key=group)
#
#
