# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import sys
import pkg_resources
from os.path import basename, isdir, splitext

from metagenomix._inputs import (sample_inputs, group_inputs,
                                 genome_key, genome_out_dir)
from metagenomix._io_utils import caller, io_update, to_do, status_update
from metagenomix.core.parameters import tech_params

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def prodigal_cmd(
        self,
        tech: str,
        fasta_fp: str,
        out: str
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fasta_fp : str
        Path to the input fasta file
    out : str
        Path to the output folder

    Returns
    --------
    cmd : str
        Prodigal command
    """
    params = tech_params(self, tech)
    cmd, cmd_rm = '', ''

    contig_fa = fasta_fp
    if fasta_fp.endswith('.fa.gz') or fasta_fp.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta_fp, fasta_fp.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta_fp.rstrip('.gz')
        contig_fa = fasta_fp.rstrip('.gz')
    elif fasta_fp.endswith('fastq.gz'):
        contig_fa = fasta_fp.replace('fastq.gz', 'fasta')
        cmd += 'seqtk seq -a %s > %s\n' % (fasta_fp, contig_fa)

    cmd = 'prodigal'
    cmd += ' -i %s' % contig_fa
    cmd += ' -a %s/protein.translations.fasta' % out
    cmd += ' -d %s/nucleotide.sequences.fasta' % out
    cmd += ' -s %s/potential.starts.fasta' % out
    cmd += ' -o %s/gene.coords.%s' % (out, params['f'])
    for param in ['f', 'p']:
        cmd += ' -%s %s' % (param, params[param])
    for boolean in ['c', 'm', 'n', 'q']:
        if params[boolean]:
            cmd += ' -%s' % params[boolean]
    cmd += '\n' + cmd_rm
    cmd += 'for i in %s/*; do gzip $i; do\n' % out
    return cmd


def get_prodigal(
        self,
        tech: str,
        fastas: dict,
        out_dir: str,
        sam_group: str,
) -> None:
    """Get the prodigal command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta file per sample name / data source
    out_dir : str
        Paths to the output folder
    sam_group : str
        Sample name or group for the current co-assembly
    """
    for genome, fasta in fastas.items():

        to_dos = status_update(
            self, tech, [fasta[0]], group=sam_group, genome=genome)

        out = out_dir
        key = (tech, sam_group)
        if genome:
            out += '/' + genome
            key += (genome,)

        self.outputs['dirs'].append(out)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out)

        proteins = '%s/protein.translations.fasta.gz' % out
        if self.config.force or to_do(proteins):
            cmd = prodigal_cmd(self, tech, fasta[0], out)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta[0], o_d=out, key=key)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def prodigal(self) -> None:
    """Fast, reliable protein-coding gene prediction for prokaryotic genomes.

    References
    ----------
    Hyatt, Doug, et al. "Prodigal: prokaryotic gene recognition and
    translation initiation site identification." BMC bioinformatics 11.1 (
    2010): 1-11.

    Notes
    -----
    GitHub  : https://github.com/hyattpd/Prodigal
    Paper   : https://doi.org/10.1186/1471-2105-11-119

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
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
            get_prodigal(self, tech, fastas, out_dir, group)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            out_dir = '/'.join([self.dir, tech, mags])
            get_prodigal(self, tech, fastas, out_dir, mags)

    else:
        tech_fastas = sample_inputs(self, ['pacbio', 'nanopore'])
        for tech, fastas in tech_fastas.items():
            out_dir = '/'.join([self.dir, tech, self.sam_pool])
            get_prodigal(self, tech, fastas, out_dir, self.sam_pool)


def macsyfinder_cmd(
        self,
        tech: str,
        proteins: str,
        sam_group: str,
        genome: str,
        out_dir: str,
        models_dir: str
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
    tech : str
        Technology: 'illumina', 'pacbio', 'nanopore', or hybrif naming
    proteins : str
        Path to the input proteins fasta file
    sam_group : str
        Sample name or group for the current co-assembly
    genome : str
        MAGs/Genomes folder name or empty string (for assembly contigs)
    out_dir : str
        Path to the output folder
    models_dir : str
        Path to the models folder

    Returns
    -------
    cmd : str
        MacSyFinder command
    outs : dict
        Paths to the output folder per model
    """
    params = tech_params(self, tech)

    cmd_header, cmd_rm = '', ''
    if proteins.endswith('.fasta.gz'):
        cmd_header += 'gunzip -c %s > %s\n' % (proteins, proteins.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % proteins.rstrip('.gz')
        proteins = proteins.rstrip('.gz')

    proteins_fp = proteins
    if tech in ['illumina', 'pacbio', 'nanopore']:
        proteins_fp = '%s_edit.fasta' % proteins.replace('.fasta', '')
        cmd_header += '%s/header_space_replace.py -i %s -o %s --n\n' % (
            RESOURCES, proteins, proteins_fp)

    cmd = ''
    outs = {}
    for model in params['models']:
        model_dir = '%s/models' % models_dir
        if not self.config.dev and not isdir(model_dir):
            self.soft.add_status(tech, self.sam_pool, 1, group=sam_group,
                                 message='no model %s' % model, genome=genome)
            continue

        model_out_dir = '%s/%s' % (out_dir, model)
        self.outputs['dirs'].append(model_out_dir)
        outs[model] = model_out_dir

        res = '%s/macsyfinder.log' % model_out_dir
        if self.config.force or to_do(res):
            cmd += 'macsyfinder'
            if self.soft.prev in self.config.tools['assembling']:
                cmd += ' --db-type ordered_replicon'
                if self.soft.prev == 'spades_plasmid':
                    cmd += ' --replicon-topology circular'
                else:
                    cmd += ' --replicon-topology linear'
            elif '_drep_' in model_out_dir or '_metawrap' in model_out_dir:
                cmd += ' --db-type ordered_replicon'
                cmd += ' --replicon-topology circular'
            else:
                cmd += ' --db-type unordered'
            for param in [
                'e_value_search', 'i_evalue_sel', 'coverage_profile',
                'mandatory_weight', 'accessory_weight', 'exchangeable_weight',
                'redundancy_penalty', 'out_of_cluster'
            ]:
                cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
            cmd += ' --out-dir %s' % model_out_dir
            cmd += ' --res-search-suffix _hmm.tsv'
            cmd += ' --res-extract-suffix _out.tsv'
            cmd += ' --worker %s' % params['cpus']
            cmd += ' --sequence-db %s' % proteins_fp
            cmd += ' --models-dir %s' % model_dir
            cmd += ' --models %s all' % model
            cmd += ' --verbosity\n'
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)
    if cmd:
        cmd = cmd_header + cmd + cmd_rm
    return cmd, outs


def get_macsyfinder(
        self,
        tech: str,
        fastas: dict,
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample name or group for the current co-assembly
    """
    for genome, fasta in fastas.items():

        key = genome_key(tech, sam_group, genome)
        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        proteins = fasta[0]
        if self.soft.prev == 'prodigal':
            proteins = '%s/protein.translations.fasta.gz' % fasta[0]

        to_dos = status_update(
            self, tech, [proteins], group=sam_group, genome=genome)

        models_dir = self.databases.paths['macsyfinder']
        cmd, outs = macsyfinder_cmd(
            self, tech, proteins, sam_group, genome, out_dir, models_dir)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(outs)
        if cmd:
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=proteins, i_d=models_dir, o_d=out_dir, key=key)


def macsyfinder(self) -> None:
    """MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.

    References
    ----------
    Abby, Sophie S., et al. "MacSyFinder: a program to mine genomes for
    molecular systems with an application to CRISPR-Cas systems." PloS one
    9.10 (2014): e110726.

    Notes
    -----
    GitHub  : https://github.com/gem-pasteur/macsyfinder
    Docs    : https://macsyfinder.readthedocs.io/en/latest/
    Paper   : https://doi.org/10.1371/journal.pone.0110726

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
    if self.sam_pool in self.pools:
        if self.soft.prev != 'prodigal':
            sys.exit('[macsyfinder] Runs on protein data (prodigal only)')
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_macsyfinder(self, tech, fastas, group)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_macsyfinder(self, tech, fastas, mags)

    else:
        prev = self.soft.prev
        if prev not in ['plass', 'prodigal']:
            sys.exit('[macsyfinder] Runs on protein data (plass or prodigal)')
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_macsyfinder(self, tech, fastas, self.sam_pool)


def write_hmms_cmd(self) -> tuple:
    """Write a file to contain on each line, the paths to the .hmm
    files to search as part of integronfinder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to the output directory.
        .databases.hmms_dias : dict
            .hmm files per name of profile.

    Returns
    -------
    cmd = str
        Command to write the HMMs file
    hmms_fp : str
        Path to the file containing the paths to the .hmm files to search.
    """
    cmd = ''
    hmms_fp = ''
    if self.databases.hmms_dias:
        hmms_fp = '%s/hmms.txt' % self.dir
        first = True
        for i, (target, gen_hmm) in enumerate(self.databases.hmms_dias.items()):
            if self.config.dev and i == 1:
                break
            for j, (hmm, _) in enumerate(gen_hmm.values()):
                if self.config.dev and j == 1:
                    break
                if first:
                    cmd += 'echo -e "%s" > %s\n' % (hmm, hmms_fp)
                    first = False
                else:
                    cmd += 'echo -e "%s" >> %s\n' % (hmm, hmms_fp)
    return cmd, hmms_fp


def integronfinder_cmd(
        self,
        tech: str,
        fastx: str,
        out: str,
        hmms_fp: str
) -> str:
    """Collect command for integron_finder.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Path to the reference databases
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastx : str
        Path to the input fasta/fastq(.gz) file
    out : str
        Path to the output folder
    hmms_fp : str

    Returns
    -------
    cmd : str
        integron_finder command
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    fasta = fastx
    if fastx.endswith('.fa.gz') or fastx.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fastx, fastx.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fastx.rstrip('.gz')
        fasta = fastx.rstrip('.gz')
    if fastx.endswith('.fastq.gz') or fastx.endswith('.fastq'):
        fasta = '%s.fasta' % fastx.rsplit('.fastq', 1)[0]
        cmd += 'seqtk seq -a %s > %s\n' % (fastx, fasta)

    fasta_out = '%s_min%snt.fasta' % (splitext(fasta)[0], params['min_length'])
    cmd += '\n%s/filter_on_length.py' % RESOURCES
    cmd += ' -i %s' % fasta
    cmd += ' -o %s' % fasta_out
    cmd += ' -t %s\n' % params['min_length']

    cmd += 'integron_finder'
    for boolean in [
        'local_max', 'promoter_attI', 'mute', 'pdf', 'gbk', 'union_integrases'
    ]:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --verbose'
    cmd += ' --outdir %s' % out
    cmd += ' --cpu %s' % params['cpus']
    if self.databases.hmms_dias:
        cmd += ' --func-annot'
        cmd += ' --path-func-annot %s' % hmms_fp
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
        if self.soft.prev == 'plass':
            cmd += ' --topology linear'
        elif '_spades_' in out:
            cmd += ' --topology linear'
        elif '_drep_' in out or '_metawrap' in out:
            cmd += ' --topology circ'
    cmd += ' %s\n' % fasta_out
    cmd += cmd_rm
    return cmd


def get_integronfinder(
        self,
        tech: str,
        fastas: dict,
        sam_group: str,
) -> None:
    """Get the integron_finder command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Databases
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Sample name or group for the current co-assembly
    """
    hmms_cmds, hmms_fp = write_hmms_cmd(self)
    if hmms_cmds:
        hmms_sh = '%s/hmms.sh' % self.dir
        hmms_dir = self.dir.replace('${SCRATCH_FOLDER}', '')
        if not isdir(hmms_dir):
            os.makedirs(hmms_dir)
        with open(hmms_sh.replace('${SCRATCH_FOLDER}', ''), 'w') as o:
            o.write(hmms_cmds)
        hmms_cmd = '\nif [ ! -f %s ]; then sh %s; fi\n' % (hmms_fp, hmms_sh)

    for genome, fasta_ in fastas.items():

        key = genome_key(tech, sam_group, genome)
        out = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out)

        fasta = fasta_[0]
        if self.soft.prev == 'prodigal':
            fasta = '%s/nucleotide.sequences.fasta.gz' % fasta_[0]

        to_dos = status_update(
            self, tech, [fasta], group=sam_group, genome=genome)

        fpo = '%s/Results_Integron_Finder_mysequences/mysequences.summary' % out
        if self.config.force or to_do(fpo):
            cmd = hmms_cmd + integronfinder_cmd(self, tech, fasta, out, hmms_fp)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=[fasta, hmms_sh], o_d=out, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def integronfinder(self) -> None:
    """Finds integrons in DNA sequences.

    Integrons are major genetic element, notorious for their major
    implication in the spread of antibiotic resistance genes. More generally,
    integrons are gene-capturing platform, whose broader evolutionary role
    remains poorly understood. IntegronFinder is able to detect with high
    accuracy integron in DNA sequences. It is accurate because it combines
    the use of HMM profiles for the detection of the essential protein,
    the site-specific integron integrase, and the use of Covariance Models
    for the detection of the recombination site, the attC site.

    References
    ----------
    Cury, Jean, et al. "Identification and analysis of integrons and cassette
    arrays in bacterial genomes." Nucleic acids research 44.10 (2016):
    4539-4550.

    Notes
    -----
    GitHub  : https://github.com/gem-pasteur/Integron_Finder
    Docs    : https://integronfinder.readthedocs.io/en/latest/
    Paper   : https://doi.org/10.1093/nar/gkw319

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
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_integronfinder(self, tech, fastas, group)

    elif set(self.inputs) == {''}:
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_integronfinder(self, tech, fastas, mags)

    else:
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_integronfinder(self, tech, fastas, self.sam_pool)


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
#             out = '%s/%s.tsv' % (sam_dir, self.sam_pool)
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
#                     tmp = '%s/%s_tmp' % (sam_dir, self.sam_pool)
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


def diamond_cmd(
        params: dict,
        fp: str,
        db_path: str,
        out: str,
        tmp_dir: str
) -> str:
    """Collect the command line for DIAMOND.

    Notes
    -----
    GitHub  : https://github.com/bbuchfink/diamond
    Paper   : https://doi.org/10.1038/s41592-021-01101-x

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


def diamond(
        self,
        tech: str,
        sam_group: str,
        proteins: str,
        out_dir: str,
        key: tuple,
        to_dos: list
) -> dict:
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    proteins : str
        Path to the input proteins fasta file
    out_dir : str
        Path to the main output folder
    key : tuple
        Variables names for the current analytic level
    to_dos : list
        List of file that

    Returns
    -------
    outs : dict
        Paths to the output folder per database
    """
    outs = {}
    params = tech_params(self, tech)
    # tmp = '$TMPDIR/%s_%s_%s' % (self.soft.name, self.sam_pool, '_'.join(key))
    # for db, db_paths in params['databases'].items():
    #     for db_path in db_paths:
    #         db_base = splitext(basename(db_path))[0]
    #         db_out_dir = '%s/%s' % (out_dir, db_base)
    #         self.outputs['dirs'].append(db_out_dir)
    #
    #         out = '%s/%s.tsv' % (db_out_dir, db_base)
    #         outs.setdefault((db, db_base), []).append(db_out_dir)
    #
    #         if self.config.force or to_do(out):
    #             cmd = diamond_cmd(params, proteins, db_path, out, tmp)
    #             if to_dos:
    #                 self.outputs['cmds'].setdefault(key, []).append(False)
    #             else:
    #                 self.outputs['cmds'].setdefault(key, []).append(cmd)
    #             io_update(self, i_f=proteins, o_d=db_out_dir, key=key)
    #             self.soft.add_status(
    #                 tech, self.sam_pool, 1, group=sam_group, genome=db)
    #         else:
    #             self.soft.add_status(
    #                 tech, self.sam_pool, 0, group=sam_group, genome=db)
    return outs


def hmmer_cmd(
        params: dict,
        fp: str,
        db_path: str,
        out: str
) -> str:
    """Collect the command line for HMMER.

    Notes
    -----
    Docs    : http://hmmer.org/documentation.html
    Paper   : https://doi.org/10.1371/journal.pcbi.1002195

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


def hmmer(
        self,
        tech: str,
        sam_group: str,
        proteins: str,
        out_dir: str,
        key: tuple,
        to_dos: list
) -> dict:
    """Go through every sample / co-assembly group to collect the
    hmmer command lines for the different search terms or pfam models.

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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or group for the current co-assembly
    proteins : str
        Path to the input proteins fasta file
    out_dir : str
        Path to the main output folder
    key : tuple
        Variables names for the current analytic level

    Returns
    -------
    outs : dict
        Paths to the output folder per database
    """
    outs = {}
    params = tech_params(self, tech)
    # if params.get('descriptions'):
    #     print(self.databases.pfams)
    #     print(self.databases.hmms_dias.keys())
    #     print(selfdatabases)
    #     db_out_dir = '%s/%s' % (out_dir, db_base)
    #     self.outputs['dirs'].append(db_out_dir)
    #
    #     out = '%s/%s.tsv' % (db_out_dir, db_base)
    #     outs.setdefault((db, db_base), []).append(db_out_dir)
    #
    #     if self.config.force or to_do(out):
    #         cmd = hmmer_cmd(params, proteins, db_path, out)
    #         if to_dos:
    #             self.outputs['cmds'].setdefault(key, []).append(False)
    #         else:
    #             self.outputs['cmds'].setdefault(key, []).append(cmd)
    #         io_update(self, i_f=proteins, o_d=db_out_dir, key=key)
    #         self.soft.add_status(
    #             tech, self.sam_pool, 1, group=sam_group, genome=db)
    #     else:
    #         self.soft.add_status(
    #             tech, self.sam_pool, 0, group=sam_group, genome=db)
    return outs


def get_search(
        self,
        tech: str,
        fastas: dict,
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Paths to the input fasta file per sample name / data source
    sam_group : str
        Sample name or group for the current co-assembly
    """
    for genome, fasta in fastas.items():

        key = genome_key(tech, sam_group, genome)
        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        proteins = fasta[0]
        if self.soft.prev == 'prodigal':
            proteins = '%s/protein.translations.fasta' % fasta[0]

        to_dos = status_update(
            self, tech, [proteins], group=sam_group, genome=genome)

        search_ = caller(self, __name__)
        outs = search_(self, tech, sam_group, proteins, out_dir, key, to_dos)
        if outs:
            self.outputs['outs'].setdefault((tech, sam_group), {}).update(outs)


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
    name = self.soft.name
    if self.sam_pool in self.pools:
        if self.soft.prev != 'prodigal':
            sys.exit('[%s] Runs on protein data (plass, prodigal...)' % name)
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_search(self, tech, fastas, group)

    elif set(self.inputs) == {''}:
        if self.soft.prev != 'prodigal':
            sys.exit('[%s] Runs on protein data (plass, prodigal...)' % name)
        for (tech, mags), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_search(self, tech, fastas, mags)

    else:
        prev = self.soft.prev
        if prev not in ['plass', 'prodigal']:
            sys.exit('[%s] Runs on protein data (plass, prodigal...)' % name)
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_search(self, tech, fastas, self.sam_pool)


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
    cmd, cmd_rm = '', ''
    if contigs.endswith('.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % contigs.rstrip('.gz')

    cmd += 'prokka'
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
    if contigs.endswith('.gz'):
        cmd += ' %s\n' % contigs.rstrip('.gz')
    else:
        cmd += ' %s\n' % contigs

    cmd += 'gzip %s/%s.out\n' % (out_dir, pref)
    cmd += cmd_rm
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
        tech: str,
        group: str,
        inputs: list,
        out_dir: str
) -> tuple:
    """Get the Prokka commands for the different configs
    (at each taxonomic level).

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    group : str
        Name of the current co-assembly
    inputs : list
        Paths to the input files
    out_dir : str
        Path to the output folder

    Returns
    -------
    to_dos : list
    cmd : str
        Prokka command
    """
    cmd = ''
    configs, cols = prokka_configs(self)
    to_dos = []
    for config in configs:

        to_dos.extend(status_update(
            self, tech, [inputs[0]], group=group, genome=list(config)[0]))

        pref = '_'.join([config[x] for x in cols if config[x]])
        if config.get('proteins'):
            pref += '_%s' % splitext(basename(config['proteins']))[0]
        file_out = '%s/%s.out.gz' % (out_dir, pref)
        if self.config.force or to_do(file_out):
            cmd += prokka_cmd(self, inputs[0], out_dir, pref, config, cols)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=list(config)[0])
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=list(config)[0])
    return to_dos, cmd


def prokka(self) -> None:
    """Prokka: rapid prokaryotic genome annotation.
    Whole genome annotation is the process of identifying features of
    interest in a set of genomic DNA sequences, and labelling them with
    useful information. Prokka is a software tool to annotate bacterial,
    archaeal and viral genomes quickly and produce standards-compliant
    output files,

    References
    ----------
    Seemann, Torsten. "Prokka: rapid prokaryotic genome annotation."
    Bioinformatics 30.14 (2014): 2068-2069.

    Notes
    -----
    GitHub  : https://github.com/tseemann/prokka
    Paper   : https://doi.org/10.1093/bioinformatics/btu153

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
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'].setdefault(group, []).append(out_dir)

            to_dos, cmd = get_prokka(self, tech, group, inputs, out_dir)
            if cmd:
                key = (tech, group)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=inputs[0], o_d=out_dir, key=key)


def barrnap_cmd(
        self,
        tech: str,
        fastx: str,
        out: str
) -> str:
    """Collect the command line for Barrnap.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastx : str
        Path to the input fasta file
    out : str
        Path to the output file

    Returns
    -------
    cmd : str
        Barrnap command
    """
    params = tech_params(self, tech)

    fasta = fastx
    cmd, cmd_rm = '', ''
    if fastx.endswith('.fa.gz') or fastx.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fastx, fastx.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fastx.rstrip('.gz')
        fasta = fastx.rstrip('.gz')
    elif fastx.endswith('.fastq.gz') or fastx.endswith('.fastq'):
        fasta = '%s.fasta' % fastx.rsplit('.fastq', 1)[0]
        cmd += 'seqtk seq -a %s > %s\n' % (fastx, fasta)

    cmd += 'barrnap'
    cmd += ' --kingdom %s' % params['kingdom']
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --reject %s' % params['reject']
    cmd += ' --lencutoff %s' % params['lencutoff']
    cmd += ' --evalue %s' % params['evalue']
    if params['incseq']:
        cmd += ' --incseq'
    cmd += ' --outseq %s' % out
    cmd += ' %s\n' % fasta
    cmd += cmd_rm
    return cmd


def ccfind_cmd(
        self,
        tech: str,
        fastx: str,
        out_dir: str
) -> str:
    """Collect the command line for ccfind.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastx : str
        Path to the input fasta file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        ccfind command
    """
    params = tech_params(self, tech)
    cmd, cmd_rm = '', ''
    if self.config.force:
        cmd += 'if [ -d "%s" ]; then rm -rf %s; fi\n' % (out_dir, out_dir)

    fasta = fastx
    cmd, cmd_rm = '', ''
    if fastx.endswith('.fa.gz') or fastx.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fastx, fastx.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fastx.rstrip('.gz')
        fasta = fastx.rstrip('.gz')
    elif fastx.endswith('.fastq.gz') or fastx.endswith('.fastq'):
        fasta = '%s.fasta' % fastx.rsplit('.fastq', 1)[0]
        cmd += 'seqtk seq -a %s > %s\n' % (fastx, fasta)

    cmd += 'ccfind %s %s' % (fasta, out_dir)
    cmd += ' --terminal-fragment-size %s' % params['terminal_fragment_size']
    cmd += ' --min-percent-identity %s' % params['min_percent_identity']
    cmd += ' --min-aligned-length %s' % params['min_aligned_length']
    cmd += ' --ncpus %s' % params['cpus']
    if params['preserve_tmpdir']:
        cmd += ' --preserve-tmpdir'
    cmd += '\n' + cmd_rm
    return cmd


def get_ccfind(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """Get the ccfind command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for ccfind
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Path(s) to the input fasta file(s) per genome/MAGs
    sam_group : str
        Sample name or name of a co-assembly pool's group
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)

        to_dos = status_update(
            self, tech, [fasta[0]], group=sam_group, genome=genome)

        out = '%s/circ.detected.list' % out_dir
        if self.config.force or not to_do(out):
            cmd = ccfind_cmd(self, tech, fasta[0], out_dir)
            key = genome_key(tech, sam_group, genome)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta[0], o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def get_barrnap(
        self,
        tech: str,
        fastas: dict,
        sam_group: str
) -> None:
    """Get the Barrnap command and fill the pipeline data structures.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Barrnap
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas : dict
        Path(s) to the input fasta file(s) per genome/MAGs
    sam_group : str
        Sample name or name of a co-assembly pool's group
    """
    for genome, fasta in fastas.items():

        key = genome_key(tech, sam_group, genome)

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        out = '%s/rrna_seqs.fas' % out_dir
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out)

        contigs = fasta[0]
        to_dos = status_update(
            self, tech, [contigs], group=sam_group, genome=genome)

        if self.config.force or not to_do(out):
            cmd = barrnap_cmd(self, tech, contigs, out)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def dispatch(self) -> None:
    """Run either barrnap or ccfind on the per-sample fastq files,
    or on the fasta resulting from the assembly, on per co-assembly MAGs for
    the genome dereplication step.

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
    g_barrnap_ccfind = getattr(sys.modules[__name__], 'get_%s' % self.soft.name)

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            g_barrnap_ccfind(self, tech, fastas, group)
    else:
        tech_fastas = sample_inputs(self, ['pacbio', 'nanopore'])
        for tech, fastas in tech_fastas.items():
            g_barrnap_ccfind(self, tech, fastas, self.sam_pool)


def barrnap(self) -> None:
    """Dispatch the collection of barrnap commands using generic function.

    Barrnap predicts the location of ribosomal RNA genes in genomes. It
    supports bacteria (5S,23S,16S), archaea (5S,5.8S,23S,16S), metazoan
    mitochondria (12S,16S) and eukaryotes (5S,5.8S,28S,18S).
    It takes FASTA DNA sequence as input, and write GFF3 as output. It uses
    the new nhmmer tool that comes with HMMER 3.1 for HMM searching in
    RNA:DNA style. Multithreading is supported and one can expect roughly
    linear speed-ups with more CPUs.

    References
    ----------
    Seemann T
    barrnap 0.9 : rapid ribosomal RNA prediction
    https://github.com/tseemann/barrnap

    Notes
    -----
    GitHub  : https://github.com/tseemann/barrnap

    Parameters
    ----------
    self : Commands class instance
    """
    dispatch(self)


def ccfind(self) -> None:
    """Dispatch the collection of ccfind commands using generic function.

    ccfind - Circular Complete genome FINDer.
    ccfind is a general tool to detect circular complete genomes with clues
    of terminal redundancy, originally designed for identification
    of complete virus genomes from metagenome assembly. It can be used for
    any contig/genome, but cutoff values should be carefully considered.
    It should be noted that terminal redundancy (circularity) does not
    necessarily mean completion of the sequence. Partial genomes might be
    detected as circular contigs for some reasons (e.g., sequence repeats).

    References
    ----------
    Nishimura, Yosuke, et al. "Environmental viral genomes shed new light on
    virus-host interactions in the ocean." MSphere 2.2 (2017): e00359-16.

    Notes
    -----
    GitHub  : https://github.com/yosuken/ccfind
    Paper   : https://doi.org/10.1128/mSphere.00359-16

    Parameters
    ----------
    self : Commands class instance
    """
    dispatch(self)


def antismash_cmd(
        self,
        fasta: str,
        out_dir: str
) -> str:
    """Collect the command line for Antismash.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
        .soft.prev : str
            Previous software
    fasta : str
        Path to the input genome fasta file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Antismash command
    """
    cmd = '\nantismash'
    for boolean in ['rre', 'asf', 'cassis', 'pfam2go', 'tigrfam', 'fullhmmer',
                    'cc_mibig', 'cb_general', 'smcog_trees', 'clusterhmmer',
                    'cb_subclusters', 'cb_knownclusters']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --tta-threshold %s' % self.soft.params['tta_threshold']
    cmd += ' --genefinding-tool %s' % self.soft.params['genefinding_tool']
    cmd += ' --html-description after_%s' % self.soft.prev
    cmd += ' --taxon %s' % self.soft.params['taxon']
    cmd += ' --cpus %s' % self.soft.params['cpus']
    cmd += ' --output-basename %s' % splitext(basename(fasta))[0]
    cmd += ' --html-title "%s"' % splitext(basename(fasta))[0]
    cmd += ' --output-dir %s' % out_dir
    cmd += ' %s\n' % fasta
    return cmd


def get_antismash(
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
        Sample name or name of a co-assembly pool's group
    """
    for genome, dirs in fastas.items():

        fasta = dirs[0]
        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)

        key = genome_key(tech, sam_group, genome)

        to_dos = status_update(
            self, tech, [fasta], group=sam_group, genome=genome)

        if self.config.force or not glob.glob('%s/*' % out_dir):
            cmd = antismash_cmd(self, fasta, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, o_d=out_dir, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def antismash(self) -> None:
    """Antibiotics and Secondary Metabolite Analysis SHell.
    antiSMASH allows the rapid genome-wide identification, annotation and
    analysis of secondary metabolite biosynthesis gene clusters in bacterial
    and fungal genomes. It integrates and cross-links with a large number of
    in silico secondary metabolite analysis softwares that have been published
    earlier.

    References
    ----------
    Blin, Kai, et al. "antiSMASH 6.0: improving cluster detection and
    comparison capabilities." Nucleic Acids Research 49.W1 (2021): W29-W35.

    Notes
    -----
    GitHub  : https://github.com/antismash
    Paper   : https://doi.org/10.1093/nar/gkab335

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
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_antismash(self, fastas, tech, group)

    elif self.soft.prev == 'drep':
        for (tech, bin_algo), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_antismash(self, fastas, tech, bin_algo)
    else:
        tech_fastas = sample_inputs(self, ['pacbio', 'nanopore'])
        for tech, fastas in tech_fastas.items():
            get_antismash(self, fastas, tech, self.sam_pool)


def tiara_cmd(self, contig: str, out_dir: str, out_fp: str):
    """Get the tiara command line.

    Parameters
    ----------
    self
    contig : str
        Path to the input contigs file
    out_dir : str
        Path to the output folder
    out_fp : str
        Path to the output file

    Returns
    -------
    cmd : str
        tiara command line
    """
    cmd, cmd_rm = '', ''
    if contig.endswith('.fa.gz') or contig.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contig, contig.rstrip('.gz'))
        contig = contig.rstrip('.gz')
        cmd_rm += 'rm %s\n' % contig

    cmd += 'cd %s\n' % out_dir

    cmd += 'tiara'
    cmd += ' --input %s' % contig
    cmd += ' --output %s' % out_fp
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['probabilities', 'verbose', 'gzip']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    for param in ['min_len', 'first_stage_kmer', 'second_stage_kmer']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    for param in ['prob_cutoff', 'to_fasta']:
        cmd += ' --%s %s' % (param, ' '.join(self.soft.params[param]))
    cmd += '\n' + cmd_rm
    cmd += 'for i in %s/*; gzip $i; done\n' % out_dir
    return cmd


def tiara(self) -> None:
    """Deep-learning-based approach for identification of eukaryotic
    sequences in the metagenomic data powered by PyTorch. The sequences are
    classified in two stages:
    - In the first stage, the sequences are classified to classes: archaea,
      bacteria, prokarya, eukarya, organelle and unknown.
    - In the second stage, the sequences labeled as organelle in the first
      stage are classified to either mitochondria, plastid or unknown.

    Notes
    -----
    GitHub  : https://github.com/ibe-uw/tiara
    Paper   : https://doi.org/10.1093/bioinformatics/btab672

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline
        .dir : str
            Path to pipeline output folder for tiara
        .sam_pool : str
            Co-assembly name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.config.tools[self.soft.prev] != 'assembling':
        sys.exit('[tiara] can only be run on assembly output')
    for (tech, group), contigs in self.inputs[self.sam_pool].items():

        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir
        out_fp = '%s/classifications.txt.gz' % out_dir

        contig = contigs[0]
        to_dos = status_update(self, tech, [contig], group=group)
        if self.config.force or to_do(out_fp):
            cmd = tiara_cmd(self, contig, out_dir, out_fp)
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=contig, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def metaclade2(self):
    """Novel profile-based domain annotation pipeline based on the multi-source
    domain annotation strategy. It provides a domain annotation realised
    directly from reads, and reaches an improved identification of the
    catalog of functions in a microbiome. MetaCLADE can be applied to either
    metagenomic or metatranscriptomic datasets.

    References
    ----------
    Ugarte, A., Vicedomini, R., Bernardes, J. and Carbone, A., 2018. A
    multi-source domain annotation pipeline for quantitative metagenomic and
    metatranscriptomic functional profiling. Microbiome, 6(1), pp.1-27.

    Notes
    -----
    GitLab  : http://gitlab.lcqb.upmc.fr/vicedomini/metaclade2
    Paper   : https://doi.org/10.1186/s40168-018-0532-2
    Docs    : http://www.lcqb.upmc.fr/metaclade/

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def deeptmhmm(self):
    """A Deep Learning Model for Transmembrane Topology Prediction and
    Classification

    Protein structure prediction using deep learning methods have seen
    several advancements within the last years. In this project,
    we investigate deep learning for prediction of the membrane topology of
    transmembrane proteins. Transmembrane proteins are relevant for drug
    development since they make up more than 50% of all human drug targets.

    DeepTMHMM is currently the most complete and best-performing method for
    the prediction of the topology of both alpha-helical and beta-barrel
    transmembrane proteins. The model encodes the primary amino acid sequence
    by a pre-trained language model and decodes the topology by a state
    space model to produce topology and type predictions at unprecedented
    accuracy. DeepTMHMM makes it possible to scan full proteomes in order to
    detect both classes of transmembrane proteins, and we anticipate our
    method to be very valuable for the research community.

    References
    ----------
    Hallgren, J., Tsirigos, K.D., Pedersen, M.D., Armenteros, J.J.A.,
    Marcatili, P., Nielsen, H., Krogh, A. and Winther, O., 2022. DeepTMHMM
    predicts alpha and beta transmembrane proteins using deep neural
    networks. bioRxiv.

    Notes
    -----
    Paper   : https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1
    Demo    : https://colab.research.google.com/drive/1EijAjVAcNsgyI4nX-coB4UkoMzsG5SlI?usp=sharing
    Docs    : https://dtu.biolib.com/DeepTMHMM

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def itassermtd(self):
    """I-TASSER-MTD is a hierarchical protocol to predict structures and
    functions of multi-domain (MTD) proteins. It first predicts the domain
    boundaries by FUpred and ThreaDom based on the deep-learning contact-map
    prediction and multiple threading alignments. Next, the structure model
    of each individual domain is constructed independently by I-TASSER guided
    by the deep learning predicted spatial restraints. Finally,
    the individual domain models are assembled into full-length structure by
    DEMO under guidance of quaternary structural templates and deep-learning
    distance profiles. Meanwhile, the protein functions at both domain level
    and full-chain level are annotated by COFACTOR based on structures,
    sequences, and protein-protein interaction networks.

    References
    ----------
    Zhou, X., Zheng, W., Li, Y., Pearce, R., Zhang, C., Bell, E.W., Zhang,
    G. and Zhang, Y., 2022. I-TASSER-MTD: a deep-learning-based platform for
    multi-domain protein structure and function prediction. Nature Protocols,
    17(10), pp.2326-2353.

    Notes
    -----
    Paper   : https://doi.org/10.1038/s41596-022-00728-0
    Docs    : https://zhanggroup.org/I-TASSER-MTD/about.html

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def itasser(self):
    """I-TASSER Suite is a package of standalone computer programs, developed
    for high-resolution protein structure prediction, refinement, and
    structure-based function annotations. A detailed instruction on how to
    download and install the Suite can be found at README5.2.txt. Please
    report bugs and questions at I-TASSER message board and some members will
    study the problems and answer them asap. The I-TASSER Suite is free for
    academic and non-profit researchers.

    References
    ----------
    Yang, J., Yan, R., Roy, A., Xu, D., Poisson, J. and Zhang, Y., 2015. The
    I-TASSER Suite: protein structure and function prediction. Nature
    methods, 12(1), pp.7-8.

    Notes
    -----
    Paper   : https://doi.org/10.1038/nmeth.3213
    Docs    : https://zhanggroup.org/I-TASSER/

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def graspx(self):
    """Homology detection algorithm named GRASP (Guided Reference-based
    Assembly of Short Peptides) that identifies the homologs of a given
    reference protein sequence in a database of short peptide metagenomic
    sequences. GRASP was developed to implement a simultaneous alignment and
    assembly algorithm for annotation of short peptides identified on
    metagenomic reads.

    References
    ----------
    Zhong, C., Yang, Y. and Yooseph, S., 2016. GRASPx: efficient
    homolog-search of short peptide metagenome database through simultaneous
    alignment and assembly. BMC bioinformatics, 17(8), pp.611-621.

    Notes
    -----
    SourceForge : https://sourceforge.net/projects/graspx/
    Paper       : https://doi.org/10.1186/s12859-016-1119-1
    Docs        : https://cbb.ittc.ku.edu/GRASPx.html

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def rundbcan(self):
    """Run_dbcan V3, using genomes/metagenomes/proteomes of any assembled
    organisms (prokaryotes, fungi, plants, animals, viruses) to search for
    CAZymes.

    References
    ----------
    Huang, L., Zhang, H., Wu, P., Entwistle, S., Li, X., Yohe, T., Yi, H.,
    Yang, Z. and Yin, Y., 2018. dbCAN-seq: a database of carbohydrate-active
    enzyme (CAZyme) sequence and annotation. Nucleic Acids Research, 46(D1),
    pp.D516-D521.

    Notes
    -----
    GitHub  : https://github.com/linnabrown/run_dbcan
    Paper   : https://doi.org/10.1093/nar/gkx894
    Docs    : bcb.unl.edu/dbCAN2

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def ioncom(self):
    """IonCom is a ligand-specific method for small ligand (including metal
    and acid radical ions) binding site prediction. Starting from given
    sequences or structures of the query proteins, IonCom performs a
    composite binding-site prediction that combines ab initio training and
    template-based transferals. To enhance specificity and sensitivity,
    the server focuses on binding site prediction of thirteen most important
    small ligand molecules, including nine metal ions (Zn2+, Cu2+, Fe2+,
    Fe3+, Ca2+, Mg2+, Mn2+, Na+, K+) and four acid radical ions (CO32-, NO2-,
    SO42-, PO43-).

    References
    ----------
    Hu, X., Dong, Q., Yang, J. and Zhang, Y., 2016. Recognizing metal and
    acid radical ion-binding sites by integrating ab initio modeling with
    template-based transferals. Bioinformatics, 32(21), pp.3260-3269.

    Notes
    -----
    Paper   : https://doi.org/10.1093/bioinformatics/btw396
    Docs    : http://zhanglab.ccmb.med.umich.edu/IonCom

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def gmove(self):
    """Gene modelling using various evidence.

    References
    ----------


    Notes
    -----
    GitHub  : https://github.com/institut-de-genomique/Gmove
    Paper   :
    Docs    :

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def gmove(self):
    """Gene modelling using various evidence.

    References
    ----------


    Notes
    -----
    GitHub  : https://github.com/institut-de-genomique/Gmove
    Paper   :
    Docs    :

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def srst2(self):
    """

    References
    ----------


    Notes
    -----
    GitHub  : https://github.com/katholt/srst2
    Paper   :
    Docs    :

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def pirate(self):
    """

    References
    ----------


    Notes
    -----
    GitHub  : https://github.com/tseemann/PIRATE
    Paper   :
    Docs    :

    Parameters
    ----------
    self : Commands class instance
    """
    pass


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
#     if self.sam_pool in self.pools:
#         for group in self.pools[self.sam_pool]:
#             o_dir, fp = get_out_dir(self, self.sam_pool, group)
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
#             self.outputs['outs'].setdefault(self.sam_pool, []).append(o_dir)
#             io_update(self, i_f=fp, i_d=[ion_com, i_tasser, i_tasser_libs],
#                       o_d=o_dir, key=group)
#
#


# TO ADD:  bakta, genbank, microscope, patric, pgap, prokka, pseudomonasdb, rast