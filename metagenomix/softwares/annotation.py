# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
import pkg_resources
from os.path import basename, dirname, isdir, splitext

from metagenomix._inputs import (
    sample_inputs, group_inputs, genome_key, genome_out_dir, get_reads,
    get_group_reads, get_contigs, get_contigs_from_path, get_plasmids, get_args)
from metagenomix._io_utils import caller, io_update, to_do, status_update
from metagenomix.core.parameters import tech_params

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def prodigal_cmd(
        fasta_fp: str,
        out: str,
        params: dict
) -> str:
    """Create command lines for Prodigal.

    Parameters
    ----------
    fasta_fp : str
        Path to the input fasta file
    out : str
        Path to the output folder
    params : dict
        Paramters

    Returns
    --------
    cmd : str
        Prodigal command
    """
    cmd, cmd_rm = '', ''

    contig_fa = fasta_fp
    if fasta_fp.endswith('.fa.gz') or fasta_fp.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta_fp, fasta_fp.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta_fp.rstrip('.gz')
        contig_fa = fasta_fp.rstrip('.gz')
    elif fasta_fp.endswith('fastq.gz'):
        contig_fa = fasta_fp.replace('fastq.gz', 'fasta')
        cmd += 'seqtk seq -a %s > %s\n' % (fasta_fp, contig_fa)

    cmd += 'prodigal'
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
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out
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
    params = tech_params(self, tech)
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

        coords = '%s/gene.coords.%s.gz' % (out, params['f'])
        proteins = '%s/protein.translations.fasta.gz' % out
        if self.config.force or to_do(proteins) or to_do(coords):
            cmd = prodigal_cmd(fasta[0], out, params)
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
        cmd_header += 'python %s/header_space_replace.py -i %s -o %s --n\n' % (
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

        res = '%s/all_systems.tsv.gz' % model_out_dir
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
        cmd += 'tar cpf %s/hmmer_results.tar -C %s hmmer_results\n' % (
               model_out_dir, model_out_dir)
        cmd += 'rm %s/hmmer_results\n' % model_out_dir
        cmd += 'for i in %s/*; do gzip -q $i; done\n' % model_out_dir
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
        out_dir: str,
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
    out_dir : str
        Path to the output folder

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

    fasta_out = '%s/out.fasta' % dirname(fasta)
    cmd += '\n%s/filter_on_length.py' % RESOURCES
    cmd += ' -i %s' % fasta
    cmd += ' -o %s' % fasta_out
    cmd += ' -t %s\n' % params['min_length']

    if params.get('path'):
        cmd += '%s/parallel_integron_finder.nf' % params['path']
        cmd += ' --distance-threshold %s' % params['distance_threshold']
    else:
        cmd += 'integron_finder'
        cmd += ' --distance-thresh %s' % params['distance_threshold']
    for boolean in [
        'pdf', 'gbk', 'keep_tmp', 'mute', 'gembase', 'local_max',
        'promoter_attI', 'union_integrases', 'keep_palindromes',
        'circ', 'linear'
    ]:
        if params[boolean]:
            if params.get('path'):
                if boolean in ['gembase', 'mute']:
                    continue
            else:
                cmd += ' --%s' % boolean.replace('_', '-')

    cmd += ' --verbose'
    cmd += ' --outdir %s' % out_dir
    cmd += ' --cpu %s' % params['cpus']

    if params['no_proteins']:
        cmd += ' --no-proteins'
    elif params['func_annot']:
        cmd += ' --func-annot'
        if params['path_func_annot']:
            cmd += ' --path-func-annot %s' % params['path_func_annot']

    if params['attc_model']:
        cmd += ' --attc-model %s' % params['attc_model']
        cmd += ' --evalue-attc %s' % params['evalue_attc']
        cmd += ' --max-attc-size %s' % params['max_attc_size']
        cmd += ' --min-attc-size %s' % params['min_attc_size']

    if params['linear']:
        cmd += ' --linear'
    if params['circ']:
        cmd += ' --circ'
    if params['topology_file']:
        # this file is automatically created for contigs annotates as plasmids
        cmd += ' --topology-file %s' % params['topology_file']

    if params.get('path'):
        cmd += ' --replicons %s\n' % fasta_out
    else:
        cmd += ' %s\n' % fasta_out

    out = 'Results_Integron_Finder_out'
    cmd += 'tar cpfz %s/%s.tar.gz -C %s %s\n' % (out_dir, out, out_dir, out)
    cmd += 'rm -rf %s/%s\n' % (out_dir, out)

    cmd += cmd_rm
    return cmd


def get_integronfinder(
        self,
        tech: str,
        group: str,
        inputs: dict,
        contigs: str
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
    group : str
        Co-assembly group / Sample name
    inputs : dict
        Paths to the input fasta files per genome/MAG
    contigs : str
        Empty of not run after a plasmid-detection tool
    """
    for genome, fastas in inputs.items():

        key = genome_key(tech, group, genome)
        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        if contigs:
            if self.soft.prev == 'hamronization':
                selection, fasta, cmds, rms = get_args(fastas, contigs)
                i_f = [selection, contigs]
            else:
                plasmids, fasta, cmds, rms = get_plasmids(self, fastas, contigs)
                i_f = [plasmids, contigs]
        else:
            fasta, cmds, rms = fastas[0], '', ''
            i_f = [fasta]
        to_dos = status_update(self, tech, i_f, group=group, genome=genome)

        fpo = '%s/Results_Integron_Finder_out.tar.gz' % out_dir
        if self.config.force or to_do(fpo):
            # cmd = hmms_cmd + integronfinder_cmd(self, tech, fasta, out)
            cmds += integronfinder_cmd(self, tech, fasta, out_dir)
            cmds += rms
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmds)
            # io_update(self, i_f=[fasta, hmms_sh], o_d=out, key=key)
            io_update(self, i_f=i_f, o_d=out, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=genome)


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
    invalid = True
    previous = self.config.tools[self.soft.prev]
    if self.soft.prev in self.config.tools['assembling']:
        invalid = False
    elif previous in ['MAG (dereplication)', 'annotation (plasmid)']:
        invalid = False
    if invalid:
        sys.exit('[integronfinder] Only on assembly, MAGs, plasmid annotations')

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            contigs = ''
            if previous == 'annotation (plasmid)':
                contigs = get_contigs_from_path(self, tech, group)
            elif self.soft.prev == 'hamronization':
                contigs = get_contigs_from_path(self, tech, group)
            get_integronfinder(self, tech, group, fastas, contigs)


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
            proteins = '%s/protein.translations.fasta.gz' % fasta[0]

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
        prefix: str,
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
    prefix : str
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
    for boolean in [
        'metagenome',
        'notrna',
        'norrna',
        'addgenes',
        'addmrna',
        'compliant',
        'rawproduct',
        'cdsrnaolap',
        'fast',
        'noanno',
        'rnammer',
    ]:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --force'
    for col in cols:
        if col and config[col]:
            cmd += ' --%s %s' % (col, config[col])
    cmd += ' --prefix %s' % prefix
    cmd += ' --outdir %s' % out_dir
    for param in ['centre', 'gram']:
        if self.soft.params[param]:
            cmd += ' --%s %s' % (param, self.soft.params[param])
    for param in ['gffver', 'accver', 'rfam']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    if contigs.endswith('.gz'):
        cmd += ' %s\n' % contigs.rstrip('.gz')
    else:
        cmd += ' %s\n' % contigs
    cmd += 'for i in %s/%s*; do gzip -q $i; done\n' % (out_dir, prefix)
    cmd += cmd_rm
    return cmd


def prokka_config(
        self, configs: list
) -> list:
    """Parse the user-defined taxonomic config file to collect the taxonomic
    levels to annotate more precisely using Prokka.

    Parameters
    ----------
    self : Commands class instance
        .soft.params : dict
            Parameters
    configs : list
        List of field#:taxon dicts obtained by parsing the user-defined config

    Returns
    -------
    cols : list
        Columns of the taxonomic config: 'genus', 'species', 'strain', 'plasmid'
    """
    cols = ['genus', 'species', 'strain', 'plasmid', 'proteins']
    with open(self.soft.params['config']) as f:
        for ldx, line in enumerate(f):
            ls = line.strip('\n').split('\t')
            if not ldx:
                continue
            configs.append({cols[idx]: val for idx, val in enumerate(ls)})
    return cols


def get_prokka_configs(self) -> tuple:
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
    configs = []
    if self.soft.params['config']:
        cols = prokka_config(self, configs)
    if not configs:
        cols = ['']
        configs = [{'': 'no_config'}]
    return configs, cols


def prokka_configs(
        self,
        tech: str,
        group: str,
        fasta: list,
        configs: list,
        cols: list,
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    group : str
        Name of the current co-assembly
    fasta : list
        Path to one input genome/MAG fasta file
    configs : list
        Dicts of config files per taxonomic level
    cols : list
        Taxonomic levels
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Prokka command
    """
    cmd = ''
    for config in configs:

        prefix = '_'.join([config[x] for x in cols if config[x]])
        if config.get('proteins'):
            prefix += '_%s' % splitext(basename(config['proteins']))[0]
        out = '%s/%s.txt.gz' % (out_dir, prefix)
        if self.config.force or to_do(out):
            cmd += prokka_cmd(self, fasta[0], out_dir, prefix, config, cols)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=group, genome=list(config)[0])
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=group, genome=list(config)[0])
    return cmd


def get_prokka(
        self,
        tech: str,
        group: str,
        fastas: dict,
        configs: list,
        cols: list
) -> None:
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
    fastas : dict
        Paths to the input files
    configs : list
        Dicts of config files per taxonomic level
    cols : list
        Taxonomic levels
    """
    for genome, fas in fastas.items():
        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        to_dos = status_update(self, tech, fas[:1], group=group, genome=genome)
        cmd = prokka_configs(self, tech, group, fas, configs, cols, out_dir)
        if cmd:
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fas[0], o_d=out_dir, key=key)


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
        configs, cols = get_prokka_configs(self)
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_prokka(self, tech, group, fastas, configs, cols)


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

    cmd += 'rm -rf %s\n' % out_dir
    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s\n' % self.soft.params['path']
    cmd += 'ccfind'
    cmd += ' %s' % fasta
    cmd += ' %s' % out_dir
    cmd += ' --terminal-fragment-size %s' % params['terminal_fragment_size']
    cmd += ' --min-percent-identity %s' % params['min_percent_identity']
    cmd += ' --min-aligned-length %s' % params['min_aligned_length']
    if params['preserve_tmpdir']:
        cmd += ' --preserve-tmpdir'
    cmd += ' --ncpus %s\n' % params['cpus']
    cmd += cmd_rm
    cmd += 'for i in %s/result/*.list; do gzip -q $i; done\n' % out_dir
    cmd += 'for i in %s/result/*.fasta; do gzip -q $i; done\n' % out_dir
    cmd += 'for i in %s/result/intermediate/*; do gzip -q $i; done\n' % out_dir
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

        out = '%s/result/circ.detected.list.gz' % out_dir
        if self.config.force or to_do(out):
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

        if self.config.force or to_do(out):
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
    get_ = getattr(sys.modules[__name__], 'get_%s' % self.soft.name)
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_(self, tech, fastas, group)
    else:
        tech_fastas = sample_inputs(self, ['pacbio', 'nanopore'])
        for tech, fastas in tech_fastas.items():
            get_(self, tech, fastas, self.sam_pool)


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
    cmd += 'tar cpfz %s.tar.gz -C %s %s\n' % (
        out_dir, dirname(out_dir), basename(out_dir))
    cmd += 'rm -rf %s\n' % out_dir
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
        out_fp = '%s.tar.gz' % out_dir
        if self.config.force or to_do(out_fp):
            cmd = antismash_cmd(self, fasta, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=fasta, o_d=out_dir, o_f=out_fp, key=key)
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
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['probabilities', 'verbose', 'gzip']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    for param in ['min_len', 'first_stage_kmer', 'second_stage_kmer']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    for param in ['prob_cutoff', 'to_fasta']:
        cmd += ' --%s %s' % (param, ' '.join(self.soft.params[param]))
    cmd += ' --output %s\n' % out_fp
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    cmd += cmd_rm
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

        contig = contigs[0]
        to_dos = status_update(self, tech, [contig], group=group)

        out_fp = '%s/classifications.txt' % out_dir
        if self.config.force or to_do('%s.gz' % out_fp):
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


def concat_group_reads(
        group: str,
        reads: dict
) -> str:
    """

    Parameters
    ----------
    group : str
        Name of the current co-assembly group
    reads : dict
        Path(s) to reads per sample(s) of the current co-assembly group

    Returns
    -------
    cmd : str
        gzip and concat of reads fastq file(s)
    """
    cmd = ''
    fqs = {0: [], 1: []}
    for sam, tech_files in reads.items():
        for (tech, _), fastqs in tech_files.items():
            if len(fastqs) < 2:
                continue
            for r in [0, 1]:
                fastq = fastqs[r].replace('${SCRATCH_FOLDER}', '')
                if fastq.endswith('.gz'):
                    fqs[r].append(fastq)
                else:
                    fastq_gz = 'tmp/%s.gz' % basename(fastq)
                    cmd += 'gzip -c %s > %s\n' % (fastq, fastq_gz)
                    fqs[r].append(fastq_gz)

    cmd += 'cat %s > reads/%s_1.fastq.gz\n' % (' '.join(fqs[0]), group)
    cmd += 'cat %s > reads/%s_2.fastq.gz\n' % (' '.join(fqs[1]), group)
    return cmd


def diting_cmd(
        self,
        tech: str,
        contigs: dict,
        all_reads: dict,
        out: str,
        pool: str = None
) -> str:
    """Get command line for DiTing.

    Parameters
    ----------
    self
    tech : str
    contigs : dict
    all_reads : dict
    out : str
    pool : str
        Name of the co-assembly

    Returns
    -------
    cmd : str
       DiTing command
    """
    params = tech_params(self, tech)

    cmd = 'cd %s\n' % out
    cmd += 'mkdir tmp\n'
    cmd += 'mkdir reads\n'
    cmd += 'mkdir contigs\n'
    for group, contig in contigs.items():
        reads = get_group_reads(self, tech, group, all_reads, pool)
        cmd += concat_group_reads(group, reads)
        if contig.endswith('.fa.gz') or contig.endswith('.fasta.gz'):
            cmd += 'gunzip -c %s > contigs/%s.fa\n' % (contig, group)
        else:
            cmd += 'cp %s contigs/%s.fa\n' % (contig, group)

    cmd += 'diting.py'
    cmd += ' --reads reads'
    cmd += ' --assembly contigs'
    cmd += ' --outdir %s' % out
    if params['noclean']:
        cmd += ' --noclean'
    cmd += ' --threads %s\n' % params['cpus']
    # make the figures
    cmd += 'diting.py --visualization pathways_relative_abundance.tab\n'
    # remove the temporary, used files for gzip and concat reads, and contigs
    cmd += 'rm -rf tmp reads contigs KEGG_annotation/hmmout BBMap\n'
    return cmd


def diting_sample(self, all_reads: dict):
    """Run diting on a single sample.

    Parameters
    ----------
    self : Commands class instance
    all_reads : dict
    """
    for (tech, group), inputs in self.inputs[self.sam_pool].items():

        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir

        contigs = {group: inputs[0]}
        to_dos = status_update(self, tech, [inputs[0]], group=group)

        png = '%s/carbon_cycle_sketch.png' % out_dir
        if self.config.force or to_do(png):
            cmd = diting_cmd(self, tech, contigs, all_reads, out_dir)
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def diting_all(self, all_reads: dict):
    """Run diting on a single sample.

    Parameters
    ----------
    self : Commands class instance
    all_reads : dict
    """
    for pool, group_sams in self.pools.items():
        for tech in set([x[0] for x in self.inputs[pool]]):
            out_dir = '%s/%s/%s' % (self.dir, tech, pool)
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'][pool] = out_dir

            contigs = get_contigs(self, tech, pool)
            contigs_list = list(contigs.values())
            to_dos = status_update(self, tech, contigs_list)

            png = '%s/carbon_cycle_sketch.png' % out_dir
            if self.config.force or to_do(png):
                cmd = diting_cmd(self, tech, contigs, all_reads, out_dir, pool)
                key = (tech, pool)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=contigs_list, o_d=out_dir, key=key)
                self.soft.add_status(tech, pool, 1)
            else:
                self.soft.add_status(tech, pool, 0)


def diting(self):
    """DiTing is designed to determine the relative abundance of metabolic
    and biogeochemical functional pathways in a set of given metagenomic/
    metatranscriptomic data. The input is expected to be a folder containing
    a group of paired-end clean reads. These reads will be assembled, annotated,
    and parsed for producing a table of relative abundance of elemental/
    biogeochemical cycling pathways (e.g., Nitrogen, Carbon, Sulfur) in each
    sample. Sketch maps and heatmaps will also be produced accordingly for
    comparing biogeochemical functions visually.

    References
    ----------
    Xue, C.X., Lin, H., Zhu, X.Y., Liu, J., Zhang, Y., Rowley, G., Todd,
    J.D., Li, M. and Zhang, X.H., 2021. DiTing: a pipeline to infer and
    compare biogeochemical pathways from metagenomic and metatranscriptomic
    data. Frontiers in microbiology, p.2118.

    Notes
    -----
    GitHub  : https://github.com/xuechunxu/DiTing
    Paper   : https://doi.org/10.3389/fmicb.2021.698286

    Parameters
    ----------
    self : Commands class instance
    """
    if self.config.tools[self.soft.prev] != 'assembling':
        sys.exit('[tiara] can only be run on assembly output')
    all_reads = get_reads(self)

    if self.soft.params['samples'] == 'all':
        diting_all(self, all_reads)
    else:
        diting_sample(self, all_reads)


def eggnogmapper_cmd(
        self,
        tech: str,
        prefix: str,
        proteins: str,
        out_dir: str
) -> str:
    """Get EggNOG-mapper command line

    Parameters
    ----------
    self
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    prefix : str
        Prefix name for the output files
    proteins : str
        Path to the input proteins fasta file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        EggNOG-mapper command line
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    if proteins.endswith('.fa.gz') or proteins.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (proteins, proteins.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % proteins.rstrip('.gz')
        proteins = proteins.rstrip('.gz')

    if params['data_dir']:
        cmd += 'export EGGNOG_DATA_DIR=%s\n' % params['data_dir']
    cmd += 'emapper.py'
    cmd += ' -i %s' % proteins
    cmd += ' --output %s' % prefix
    cmd += ' --output_dir %s' % out_dir
    cmd += ' --temp_dir $TMPDIR'
    booleans = [
        'resume',
        'override',
        'translate',
        'outfmt_short',
        'dmnd_ignore_warnings',
        'report_no_hits',
        'cut_ga',
        'no_annot',
        'dbmem',
        'report_orthologs',
        'md5',
        'no_file_comments',
        'excel'
    ]
    for boolean in booleans:
        if params[boolean]:
            cmd += ' --%s' % boolean

    parameters = [
        'evalue',
        'start_sens',
        'sens_steps',
        'final_sens',
        'port',
        'end_port',
        'num_servers',
        'num_workers',
        'hmm_maxhits',
        'hmm_maxseqlen',
        'Z',
        'seed_ortholog_evalue',
        'mp_start_method',
        'itype',
        'qtype',
        'dbtype',
        'clean_overlaps',
        'tax_scope',
        'tax_scope_mode',
        'target_orthologs',
        'go_evidence',
        'pfam_realign',
        'decorate_gff',
        'decorate_gff_ID_field'
    ]
    for param in parameters:
        cmd += ' --%s %s' % (param, params[param])
        if param == 'itype' and params[param] in ['genome', 'metagenome']:
            cmd += ' --genepred %s' % params['genepred']
            if params['genepred'] == 'prodigal':
                if params['training_genome']:
                    cmd += ' --training_genome %s' % params['training_genome']
                if params['training_file']:
                    cmd += ' --training_file %s' % params['training_file']
            if params['genepred'] == 'search':
                if params['allow_overlaps']:
                    cmd += ' --allow_overlaps %s' % params['allow_overlaps']
                    cmd += ' --overlap_tol %s' % params['overlap_tol']
        if param == 'qtype':
            if params['usemem'] and not params['dbmem']:
                cmd += ' --usemem'

    parameters = [
        'annotate_hits_table',
        'cache',
        'data_dir',
        'trans_table',
        'pident',
        'query_cover',
        'subject_cover',
        'score',
    ]
    no_m = True
    for param in parameters:
        if params[param]:
            if param in ['pident', 'query_cover', 'subject_cover']:
                if params['m'] == 'hmmer':
                    continue
            cmd += ' --%s %s' % (param, params[param])
            if param == 'annotate_hits_table':
                cmd += ' -m no_search'
                no_m = False
            if param == 'cache':
                cmd += ' -m cache'
                no_m = False

    parameters = [
        'dmnd_algo',
        'dmnd_db',
        'dmnd_iterate',
        'sensmode',
        'matrix',
        'dmnd_frameshift',
        'gapopen',
        'gapextend',
        'block_size',
        'index_chunks',
        'mmseqs_db',
        'mmseqs_sub_mat',
        'database',
        'servers_list',
        'seed_ortholog_score',
        'target_taxa',
        'excluded_taxa',
    ]
    for param in parameters:
        if params[param]:
            cmd += ' --%s %s' % (param, params[param])

    if no_m:
        cmd += ' -m %s' % params['m']

    cmd += ' --cpu %s\n' % params['cpus']
    cmd += cmd_rm
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    return cmd


def get_eggnogmapper(
        self,
        tech: str,
        fastas: dict,
        sam_group: str,
) -> None:
    """Get the EggNOG-mapper command and fill the pipeline data structures.

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

        prefix = '-'.join([x for x in [sam_group, genome] if x])
        out = '%s/%s.emapper.hits.gz' % (out_dir, prefix)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)

        if self.config.force or to_do(out):
            cmd = eggnogmapper_cmd(self, tech, prefix, proteins, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=proteins, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1)
        else:
            self.soft.add_status(tech, self.sam_pool, 0)


def eggnogmapper(self):
    """EggNOG-mapper is a tool for fast functional annotation of novel
    sequences. It uses precomputed orthologous groups and phylogenies from
    the eggNOG database (http://eggnog5.embl.de) to transfer functional
    information from fine-grained orthologs only.

    References
    ----------
    Cantalapiedra, C.P., Hernández-Plaza, A., Letunic, I., Bork, P. and
    Huerta-Cepas, J., 2021. eggNOG-mapper v2: functional annotation,
    orthology assignments, and domain prediction at the metagenomic scale.
    Molecular biology and evolution, 38(12), pp.5825-5829.

    Notes
    -----
    GitHub  : https://github.com/eggnogdb/eggnog-mapper
    Docs    : http://eggnog-mapper.embl.de/
    Paper   : https://doi.org/10.1093/molbev/msab293

    Parameters
    ----------
    self
    """
    if self.sam_pool in self.pools:
        if self.soft.prev != 'prodigal':
            sys.exit('[eggnogmapper] Runs on protein data (prodigal only)')
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            fastas = group_inputs(self, inputs)
            get_eggnogmapper(self, tech, fastas, group)
    else:
        prev = self.soft.prev
        if prev not in ['plass', 'prodigal']:
            sys.exit('[eggnogmapper] Runs on protein data (plass or prodigal)')
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_eggnogmapper(self, tech, fastas, self.sam_pool)


def keggcharter_cmd(
        self,
        tech: str,
        tables: list,
        out_dir: str
) -> str:
    """Get the KEGGCharter command.

    Parameters
    ----------
    self
    tech : str
        Technology
    tables : list
        Woltka's ko.tsv data tables
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str

    """
    params = tech_params(self, tech)

    cmd = ''
    for tb in tables:

        cols = []
        with open(tb.replace('${SCRATCH_FOLDER}', '')) as f:
            for line in f:
                cols = line.strip().split('\t')[1:]
                break

        tx = splitext(tb)[0].split('/ko_')[-1]
        if tx not in params['taxonomic_levels']:
            continue
        o = None
        if self.soft.prev == 'woltka':
            o = tb.replace('.tsv', '_charter.tsv')
            cmd += 'head -1 %s | sed "s/#OTU ID/%s\\tko/" > %s\n' % (tb, tx, o)
            cmd += 'tail -n +2 %s | sed "s/|/\\t/" >> %s\n' % (tb, o)

        cmd += 'keggcharter.py'
        cmd += ' --output %s/%s' % (out_dir, tx)
        cmd += ' --taxa-column %s' % tx
        cmd += ' --number-of-taxa 10'
        cmd += ' --metabolic-maps %s' % ','.join(params['metabolic_maps'])
        for t in ['taxa_list']:
            if not params[t]:
                continue
            cmd += ' --%s %s' % (t.replace('_', '-'), ','.join(params[t]))
        if self.soft.prev == 'woltka':
            cmd += ' --ko-column ko'
            cmd += ' --genomic-columns %s' % ','.join(cols)
        elif params['input_quantification'] or not params['genomic_columns']:
            cmd += ' --input-quantification'
        else:
            cmd += ' --genomic-columns %s' % ','.join(params['genomic_columns'])
        if params['resume']:
            cmd += ' --resume'
        cmd += ' --file %s\n' % o
        cmd += 'rm %s\n' % o

    return cmd


def get_keggcharter(
        self,
        tech: str,
        ali_group: str,
        inputs: list
) -> None:
    """

    Parameters
    ----------
    self
    tech : str
    ali_group
    inputs
    """
    for inp in inputs:

        key = genome_key(tech, ali_group)
        out_dir = '/'.join([self.dir, tech, ali_group])
        self.outputs['dirs'].append(out_dir)

        tabs = []
        if self.soft.prev == 'woltka':
            tax = ['none', 'species', 'genus', 'family', 'order', 'class',
                   'phylum']
            tabs = ['%s/kegg/ko_%s.tsv' % (inp, t) for t in tax]
            tabs = [tab for tab in tabs if not to_do(tab)]

        to_dos = status_update(self, tech, tabs, group=ali_group)

        out = '%s/out' % out_dir
        self.outputs['outs'].setdefault((tech, ali_group), []).append(out_dir)

        if self.config.force or to_do(out):
            cmd = keggcharter_cmd(self, tech, tabs, out_dir)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=tabs, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1)
        else:
            self.soft.add_status(tech, self.sam_pool, 0)


def keggcharter(self):
    """

    References
    ----------
    Sequeira, J.C., Rocha, M., Alves, M.M. and Salvador, A.F., 2022. UPIMAPI,
    reCOGnizer and KEGGCharter: Bioinformatics tools for functional
    annotation and visualization of (meta)-omics datasets. Computational and
    Structural Biotechnology Journal, 20, pp.1798-1810.

    Notes
    -----
    GitHub  : https://github.com/iquasere/KEGGCharter
    Paper   : https://doi.org/10.1016/j.csbj.2022.03.042

    Parameters
    ----------
    self : Commands class instance
    """
    prevs = ['woltka', 'eggnogmapper']
    if self.soft.prev not in prevs:
        sys.exit('[keggcharter] Only possible after "%s"' % '", "'.join(prevs))
    elif self.soft.prev == 'woltka':
        for (tech, aligner), inputs in self.inputs.items():
            get_keggcharter(self, tech, aligner, inputs)

    elif self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            print(tech)
            print(group)
            print(inputs)
            print(inputscsa)
            fastas = group_inputs(self, inputs)
            get_keggcharter(self, tech, fastas, group)


def trf_cmd(
        fasta: str,
        params: dict,
        out_dir: str,
        out: str
) -> str:
    """Collect trf command.

    Parameters
    ----------
    fasta : str
        Path to the input fasta file
    params : dict
        Parameters
    out_dir : str
        Paths to the output folder
    out : str
        Path to the output file

    Returns
    -------
    cmd : str
        trf command
    """
    cmd_rm = ''
    if params['path']:
        cmd = 'export PATH=$PATH:%s\n' % params['path']
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta.rstrip('.gz')
        fasta = fasta.rstrip('.gz')

    cmd += 'mkdir -p %s\n' % out_dir
    cmd += 'cd %s\n' % out_dir
    cmd += 'trf %s' % fasta
    for param in ['match', 'mismatch', 'delta', 'match_probability',
                  'indel_probability', 'min_score', 'max_period']:
        cmd += ' %s' % params[param]
    for boolean in ['m', 'f', 'd', 'h', 'r']:
        if params[boolean]:
            cmd += ' -%s' % boolean
    cmd += ' -l %s' % params['l']
    cmd += ' -ngs > %s\n' % out
    cmd += cmd_rm
    return cmd


def get_trf(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .soft.status
            Current status of the pipeline in terms of available outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    folders : dict
        Paths to the input fasta files per genome/MAG
    group : str
        Name of the current co-assembly group
    """
    for genome, inputs in folders.items():

        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        to_dos = status_update(self, tech, inputs, group=group, genome=genome)

        params = tech_params(self, tech)
        for fasta in inputs:
            base = basename(fasta).rstrip('.gz')
            nums = '.'.join([str(params[x]) for x in [
                'match', 'mismatch', 'delta', 'match_probability',
                'indel_probability', 'min_score', 'max_period']])
            out = '%s/%s.%s.out' % (out_dir, base, nums)
            if self.config.force or to_do(out):
                cmd = trf_cmd(fasta, params, out_dir, out)
                key = genome_key(tech, group, genome)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta, o_d=out_dir, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=group, genome=genome)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=group, genome=genome)


def trf(self):
    """A tandem repeat in DNA is two or more adjacent, approximate copies of
    a pattern of nucleotides. Tandem Repeats Finder is a program to locate
    and display tandem repeats in DNA sequences. In order to use the program,
    the user submits a sequence in FASTA format. There is no need to specify
    the pattern, the size of the pattern or any other parameter.

    References
    ----------
    Benson G. Tandem repeats finder: a program to analyze DNA sequences.
    Nucleic Acids Res. 1999; 27(2):573–580.

    Notes
    -----
    GitHub  : https://github.com/Benson-Genomics-Lab/TRF
    Paper   : https://doi.org/10.1093/nar/27.2.573
    Docs    : https://tandem.bu.edu/trf/trf.html

    Parameters
    ----------
    self : Commands class instance
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_trf(self, tech, folders, group)


def kmerssr_cmd(
        fasta: str,
        params: dict,
        out: str
) -> str:
    """Collect kmerssr command.

    Parameters
    ----------
    fasta : str
        Path to the input fasta file
    params : dict
        Parameters
    out : str
        Paths to the output file

    Returns
    -------
    cmd : str
        kmerssr command
    """
    cmd_rm = ''
    if params['path']:
        cmd = 'export PATH=$PATH:%s\n' % params['path']
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta.rstrip('.gz')
        fasta = fasta.rstrip('.gz')

    cmd += 'kmer-ssr'
    cmd += ' -i %s' % fasta
    cmd += ' -o %s' % out
    for param in ['a', 'p', 'l', 'L', 'n', 'N', 'r', 'R', 'Q']:
        cmd += ' -%s %s' % (param, params[param])
    if params['s']:
        cmd += ' -s %s' % params['s']
    for boolean in ['A', 'e', 'd']:
        if params[boolean]:
            cmd += ' -%s' % boolean
    cmd += ' -t %s\n' % params['cpus']
    cmd += cmd_rm
    return cmd


def get_kmerssr(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .soft.status
            Current status of the pipeline in terms of available outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    folders : dict
        Paths to the input fasta files per genome/MAG
    group : str
        Name of the current co-assembly group
    """
    for genome, inputs in folders.items():

        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        to_dos = status_update(self, tech, inputs, group=group, genome=genome)

        params = tech_params(self, tech)
        for fasta in inputs:
            out = '%s/output.tsv' % out_dir
            if self.config.force or to_do(out):
                cmd = kmerssr_cmd(fasta, params, out)
                key = genome_key(tech, group, genome)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta, o_d=out_dir, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=group, genome=genome)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=group, genome=genome)


def kmerssr(self):
    """Kmer-SSR is a software tool developed to find Simple Sequence Repeats
    (SSRs) in a sequence (presumably of DNA or RNA). SSRs are sometimes
    referred to as Short Tandem Repeats (STRs) or microsatellites. SSRs are
    genetic markers with several interesting and meaningful biological
    implications. For example, SSRs can play significant roles in genome
    alignment against a reference and species identification.

    References
    ----------
    Pickett, B.D., Miller, J.B. and Ridge, P.G., 2017. Kmer-SSR: a fast and
    exhaustive SSR search algorithm. Bioinformatics, 33(24), pp.3922-3928.

    Notes
    -----
    GitLab  : https://github.com/ridgelab/Kmer-SSR
    Paper   : https://doi.org/10.1093/bioinformatics/btx538

    Parameters
    ----------
    self : Commands class instance
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_kmerssr(self, tech, folders, group)


def divissr_cmd(
        fasta: str,
        params: dict,
        gff: str,
        out: str
) -> str:
    """Collect divissr command.

    Parameters
    ----------
    fasta : str
        Path to the input fasta file
    params : dict
        Parameters
    gff : str
        Path to the prodigal GFF file
    out : str
        Path to the output file

    Returns
    -------
    cmd : str
        divissr command
    """
    cmd, cmd_rm = '', ''
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta.rstrip('.gz')
        fasta = fasta.rstrip('.gz')

    cmd += 'divissr'
    cmd += ' --input %s' % fasta
    cmd += ' --output %s' % out
    for param in ['min_motif_size', 'max_motif_size', 'min_length', 'gene_key',
                  'comp_dist', 'up_promoter', 'down_promoter', 'anno_format']:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    for boolean in ['filter_reads', 'analyse', 'compound']:
        if params[boolean]:
            cmd += ' -%s' % boolean
    if gff:
        cmd += ' --annotate %s' % gff
    cmd += ' --format fasta\n'
    cmd += cmd_rm
    return cmd


def get_divissr(
        self,
        tech: str,
        folders: dict,
        group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .soft.status
            Current status of the pipeline in terms of available outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    folders : dict
        Paths to the input fasta files per genome/MAG
    group : str
        Name of the current co-assembly group
    """
    for genome, inputs in folders.items():

        gff = None

        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        to_dos = status_update(self, tech, inputs, group=group, genome=genome)

        params = tech_params(self, tech)
        for fasta in inputs:
            out = '%s/output.tsv' % out_dir
            if self.config.force or to_do(out):
                cmd = divissr_cmd(fasta, params, gff, out)
                key = genome_key(tech, group, genome)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta, o_d=out_dir, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=group, genome=genome)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=group, genome=genome)


def divissr(self):
    """DiviSSR is a DNA tandem repeat identification tool. Tandem repeats
    are important genomic sequences which have functional and evolutionary
    significance.

    References
    ----------
    Avvaru, A.K., Mishra, R.K. and Sowpati, D.T., 2021. DiviSSR: Simple
    arithmetic for efficient identification of tandem repeats. bioRxiv,
    pp.2021-10.

    Notes
    -----
    GitLab  : https://github.com/avvaruakshay/divissr
    Paper   : https://doi.org/10.1101/2021.10.05.462997

    Parameters
    ----------
    self : Commands class instance
    """
    if self.sam_pool in self.pools:
        # if 'prodigal' in self.softs and self.softs['prodigal'].prev == :
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_divissr(self, tech, folders, group)


def size_cmd(fasta, out):
    if fasta.endswith('.gz'):
        cmd = 'zcat %s | grep -c ">" > %s\n' % (fasta, out)
    else:
        cmd = 'grep -c ">" %s > %s\n' % (fasta, out)
    cmd += 'gzip %s\n' % out
    return cmd


def get_size(self, tech, folders, group):
    for genome, inputs_ in folders.items():

        out_dir = genome_out_dir(self, tech, group, genome)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        if self.soft.prev == 'prodigal':
            fps = ['nucleotide.sequences', 'potential.starts',
                   'protein.translations']
            inputs = ['%s/%s.fasta.gz' % (inputs_[0], fp) for fp in fps]
        else:
            inputs = inputs_
        to_dos = status_update(self, tech, inputs, group=group, genome=genome)

        for fasta in inputs:
            base = splitext(basename(fasta.replace('.gz', '')))[0]
            out = '%s/%s.size.tsv' % (out_dir, base)
            if self.config.force or to_do('%s.gz' % out):
                cmd = size_cmd(fasta, out)
                key = genome_key(tech, group, genome)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta, o_d=out_dir, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=group, genome=genome)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=group, genome=genome)


def size(self):
    """Measure the size of each fasta entry.

    Parameters
    ----------
    self : Commands class instance
    """
    predicting = self.config.tools['annotation (protein prediction)']
    assembling = self.config.tools['assembling']
    if self.soft.prev not in (predicting + assembling):
        sys.exit("[size] Calculation only after an assembly/protein prediction")

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_size(self, tech, folders, group)


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


def circlator(self):
    """Circlator will attempt to identify each circular sequence and
    output a linearised version of it. It does this by assembling all reads
    that map to contig ends and comparing the resulting contigs with the
    input assembly.

    References
    ----------
    Hunt, M., Silva, N.D., Otto, T.D., Parkhill, J., Keane, J.A. and Harris,
    S.R., 2015. Circlator: automated circularization of genome assemblies
    using long sequencing reads. Genome biology, 16(1), pp.1-10.

    Notes
    -----
    GitLab  : https://github.com/sanger-pathogens/circlator
    Paper   : https://doi.org/10.1186/s13059-015-0849-0
    Docs    : https://sanger-pathogens.github.io/circlator/

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


def srst2(self):
    """Short Read Sequence Typing for Bacterial Pathogens.
    This program is designed to take Illumina sequence data, a MLST database
    and/or a database of gene sequences (e.g. resistance genes, virulence
    genes, etc) and report the presence of STs and/or reference genes.

    References
    ----------
    Inouye, M., Dashnow, H., Raven, L.A., Schultz, M.B., Pope, B.J., Tomita,
    T., Zobel, J. and Holt, K.E., 2014. SRST2: Rapid genomic surveillance for
    public health and hospital microbiology labs. Genome medicine, 6(11),
    pp.1-16.

    Notes
    -----
    GitHub  : https://github.com/katholt/srst2
    Paper   : https://doi.org/10.1186/s13073-014-0090-6

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def pirate(self):
    """Cataloguing genes and their distributions within natural bacterial
    populations is essential for understanding evolutionary processes and the
    genetic bases of adaptation. genes that are shared between different
    bacterial strains and species is essential for understanding the genomic
    variation that underlies the enormous phenotypic variation observed in
    the microbial world. Here we present a pangenomics toolbox, PIRATE,
    which identifies and classifies orthologous gene families in bacterial
    pangenomes over a wide range of sequence similarity thresholds. PIRATE
    builds upon recent scalable software developments for the rapid
    interrogation of pangenomes from large dat thousands of genomes. PIRATE
    clusters genes (or other annotated features) over a wide range of
    amino-acid or nucleotide identity thresholds, and classifies paralogous
    genes families into either putative gene fission/fusion events or gene
    duplications. Furthermore, PIRATE provides a measure of allelic variance
    and cluster homology, and orders the resulting pangenome on a pangenome
    graph. Additional scripts are provided for comparison and visualization.
    PIRATE provides a robust framework for analysing the pangenomes of
    bacteria, from largely clonal to panmictic species.

    References
    ----------
    Bayliss, S.C., Thorpe, H.A., Coyle, N.M., Sheppard, S.K. and Feil, E.J.,
    2019. PIRATE: A fast and scalable pangenomics toolbox for clustering
    diverged orthologues in bacteria. Gigascience, 8(10), p.giz119.

    Notes
    -----
    GitHub  : https://github.com/SionBayliss/PIRATE
    Paper   : https://doi.org/10.1093/gigascience/giz119

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def mummer2circos(self):
    """Generate circular bacterial genome plots based on BLAST or
    NUCMER/PROMER alignments. Generate SVG and PNG images with circos.

    Notes
    -----
    GitHub  : https://github.com/metagenlab/mummer2circos

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