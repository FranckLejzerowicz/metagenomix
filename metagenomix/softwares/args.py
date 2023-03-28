# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename, dirname, splitext
from metagenomix._inputs import (
    sample_inputs, group_inputs, genome_key, genome_out_dir, get_contigs)
from metagenomix._io_utils import (
    caller, io_update, to_do, tech_specificity, not_paired,
    status_update, get_assembly)
from metagenomix.core.parameters import tech_params


def predict_cmd(
        self,
        fasta: str,
        prefix: str,
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
    prefix : str
        Paths to the output file's prefix
    typ : str
        Type of input data: 'nucl' (DNA) or 'prot' (protein)

    Returns
    -------
    cmd : str
        deeparg predict commands
    """
    cmd, cmd_rm = '', ''
    if fasta.endswith('.fa.gz') or fasta.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (fasta, fasta.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % fasta.rstrip('.gz')
        fasta = fasta.rstrip('.gz')

    cmd += 'deeparg predict'
    cmd += ' --input-file %s' % fasta
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
    cmd += ' --output-file %s\n' % prefix
    cmd += 'for i in %s*; do gzip -q $i; done\n' % prefix
    cmd += cmd_rm
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
        typ_seqs = {'prot': '%s/protein.translations.fasta.gz' % fasta,
                    'nucl': '%s/nucleotide.sequences.fasta.gz' % fasta}
    elif self.soft.prev == 'plass':
        typ_seqs = {'prot': '%s/prot_contigs.fasta.gz' % dirname(fasta),
                    'nucl': '%s/nucl_contigs.fasta.gz' % dirname(fasta)}
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
        Sample name or name of a co-assembly pool's group
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        typ_seqs = predict_inputs(self, fasta[0])
        for typ, seq in typ_seqs.items():

            key = genome_key(tech, sam_group, genome)
            to_dos = status_update(
                self, tech, [seq], group=sam_group, genome=genome)

            base = splitext(basename(seq.rstrip('.gz')))[0]
            prefix = '%s/%s' % (out_dir, base)
            arg = '%s.mapping.ARG.gz' % prefix
            pot_arg = '%s.mapping.potential.ARG.gz' % prefix
            outs = [arg, pot_arg]
            self.outputs['outs'].setdefault((tech, sam_group), []).append(outs)

            # check if tool already run (or if --force) to allow getting command
            if self.config.force or to_do(arg):
                # collect the command line
                cmd = predict_cmd(self, seq, prefix, typ)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=seq, o_d=out_dir, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=sam_group, genome=genome)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=sam_group, genome=genome)


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
            Co-assembly pools and sample per group
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
    else:
        tech_fastas = sample_inputs(self)
        for tech, fastas in tech_fastas.items():
            get_predict(self, fastas, tech, self.sam_pool)


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
    for param in [
        'deeparg_data_path', 'deeparg_identity', 'deeparg_probability',
        'deeparg_evalue', 'gene_coverage', 'bowtie_16s_identity'
    ]:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])
    cmd += ' --output_file %s\n' % prefix
    cmd += 'for i in %s*; do gzip -q $i; done\n' % prefix
    return cmd


def short(self) -> None:
    """Run deeparg short_reads_pipeline module.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        to_dos = status_update(self, tech, fastqs)

        out = '%s/%s/%s' % (self.dir, tech, self.sam_pool)
        self.outputs['dirs'].append(out)

        prefix = out + '/' + self.sam_pool
        arg = '%s.mapping.ARG.gz' % prefix
        pot_arg = '%s.potential.ARG.gz' % prefix
        outs = [arg, pot_arg]
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(outs)

        if self.config.force or to_do(arg):
            cmd = short_cmd(self, fastqs, prefix)
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = [cmd]
            io_update(self, i_f=fastqs, o_d=out, key=tech)
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def deeparg(self) -> None:
    """A deep learning based approach to predict Antibiotic Resistance Genes
    (ARGs) from metagenomes. It provides two models,deepARG-SS and deepARG-LS.

    References
    ----------
    Arango-Argoty, Gustavo, et al. "DeepARG: a deep learning approach for
    predicting antibiotic resistance genes from metagenomic data." Microbiome
    6.1 (2018): 1-15.

    Notes
    -----
    BitBucket   : https://bitbucket.org/gusphdproj/deeparg-ss/src/master/
    Website     : https://bench.cs.vt.edu/deeparg
    Docs        : https://readthedocs.org/projects/deeparg/
    Paper       : https://doi.org/10.1186/s40168-018-0401-z

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    # This function splits the name of the software and calls as function
    # the last underscore-separated field (which is in this module)
    module_call = caller(self, __name__)
    module_call(self)


def mmarc_cmd(
        self,
        tech: str,
        inputs: list,
        out: str,
        genome: str,
        level: str
) -> str:
    """Collect mmarc (meta-marc) command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Paths to the input files
    out : str
        Paths to the output folder
    genome : str
        MAGs/Genomes folder name or empty string (for assembly contigs)
    level : str
        Model level

    Returns
    -------
    cmd : str
        mmarc (meta-marc) command
    """
    cmd, cmd_rm = '', ''
    if genome or len(inputs) == 1:
        if inputs[0].endswith('.fa.gz') or inputs[0].endswith('.fasta.gz'):
            cmd += 'gunzip -c %s > %s\n' % (out, out.rstrip('.gz'))
            cmd_rm += 'rm %s\n' % out.rstrip('.gz')
            inp = inputs[0].rstrip('.gz')
        cmd += '\n./mmarc --input %s' % inp
    else:
        cmd += '\n./mmarc --i1 %s --i2 %s' % tuple(inputs)
        if tech == 'illumina':
            cmd += ' --dedup'
        cmd += ' --multicorrect'
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --filename output'
    cmd += ' --skewness %s/skewness.txt' % out
    cmd += ' --graphs %s/graphs' % out
    cmd += ' --level %s' % level
    for boolean in ['dedup', 'multicorrect']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    for param in ['coverage', 'evalue', 'kmer']:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])
    cmd += ' --output %s\n' % out
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out
    cmd += 'gzip -q %s/duplicate_tables/output_dupcounts.txt\n' % out
    cmd += cmd_rm
    return cmd


def get_metamarc(
        self,
        tech: str,
        folders: dict,
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
        .soft.status
            Current status of the pipeline in terms of available outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    folders : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Name of the current sample or co-assembly group
    """
    for genome, inputs in folders.items():

        out_dir = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)
        to_dos = status_update(
            self, tech, inputs, self.sam_pool, group=sam_group, genome=genome)

        outs, cmd = [], ''
        key = genome_key(tech, sam_group, genome)
        for level in self.soft.params['level']:
            out = '%s/model_level_%s' % (out_dir, level)
            self.outputs['dirs'].append(out_dir)
            outs.append(out)
            out_fp = '%s/output_dedup.fasta.gz' % out
            if self.config.force or to_do(out_fp):
                # collect the command line
                cmd += mmarc_cmd(self, tech, inputs, out, genome, level)

        if cmd:
            path = self.soft.params['path']
            full_cmd = 'cd %s\n' % out_dir
            full_cmd += 'cp -r %s/bin %s/.\n' % (path, out_dir)
            full_cmd += 'cp -r %s/src %s/.\n' % (path, out_dir)
            full_cmd += 'cp -r %s/hmmer-3.1b2 %s/.\n' % (path, out_dir)
            full_cmd += 'cd bin\n'
            full_cmd += cmd
            full_cmd += '\ncd ..\n'
            full_cmd += 'rm -rf bin\n'
            full_cmd += 'rm -rf src\n'
            full_cmd += 'rm -rf hmmer-3.1b2\n'
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(full_cmd)
            io_update(self, i_f=inputs, o_d=outs, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def metamarc(self) -> None:
    """Meta-MARC is a set of profile Hidden Markov Models developed for the
    purpose of screening and profiling resistance genes in DNA-based
    metagenomic data. This tool was developed for the characterization of
    various resistance classes, mechanisms, and gene/operon groups from raw
    sequencing data much in the way that microbiome softwares profile the
    bacterial taxonomy of metagenomic samples. Meta-MARC is not intended to
    be used as a final annotation or all-in-one tool; this software simply
    offers the user another view of complex metagenomic data. Users are
    encouraged to perform further statistical testing on their data to
    validate results obtained from Meta-MARC.

    References
    ----------
    Lakin, Steven M., et al. "Hierarchical Hidden Markov models enable
    accurate and diverse detection of antimicrobial resistance sequences."
    Communications biology 2.1 (2019): 1-11.

    Notes
    -----
    GitHub  : https://github.com/lakinsm/meta-marc
    Paper   : https://doi.org/10.1038/s42003-019-0545-9

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_metamarc(self, tech, folders, group)
    else:
        tech_fastas = sample_inputs(self, raw=True)
        for tech, fastas in tech_fastas.items():
            get_metamarc(self, tech, fastas, self.sam_pool)


def karga_kargva_cmd(
        self,
        tech: str,
        merged: str,
        out_dir: str,
        genes: str,
        reads: str,
        db: str,
        db_path: str
) -> str:
    """Collect KARGA or KARGVA command.

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .soft.params
            Parameters
    tech : str
        Technology
    merged : str
        Path to the input file
    out_dir : str
        Path to the output folder
    genes : str
        Path to the output mappedGenes.csv
    reads : str
        Path to the output mappedReads.csv
    db : str
        Name of the database
    db_path : str
        Path to ARG/MGE fasta file with resistance annotation/mutation in header

    Returns
    -------
    cmd : str
        KARGA or KARGVA command
    """
    params = tech_params(self, tech)
    cmd = 'java %s' % self.soft.name.upper()
    cmd += ' f:%s' % merged
    cmd += ' d:%s' % db_path
    for param in ['k', 'i', 's', 'm', 'r']:
        if param in ['s', 'r'] and self.soft.name != 'karga':
            continue
        if self.soft.params[param]:
            cmd += ' %s:%s' % (param, params[param])
    cmd += ' -Xmx%s%s\n' % (params['mem'], params['mem_dim'])

    rad = merged.rstrip('.gz').rstrip('.fastq')
    reads_ = '%s_%s_mappedReads.csv' % (rad, self.soft.name.upper())
    genes_ = '%s_%s_mappedGenes.csv' % (rad, self.soft.name.upper())
    cmd += 'mkdir -p %s/db_%s\n' % (out_dir, db)
    cmd += 'mv %s %s\n' % (reads_, reads)
    cmd += 'mv %s %s\n' % (genes_, genes)
    cmd += 'grep -v ",?" %s > %s.tmp\n' % (reads, reads)
    cmd += 'mv %s.tmp %s\n' % (reads, reads)
    cmd += 'gzip %s %s\n' % (reads, genes)
    return cmd


def get_karga_kargva(
        self,
        tech: str,
        fastas_folders: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .outputs : dict
            All outputs
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastas_folders : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Name of the current sample or co-assembly group
    """
    for genome, inputs in fastas_folders.items():

        fastx = inputs[0]
        out_dir = genome_out_dir(self, tech, sam_group)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)
        to_dos = status_update(
            self, tech, inputs, group=sam_group, genome=genome)

        merge_cmd, rm_cmd = '', ''
        if len(inputs) == 1:
            merged = fastx
        elif len(inputs) > 1:
            merged = '%s/merged_%s' % (out_dir, basename(fastx))
            merge_cmd += 'cat %s > %s\n' % (' '.join(inputs), merged)
            rm_cmd += 'rm %s\n' % merged

        outs, cmd = [], ''
        bas = basename(fastx).rstrip('.gz').rstrip('.fastq')
        for db, db_path in self.soft.params['databases'].items():
            self.outputs['dirs'].append(out_dir)
            out_db = '%s/db_%s' % (out_dir, db)
            outs.append(out_db)
            rad = '%s/%s_%s_mapped' % (out_db, bas, self.soft.name.upper())
            genes, reads = '%sGenes.csv' % rad, '%sReads.csv' % rad
            if self.config.force or to_do(genes+'.gz') or to_do(reads+'.gz'):
                # collect the command line
                cmd += karga_kargva_cmd(
                    self, tech, merged, out_dir, genes, reads, db, db_path)
        if cmd:
            if self.soft.name == 'karga':
                ps = ['openjdk-8', 'AMRGene.class']
            elif self.soft.name == 'kargva':
                ps = ['AMRVariantGene.class', 'kargva_db_v5.fasta']
            path = self.soft.params['path']
            full_cmd = 'cd %s\n' % out_dir
            full_cmd += 'cp -r %s/KARG* %s/.\n' % (path, out_dir)
            for p in ps:
                full_cmd += 'cp -r %s/%s %s/.\n' % (path, p, out_dir)
            full_cmd += merge_cmd + cmd
            for p in ps:
                full_cmd += 'rm -rf %s\n' % p
            full_cmd += 'rm -rf KARG*\n'
            full_cmd += rm_cmd

            key = genome_key(tech, sam_group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(full_cmd)
            io_update(self, i_f=inputs, o_d=outs, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def karga_kargva(self) -> None:
    """Get the input for KARGA or KARGVA

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if not self.sam_pool or self.sam_pool in self.pools:
        sys.exit('[%s] Only run on non-assembled reads' % self.soft.name)
    else:
        tech_fastas = sample_inputs(self, raw=True)
        for tech, fastas in tech_fastas.items():
            get_karga_kargva(self, tech, fastas, self.sam_pool)


def karga(self) -> None:
    """K-mer-based antibiotic gene resistance analyzer, a multi-platform Java
    toolkit for identifying ARGs from metagenomic short read data. KARGA does
    not perform alignment; it uses an efficient double-lookup strategy,
    statistical filtering on false positives, and provides individual read
    classification as well as covering of the database resistome. On simulated
    data, KARGA’s antibiotic resistance class recall is 99.89% for error/
    mutation rates within 10%, and of 83.37% for error/mutation rates between
    10% and 25%, while it is 99.92% on ARGs with rearrangements. On empirical
    data, KARGA provides higher hit score (≥1.5-fold) than AMRPlusPlus, DeepARG,
    and MetaMARC. KARGA has also faster runtimes than all other softwares (2x faster
    than AMRPlusPlus, 7x than DeepARG, and over 100x than MetaMARC).

    References
    ----------
    M. Prosperi and S. Marini, "KARGA: Multi-platform Toolkit for k-mer-based
    Antibiotic Resistance Gene Analysis of High-throughput Sequencing Data,"
    2021 IEEE EMBS International Conference on Biomedical and Health Informatics
    (BHI), 2021, pp. 1-4, doi: 10.1109/BHI50953.2021.9508479.

    Notes
    -----
    GitHub  : https://github.com/DataIntellSystLab/KARGA
    Paper   : https://ieeexplore.ieee.org/document/9508479

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    karga_kargva(self)


def kargva(self) -> None:
    """K-mer-based Antibiotic Resistance Gene Variant Analyzer (KARGVA).
    KARGVA is a Multi-platform Toolkit for Identification of Antibiotic
    Resistance from Sequencing Data Conferred by Point Mutations in
    Bacterial Genes.

    References
    ----------
    M. Prosperi and S. Marini, "KARGA: Multi-platform Toolkit for k-mer-based
    Antibiotic Resistance Gene Analysis of High-throughput Sequencing Data,"
    2021 IEEE EMBS International Conference on Biomedical and Health Informatics
    (BHI), 2021, pp. 1-4, doi: 10.1109/BHI50953.2021.9508479.

    Notes
    -----
    GitHub  : https://github.com/DataIntellSystLab/KARGVA
    Paper   : https://ieeexplore.ieee.org/document/9508479

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    karga_kargva(self)


def abricate_cmd(
        self,
        tech: str,
        inputs: list,
        out_d: str,
) -> str:
    """Collect abricate command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech : str
        Technology
    inputs : list
        Paths to the input files
    out_d : str
        Paths to the output folder

    Returns
    -------
    cmd : str
        abricate command
    """
    cmd, seqtk_cmd, cmd_rm = '', '', ''
    if tech in ['nanopore', 'pacbio']:
        outs = []
        for idx, inp in enumerate(inputs):
            out = '%s/%s.fa' % (out_d, idx)
            outs.append(out)
            seqtk_cmd += 'seqtk seq -A %s > %s\n' % (inp, out)
        out = '%s/%s.fasta' % (out_d, splitext(basename(inputs[0]))[0])
        seqtk_cmd += 'cat %s > %s\n' % (' '.join(outs), out)
        cmd_rm += 'rm %s\n' % ' '.join(outs)
    else:
        out = inputs[0]
        if out.endswith('.fa.gz') or out.endswith('.fasta.gz'):
            cmd += 'gunzip -c %s > %s\n' % (out, out.rstrip('.gz'))
            cmd_rm += 'rm %s\n' % out.rstrip('.gz')
            out = out.rstrip('.gz')

    for db in self.soft.params['databases']:
        cmd += '\nabricate'
        cmd += ' %s' % out
        cmd += ' --threads %s' % self.soft.params['cpus']
        cmd += ' --db %s' % db
        for boolean in ['setupdb', 'noheader', 'csv', 'nopath']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean
        for param in ['minid', 'mincov']:
            cmd += ' --%s %s' % (param, self.soft.params[param])
        cmd += ' > %s/%s.txt\n' % (out_d, db)
        cmd += 'gzip -q %s/%s.txt\n' % (out_d, db)

    if seqtk_cmd:
        cmd = seqtk_cmd + cmd + cmd_rm
    return cmd


def get_abricate(
        self,
        tech: str,
        folders: dict,
        sam_group: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    folders : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Name of the current sample or co-assembly group
    """
    for genome, inputs in folders.items():

        out_d = genome_out_dir(self, tech, sam_group, genome)
        self.outputs['dirs'].append(out_d)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_d)
        to_dos = status_update(
            self, tech, inputs, group=sam_group, genome=genome)

        o = ['%s/%s.txt.gz' % (out_d, x) for x in self.soft.params['databases']]
        # check if tool already run (or if --force) to allow getting command
        if self.config.force or sum([to_do(x) for x in o]):
            # collect the command line
            cmd = abricate_cmd(self, tech, inputs, out_d)
            key = genome_key(tech, sam_group, genome)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=inputs, o_d=out_d, key=key)
            self.soft.add_status(
                tech, self.sam_pool, 1, group=sam_group, genome=genome)
        else:
            self.soft.add_status(
                tech, self.sam_pool, 0, group=sam_group, genome=genome)


def abricate(self) -> None:
    """Mass screening of contigs for antimicrobial resistance or virulence
    genes. It comes bundled with multiple databases: NCBI, CARD, ARG-ANNOT,
    Resfinder, MEGARES, EcOH, PlasmidFinder, Ecoli_VF and VFDB.

    Notes
    -----
    GitHub  : https://github.com/tseemann/abricate

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_abricate(self, tech, folders, group)
    else:
        tech_fastas = sample_inputs(self, ['nanopore', 'pacbio'], raw=True)
        for tech, fastas in tech_fastas.items():
            get_abricate(self, tech, fastas, self.sam_pool)


def amrplusplus2_cmd(
        self,
        fastqs: list,
        out_dir: str,
) -> str:
    """Collect amrplusplus2 command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastqs : list
        Paths to the input files
    out_dir : str
        Paths to the output folder

    Returns
    -------
    cmd : str
        amrplusplus2 command
    """
    cmd = 'mkdir %s/copy\n' % out_dir
    cmd += 'cp -r %s/* %s/copy/.\n' % (self.soft.params['path'], out_dir)
    cmd += 'cd %s/copy\n' % out_dir

    cmd += 'nextflow run main_AmrPlusPlus_v2.nf'
    # input argument
    if len(fastqs) > 1:
        cmd += ' --reads "%s"' % fastqs[0].replace('R1', '{R1,R2}')
    else:
        cmd += ' --reads %s' % fastqs[0]

    # optional arguments - "valued"
    for param in ['adapters', 'fqc_adapters', 'host_index', 'host',
                  'kraken_db', 'amr_index', 'amr', 'annotation',
                  'snp_annotation', 'snp_confirmation']:
        if self.soft.params[param]:
            cmd += ' --%s "%s"' % (param, self.soft.params[param])

    for param in ['leading', 'trailing', 'slidingwindow', 'minlen',
                  'threshold', 'min', 'max', 'skip', 'samples']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --output "%s"' % out_dir
    cmd += ' -work-dir "%s/work_dir"\n' % out_dir
    cmd += 'cd %s\n' % out_dir
    cmd += 'rm -rf copy\n'
    cmd += 'rm -rf work_dir\n'
    cmd += 'rm -rf RunQC\n'
    cmd += 'rm -rf BuildHostIndex\n'
    cmd += 'rm -rf AlignReadsToHost\n'
    cmd += 'rm -rf NonHostReads\n'
    cmd += 'rm -rf RemoveHostDNA\n'
    cmd += 'rm -rf RunKraken\n'
    cmd += 'rm -rf KrakenResults\n'
    cmd += 'rm -rf FilteredKrakenResults\n'
    return cmd


def get_amrplusplus2(
        self,
        tech: str,
        fastqs: list,
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fastqs : list
        Paths to the input fastqs files
    """
    out_dir = genome_out_dir(self, tech, self.sam_pool)
    self.outputs['dirs'].append(out_dir)
    self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(out_dir)

    to_dos = status_update(self, tech, fastqs)

    # check if tool already run (or if --force) to allow getting command
    out = '%s/ResistomeResults/AMR_analytic_matrix.csv' % out_dir
    if self.config.force or to_do(out):
        # collect the command line
        cmd = amrplusplus2_cmd(self, fastqs, out_dir)
        key = genome_key(tech, self.sam_pool)
        if to_dos:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_f=fastqs, o_d=out_dir, key=key)
        self.soft.add_status(tech, self.sam_pool, 1)
    else:
        self.soft.add_status(tech, self.sam_pool, 0)


def amrplusplus2(self) -> None:
    """AmrPlusPlus is an easy to use app that identifies and characterizes
    resistance genes within sequence data.

    References
    ----------
    Lakin, Steven M., et al. "MEGARes: an antimicrobial resistance database
    for high throughput sequencing." Nucleic acids research 45.D1 (2017):
    D574-D580.

    Notes
    -----
    GitHub  : https://github.com/meglab-metagenomics/amrplusplus_v2
    Docs    : http://megares.meglab.org/amrplusplus/latest/html/index.html
    Paper   : https://doi.org/10.1093/nar/gkw1009

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if not self.sam_pool or self.sam_pool in self.pools:
        sys.exit('[%s] Only run on non-assembled reads' % self.soft.name)
    else:
        tech_fastqs = sample_inputs(self, raw=True)
        for tech, fastqs in tech_fastqs.items():
            get_amrplusplus2(self, tech, fastqs[''])


def metacompare_cmd(
        self,
        contigs: str,
        genes: str,
        out: str
) -> str:
    """Collect the command line for MetaCompare.

    Parameters
    ----------
    self
    contigs : str
        Path to the contigs file
    genes : str
        Path to the protein sequence predictions file
    out : str
        Path to the output file

    Returns
    -------
    cmd : str
        MetaCompare command line
    """
    cmd, cmd_rm = '', ''
    if contigs.endswith('.fa.gz') or contigs.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % contigs.rstrip('.gz')
        contigs = contigs.rstrip('.gz')

    if genes.endswith('.gz'):
        cmd += 'gunzip -c %s > %s\n' % (genes, genes.rstrip('.gz'))
        cmd_rm += 'rm %s\n' % genes.rstrip('.gz')
        genes = genes.rstrip('.gz')

    cmd += 'mkdir -p %s\n' % dirname(out)

    cmd += '%s/metacmp.py' % self.soft.params['path']
    cmd += ' -c %s' % contigs
    cmd += ' -g %s' % genes
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -v 1 > %s\n' % out.rstrip('.gz')

    cmd += cmd_rm
    cmd += 'gzip %s\n' % out.rstrip('.gz')
    return cmd


def get_metacompare(
        self,
        tech: str,
        proteins: dict,
        contigs: str,
        group: str) -> None:
    """

    Parameters
    ----------
    self
    tech : str
    contigs : str
        Path to the contigs file
    proteins : dict
    group : str
    """
    for genome, prodigal_dir in proteins.items():
        out_dir = genome_out_dir(self, tech, self.sam_pool)
        out_dir = out_dir + '/' + group
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)

        genes = '%s/nucleotide.sequences.fasta.gz' % prodigal_dir[0]
        to_dos = status_update(self, tech, [genes])
        to_dos.extend(status_update(self, tech, [contigs]))

        # check if tool already run (or if --force) to allow getting command
        out = '%s/output.txt.gz' % out_dir
        if self.config.force or to_do(out):
            cmd = metacompare_cmd(self, contigs, genes, out)
            key = genome_key(tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=[contigs, genes], o_f=out, key=key)
            self.soft.add_status(tech, self.sam_pool, 1)
        else:
            self.soft.add_status(tech, self.sam_pool, 0)


def metacompare(self) -> None:
    """MetaCompare is a computational pipeline for prioritizing resistome
    risk by estimating the potential for ARGs to be disseminated into human
    pathogens from a given environmental sample based on metagenomic
    sequencing data.

    References
    ----------
    Oh, M., Pruden, A., Chen, C., Heath, L.S., Xia, K. and Zhang, L.,
    2018. MetaCompare: a computational pipeline for prioritizing
    environmental resistome risk. FEMS microbiology ecology, 94(7), p.fiy079.

    Notes
    -----
    GitHub  : https://github.com/minoh0201/MetaCompare
    Paper   : https://doi.org/10.1093/femsec/fiy079

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MetaCompare
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.soft.prev not in ['prodigal']:
        sys.exit('[metacompare] Only after protein predictions')

    assembler, assembly = get_assembly(self)

    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            proteins = group_inputs(self, inputs)
            contigs = assembly[self.sam_pool][(tech, group)][0]
            get_metacompare(self, tech, proteins, contigs, group)


def abritamr_cmd(
        self,
        contigs: dict,
        out_dir: str,
) -> str:
    """Collect the command line for abriTAMR.

    Parameters
    ----------
    self
    contigs : dict
        Path(s) to input fasta(s) per co-assembly group
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        abriTAMR command
    """
    cmd, cmd_rm = '', ''
    contig_paths = {}
    for group, contig in contigs.items():
        if contig.endswith('.fa.gz') or contig.endswith('.fasta.gz'):
            cmd += 'gunzip -c %s > %s\n' % (contig, contig.rstrip('.gz'))
            contig = contig.rstrip('.gz')
            contig_paths[group] = contig
            cmd_rm += 'rm %s\n' % contig

    if self.soft.params['samples'] == 'all':
        txt = '%s/contigs.txt' % out_dir
        for cdx, (group, contig_path) in enumerate(contig_paths.items()):
            if cdx:
                cmd += 'echo -e "%s\\t%s" >> %s\n' % (group, contig_path, txt)
            else:
                cmd += 'echo -e "%s\\t%s" > %s\n' % (group, contig_path, txt)
        if cmd:
            cmd += 'envsubst < %s > %s.tmp\n' % (txt, txt)
            cmd += 'mv %s.tmp %s\n' % (txt, txt)

    cmd += 'cd %s\n' % out_dir
    cmd += 'abritamr run'
    if self.soft.params['samples'] == 'all':
        cmd += ' --contigs %s' % txt
    else:
        cmd += ' --contigs %s' % contig
    cmd += ' --identity %s' % self.soft.params['identity']
    if self.soft.params['species']:
        cmd += ' --species %s' % self.soft.params['species']
    cmd += ' --jobs %s\n' % self.soft.params['cpus']

    # cmd += 'abritamr report'
    # cmd += ' --qc %s'
    # cmd += ' --runid %s'
    # cmd += ' --matches summary_matches.txt'
    # cmd += ' --partials summary_partials.txt'
    # cmd += ' --sop %s' % self.soft.params['sop']
    # cmd += ' --sop_name %s\n'

    cmd += 'gzip -q *.txt\n'
    cmd += 'gzip -q */amrfinder.out\n'
    cmd += cmd_rm
    return cmd


def abritamr_all(self):
    """Run abritamr on all samples contigs.

    Parameters
    ----------
    self : Commands class instance
    """
    for pool, group_sams in self.pools.items():
        for tech in set([x[0] for x in self.inputs[pool]]):
            out_dir = '%s/%s/%s' % (self.dir, tech, pool)
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'][pool] = out_dir

            contigs = get_contigs(self, tech, pool)
            contigs_list = list(contigs.values())
            to_dos = status_update(self, tech, contigs_list)

            out = '%s/summary_matches.txt.gz' % out_dir
            if self.config.force or to_do(out):
                cmd = abritamr_cmd(self, contigs, out_dir)
                key = (tech, pool)
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=contigs_list, o_d=out_dir, key=key)
                self.soft.add_status(tech, pool, 1)
            else:
                self.soft.add_status(tech, pool, 0)


def abritamr_sample(self):
    """Run abritamr on a single sample.

    Parameters
    ----------
    self : Commands class instance
    """
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir

        contigs = {group: inputs[0]}
        to_dos = status_update(self, tech, [inputs[0]], group=group)

        out = '%s/summary_matches.txt.gz' % out_dir
        if self.config.force or to_do(out):
            cmd = abritamr_cmd(self, contigs, out_dir)
            key = (tech, group)
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=inputs[0], o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def abritamr(self) -> None:
    """abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on
    a single (or list ) of given isolates and collates the results into a
    table, separating genes identified into functionally relevant groups.

    References
    ----------
    Sherry, N.L., Horan, K.A., Ballard, S.A., Gonҫalves da Silva, A., Gorrie,
    C.L., Schultz, M.B., Stevens, K., Valcanis, M., Sait, M.L., Stinear,
    T.P. and Howden, B.P., 2023. An ISO-certified genomics workflow for
    identification and surveillance of antimicrobial resistance. Nature
    Communications, 14(1), pp.1-12.

    Notes
    -----
    GitHub  : https://github.com/MDU-PHL/abritamr
    Paper   : https://doi.org/10.1038/s41467-022-35713-4
    Docs

    Parameters
    ----------
    self : Commands class instance
    """
    if self.config.tools[self.soft.prev] != 'assembling':
        sys.exit('[tiara] can only be run on assembly output')

    if self.soft.params['samples'] == 'all':
        abritamr_all(self)
    else:
        abritamr_sample(self)


def staramr_cmd(
        self,
        tech,
        contigs,
        out_dir,
) -> str:
    """Get the staramr command line

    Parameters
    ----------
    self
    tech : str
        Tachnology
    contigs : str
        Path to the input contigs file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        staramr command line
    """
    params = tech_params(self, tech)

    cmd, cmd_rm = '', ''
    if contigs.endswith('.fa.gz') or contigs.endswith('.fasta.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        contigs = contigs.rstrip('.gz')
        cmd_rm += 'rm %s\n' % contigs

    cmd += 'rm -rf %s\n' % out_dir
    cmd += 'staramr search'
    for param in [ 'pointfinder_organism', 'plasmidfinder_database_type',
                   'mlst_scheme', 'exclude_genes_file']:
        if params[param]:
            cmd += ' --%s %s' % (param.replace('_', '-'), params[param])

    for boolean in ['ignore_invalid_files', 'no_exclude_gene',
                    'exclude_negatives', 'exclude_resistance_phenotypes',
                    'report_all_blast']:
        if params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')

    for param in ['genome_size_lower_bound', 'genome_size_upper_bound',
                  'minimum_N50_value', 'minimum_contig_length',
                  'unacceptable_number_contigs', 'pid_threshold',
                  'percent_length_overlap_resfinder',
                  'percent_length_overlap_pointfinder',
                  'percent_length_overlap_plasmidfinder']:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    cmd += ' --nprocs %s' % params['cpus']
    cmd += ' --output-dir %s' % out_dir
    cmd += ' %s\n' % contigs

    cmd += cmd_rm

    return cmd


def get_staramr(self, tech, folders, group):
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
        Name of the current sample or co-assembly group
    """
    for genome, inputs in folders.items():
        out_dir = genome_out_dir(self, tech, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'].setdefault((tech, group), []).append(out_dir)
        to_dos = status_update(
            self, tech, [inputs[0]], self.sam_pool, group=group)

        key = genome_key(tech, group)

        out_fp = '%s/summary.tsv' % out_dir
        if self.config.force or to_do(out_fp):
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                cmd = staramr_cmd(self, tech, inputs[0], out_dir)
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=inputs[0], o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def staramr(self) -> None:
    """staramr (*AMR) scans bacterial genome contigs against the ResFinder,
    PointFinder, and PlasmidFinder databases (used by the ResFinder
    webservice and other webservices offered by the Center for Genomic
    Epidemiology) and compiles a summary report of detected antimicrobial
    resistance genes. The star|* in staramr indicates that it can handle all
    of the ResFinder, PointFinder, and PlasmidFinder databases.

    References
    ----------
    Bharat, A., Petkau, A., Avery, B.P., Chen, J.C., Folster, J.P., Carson,
    C.A., Kearney, A., Nadon, C., Mabon, P., Thiessen, J. and Alexander,
    D.C., 2022. Correlation between Phenotypic and In Silico Detection of
    Antimicrobial Resistance in Salmonella enterica in Canada Using Staramr.
    Microorganisms, 10(2), p.292.

    Notes
    -----
    GitHub  : https://github.com/phac-nml/staramr
    Paper   : https://doi.org/10.3390/microorganisms10020292

    Parameters
    ----------
    self : Commands class instance
    """
    if self.sam_pool in self.pools:
        for (tech, group), inputs in self.inputs[self.sam_pool].items():
            folders = group_inputs(self, inputs)
            get_staramr(self, tech, folders, group)
    else:
        sys.exit('[staramr] Only work after assembly')


def amrfinderplus(self) -> None:
    """This software and the accompanying database are designed to find
    acquired antimicrobial resistance genes and some point mutations in
    protein or assembled nucleotide sequences. We have also added "plus"
    stress, head, and biocide resistance as well as some virulence factors
    and E. coli antigens.

    References
    ----------
    Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J.G., Haendiges,
    J., Haft, D.H., Hoffmann, M., Pettengill, J.B., Prasad, A.B., Tillman,
    G.E. and Tyson, G.H., 2021. AMRFinderPlus and the Reference Gene Catalog
    facilitate examination of the genomic links among antimicrobial
    resistance, stress response, and virulence. Scientific reports, 11(1),
    pp.1-9.

    Notes
    -----
    GitHub  : https://github.com/ncbi/amr
    Docs    : https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
    Wiki    : https://github.com/ncbi/amr/wiki
    Paper   : https://doi.org/10.1038/s41598-021-91456-0

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for AMRFinderPlus
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass


def ariba(self) -> None:
    """

    References
    ----------

    Notes
    -----
    GitHub  : https://github.com/sanger-pathogens/ariba

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def argsoap(self) -> None:
    """Fast annotation and classification of antibiotic resistance gene-like
    sequences from metagenomic data.

    References
    ----------
    Yin, X., Jiang, X.T., Chai, B., Li, L., Yang, Y., Cole, J.R., Tiedje,
    J.M. and Zhang, T., 2018. ARGs-OAP v2. 0 with an expanded SARG database
    and Hidden Markov Models for enhancement characterization and
    quantification of antibiotic resistance genes in environmental
    metagenomes. Bioinformatics, 34(13), pp.2263-2270.

    Notes
    -----
    GitHub  : https://github.com/biofuture/Ublastx_stageone
    Docs    : https://galaxyproject.org/use/args-oap/
    Paper   : https://doi.org/10.1093/bioinformatics/bty053

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for ARGs-OAP
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass


def resfinder(self) -> None:
    """ResFinder identifies acquired antimicrobial resistance genes in total
    or partial sequenced isolates of bacteria.

    References
    ----------
    Bortolaia, V., Kaas, R.S., Ruppe, E., Roberts, M.C., Schwarz, S.,
    Cattoir, V., Philippon, A., Allesoe, R.L., Rebelo, A.R., Florensa,
    A.F. and Fagelhauer, L., 2020. ResFinder 4.0 for predictions of
    phenotypes from genotypes. Journal of Antimicrobial Chemotherapy, 75(12),
    pp.3491-3500.

    Notes
    -----
    BitBucket : https://bitbucket.org/genomicepidemiology/resfinder/src/master/
    Paper     : https://doi.org/10.1093/jac/dkaa345

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for ResFinder
        .sam_pool : str
            Sample of co-assembly group name
        .pools : dict
            Co-assembly pools and sample per group
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.name : str
            Name of the current software in the pipeline
        .soft.prev : str
            Name of the previous software in the pipeline
        .soft.status
            Current status of the pipeline in terms of available outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    pass
