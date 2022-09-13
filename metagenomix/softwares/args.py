# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename, dirname, splitext
from metagenomix._inputs import (sample_inputs, group_inputs, genome_key,
                                 genome_out_dir)
from metagenomix._io_utils import (caller, io_update, to_do, tech_specificity,
                                   not_paired, status_update)


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
        Sample name or name of a co-assembly pool's group
    """
    for genome, fasta in fastas.items():

        out_dir = genome_out_dir(self, tech, fasta[0], sam_group, genome)
        self.outputs['dirs'].append(out_dir)

        typ_seqs = predict_inputs(self, fasta[0])
        for typ, seq in typ_seqs.items():

            key = genome_key(tech, sam_group, genome)
            status_update(self, tech, [seq], group=sam_group, genome=genome)

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
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=fasta[0], o_d=out_dir, key=key)
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

    elif self.soft.prev == 'drep':
        for (tech, bin_algo), inputs in self.inputs[''].items():
            fastas = group_inputs(self, inputs)
            get_predict(self, fastas, tech, bin_algo)
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
    # iterate over the inputs
    for (tech, sam), fastqs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, fastqs, tech, sam, ['illumina']):
            continue
        if not_paired(self, tech, sam, fastqs):
            continue
        status_update(self, tech, fastqs)

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
    cmd = '\n./mmarc'
    if genome:
        cmd += ' --input %s' % inputs[0]
    else:
        if len(inputs) == 1:
            cmd += ' --input %s' % inputs[0]
        if len(inputs) == 2:
            cmd += ' --i1 %s --i2 %s' % tuple(inputs)
        if tech == 'illumina':
            cmd += ' --dedup'
        cmd += ' --multicorrect'
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --output %s' % out
    cmd += ' --filename output'
    cmd += ' --skewness %s/skewness.txt' % out
    cmd += ' --graphs %s/graphs' % out
    cmd += ' --level %s' % level

    for boolean in ['dedup', 'multicorrect']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean

    for param in ['coverage', 'evalue', 'kmer']:
        cmd += ' --%s %s' % (param.replace('_', '-'), self.soft.params[param])

    return cmd


def get_metamarc(
        self,
        tech: str,
        fastas_folders: dict,
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
    fastas_folders : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Name of the current sample or co-assembly group
    """
    for genome, inputs in fastas_folders.items():

        fasta = inputs[0]
        out_dir = genome_out_dir(self, tech, fasta, sam_group, genome)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)
        status_update(
            self, tech, self.sam_pool, inputs, group=sam_group, genome=genome)

        outs, cmd = [], ''
        key = genome_key(tech, sam_group, genome)
        for level in self.soft.params['level']:
            out = '%s/model_level_%s' % (out_dir, level)
            self.outputs['dirs'].append(out_dir)
            outs.append(out)
            skewness = '%s/skewness.txt' % out
            # check if tool already run (or if --force) to allow getting command
            if self.config.force or to_do(skewness):
                # collect the command line
                cmd += mmarc_cmd(self, tech, inputs, out, genome, level)

        if cmd:
            path = self.soft.params['path']
            full_cmd = 'cd %s\n' % out_dir
            full_cmd += 'cp -r %s %s/.\n' % (path, out_dir)
            full_cmd += 'cp -r %s/../hmmer-3.1b2 %s/.\n' % (path, out_dir)
            full_cmd += 'cp -r %s/../hmmer-3.1b2 %s/.\n' % (path, out_dir)
            full_cmd += 'cd bin\n'
            full_cmd += cmd
            full_cmd += '\nrm -rf bin\n'
            full_cmd += 'rm -rf hmmer-3.1b2\n'
            # add is to the 'cmds'
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
            fastas_folders = group_inputs(self, inputs)
            get_metamarc(self, tech, fastas_folders, group)

    elif self.soft.prev == 'drep':
        for (tech, bin_algo), inputs in self.inputs[''].items():
            folders = group_inputs(self, inputs)
            get_metamarc(self, tech, folders, bin_algo)

    else:
        tech_fastas = sample_inputs(self, raw=True)
        for tech, fastas in tech_fastas.items():
            get_metamarc(self, tech, fastas, self.sam_pool)


def karga_kargva_cmd(
        self,
        merged: str,
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
    merged : str
        Path to the input file
    db_path : str
        Path to ARG/MGE fasta file with resistance annotation/mutation in header

    Returns
    -------
    cmd : str
        KARGA or KARGVA command
    """
    cmd = 'java %s' % self.soft.name.upper()
    cmd += ' f:%s' % merged
    cmd += ' d:%s' % db_path
    for param in ['k', 'i', 's']:
        if param == 's' and self.soft.name != 'karga':
            continue
        if self.soft.params[param]:
            cmd += ' %s:%s' % (param, self.soft.params[param])
    cmd += ' -Xmx16GB\n'
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
        out_dir = genome_out_dir(self, tech, fastx, sam_group)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_dir)
        status_update(self, tech, inputs, group=sam_group, genome=genome)

        merge_cmd = ''
        if len(inputs) == 1:
            merged = fastx
        elif len(inputs) > 1:
            merged = '%s/%s' % (out_dir, basename(fastx))
            merge_cmd += 'cat %s > %s\n' % (' '.join(inputs), merged)

        outs, cmd = [], ''
        for db, db_path in self.soft.params['databases'].items():
            out = '%s/db_%s' % (out_dir, db)
            self.outputs['dirs'].append(out_dir)
            outs.append(out)
            bas = splitext(basename(merged))[0]
            fp = '%s/%s_%s_mappedReads.csv' % (out, self.soft.name.upper(), bas)
            # check if tool already run (or if --force) to allow getting command
            if self.config.force or to_do(fp):
                # collect the command line
                cmd += karga_kargva_cmd(self, merged, db_path)

        if cmd:
            if self.soft.name == 'karga':
                ps = ['openjdk-8', 'AMRGene.class', 'KARGA$1.class',
                      'KARGA.class', 'KARGA.java']
            elif self.soft.name == 'kargva':
                ps = ['AMRVariantGene.class', 'KARGVA$1.class', 'KARGVA.class',
                      'KARGVA.java', 'kargva_db_v5.fasta']
            path = self.soft.params['path']
            full_cmd = 'cd %s\n' % out_dir
            for p in ps:
                full_cmd += 'cp -r %s/%s %s/.\n' % (path, p, out_dir)
            full_cmd += merge_cmd + cmd
            for p in ps:
                full_cmd += 'rm -rf %s\n' % p

            # add is to the 'cmds'
            key = genome_key(tech, sam_group)
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
    inputs : list
        Paths to the input files
    out_d : str
        Paths to the output folder

    Returns
    -------
    cmd : str
        abricate command
    """
    cmd, seqtk_cmd = '', ''
    for db in self.soft.params['databases']:

        cmd += '\nabricate'

        if tech in ['nanopore', 'pacbio']:
            outs = []
            for idx, inp in enumerate(inputs):
                out = '%s/%s.fa' % (out_d, idx)
                outs.append(out)
                seqtk_cmd += 'seqtk seq -A %s > %s/%s\n' % (inp, out_d, out)
            out = '%.fasta' % splitext(basename(inputs[0]))[0]
            seqtk_cmd += 'cat %s > %s/%s\n' % (' '.join(outs), out_d, out)
        else:
            out = inputs[0]

        cmd += ' %s' % out
        cmd += ' --threads %s' % self.soft.params['cpus']
        cmd += ' --db %s' % db
        for boolean in ['setupdb', 'noheader', 'csv', 'nopath']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean
        for param in ['minid', 'mincov']:
            cmd += ' --%s %s' % (param, self.soft.params[param])
        cmd += ' > %s/%s.txt' % (out_d, db)

    if seqtk_cmd:
        cmd = seqtk_cmd + cmd

    return cmd


def get_abricate(
        self,
        tech: str,
        fastas_folders: dict,
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
    fastas_folders : dict
        Paths to the input fasta files per genome/MAG
    sam_group : str
        Name of the current sample or co-assembly group
    """
    for genome, inputs in fastas_folders.items():

        fasta = inputs[0]
        out_d = genome_out_dir(self, tech, fasta, sam_group, genome)
        self.outputs['dirs'].append(out_d)
        self.outputs['outs'].setdefault((tech, sam_group), []).append(out_d)

        status_update(self, tech, [fasta], group=sam_group, genome=genome)

        outs = ['%s/%s.txt' % (out_d, x) for x in self.soft.params['databases']]
        # check if tool already run (or if --force) to allow getting command
        if self.config.force or sum([to_do(x) for x in outs]):
            # collect the command line
            cmd = abricate_cmd(self, tech, inputs, out_d)
            # add is to the 'cmds'
            key = genome_key(tech, sam_group, genome)
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
            fastas_folders = group_inputs(self, inputs)
            get_abricate(self, tech, fastas_folders, group)

    elif self.soft.prev == 'drep':
        for (tech, bin_algo), inputs in self.inputs[''].items():
            folders = group_inputs(self, inputs)
            get_abricate(self, tech, folders, bin_algo)

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
    if len(fastqs) > 1:
        cmd += ' --reads "%s"' % fastqs[0].replace('R1', '{R1,R2}')
    else:
        cmd += ' --reads %s' % fastqs[0]
    for path in ['adapters', 'fqc_adapters', 'host_index', 'host',
                 'kraken_db', 'amr_index', 'amr', 'annotation',
                 'snp_annotation', 'snp_confirmation']:
        if self.soft.params[path]:
            cmd += ' --%s "%s"' % (path, self.soft.params[path])
    for param in ['leading', 'trailing', 'slidingwindow', 'minlen',
                  'threshold', 'min', 'max', 'skip', 'samples']:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --output "%s"' % out_dir
    cmd += ' -work "%s/work_dir"\n' % out_dir
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
        sam: str
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
    sam : str
        Name of the current sample
    """
    out_dir = genome_out_dir(self, tech, fastqs[0], sam)
    self.outputs['dirs'].append(out_dir)
    self.outputs['outs'].setdefault((tech, sam), []).append(out_dir)

    status_update(self, tech, fastqs)

    # check if tool already run (or if --force) to allow getting command
    out = '%s/ResistomeResults/AMR_analytic_matrix.csv' % out_dir
    if self.config.force or to_do(out):
        # collect the command line
        cmd = amrplusplus2_cmd(self, fastqs, out_dir)
        # add is to the 'cmds'
        key = genome_key(tech, sam)
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_f=fastqs, o_d=out_dir, o_f=out, key=key)
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
            get_amrplusplus2(self, tech, fastqs[''], self.sam_pool)
