# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob
import pkg_resources
from os.path import basename, dirname, isdir, isfile, splitext
from metagenomix._io_utils import (
    min_nlines, io_update, to_do, tech_specificity, status_update)
from metagenomix.softwares.alignment import get_alignment_cmd
from metagenomix.core.parameters import tech_params

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources/scripts')


def shogun_append_cmd(
        self,
        tech: str,
        cmds: list,
        tab: str
) -> None:
    """Add or not the current SHOGUN command to the list of commands to run
    depending on whether the file already exists or exists but is only a header.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    cmds : list
        List of current SHOGUN command lines
    tab : str
        Path to the output table of the current command
    """
    if self.config.force:
        self.outputs['cmds'].setdefault((tech,), []).extend(cmds)
    elif to_do(tab) or not min_nlines(tab):
        cmd = 'file="%s"\n' % tab
        cmd += 'if [ -f "$file" ]\n'
        cmd += 'then\n'
        cmd += "FILESIZE=`ls -l $file | cut -d ' ' -f5`\n"
        cmd += "if [[ $FILESIZE -lt 1000 ]]\n"
        cmd += 'then\n'
        cmd += '\n'.join(cmds)
        cmd += 'fi\n'
        self.outputs['cmds'].setdefault((tech,), []).append(cmd)


def shogun_redistribute(
        self,
        tech: str,
        db: str,
        aligner: str,
        tax_norm: str,
        db_path: str
) -> None:
    """Get the SHOGUN command to redistribute a taxonomic profile
    at various taxonomic levels.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    db : str
        Database name
    aligner : str
        Name of the aligner
    tax_norm : str
        Path to the input taxonomic table file
    db_path : str
        Path to the taxonomic database
    """
    redist_cmds = []
    for level in ['phylum', 'genus', 'strain']:
        redist = '%s.redist.%s.tsv' % (splitext(tax_norm)[0], level)
        cmd = 'shogun redistribute'
        cmd += ' -i %s' % tax_norm
        cmd += ' -d %s' % db_path
        cmd += ' -l %s' % level
        cmd += ' -o %s' % redist
        redist_cmds.append(cmd)
        self.outputs['outs'][(tech, self.sam_pool)].setdefault(
            (db, aligner, tech), []).append(redist)
        io_update(self, o_f=redist, key=tech)
    redist_out = '%s.redist.strain.tsv' % splitext(tax_norm)[0]
    shogun_append_cmd(self, tech, redist_cmds, redist_out)


def shogun_assign_taxonomy(
        aligner: str,
        ali: str,
        tax: str,
        db_path: str,
        sub_db: str
) -> str:
    """Get the assign taxonomy command.

    Parameters
    ----------
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    tax : str
        Path to the output table
    db_path : str
        Path to the database
    sub_db : str
        Suffix to the database

    Returns
    -------
    cmd : str
        assign taxonomy command
    """
    cmd = 'shogun assign_taxonomy'
    cmd += ' -a %s' % aligner
    cmd += ' -i %s' % ali
    cmd += ' -d %s%s' % (db_path, sub_db)
    cmd += ' -o %s\n' % tax
    return cmd


def shogun_normalize(
        input_fp: str,
        output_fp: str
) -> str:
    """Get the SHOGUN command to normalize table.

    Parameters
    ----------
    input_fp : str
        Path to the input table to normalize.
    output_fp : str
        Path to the output normalized table.

    Returns
    -------
    cmd : str
        Command line for table normalization.
    """
    cmd = 'shogun normalize -i %s -o %s\n' % (input_fp, output_fp)
    return cmd


def align_cmd(
        params: dict,
        fasta: str,
        out: str,
        db_path: str,
        aligner: str
) -> str:
    """Get the SHOGUN alignment command.

    Parameters
    ----------
    params : dict
        Parameters for the current technology
    fasta : str
        Path to the input fasta file
    db_path : str
        Path to the shogun database
    out : str
        Path to the output folder
    aligner : str
        Name of the aligner

    Returns
    -------
    cmd : str
        Alignment command using SHOGUN aligner
    """
    cmd = 'shogun align -a %s' % aligner
    cmd += ' -i %s' % fasta
    cmd += ' -d %s' % db_path
    cmd += ' -t %s' % params['cpus']
    cmd += ' -o %s' % out
    return cmd


def get_alignment_basename(
        aligner: str
) -> str:
    """Get the default basename for SHOGUN-output alignments.

    Parameters
    ----------
    aligner : str
        Name of the Aligner

    Returns
    -------
    ali_base : str
        Basename for the output alignment file
    """
    if aligner == 'bowtie2':
        ali_base = 'alignment.bowtie2.sam'
    elif aligner == 'burst':
        ali_base = 'alignment.burst.b6'
    else:
        raise ValueError('No output basename for aligner "%s"' % aligner)
    return ali_base


def get_orients(
        inputs: list
) -> list:
    """Return the reads orientations.

    Parameters
    ----------
    inputs : list
        Path to the input fastx files

    Returns
    -------
    orients : list
        Reads orientations
    """
    if len(inputs) == 1:
        orients = ['']
    elif len(inputs) == 2:
        orients = ['1', '2']
    elif len(inputs) == 3:
        orients = ['', '1', '2']
    else:
        raise IOError('Input to shogun must be 1, 2 or 3 fasta files')
    return orients


def get_combine_cmd(
        self,
        tech: str,
        inputs: list,
        out: str,
        combine_cmds: list
) -> list:
    """

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Path to the input files
    out : str
        Path to pipeline output folder for SHOGUN
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files

    Returns
    -------
    fastas : list
        Paths to fasta files to combine for SHOGUN
    """
    fastas = []
    # get the fastq versions of input file
    orients = list(get_orients(inputs))
    for pdx, path_ in enumerate(inputs):
        path = path_
        if path_.endswith('fastq') or path_.endswith('fastq.gz'):
            # replace non-fasta by fasta extensions
            path = path_.replace('.fastq', '.fasta').replace('.gz', '')
            # Prepare the extraction command and collect it
            if to_do(path):
                io_update(self, i_f=path_, o_f=path, key=tech)
                to_fasta_cmd = 'seqtk seq -A %s > %s' % (path_, path)
                combine_cmds.append(to_fasta_cmd)
                self.soft.add_status(tech, self.sam_pool, 1)
            else:
                self.soft.add_status(tech, self.sam_pool, 0)

        path_out = '%s/%s' % (out, basename(path))
        edit_fasta = '%s/fasta4shogun.py -i %s -o %s -s %s' % (
            RESOURCES, path, path_out, self.sam_pool)
        orient = orients[pdx]
        if orient:
            edit_fasta += ' -r %s' % orient
        combine_cmds.append(edit_fasta)
        fastas.append(path_out)

    return fastas


def combine_inputs(
        self,
        tech: str,
        inputs: list,
        out: str,
        combine_cmds: list
) -> str:
    """Combine the fastq/a of the sample to run shogun.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Path to the input files
    out : str
        Path to pipeline output folder for SHOGUN
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files

    Returns
    -------
    fasta : str
        Path to the combined sequences fasta file
    """
    fastas = get_combine_cmd(self, tech, inputs, out, combine_cmds)
    fasta = '%s/combined.fasta' % out
    combine_cmds.extend([
        'cat %s > %s' % (' '.join(fastas), fasta), 'rm %s' % ' '.join(fastas)])
    return fasta


def align(
        self,
        tech: str,
        fasta: str,
        out: str,
        db: str,
        aligner: str,
        ali_cmds: list
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fasta : str
        Path to the combined sequences fasta file
    out : str
        Path to pipeline output folder for SHOGUN
    db : str
        Database name
    aligner : str
    ali_cmds : list

    Returns
    -------
    ali : str
        Alignment output file
    """
    out_dir = '%s/%s/%s' % (out, db, aligner)
    self.outputs['dirs'].append(out_dir)

    ali = '%s/%s' % (out_dir, get_alignment_basename(aligner))
    params = tech_params(self, tech)
    db_path = '%s/shogun' % self.databases.paths[db]
    cmd = align_cmd(params, fasta, out_dir, db_path, aligner)
    cmd = get_alignment_cmd([fasta], cmd, ali)
    ali_cmds.append(cmd)
    if not cmd.startswith('shogun'):
        io_update(self, i_f=ali, key=tech)
    self.outputs['outs'][(tech, self.sam_pool)].setdefault(
        (db, aligner), []).append(ali)
    io_update(self, o_d=out_dir, key=tech)

    return ali


def format_sam(
        sam_: str,
        ali_cmds: list,
        sample: str
) -> str:
    sam = '%s_formatted.sam' % splitext(sam_)[0]
    cmd = '%s/sam4shogun.py -i %s -o %s -s %s' % (
        RESOURCES, sam_, sam, sample)
    ali_cmds.append(cmd)
    return sam


def get_dir(
        out_dir: str,
        db: str,
        aligner: str
) -> str:
    """

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for SHOGUN
    db : str
        Database name
    aligner : str
        Name of the aligner

    Returns
    -------
    out : str
        Path to pipeline output folder for SHOGUN
    """
    out = '%s/%s' % (out_dir, db)
    if aligner:
        # i.e., if the previous software is not bowtie2
        out += '/%s' % aligner
    return out


def shogun_assign_normalize(
        self,
        tech: str,
        tab: str,
        aligner: str,
        ali: str,
        db_path: str,
        sub_db: str
) -> None:
    """Get the assignment and the normalization SHOGUN commands for
    the taxonomy or functional classifications.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    tab : str
        Path to the output table
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    db_path : str
        Path to the database
    sub_db : str
        Suffix to the database
    """
    cmd = ''
    tab_norm = '%s_norm.tsv' % splitext(tab)[0]
    if self.config.force or to_do(tab):
        cmd += shogun_assign_taxonomy(aligner, ali, tab, db_path, sub_db)
    if self.config.force or to_do(tab_norm):
        cmd += shogun_normalize(tab, tab_norm)
    if cmd:
        shogun_append_cmd(self, tech, list([cmd]), tab_norm)


def get_paths(
        self,
        tech: str,
        out: str,
        aligner: str,
        db: str,
        tax_fun: str
) -> tuple:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    out : str
        Path to pipeline output folder for SHOGUN
    aligner : str
        Name of the aligner
    db : str
        Database name
    tax_fun : str
        "taxonomy" or "function"

    Returns
    -------
    tax : str
        Path to the output table
    norm : str
        Path to the output table (relative abundances)
    """
    out_dir = get_dir(out, db, aligner)
    self.outputs['dirs'].append(out_dir)

    tab = '%s/%s.tsv' % (out_dir, tax_fun)
    norm = '%s/%s_norm.tsv' % (out_dir, tax_fun)
    self.outputs['outs'][(tech, self.sam_pool)].setdefault(
        (db, aligner), []).extend([tab, norm])
    io_update(self, i_d=out_dir, o_f=[tab, norm], key=tech)
    return tab, norm


def shogun_coverage_cmd(
        self,
        tech: str,
        ali: str,
        db_path: str,
        level: str,
        cov_tab: str
) -> None:
    """Get the SHOGUN command that calculates the coverage per taxon.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    ali : str
        Path to the alignment
    db_path : str
        Path to the taxonomic database
    level : str
        Taxonomic level for coverage calculation
    cov_tab : str
        Path to the output coverage table file
    """
    cmd = 'shogun coverage'
    cmd += ' -i %s' % ali
    cmd += ' -d %s' % db_path
    cmd += ' -l %s' % level
    cmd += ' -o %s' % cov_tab
    shogun_append_cmd(self, tech, list([cmd]), cov_tab)


def shogun_coverage(
        self,
        tech: str,
        out: str,
        aligner: str,
        ali: str,
        db: str
) -> None:
    db_path = '%s/shogun' % self.databases.paths[db]
    if aligner == 'burst':
        cov_tab = '%s/coverage.tsv' % out
        self.outputs['outs'][(tech, self.sam_pool)].setdefault(
            (db, aligner), []).append(cov_tab)
        if self.config.force or to_do(cov_tab):
            shogun_coverage_cmd(self, tech, ali, db_path, 'strain', cov_tab)
            io_update(self, o_f=cov_tab, key=tech)


def shogun_taxonomy(
        self,
        tech: str,
        out: str,
        aligner: str,
        ali: str,
        db: str
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    out : str
        Path to pipeline output folder for SHOGUN
    aligner : str
        Name of the aligner
    ali : str
        Path to the alignment
    db : str
        Database name

    Returns
    -------
    norm : str
        Path to the normalized output table
    """
    db_path = '%s/shogun' % self.databases.paths[db]
    tab, norm = get_paths(self, tech, out, aligner, db, 'taxonomy')
    shogun_assign_normalize(self, tech, tab, aligner, ali, db_path, '')
    shogun_redistribute(self, tech, db, aligner, norm, db_path)
    return norm


def shogun_function(
        self,
        tech: str,
        out: str,
        aligner: str,
        ali: str,
        db: str,
        norm: str
) -> None:
    db_path = '%s/shogun' % self.databases.paths[db]
    for sub in ['kegg', 'refseq', 'uniprot']:
        fun, f_norm = get_paths(
            self, tech, out, aligner, db, 'function-%s' % sub)
        shogun_assign_normalize(
            self, tech, fun, aligner, ali, db_path, '/functions-%s' % sub)
        shogun_redistribute(self, tech, db, aligner, f_norm, db_path)

    out_dir = out + '/functional'
    for level in ['genus', 'species']:
        shogun_functional(self, tech, norm, fun, out_dir, level)


def shogun_functional(
        self,
        tech: str,
        norm: str,
        db_path: str,
        out_dir: str,
        level: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    norm : str
        Path to the input table file
    db_path : str
        Path to the functional database
    out_dir : str
        Path to the output functional table
    level : str
        Taxonomic level
    """
    cmd = 'shogun functional'
    cmd += ' -i %s' % norm
    cmd += ' -o %s' % out_dir
    cmd += ' -d %s' % db_path
    cmd += ' -l %s' % level
    base = splitext(basename(norm))[0]
    out = '%s/%s.%s.normalized.txt' % (out_dir, base, level)
    shogun_append_cmd(self, tech, list([cmd]), out)


def shogun(self) -> None:
    """Shallow seq pipeline for optimal shotgun data usage.
    Get the SHOGUN commands that consist of the classifications using the
    databases built for shogun and the alignment in shogun, or not as this
    may be run after bowtie2, in which case the classifications are performed
    on the alignment obtained before, provided that it was done on a database
    formatted for SHOGUN.

    References
    ----------
    Hillmann, Benjamin, et al. "SHOGUN: a modular, accurate and scalable
    framework for microbiome quantification." Bioinformatics 36.13 (2020):
    4088-4090.

    Notes
    -----
    GitHub  : https://github.com/knights-lab/SHOGUN
    Paper   : https://doi.org/10.1093/bioinformatics/btaa277

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
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam):
            continue

        self.outputs['outs'][(tech, self.sam_pool)] = {}
        params = tech_params(self, tech)
        combine_cmds, ali_cmds = [], []

        out = '/'.join([self.dir, tech, self.sam_pool])
        self.outputs['dirs'].append(out)
        io_update(self, o_d=out, key=tech)

        if self.soft.prev == 'bowtie2':
            to_dos = status_update(self, tech, list(inputs.values()))
            for (db, aligner), sam in inputs.items():
                io_update(self, i_f=sam, key=tech)
                ali = format_sam(sam, ali_cmds, self.sam_pool)
                self.outputs['outs'][(tech, self.sam_pool)].setdefault(
                    (db, 'bowtie2'), []).append(ali)
                shogun_taxonomy(self, tech, out, 'bowtie2', ali, db)
        elif params['databases']:
            io_update(self, i_f=inputs, key=tech)
            to_dos = status_update(self, tech, inputs)
            fasta = combine_inputs(self, tech, inputs, out, combine_cmds)
            for db, aligners in params['databases'].items():
                for aligner in aligners:
                    ali = align(self, tech, fasta, out, db, aligner, ali_cmds)
                    norm = shogun_taxonomy(self, tech, out, aligner, ali, db)
                    shogun_coverage(self, tech, out, aligner, ali, db)
                    shogun_function(self, tech, out, aligner, ali, db, norm)
        if self.outputs['cmds']:
            cmd = combine_cmds + ali_cmds + self.outputs['cmds'][(tech,)]
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = cmd
            self.soft.add_status(tech, sam, 1)
        else:
            self.soft.add_status(tech, sam, 0)


def woltka_aligments(
        self,
        tech: str
) -> dict:
    """Get the .bam alignment paths per sample.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'

    Returns
    -------
    alignments : dict
        Alignments per sample and per pairing
    """
    alignments = {}
    for sample, sam_inputs in self.inputs.items():
        if sam_inputs.get((tech, sample)):
            for (db, aligner), bam in sam_inputs[(tech, sample)].items():
                if bam and db == 'wol':
                    if aligner not in alignments:
                        alignments[aligner] = {}
                    alignments[aligner][sample] = bam
    return alignments


def woltka_map(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        alis: dict
) -> tuple:
    """Write the mapping file that servers as input to Woltka.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to the output folder
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pairing : str
        Type of reads pairing during alignment
    aligner : str
        Aligner used to create the input alignments
    alis : dict
        Alignments

    Returns
    -------
    map_fp : str
        Path to the output woltka samples file
    cmds : list
        Commands to create the samples.map file
    rms : list
        Commands to remove the .sam files
    """
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    self.outputs['dirs'].append(out_dir)

    sam_cmds, rms = '', ''
    for idx, sample in enumerate(alis.keys()):
        bam = alis[sample]
        sam_cmds += 'samtools view %s > %s.sam\n' % (bam, bam.rstrip('.bam'))
        rms += 'rm %s.sam\n' % bam.rstrip('.bam')

    cmds = ''
    map_fp = '%s/samples.map' % out_dir
    for idx, sample in enumerate(alis.keys()):
        echo = 'echo -e "%s\\t%s.sam"' % (sample, alis[sample].rstrip('.bam'))
        if idx:
            cmds += '%s >> %s\n' % (echo, map_fp)
        else:
            cmds += '%s > %s\n' % (echo, map_fp)

    if cmds:
        cmds += sam_cmds
        cmds += 'envsubst < %s > %s.tmp\n' % (map_fp, map_fp)
        cmds += 'mv %s.tmp %s\n' % (map_fp, map_fp)

    return map_fp, cmds, rms


def woltka_tax(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        map: str,
        database: str,
        params: dict
) -> tuple:
    """Get the taxonomic classification outputs and prepare the
    Woltka commands for this classification.

    Parameters
    ----------
    self : Commands class instance
        .out_dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pairing : str
        Type of reads pairing during alignment
    aligner : str
        Aligner used to create the input alignments
    map : str
        Path to the Woltka input file
    database : str
        Path to the WOL database
    params : dict
        Run parameters

    Returns
    -------
    tax_outmap : str
        Path to the folder containing the taxonomic maps.
    tax_to_do : list
        Empty if nothing is to be done
    """
    key = (tech, aligner)
    out = '/'.join([self.dir, tech, aligner, pairing])
    tax_out, tax_outmap = '%s/taxa' % out, '%s/taxmap' % out
    if params['outmap']:
        tax_outmap = params['outmap']
    tax_to_do = []
    taxa = params['taxa']
    for tdx, taxon in enumerate(taxa):
        out_dir = '%s/%s.tsv' % (tax_out, taxon)
        # self.outputs['outs'].setdefault(key, []).append(out_dir)
        if to_do(out_dir):
            tax_to_do.append(taxon)
            io_update(self, o_f=out_dir, key=key)
            self.soft.add_status(tech, 'all samples', 1, group='taxonomy',
                                 genome=taxon)
        else:
            self.soft.add_status(tech, 'all samples', 0, group='taxonomy',
                                 genome=taxon)

    taxid = '%s/taxonomy/taxid.map' % database
    nodes = '%s/taxonomy/nodes.dmp' % database
    names = '%s/taxonomy/names.dmp' % database
    if len(tax_to_do):
        cur_cmd = '\n# taxonomic\n'
        cur_cmd += 'woltka classify'
        cur_cmd += ' -i %s' % map
        cur_cmd += ' --map %s' % taxid

        if params['nodes_dmp']:
            cur_cmd += ' --nodes %s' % params['nodes_dmp']
        else:
            cur_cmd += ' --nodes %s' % nodes

        if params['names_dmp']:
            cur_cmd += ' --names %s' % params['names_dmp']
        else:
            cur_cmd += ' --names %s' % names

        if params['rank']:
            cur_cmd += ' --rank %s' % ','.join(taxa)
        else:
            cur_cmd += ' --rank %s' % ','.join(taxa)

        if params['add_rank']:
            cur_cmd += ' --add-rank'

        if params['add_lineage']:
            cur_cmd += ' --add-lineage'

        if params['name_as_id']:
            cur_cmd += ' --name-as-id'

        if params['to_biom']:
            cur_cmd += ' --to-biom'
        else:
            cur_cmd += ' --to-tsv'

        cur_cmd += ' --outmap %s' % tax_outmap
        cur_cmd += ' -o %s' % tax_out
        self.outputs['cmds'].setdefault(key, []).append(cur_cmd)
        io_update(self, o_d=tax_outmap, key=key)
    else:
        io_update(self, i_d=tax_outmap, key=key)
    return tax_outmap, tax_to_do


def woltka_go(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        map: str,
        taxmap: str,
        database: str,
        params: dict
) -> list:
    """Get the taxonomic classification outputs and prepare the Woltka
    commands for this classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pairing : str
        Type of reads pairing during alignment
    aligner : str
        Aligner used to create the input alignments
    map : str
        Path to the Woltka input file.
    taxmap : str
        Path to the folder containing the taxonomic maps.
    database : str
        WOL database path
    params : dict
        Run parameters

    Returns
    -------
    go_to_do : list
        Empty if nothing is to done
    """
    coords = '%s/proteins/coords.txt.xz' % database
    uniref_map = '%s/function/uniref/uniref.map.xz' % database
    go_rt = '%s/function/go' % database

    go_to_do = []
    gos = params['go']
    taxa = params['taxa']
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    go_dir = '%s/go' % out_dir
    key = (tech, aligner)
    io_update(self, o_d=go_dir, key=key)
    for go in gos:
        cur_out = '%s/%s.tsv' % (go_dir, go)
        cur_map = '%s/%s.map.xz' % (go_rt, go)
        cmd = '\n# go: %s [no stratification]\n' % go
        cmd += 'woltka classify'
        cmd += ' -i %s' % map
        cmd += ' --coords %s' % coords
        if params['map_as_rank']:
            cmd += ' --map-as-rank'
        cmd += ' --rank %s' % go
        cmd += ' --map %s' % uniref_map
        cmd += ' --to-tsv'
        cmd += ' --map %s' % cur_map
        cmd += ' -o %s' % cur_out
        if to_do(cur_out):
            go_to_do.append(cur_out)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=cur_out, key=key)
            self.soft.add_status(tech, 'all samples', 1, group='go',
                                 genome=go)
        else:
            self.soft.add_status(tech, 'all samples', 0, group='go',
                                 genome=go)
        # self.outputs['outs'].setdefault(key, []).append(cur_out)

    for stratif in taxa:
        for go in gos:
            go_out = '%s/%s_%s.tsv' % (go_dir, go, stratif)
            cur_map = '%s/%s.map.xz' % (go_rt, go)
            cmd = '\n# go: %s [%s]\n' % (go, stratif)
            cmd += 'woltka classify'
            cmd += ' -i %s' % map
            cmd += ' --coords %s' % coords
            if params['map_as_rank']:
                cmd += ' --map-as-rank'
            cmd += ' --rank %s' % go
            cmd += ' --stratify %s/%s' % (taxmap, stratif)
            cmd += ' --map %s' % uniref_map
            cmd += ' --map %s' % cur_map
            cmd += ' --to-tsv'
            cmd += ' -o %s' % go_out
            if to_do(go_out):
                go_to_do.append(go_out)
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                self.soft.add_status(tech, 'all samples', 1,
                                     group='go (%s)' % stratif, genome=go)
            else:
                self.soft.add_status(tech, 'all samples', 0,
                                     group='go (%s)' % stratif, genome=go)
            # self.outputs['outs'].setdefault(key, []).append(go_out)
    return go_to_do


def woltka_genes(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        woltka_map: str,
        taxmap: str,
        database: str,
        params: dict
) -> tuple:
    """Get the Woltka commands for the gene-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pairing : str
        Type of reads pairing during alignment
    aligner : str
        Aligner used to create the input alignments
    woltka_map : str
        Path to the Woltka input file.
    taxmap : str
        Path to the folder containing the taxonomic maps.
    database : str
        WOL database path
    params : dict
        Run parameters

    Returns
    -------
    genes_tax : dict
        Paths to the genes classification outputs
    genes_to_do : list
        Empty if nothing to be done
    """
    key = (tech, aligner)
    coords = '%s/proteins/coords.txt.xz' % database
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    genes_out = '%s/genes' % out_dir
    io_update(self, o_d=genes_out, key=key)
    genes = '%s/genes.biom' % genes_out
    genes_tax = {'': genes}
    genes_to_do = []
    if to_do(genes):
        genes_to_do.append(genes)
        cmd = '\n# per gene\n'
        cmd += 'woltka classify'
        cmd += ' -i %s' % woltka_map
        cmd += ' --coords %s' % coords
        cmd += ' -o %s' % genes
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=genes, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='genes')
    else:
        io_update(self, i_f=genes, key=key)
        self.soft.add_status(tech, 'all samples', 0, group='genes')
    # self.outputs['outs'].setdefault(key, []).append(genes)

    taxa = params['taxa']
    for stratif in taxa:
        genes = '%s/genes_%s.biom' % (genes_out, stratif)
        genes_tax[stratif] = genes
        if to_do(genes):
            genes_to_do.append(genes)
            cmd = '\n# per gene [%s]\n' % stratif
            cmd += 'woltka classify'
            cmd += ' -i %s' % woltka_map
            cmd += ' --coords %s' % coords
            cmd += ' --stratify %s/%s' % (taxmap, stratif)
            cmd += ' -o %s' % genes
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=genes, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='genes (%s)' % stratif)
        else:
            io_update(self, i_f=genes, key=key)
            self.soft.add_status(
                tech, 'all samples', 0, group='genes (%s)' % stratif)

        # self.outputs['outs'].setdefault(key, []).append(genes)
    return genes_tax, genes_to_do


def woltka_uniref(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        genes_tax: dict,
        database: str,
        params: dict
) -> dict:
    """Get the Woltka commands for the uniref-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    aligner : str
        Aligner used to create the input alignments
    pairing : str
        Type of reads pairing during alignment
    genes_tax : dict
        Paths to the genes classification outputs
    database : str
        WOL database path
    params : dict
        Run parameters

    Returns
    -------
    uniref_tax : dict
        Paths to the uniref classifications
    """
    key = (tech, aligner)
    uniref_map = '%s/function/uniref/uniref.map.xz' % database
    uniref_names = '%s/function/uniref/uniref.name.xz' % database
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    uniref_out = '%s/uniref' % out_dir
    io_update(self, o_d=uniref_out, key=key)
    uniref = '%s/uniref.biom' % uniref_out
    uniref_tax = {'': uniref}
    if to_do(uniref):
        cmd = '\n# uniref\n'
        cmd += 'echo "%s -> %s"\n' % (basename(genes_tax['']), basename(uniref))
        cmd += 'woltka tools collapse'
        cmd += ' --input %s' % genes_tax['']
        cmd += ' --map %s' % uniref_map
        cmd += ' --names %s' % uniref_names
        cmd += ' --output %s' % uniref
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=uniref, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='uniref')
    else:
        io_update(self, i_f=uniref, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='uniref')
    # self.outputs['outs'].setdefault(key, []).append(uniref)

    taxa = params['taxa']
    for stratif in taxa:
        uniref = '%s/uniref_%s.biom' % (uniref_out, stratif)
        uniref_tax[stratif] = uniref
        if to_do(uniref):
            cmd = '\n# uniref [%s]\n' % stratif
            cmd += 'echo "%s -> %s"\n' % (basename(genes_tax[stratif]),
                                          basename(uniref))
            cmd += 'woltka tools collapse'
            cmd += ' --input %s' % genes_tax[stratif]
            cmd += ' --map %s' % uniref_map
            cmd += ' --names %s' % uniref_names
            cmd += ' --field %s' % params['field']
            cmd += ' --output %s' % uniref
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=uniref, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='uniref (%s)' % stratif)
        else:
            io_update(self, i_f=uniref, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='uniref (%s)' % stratif)
        # self.outputs['outs'].setdefault(key, []).append(uniref)

    return uniref_tax


def woltka_eggnog(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        uniref_tax: dict,
        database: str,
        params: dict
) -> None:
    """Get the Woltka commands for the eggnog-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    aligner : str
        Aligner used to create the input alignments
    pairing : str
        Type of reads pairing during alignment
    uniref_tax : dict
        Paths to the uniref classifications
    database : str
        WOL database path
    params : dict
        Run parameters
    """
    key = (tech, aligner)
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    eggnog_out = '%s/eggnog' % out_dir
    io_update(self, o_d=eggnog_out, key=key)
    biom = '%s/eggnog.biom' % eggnog_out
    if to_do(biom):
        cmd = '\n# eggnog [no stratification]\n'
        cmd += 'echo "%s -> %s"\n' % (basename(uniref_tax['']),
                                      basename(biom))
        cmd += 'woltka tools collapse'
        cmd += ' --input %s' % uniref_tax['']
        cmd += ' --map %s/function/eggnog/eggnog.map.xz' % database
        cmd += ' --output %s\n\n' % biom
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=biom, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='eggnog')
    else:
        io_update(self, i_f=biom, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='eggnog')
    tsv = '%s.tsv' % splitext(biom)[0]
    if to_do(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=tsv, key=key)
    # self.outputs['outs'].setdefault(key, []).extend([tsv, biom])

    taxa = params['taxa']
    for stratif in taxa:
        biom = '%s/eggnog_%s.biom' % (eggnog_out, stratif)
        if to_do(biom):
            cmd = '\n# eggnog [%s]\n' % stratif
            cmd += 'echo "%s -> %s"\n' % (basename(uniref_tax[stratif]),
                                          basename(biom))
            cmd += 'woltka tools collapse'
            cmd += ' --input %s' % uniref_tax[stratif]
            cmd += ' --map %s/function/eggnog/eggnog.map.xz' % database
            cmd += ' --field %s' % params['field']
            cmd += ' --output %s\n\n' % biom
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=biom, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='eggnog (%s)' % stratif)
        else:
            io_update(self, i_f=biom, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='eggnog (%s)' % stratif)

        tsv = '%s.tsv' % splitext(biom)[0]
        if to_do(tsv):
            cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
            cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
            cmd += 'rm %s.tmp\n' % tsv
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=tsv, key=key)
        # self.outputs['outs'].setdefault(key, []).extend([tsv, biom])


def woltka_cazy(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        genes: str,
        genes_tax: dict,
        database: str,
        params: dict
) -> None:
    """Get the Woltka commands for the cazy-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .databases
            All databases class instance
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    aligner : str
        Aligner used to create the input alignments
    pairing : str
        Type of reads pairing during alignment
    genes : str
        Path to the genes classification.
    genes_tax : dict
        Path to the stratified genes classification output.
    database : str
        WOL database path
    params : dict
        Run parameters
    """
    key = (tech, aligner)
    cazy_map = '%s/function/cazy/3tools.txt' % database
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    cazy_out = '%s/cazy' % out_dir
    io_update(self, o_d=cazy_out, key=key)
    biom = '%s/cazy.biom' % cazy_out
    if to_do(biom):
        cmd = '\n# cazy\n'
        cmd += 'echo "%s -> %s"\n' % (basename(genes), basename(biom))
        cmd += 'woltka tools collapse'
        cmd += ' --input %s' % genes
        cmd += ' --map %s' % cazy_map
        cmd += ' --output %s\n\n' % biom
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=biom, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='cazy')
    else:
        io_update(self, i_f=biom, key=key)
        self.soft.add_status(tech, 'all samples', 1, group='cazy')
    tsv = '%s.tsv' % splitext(biom)[0]
    if to_do(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=tsv, key=key)
    # self.outputs['outs'].setdefault(key, []).extend([biom, tsv])

    taxa = params['taxa']
    for stratif in taxa:
        cazy_map = '%s/function/cazy/3tools.txt' % database
        biom = '%s/cazy_%s.biom' % (cazy_out, stratif)
        if to_do(biom):
            cmd = '\n# cazy [%s]\n' % stratif
            cmd += 'echo "%s -> %s"\n' % (basename(genes_tax[stratif]),
                                          basename(biom))
            cmd += 'woltka tools collapse'
            cmd += ' --input %s' % genes_tax[stratif]
            cmd += ' --map %s' % cazy_map
            cmd += ' --field %s' % params['field']
            cmd += ' --output %s\n\n' % biom
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=biom, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='cazy (%s)' % stratif)
        else:
            io_update(self, i_f=biom, key=key)
            self.soft.add_status(
                tech, 'all samples', 1, group='cazy (%s)' % stratif)
        tsv = '%s.tsv' % splitext(biom)[0]
        if to_do(tsv):
            cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
            cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
            cmd += 'rm %s.tmp\n' % tsv
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, o_f=tsv, key=key)
        # self.outputs['outs'].setdefault(key, []).extend([biom, tsv])


def woltka_metacyc(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        genes_tax: dict,
        database: str,
        params: dict
) -> None:
    """Get the Woltka commands for the metacyc-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    aligner : str
        Aligner used to create the input alignments
    pairing : str
        Type of reads pairing during alignment
    genes_tax : dict
        Paths to the genes classification outputs
    database : str
        WOL database path
    params : dict
        Run parameters
    """
    metacyc_dir = '%s/function/metacyc' % database
    metacyc = [('protein', 'protein_name.txt', 'protein.map.xz'),
               ('enzrxn', 'enzrxn_name.txt', 'protein-to-enzrxn.txt'),
               ('reaction', 'reaction_name.txt', 'enzrxn-to-reaction.txt'),
               ('pathway', 'pathway_name.txt', 'reaction-to-pathway.txt'),
               ('super_pathway', '', 'pathway-to-super_pathway.txt'),
               ('regulation', '', 'enzrxn-to-regulation.txt'),
               ('regulator', '', 'regulation-to-regulator.txt'),
               ('ec', '', 'reaction-to-ec.txt')]
    files = {}
    files_tax = {}
    key = (tech, aligner)
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    woltka_fun_out = '%s/metacyc' % out_dir
    io_update(self, o_d=woltka_fun_out, key=key)
    cmd = ''
    for idx, (level, names, maps) in enumerate(metacyc):
        if '-to-' in maps:
            input_biom = files[maps.split('-to-')[0]]
        else:
            input_biom = genes_tax['']
        biom = '%s/%s.biom' % (woltka_fun_out, level)
        tsv = '%s.tsv' % splitext(biom)[0]
        if to_do(biom):
            cmd += '\n# metacyc: %s [no stratification]\n' % level
            cmd += 'echo "%s -> %s"\n' % (basename(input_biom),
                                          basename(biom))
            cmd += 'woltka tools collapse'
            cmd += ' --input %s' % input_biom
            if names:
                cmd += ' --names %s/%s' % (metacyc_dir, names)
            cmd += ' --map %s/%s' % (metacyc_dir, maps)
            cmd += ' --output %s\n' % biom
            io_update(self, o_f=biom, key=key)
            self.soft.add_status(tech, 'all samples', 1,
                                 group='metacyc', genome=level)
        else:
            io_update(self, i_f=biom, key=key)
            self.soft.add_status(tech, 'all samples', 0,
                                 group='metacyc', genome=level)

        if to_do(tsv):
            cmd += 'biom convert'
            cmd += ' -i %s' % biom
            cmd += ' -o %s.tmp' % tsv
            cmd += ' --to-tsv\n'
            cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
            cmd += 'rm %s.tmp\n' % tsv
            io_update(self, o_f=tsv, key=key)
        else:
            io_update(self, i_f=tsv, key=key)
        files[level] = biom

        taxa = params['taxa']
        for stratif in taxa:
            if stratif not in files_tax:
                files_tax[stratif] = {}
            if '-to-' in maps:
                input_biom = files_tax[stratif][maps.split('-to-')[0]]
            else:
                input_biom = genes_tax[stratif]
            biom = '%s/%s_%s.biom' % (woltka_fun_out, level, stratif)
            tsv = '%s.tsv' % splitext(biom)[0]
            if to_do(biom):
                cmd += '\n# metacyc: %s [%s]\n' % (level, stratif)
                cmd += 'echo "%s -> %s"\n' % (basename(input_biom),
                                              basename(biom))
                cmd += 'woltka tools collapse'
                cmd += ' --input %s' % input_biom
                if names:
                    cmd += ' --names %s/%s' % (metacyc_dir, names)
                cmd += ' --field %s' % params['field']
                cmd += ' --map %s/%s' % (metacyc_dir, maps)
                cmd += ' --output %s\n' % biom
                io_update(self, o_f=biom, key=key)
                self.soft.add_status(tech, 'all samples', 1,
                                     group='metacyc (%s)' % stratif,
                                     genome=level)
            else:
                io_update(self, i_f=biom, key=key)
                self.soft.add_status(tech, 'all samples', 0,
                                     group='metacyc (%s)' % stratif,
                                     genome=level)
            if to_do(tsv):
                cmd += 'biom convert'
                cmd += ' -i %s' % biom
                cmd += ' -o %s.tmp' % tsv
                cmd += ' --to-tsv\n'
                cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                cmd += 'rm %s.tmp\n' % tsv
                io_update(self, o_f=tsv, key=key)
            else:
                io_update(self, i_f=tsv, key=key)
            files_tax[stratif][level] = biom

    if cmd:
        self.outputs['cmds'].setdefault(key, []).append(cmd)
    # self.outputs['outs'].setdefault(key, []).extend(files)


def woltka_kegg(
        self,
        tech: str,
        pairing: str,
        aligner: str,
        uniref_tax: dict,
        database: str,
        params: dict
) -> None:
    """Get the Woltka commands for the kegg-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    aligner : str
        Aligner used to create the input alignments
    pairing : str
        Type of reads pairing during alignment
    uniref_tax : dict
        Paths to the uniref classifications
    database : str
        WOL database path
    params : dict
        Run parameters
    """
    ko_names_maps = [
        ('ko', 'ko.name', 'ko.map.xz', ''),
        ('query', '', '', ''),
        ('ko-cog', '', 'ko-to-cog.txt', 'ko'),
        ('ko-disease', 'disease_name.txt', 'ko-to-disease.txt', 'ko'),
        ('ko-ec', '', 'ko-to-ec.txt', 'ko'),
        ('ko-go', '', 'ko-to-go.txt', 'ko'),
        ('ko-module', 'module_name.txt', 'ko-to-module.txt', 'ko'),
        ('ko-pathway', 'pathway_name.txt', 'ko-to-pathway.txt', 'ko'),
        ('ko-reaction', 'reaction_name.txt', 'ko-to-reaction.txt', 'ko'),
        ('ko-module-class', '', 'module-to-class.txt', 'ko-module'),
        ('ko-module-compound', 'compound_name.txt',
         'module-to-compound.txt', 'ko-module'),
        ('ko-module-ko', 'ko_name.txt', 'module-to-ko.txt', 'ko-module'),
        ('ko-module-pathway', 'pathway_name.txt',
         'module-to-pathway.txt', 'ko-module'),
        ('ko-module-reaction', 'reaction_name.txt',
         'module-to-reaction.txt', 'ko-module'),
        ('ko-pathway-class', '', 'pathway-to-class.txt', 'ko-pathway'),
        # ('ko-pathway-compound', 'compound_name.txt',
        #  'pathway-to-compound.txt', 'ko-pathway'),
        ('ko-pathway-disease', 'disease_name.txt',
         'pathway-to-disease.txt', 'ko-pathway'),
        # ('ko-pathway-ko', 'ko_name.txt', 'pathway-to-ko.txt', 'ko-pathway'),
        ('ko-pathway-module', 'module_name.txt',
         'pathway-to-module.txt', 'ko-pathway'),
        ('ko-reaction-compound', 'compound_name.txt',
         'reaction-to-compound.txt', 'ko-reaction'),
        ('ko-reaction-enzyme', '', 'reaction-to-enzyme.txt', 'ko-reaction'),
        ('ko-reaction-ko', 'ko_name.txt', 'reaction-to-ko.txt', 'ko-reaction'),
        ('ko-reaction-left_compound', 'compound_name.txt',
         'reaction-to-left_compound.txt', 'ko-reaction'),
        ('ko-reaction-module', 'module_name.txt',
         'reaction-to-module.txt', 'ko-reaction'),
        ('ko-reaction-pathway', 'pathway_name.txt',
         'reaction-to-pathway.txt', 'ko-reaction'),
        ('ko-reaction-rclass', 'rclass_name.txt',
         'reaction-to-rclass.txt', 'ko-reaction'),
        ('ko-reaction-right_compound', 'compound_name.txt',
         'reaction-to-right_compound.txt', 'ko-reaction')]
    files = []
    key = (tech, aligner)
    out_dir = '/'.join([self.dir, tech, aligner, pairing])
    kegg_out = '%s/kegg' % out_dir
    kegg_maps = '%s/kegg_queried' % kegg_out
    io_update(self, o_d=kegg_out, key=key)
    kegg_maps_local = kegg_maps.replace('${SCRATCH_FOLDER}', '')
    queried = len(glob.glob('%s/*' % kegg_maps_local))
    for (level, name, maps, prev) in ko_names_maps:
        if maps:
            biom = '%s/%s.biom' % (kegg_out, level)
            tsv = '%s.tsv' % splitext(biom)[0]
            tsv_ = str(tsv)
            if not prev:

                if to_do(biom):
                    cmd = '\n# kegg: %s [no stratification]\n' % name
                    cmd += 'echo "%s -> %s"\n' % (basename(uniref_tax['']),
                                                  basename(biom))
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % uniref_tax['']
                    cmd += ' --names %s/function/kegg/%s' % (database, name)
                    cmd += ' --map %s/function/kegg/%s' % (database, maps)
                    cmd += ' --output %s\n' % biom
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, o_f=biom, key=key)
                    self.soft.add_status(tech, 'all samples', 1,
                                         group='kegg', genome=level)
                else:
                    io_update(self, i_f=biom, key=key)
                    self.soft.add_status(tech, 'all samples', 0,
                                         group='kegg', genome=level)
                if to_do(tsv):
                    cmd = 'biom convert -i %s' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                    cmd += 'rm %s.tmp\n' % tsv
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, o_f=tsv, key=key)
                else:
                    io_update(self, i_f=biom, key=key)

            elif queried > 1:
                io_update(self, i_d=kegg_maps, key=key)
                input_biom = '%s/%s.biom' % (kegg_out, prev)
                if to_do(biom):
                    cmd = '\n# kegg: %s [no stratification]\n' % name
                    cmd += 'echo "%s -> %s"\n' % (basename(input_biom),
                                                  basename(biom))
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % input_biom
                    if name:
                        cmd += ' --names %s/%s' % (kegg_maps, name)
                    cmd += ' --map %s/%s' % (kegg_maps, maps)
                    cmd += ' --output %s\n' % biom
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, o_f=biom, key=key)
                    self.soft.add_status(tech, 'all samples', 1,
                                         group='kegg', genome=level)
                else:
                    io_update(self, i_f=biom, key=key)
                    self.soft.add_status(tech, 'all samples', 0,
                                         group='kegg', genome=level)
                if to_do(tsv):
                    cmd = 'biom convert -i %s' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                    cmd += 'rm %s.tmp\n\n' % tsv
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, o_f=tsv, key=key)
                else:
                    io_update(self, i_f=tsv, key=key)

            taxa = params['taxa']
            for stratif in taxa:
                biom = '%s/%s_%s.biom' % (kegg_out, level, stratif)
                tsv = '%s.tsv' % splitext(biom)[0]
                if not prev:
                    if to_do(biom):
                        cmd = '\n# kegg: %s [%s]\n' % (name, stratif)
                        cmd += 'echo "%s -> %s"\n' % (
                            basename(uniref_tax[stratif]), basename(biom))
                        cmd += 'woltka tools collapse'
                        cmd += ' --input %s' % uniref_tax[stratif]
                        cmd += ' --names %s/function/kegg/%s' % (database, name)
                        cmd += ' --map %s/function/kegg/%s' % (database, maps)
                        cmd += ' --field %s' % params['field']
                        cmd += ' --output %s\n' % biom
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                        io_update(self, o_f=biom, key=key)
                        self.soft.add_status(tech, 'all samples', 1,
                                             group='kegg (%s)' % stratif,
                                             genome=level)
                    else:
                        io_update(self, i_f=biom, key=key)
                        self.soft.add_status(tech, 'all samples', 0,
                                             group='kegg (%s)' % stratif,
                                             genome=level)
                    if to_do(tsv):
                        cmd = 'biom convert -i %s' % biom
                        cmd += ' -o %s.tmp --to-tsv\n' % tsv
                        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                        cmd += ' rm %s.tmp\n' % tsv
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                        io_update(self, o_f=tsv, key=key)
                    else:
                        io_update(self, i_f=biom, key=key)

                elif queried > 1:
                    io_update(self, i_d=kegg_maps, key=key)
                    input_biom = '%s/%s_%s.biom' % (kegg_out, prev, stratif)
                    if to_do(biom):
                        cmd = '\n# kegg: %s [%s]\n' % (name, stratif)
                        cmd += 'echo "%s -> %s"\n' % (basename(input_biom),
                                                      basename(biom))
                        cmd += 'woltka tools collapse'
                        cmd += ' --input %s' % input_biom
                        if name:
                            cmd += ' --names %s/%s' % (kegg_maps, name)
                        cmd += ' --map %s/%s' % (kegg_maps, maps)
                        cmd += ' --field %s' % params['field']
                        cmd += ' --output %s\n' % biom
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                        io_update(self, o_f=biom, key=key)
                        self.soft.add_status(tech, 'all samples', 1,
                                             group='kegg (%s)' % stratif,
                                             genome=level)
                    else:
                        io_update(self, i_f=biom, key=key)
                        self.soft.add_status(tech, 'all samples', 0,
                                             group='kegg (%s)' % stratif,
                                             genome=level)
                    if to_do(tsv):
                        cmd = 'biom convert -i %s' % biom
                        cmd += ' -o %s.tmp --to-tsv\n' % tsv
                        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                        cmd += 'rm %s.tmp\n\n' % tsv
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                        io_update(self, o_f=tsv, key=key)
                    else:
                        io_update(self, i_f=tsv, key=key)
        else:
            kinfo = '%s/kegg_info.txt' % kegg_maps_local
            nkegg = glob.glob('%s/*' % kegg_maps_local)
            if (to_do(kinfo) or len(nkegg) < 2) and not to_do(tsv_):
                cmd = 'mkdir -p %s\n' % kegg_maps_local
                cmd += 'cd %s\n' % kegg_maps_local
                cmd += 'cp %s %s/%s\n' % (tsv_.replace('${SCRATCH_FOLDER}', ''),
                                          kegg_maps_local, basename(tsv_))
                cmd += 'python3 %s/kegg_query.py %s/%s\n' % (
                    RESOURCES, kegg_maps_local, basename(tsv_))
                cmd += 'rm %s/%s\n' % (kegg_maps_local, basename(tsv_))
                self.outputs['bash'].append(cmd)
    # self.outputs['outs'].setdefault((tech, aligner), []).extend(files)


def woltka_pairing(
        self,
) -> str:
    """Get the type of reads pairing that was defined by the
    user if the previous software was bowtie2.

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software
        .soft.params
            Parameters

    Returns
    -------
    pairing : str
        Pairing parameter used for Bowtie2
    """
    pairing = 'single'
    hash_value = self.hashes[tuple(self.soft.path[:-1])]
    if self.soft.prev == 'bowtie2':
        if self.softs['bowtie2'][hash_value].params['paired']:
            pairing = 'paired'
    return pairing


def woltka(self) -> None:
    """Woltka is a versatile program for determining the structure and
    functional capacity of microbiomes. It mainly works with shotgun
    metagenomic data. It bridges first-pass sequence aligners with advanced
    analytical platforms (such as QIIME 2). It takes full advantage of,
    and is not limited by, the WoL reference database:
    https://biocore.github.io/wol

    References
    ----------
    Zhu, Qiyun, et al. "Phylogeny-Aware Analysis of Metagenome Community
    Ecology Based on Matched Reference Genomes while Bypassing Taxonomy."
    Msystems 7.2 (2022): e00167-22.

    Notes
    -----
    GitHub  : https://github.com/qiyunzhu/woltka
    Paper   : https://doi.org/10.1128/msystems.00167-22

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .inputs : dict
            Input files
        .prev : str
            Previous software
        .soft.params
            Parameters
        .outputs : dict
            All outputs
        .databases
            All databases class instance
    """
    db = self.databases.paths['wol']
    for tech in self.config.techs:
        params = tech_params(self, tech)
        pairing = woltka_pairing(self)
        alignments = woltka_aligments(self, tech)
        if not alignments:
            continue
        classifs = params['classifications']
        for aligner, alis in alignments.items():
            key = (tech, aligner)
            to_dos = status_update(
                self, tech, list(alis.values()), group=aligner)
            map_fp, cmds, rms = woltka_map(self, tech, pairing, aligner, alis)

            out = '/'.join([self.dir, tech, aligner, pairing])
            self.outputs['outs'].setdefault(key, []).append(out)

            taxmap, tdo = woltka_tax(self, tech, pairing, aligner, map_fp,
                                     db, params)
            if 'go' in classifs:
                gdo = woltka_go(
                    self, tech, pairing, aligner, map_fp, taxmap, db, params)
            else:
                gdo = []
            if set(classifs) & {'eggnog', 'metacyc', 'kegg'}:
                genes_tax, edo = woltka_genes(self, tech, pairing, aligner,
                                              map_fp, taxmap, db, params)
            if tdo or gdo or edo:
                io_update(self, i_f=list(alis.values()), key=key)
                if key in self.outputs['cmds']:
                    cmd = self.outputs['cmds'][key]
                    if to_dos:
                        self.outputs['cmds'][key] = [False]
                    else:
                        self.outputs['cmds'][key] = [cmds] + cmd + [rms]
            if not to_dos:
                if set(classifs) & {'eggnog', 'metacyc', 'kegg'}:
                    uniref_tax = woltka_uniref(self, tech, pairing, aligner,
                                               genes_tax, db, params)
                if 'eggnog' in classifs:
                    woltka_eggnog(self, tech, pairing, aligner, uniref_tax,
                                  db, params)
                # if 'cazy' in classifs:
                #     woltka_cazy(
                #         self, tech, pairing, aligner, genes, genes_tax, db,
                #         params)
                if 'metacyc' in classifs:
                    woltka_metacyc(self, tech, pairing, aligner, genes_tax,
                                   db, params)
                if 'kegg' in classifs:
                    woltka_kegg(self, tech, pairing, aligner, uniref_tax, db,
                                params)


def get_midas_cmd(
        self,
        tech: str,
        inputs: list,
        focus_dir: str,
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
    inputs : list
        Paths to the input files
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    focus_dir : str
        Path to the output folder.
    analysis : str
        MIDAS analysis (any of "species", "genes" or "snps").
    select : set
        Species names for which there is a reference in the database.

    Returns
    -------
    cmd : str
        Midas command line for the species level.
    """
    params = tech_params(self, tech)
    cmd = 'export PATH=$PATH:%s/scripts\n' % params['path']
    cmd += 'export PYTHONPATH=$PYTHONPATH:%s\n' % params['path']
    cmd += '\nrun_midas.py %s' % analysis
    cmd += ' %s' % focus_dir
    cmd += ' -1 %s' % inputs[0]
    if len(inputs) > 1:
        cmd += ' -2 %s' % inputs[1]
    cmd += ' -d %s' % self.databases.paths['midas']
    cmd += ' -t %s' % params['cpus']
    cmd += ' --remove_temp'
    for param in ['n', 'mapid', 'aln_cov']:
        if params[param]:
            if len(param) == 1:
                cmd += ' -%s %s' % (param, params[param])
            else:
                cmd += ' --%s %s' % (param, params[param])
    if analysis != 'species':
        for param in ['m', 's', 'species_cov', 'species_topn', 'readq', 'trim']:
            if params[param]:
                if len(param) == 1:
                    cmd += ' -%s %s' % (param, params[param])
                else:
                    cmd += ' --%s %s' % (param, params[param])
    else:
        cmd += ' --word_size %s' % params['word_size']
    if select:
        cmd += ' --species_id %s' % ','.join(list(select))
    return cmd


def midas_species(
        self,
        tech: str,
        inputs: list,
        focus_dir: str,
        analysis: str,
        to_dos: list
) -> None:
    """Collect the MIDAS commands for the species level.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Paths to the input files
    focus_dir : str
        Path to the species output folder
    analysis : str
        Current analysis: 'species', 'genes', or 'snps'
    to_dos : list
        Paths to input file that need to be created (or whether they are links)
    """
    spcs_out = '%s/species' % focus_dir
    spcs_profile = '%s/species_profile.txt' % spcs_out
    if not self.config.force and not to_do(spcs_profile):
        io_update(self, i_d=spcs_out, key=tech)
        self.soft.add_status(tech, self.sam_pool, 1, group='species')
    else:
        cmd = get_midas_cmd(self, tech, inputs, focus_dir, analysis)
        if to_dos:
            self.outputs['cmds'].setdefault((tech,), []).append(False)
        else:
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)
        io_update(self, i_f=inputs, o_d=spcs_out, key=tech)
        self.soft.add_status(tech, self.sam_pool, 0, group='species')
    self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(spcs_out)
    self.outputs['dirs'].append(spcs_out)


def midas_genes(
        self,
        tech: str,
        inputs: list,
        focus_dir: str,
        genes_out: str,
        analysis: str,
        select: set,
        to_dos: list
) -> None:
    """Collect the MIDAS commands for the genes level.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Paths to the input files
    focus_dir : str
        Path to the species output folder
    genes_out : str
        Path to the output folder
    analysis : str
        Current analysis: 'species', 'genes', or 'snps'
    select : set
        Set of taxa to focus
    to_dos : list
        Paths to input file that need to be created (or whether they are links)
    """
    if self.config.force or to_do('%s/readme.txt' % genes_out):
        cmd = get_midas_cmd(self, tech, inputs, focus_dir, analysis, select)
        if to_dos:
            self.outputs['cmds'].setdefault((tech,), []).append(False)
        else:
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)
        io_update(self, o_d=genes_out, key=tech)
        self.soft.add_status(tech, self.sam_pool, 1, group='genes')
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group='genes')
    self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(genes_out)
    self.outputs['dirs'].append(genes_out)


def midas_snps(
        self,
        tech: str,
        inputs: list,
        focus_dir: str,
        snps_out: str,
        analysis: str,
        select: set,
        to_dos: list
) -> None:
    """Collect the MIDAS commands for the species level.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Paths to the input files
    focus_dir : str
        Path to the species output folder
    snps_out : str
        Path to the output folder
    analysis : str
        Current analysis: 'species', 'genes', or 'snps'
    select : set
        Set of taxa to focus
    to_dos : list
        Paths to input file that need to be created (or whether they are links)
    """
    if self.config.force or to_do('%s/readme.txt' % snps_out):
        cmd = get_midas_cmd(self, tech, inputs, focus_dir, analysis, select)
        if to_dos:
            self.outputs['cmds'].setdefault((tech,), []).append(False)
        else:
            self.outputs['cmds'].setdefault((tech,), []).append(cmd)
        io_update(self, o_d=snps_out, key=tech)
        self.soft.add_status(tech, self.sam_pool, 1, group='snps')
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group='snps')
    self.outputs['outs'].setdefault((tech, self.sam_pool), []).append(snps_out)
    self.outputs['dirs'].append(snps_out)


def get_species_select(
        self,
        species_list: str
) -> set:
    """Get the species names for which there is a reference in the database.

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    species_list : str
        Path to file containing list of species to focus on.

    Returns
    -------
    select : set
        Species names for which there is a reference in the database.
    """
    select = set()
    if species_list:
        if self.config.dev and to_do(species_list):
            select.add('Escherichia')
            return select
        else:
            species = [x.strip() for x in open(species_list).readlines()]
        with open('%s/species_info.txt' % self.databases.paths['midas']) as f:
            for line in f:
                for genus_species in species:
                    if genus_species.replace(' ', '_') in line:
                        select.add(line.split('\t')[0])
    return select


def midas(self) -> None:
    """Metagenomic Intra-Species Diversity Analysis System (MIDAS).
    MIDAS is an integrated pipeline that leverages >30,000 reference genomes
    to estimate bacterial species abundance and strain-level genomic
    variation, including gene content and SNPs, from shotgun metagnomes.

    References
    ----------
    Nayfach, Stephen, et al. "An integrated metagenomics pipeline for strain
    profiling reveals novel patterns of bacterial transmission and
    biogeography." Genome research 26.11 (2016): 1612-1625.

    Notes
    -----
    GitHub  : https://github.com/snayfach/MIDAS
    Paper   : https://doi.org/10.1101/gr.201863.115

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
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam):
            continue
        to_dos = status_update(self, tech, inputs)

        params = tech_params(self, tech)
        for focus, species_list in params['focus'].items():
            focus_dir = '/'.join([self.dir, tech, focus, self.sam_pool])
            midas_species(
                self, tech, inputs, focus_dir, 'species', to_dos)
            select = set(get_species_select(self, species_list))
            genes = '%s/genes' % focus_dir
            midas_genes(
                self, tech, inputs, focus_dir, genes, 'genes', select, to_dos)
            snps = '%s/snps' % focus_dir
            midas_snps(
                self, tech, inputs, focus_dir, snps, 'snps', select, to_dos)


def get_kraken2_db(
        self,
        db: str
) -> str:
    """Get the path the Kraken2 database passed by the user.

    Notes
    -----
    If the database is not valid, the pipeline will stop with a useful hint.

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
    db : str
        Name of the Kraken2 database

    Returns
    -------
    str
        Path to the Kraken2 database
    """
    if db in self.databases.paths:
        for subdir in ['', '/kraken2', '/databases', '/databases/kraken2']:
            db_path = '%s%s' % (self.databases.paths[db], subdir)
            if self.config.dev:
                return db_path
            if isfile('%s/hash.k2d' % db_path):
                return db_path
            nested = glob.glob('%s/*/hash.k2d' % db_path)
            if nested:
                return dirname(nested[0])
    sys.exit('[kraken2] Database not found: %s' % db)


def get_kraken2_cmd(
        self,
        tech: str,
        inputs: list,
        out: str,
        db_path: str,
) -> str:
    """Collect the Kraken2 command.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Kraken2
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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Path to the input files
    out : str
        Path to the output folder
    db_path : str
        Path to the Kraken2 database

    Returns
    -------
    cmd : str
        Kraken2 command
    """
    params = tech_params(self, tech)

    cmd = ''
    result = '%s/result.tsv' % out
    report = '%s/report.tsv' % out
    if to_do(result):
        cmd += 'kraken2 '
        cmd += ' -db %s' % db_path
        cmd += ' --report %s' % report
        cmd += ' --threads %s' % params['cpus']
        cmd += ' --confidence %s' % params['confidence']
        if len(inputs) > 1:
            cmd += ' --paired'
        cmd += ' %s > %s\n' % (' '.join(inputs), result)
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out
    return cmd


def kraken2(self) -> None:
    """Kraken is a taxonomic sequence classifier that assigns taxonomic
    labels to DNA sequences. Kraken examines the -mers within a query
    sequence and uses the information within those -mers to query a database.
    That database maps -mers to the lowest common ancestor (LCA) of all
    genomes known to contain a given -mer.
    The first version of Kraken used a large indexed and sorted list of
    -mer/LCA pairs as its database. While fast, the large memory requirements
    posed some problems for users, and so Kraken 2 was created to provide a
    solution to those problems.
    Kraken 2 differs from Kraken 1 in several important ways:
    see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

    References
    ----------
    Wood, Derrick E., Jennifer Lu, and Ben Langmead. "Improved metagenomic
    analysis with Kraken 2." Genome biology 20.1 (2019): 1-13.

    Notes
    -----
    GitHub  : https://github.com/DerrickWood/kraken2
    Docs    : https://ccb.jhu.edu/software/kraken2
    Paper   : https://doi.org/10.1186/s13059-019-1891-0

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Kraken2
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
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam):
            continue
        to_dos = status_update(self, tech, inputs)

        params = tech_params(self, tech)
        for db in params['databases']:
            out = '/'.join([self.dir, tech, self.sam_pool, db])
            self.outputs['dirs'].append(out)
            self.outputs['outs'].setdefault((tech, sam), []).append((db, out))
            if self.config.force or to_do('%s/result.tsv.gz' % out):
                db_path = get_kraken2_db(self, db)
                cmd = get_kraken2_cmd(self, tech, inputs, out, db_path)
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                io_update(self, i_f=inputs, o_d=out, key=tech)
                self.soft.add_status(tech, sam, 1, group=db)
            else:
                self.soft.add_status(tech, sam, 0, group=db)


def get_bracken_db(
        self,
        db: str,
) -> str:
    """Get the path the Bracken database passed by the user.

    Notes
    -----
    If the database is not valid, the pipeline will stop with a useful hint.

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
    db : str
        Name of the Bracken database

    Returns
    -------
    path : str
        Path tho the Bracken database
    """
    if self.config.dev:
        return 'dummy/bracken/path'
    if db == 'default':
        return self.databases.paths['kraken2']
    elif 'bracken' in self.databases.builds[db]:
        return self.databases.builds[db]['bracken']
    elif db in self.databases.paths:
        for subdir in ['', '/bracken', '/databases/bracken', '/kraken2',
                       '/databases', '/databases/kraken2']:
            db_path = '%s%s' % (self.databases.paths[db], subdir)
            db_paths = glob.glob('%s/*kmer_distrib' % db_path)
            if db_paths:
                return db_path
            nested = glob.glob('%s/*/*kmer_distrib' % db_path)
            if nested:
                return dirname(nested[0])
    sys.exit('[bracken] Database not found: %s' % db)


def bracken_cmd(
        self,
        tech: str,
        db_path: str,
        report: str,
        out_dir: str
) -> str:
    """Collect command for Bracken.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    db_path : str
        Path tho the Bracken database
    report : str
        Path to the kraken2 input file
    out_dir : str
        Path to the output folder

    Returns
    -------
    cmd : str
        Bracken command
    """
    params = tech_params(self, tech)

    cmd = ''
    result = '%s/results.tsv' % out_dir
    if to_do(result):
        cmd += 'PATH=$PATH:%s:%s/src\n' % (params['path'], params['path'])
        cmd += 'export PATH\n'
        cmd += 'gunzip %s\n' % report
        cmd += 'bracken'
        cmd += ' -d %s' % db_path
        cmd += ' -i %s' % report.rstrip('.gz')
        cmd += ' -o %s' % result
        cmd += ' -w %s/report.tsv' % out_dir
        cmd += ' -r %s' % params['read_len']
        cmd += ' -l %s' % params['level']
        cmd += ' -t %s\n' % params['threshold']

    cmd += 'gzip %s\n' % result
    cmd += 'gzip %s/report.tsv\n' % out_dir
    return cmd


def bracken(self) -> None:
    """Bracken (Bayesian Reestimation of Abundance with KrakEN) is a highly
    accurate statistical method that computes the abundance of species in DNA
    sequences from a metagenomics sample.
    Braken uses the taxonomy labels
    assigned by Kraken, a highly accurate metagenomics classification
    algorithm, to estimate the number of reads originating from each species
    present in a sample. Kraken classifies reads to the best matching
    location in the taxonomic tree, but does not estimate abundances of
    species. We use the Kraken database itself to derive probabilities that
    describe how much sequence from each genome is identical to other genomes
    in the database, and combine this information with the assignments for a
    particular sample to estimate abundance at the species level, the genus
    level, or above. Combined with the Kraken classifier, Bracken produces
    accurate species- and genus-level abundance estimates even when a sample
    contains two or more near-identical species.

    References
    ----------
    Lu, Jennifer, et al. "Bracken: estimating species abundance in
    metagenomics data." PeerJ Computer Science 3 (2017): e104.

    Notes
    -----
    GitHub  : https://github.com/jenniferlu717/Bracken
    Docs    : http://ccb.jhu.edu/software/bracken/index.shtml
    Paper   : https://doi.org/10.7717/peerj-cs.104

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Bracken
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    if self.soft.prev != 'kraken2':
        sys.exit('[bracken] Can only be run after kraken2')
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam):
            continue
        for (db, k2) in inputs:
            report = '%s/report.tsv.gz' % k2
            to_dos = status_update(self, tech, [report])

            db_path = get_bracken_db(self, db)
            out = '/'.join([self.dir, tech, self.sam_pool, db])
            self.outputs['dirs'].append(out)
            self.outputs['outs'].setdefault((tech, sam), []).extend(out)
            if self.config.force or to_do('%s/results.tsv.gz' % out):
                cmd = bracken_cmd(self, tech, db_path, report, out)
                if to_dos:
                    self.outputs['cmds'].setdefault((tech,), []).append(False)
                else:
                    self.outputs['cmds'].setdefault((tech,), []).append(cmd)
                io_update(self, i_f=report, o_d=out, key=tech)
                self.soft.add_status(tech, sam, 1, group=(db + ' / ' + k2))
            else:
                self.soft.add_status(tech, sam, 0, group=(db + ' / ' + k2))


def get_metaxa2_db(
        self,
        db: str
) -> None:
    """Get the path the Metaxa2 database passed by the user.

    Notes
    -----
    If the database is not valid, the pipeline will stop with a useful hint.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Metaxa2
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
    db : str
        Name of the Metaxa2 database

    Returns
    -------
    str
        Path to the Metaxa2 database
    """
    if db == 'default':
        return self.databases.paths.get('metaxa2', {})
    elif db in self.databases.paths:
        for path in ['', '/databases']:
            db_path = '%s%s/metaxa2' % (self.databases.paths[db], path)
            if self.config.dev:
                return db_path
            if not isfile('%s/hash.k2d' % db_path):
                nested = glob.glob('%s/*/hash.k2d' % db_path)
                if nested:
                    return dirname(nested[0])
    sys.exit('[metaxa2] Database not found: %s' % db)


def metaxa2(self) -> None:
    """Improved Identification and Taxonomic Classification of Small and
    Large Subunit rRNA in Metagenomic Data.

    References
    ----------
    Bengtsson‐Palme, Johan, et al. "METAXA2: improved identification and
    taxonomic classification of small and large subunit rRNA in metagenomic
    data." Molecular ecology resources 15.6 (2015): 1403-1414.

    Notes
    -----
    Docs    : https://microbiology.se/software/metaxa2
    Paper   : https://doi.org/10.1111/1755-0998.12399

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for Metaxa2
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
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam):
            continue
        if self.soft.prev == 'kneaddata':
            split_term = '_1.fastq'
        else:
            split_term = '_R1.fastq'

        params = tech_params(self, tech)
        for db in params['databases']:
            db_path = get_metaxa2_db(self, db)

    for db, db_path in self.databases.paths['metaxa2'].items():
        dir_db = self.dir + '/' + db
        self.outputs['dirs'].append(dir_db)
        rad = dir_db + '/' + fas
        summary = '%s.summary.txt' % rad
        taxonomy = '%s.taxonomy.txt' % rad
        reltax = '%s.reltax.txt' % rad
        cmd = ''
        if self.config.force or to_do(summary):
            cmd += 'metaxa2'
            for idx, fastq in enumerate(self.inputs[self.sam_pool]):
                cmd += ' -%s %s' % ((idx + 1), fastq)
            cmd += ' -o %s' % rad
            if db_path:
                cmd += ' -p %s/HMMs' % db_path
                cmd += ' -d %s/blast' % db_path
            if idx:
                if fastq.endswith('fasta'):
                    cmd += ' -f paired-fasta'
                else:
                    cmd += ' -f paired-end'
            else:
                if fastq.endswith('fasta'):
                    cmd += ' -f fasta'
                else:
                    cmd += ' -f fastq'
            if fastq.endswith('gz'):
                cmd += ' -z gzip'
            cmd += ' --megablast %s' % self.soft.params['megablast']
            cmd += ' --align %s' % self.soft.params['align']
            cmd += ' --mode %s' % self.soft.params['mode']
            cmd += ' --plus %s' % self.soft.params['plus']
            cmd += ' --cpu %s' % self.soft.params['cpus']
            if self.soft.params['graphical']:
                cmd += ' --graphical'
            if self.soft.params['reltax']:
                cmd += ' --reltax'
            if db == 'wol_ssu_g':
                cmd += ' --taxlevel 8'
            else:
                cmd += ' --taxlevel 7'
            cmd += '\n'
            io_update(self, o_d=dir_db)
        else:
            io_update(self, i_f=taxonomy)
        redist_rad = '%s.redist' % rad
        redist_tax = '%s.taxonomy.summary.txt' % redist_rad
        if self.config.force or to_do(redist_tax):
            cmd += '\nmetaxa2_ttt'
            cmd += ' -i %s' % taxonomy
            cmd += ' -o %s' % redist_rad
            cmd += ' -r 0.8'
            cmd += ' -d 0.7'
            io_update(self, o_f=redist_tax)
        if cmd:
            self.outputs['cmds'].setdfault((tech,), []).append(cmd)
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(
            [summary, taxonomy, reltax, redist_tax])


def phyloflash(self):
    """phyloFlash is a pipeline to rapidly reconstruct the SSU rRNAs and
    explore phylogenetic composition of an Illumina (meta)genomic or
    transcriptomic dataset.

    References
    ----------
    Gruber-Vodicka, H.R., Seah, B.K. and Pruesse, E., 2020. phyloFlash:
    rapid small-subunit rRNA profiling and targeted assembly from metagenomes.
    Msystems, 5(5), pp.e00920-20.

    Notes
    -----
    Paper   : https://doi.org/10.1128/mSystems.00920-20
    GitHub  : https://github.com/HRGV/phyloFlash
    Docs    : http://hrgv.github.io/phyloFlash

    Parameters
    ----------
    self
    """
    pass


def kaiju(self):
    """Kaiju is a program for sensitive taxonomic classification of
    high-throughput sequencing reads from metagenomic whole genome
    sequencing or metatranscriptomics experiments.

    References
    ----------
    Menzel, P., Ng, K.L. and Krogh, A., 2016. Fast and sensitive taxonomic
    classification for metagenomics with Kaiju. Nature communications, 7(1),
    p.11257.

    Notes
    -----
    Paper   : https://doi.org/10.1038/ncomms11257
    GitHub  : https://github.com/bioinformatics-centre/kaiju
    Docs    : https://kaiju.binf.ku.dk/

    Parameters
    ----------
    self
    """
    pass


def closedref(self):
    """

    Parameters
    ----------
    self

    Returns
    -------

    """
    pass


def ngless(self):
    """Ngless is a domain-specific language for NGS (
    next-generation sequencing data) processing.

    References
    ----------
    Coelho, L.P., Alves, R., Monteiro, P., Huerta-Cepas, J., Freitas,
    A.T. and Bork, P., 2019. NG-meta-profiler: fast processing of metagenomes
    using NGLess, a domain-specific language. Microbiome, 7(1), pp.1-10.

    Notes
    -----
    GitHub  : https://github.com/ngless-toolkit/ngless
    Docs    : https://ngless.embl.de/
    Paper   : https://doi.org/10.1186/s40168-019-0684-8

    Parameters
    ----------
    self
    """
    pass


def ngmetaprofiler(self):
    """Collection of predefined ngless pipelines for processing shotgun
    metagenomes, incl:
    - human-gut.ngl for human gut samples
    - marine.ngl for marine samples
    - mouse-gut.ngl for mouse gut samples
    - dog-gut.ngl for dog gut samples
    - pig-gut.ngl for pig gut samples

    References
    ----------
    Coelho, L.P., Alves, R., Monteiro, P., Huerta-Cepas, J., Freitas,
    A.T. and Bork, P., 2019. NG-meta-profiler: fast processing of metagenomes
    using NGLess, a domain-specific language. Microbiome, 7(1), pp.1-10.

    Notes
    -----
    GitHub  : https://github.com/ngless-toolkit/ng-meta-profiler
    Paper   : https://doi.org/10.1186/s40168-019-0684-8

    Parameters
    ----------
    self
    """
    pass


def motus(self):
    """The mOTU profiler is a computational tool that estimates relative
    taxonomic abundance of known and currently unknown microbial community
    members using metagenomic shotgun sequencing data.

    References
    ----------
    Milanese, A., Mende, D.R., Paoli, L., Salazar, G., Ruscheweyh, H.J.,
    Cuenca, M., Hingamp, P., Alves, R., Costea, P.I., Coelho, L.P. and
    Schmidt, T.S., 2019. Microbial abundance, activity and population genomic
    profiling with mOTUs2. Nature communications, 10(1), pp.1-11.

    Notes
    -----
    GitHub  : https://github.com/motu-tool/mOTUs
    Paper   : https://doi.org/10.1038/s41467-019-08844-4
    Docs    : https://motu-tool.org/index.html
    Wiki    : https://github.com/motu-tool/mOTUs/wiki

    Parameters
    ----------
    self
    """
    pass