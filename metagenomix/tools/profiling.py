# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import sys
import pkg_resources
from os.path import basename, isdir, isfile, splitext

from metagenomix._io_utils import check_min_lines_count, io_update
from metagenomix.tools.alignment import get_alignment_cmd

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')


def shogun_append_cmd(self, cmds: list, tab: str):
    """Add or not the current SHOGUN command to the list of commands to run
    depending on whether the file already exists or exists but is only a header.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
    cmds : list
        List of current SHOGUN command lines
    tab : str
        Path to the output table of the current command
    """
    if self.config.force:
        self.outputs['cmds'].extend(cmds)
    elif not isfile(tab) or not check_min_lines_count(tab):
        cmd = 'file="%s"\n' % tab
        cmd += 'if [ -f "$file" ]\n'
        cmd += 'then\n'
        cmd += "FILESIZE=`ls -l $file | cut -d ' ' -f5`\n"
        cmd += "if [[ $FILESIZE -lt 1000 ]]\n"
        cmd += 'then\n'
        cmd += '\n'.join(cmds)
        cmd += 'fi\n'
        self.outputs['cmds'].append(cmd)


def shogun_redistribute(self, db: str, aligner: str, tax_norm: str,
                        db_path: str) -> None:
    """Get the SHOGUN command to redistribute a taxonomic profile
    at various taxonomic levels.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
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
        self.outputs['outs'].setdefault((db, aligner), []).append(redist)
        io_update(self, o_f=redist)
    redist_out = '%s.redist.strain.tsv' % splitext(tax_norm)[0]
    shogun_append_cmd(self, redist_cmds, redist_out)


def shogun_assign_taxonomy(aligner: str, ali: str, tax: str, db_path: str,
                           sub_db: str) -> str:
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


def shogun_normalize(input_fp: str, output_fp: str) -> str:
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


def get_ali_cmd(self, fasta: str, out: str, db_path: str, aligner: str) -> str:
    """Get the SHOGUN alignment command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
        .databases
            All databases class instance
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
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -o %s' % out
    return cmd


def get_alignment_basename(aligner):
    if aligner == 'bowtie2':
        ali_base = 'alignment.bowtie2.sam'
    elif aligner == 'burst':
        ali_base = 'alignment.burst.b6'
    else:
        raise ValueError('No output basename for aligner "%s"' % aligner)
    return ali_base


def get_orients(inputs: list):
    if len(inputs) == 1:
        orients = ['']
    elif len(inputs) == 2:
        orients = ['1', '2']
    elif len(inputs) == 3:
        orients = ['', '1', '2']
    else:
        raise IOError('Input to shogun must be 1, 2 or 3 fasta files')
    return orients


def get_combine_cmd(self, out_dir: str, combine_cmds: list) -> list:
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
    out_dir : str
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
    orients = get_orients(self.inputs[self.sam])
    for pdx, path_ in enumerate(sorted(self.inputs[self.sam])):
        io_update(self, i_f=path_)
        path = path_
        if path_.endswith('fastq') or path_.endswith('fastq.gz'):
            # replace non-fasta by fasta extensions
            path = path_.replace('.fastq', '.fasta').replace('.gz', '')
            # Prepare the extraction command and collect it
            if isfile(path_) and not isfile(path):
                to_fasta_cmd = 'seqtk seq -A %s > %s' % (path_, path)
                combine_cmds.append(to_fasta_cmd)

        path_out = '%s/%s' % (out_dir, basename(path))
        edit_fasta = 'python3 %s/scripts/fasta4shogun.py -i %s -o %s -s %s' % (
            RESOURCES, path, path_out, self.sam)
        orient = orients[pdx]
        if orient:
            edit_fasta += ' -r %s' % orient
        combine_cmds.append(edit_fasta)
        fastas.append(path_out)

    return fastas


def combine_inputs(self, out_dir: str, combine_cmds: list) -> str:
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
    out_dir : str
        Path to pipeline output folder for SHOGUN
    combine_cmds : list
        Commands to turn .fastq or .fastq.gz files into .fasta files

    Returns
    -------
    fasta : str
        Path to the combined sequences fasta file
    """
    fastas = get_combine_cmd(self, out_dir, combine_cmds)
    fasta = '%s/combined.fasta' % out_dir
    combine_cmds.extend([
        'cat %s > %s' % (' '.join(fastas), fasta), 'rm %s' % ' '.join(fastas)])
    return fasta


def align_cmds(self, fasta: str, out: str, db: str, aligner: str,
               ali_cmds: list) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
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
    ali = '%s/%s' % (out_dir, get_alignment_basename(aligner))

    db_path = '%s/shogun' % self.databases.paths[db]
    cmd = get_ali_cmd(self, fasta, out_dir, db_path, aligner)
    cmd = get_alignment_cmd([fasta], cmd, ali)
    ali_cmds.append(cmd)
    if cmd.startswith('shogun'):
        io_update(self, i_f=ali)

    self.outputs['dirs'].append(out_dir)
    self.outputs['outs'].setdefault((db, aligner), []).append(ali)
    io_update(self, o_d=out_dir)

    return ali


def format_sam(sam_: str, ali_cmds: list, sample: str) -> str:
    sam = '%s_formatted.sam' % splitext(sam_)[0]
    cmd = '%s/scripts/sam4shogun.py -i %s -o %s -s %s' % (
        RESOURCES, sam_, sam, sample)
    ali_cmds.append(cmd)
    return sam


def get_dir(out_dir: str, db: str, aligner: str) -> str:
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


def shogun_assign_normalize(self, tab: str, aligner: str, ali: str,
                            db_path: str, sub_db: str):
    """Get the assignment and the normalization SHOGUN commands for
    the taxonomy or functional classifications.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
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
    if self.config.force or not isfile(tab):
        cmd += shogun_assign_taxonomy(aligner, ali, tab, db_path, sub_db)
    if self.config.force or not isfile(tab_norm):
        cmd += shogun_normalize(tab, tab_norm)
    if cmd:
        shogun_append_cmd(self, list([cmd]), tab_norm)


def get_paths(self, out: str, aligner: str, db: str, tax_fun: str) -> tuple:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
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
    tax_norm : str
        Path to the output table (relative abundances)
    """
    out_dir = get_dir(out, db, aligner)
    self.outputs['dirs'].append(out_dir)

    tab = '%s/%s.tsv' % (out_dir, tax_fun)
    tab_norm = '%s/%s_norm.tsv' % (out_dir, tax_fun)
    self.outputs['outs'].setdefault((db, aligner), []).extend([tab, tab_norm])
    io_update(self, i_d=out_dir, o_f=[tab, tab_norm])
    return tab, tab_norm


def shogun_coverage_cmd(self, ali: str, db_path: str, level: str,
                        cov_tab: str) -> None:
    """Get the SHOGUN command that calculates the coverage per taxon.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
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
    shogun_append_cmd(self, list([cmd]), cov_tab)


def shogun_coverage(self, out: str, aligner: str, ali: str, db: str) -> None:
    db_path = '%s/shogun' % self.databases.paths[db]
    if aligner == 'burst':
        cov_tab = '%s/coverage.tsv' % out
        self.outputs['outs'].setdefault((db, aligner), []).append(cov_tab)
        if self.config.force or not isfile(cov_tab):
            shogun_coverage_cmd(self, ali, db_path, 'strain', cov_tab)
            io_update(self, o_f=cov_tab)


def shogun_taxonomy(self, out: str, aligner: str, ali: str, db: str) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
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
    tab, norm = get_paths(self, out, aligner, db, 'taxonomy')
    shogun_assign_normalize(self, tab, aligner, ali, db_path, '')
    shogun_redistribute(self, db, aligner, norm, db_path)
    return norm


def shogun_function(self, out: str, aligner: str, ali: str, db: str,
                    norm: str) -> None:
    db_path = '%s/shogun' % self.databases.paths[db]
    for sub in ['kegg', 'refseq', 'uniprot']:
        fun, f_norm = get_paths(self, out, aligner, db, 'function-%s' % sub)
        shogun_assign_normalize(self, fun, aligner, ali, db_path,
                                '/functions-%s' % sub)
        shogun_redistribute(self, db, aligner, f_norm, db_path)

    out_dir = out + '/functional'
    for level in ['genus', 'species']:
        shogun_functional(self, norm, fun, out_dir, level)


def shogun_functional(self, norm: str, db_path: str, out_dir: str,
                      level: str) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .config
            Configurations
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
    shogun_append_cmd(self, list([cmd]), out)


def shogun(self) -> None:
    """Get the SHOGUN commands that consist of the classifications using the
    databases built for shogun and the alignment in shogun, or not as this
    may be run after bowtie2, in which case the classifications are performed
    on the alignment obtained before, provided that it was done on a database
    formatted for SHOGUN.

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
    self.outputs['outs'] = dict({})
    combine_cmds, ali_cmds = [], []

    out = self.dir + '/' + self.sam
    self.outputs['dirs'].append(out)
    io_update(self, o_d=out)

    if self.soft.prev == 'bowtie2':
        for (db, pairing), sam_ in self.inputs[self.sam].items():
            ali = format_sam(sam_, ali_cmds, self.sam)
            shogun_taxonomy(self, out, '', ali, db)
    elif self.soft.params['databases']:
        fasta = combine_inputs(self, out, combine_cmds)
        for db, aligners in self.soft.params['databases'].items():
            for aligner in aligners:
                ali = align_cmds(self, fasta, out, db, aligner, ali_cmds)
                norm = shogun_taxonomy(self, out, aligner, ali, db)
                shogun_coverage(self, out, aligner, ali, db)
                shogun_function(self, out, aligner, ali, db, norm)

    if self.outputs['cmds']:
        self.outputs['cmds'] = combine_cmds + ali_cmds + self.outputs['cmds']


def woltka_aligments(self) -> dict:
    """Get the alignment paths per sample.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
        .outputs : dict
            All outputs

    Returns
    -------
    alis : dict
        Alignments
    """
    alis = {}
    for sample, sam_inputs in self.inputs.items():
        for db_pairing, sam in sam_inputs.items():
            if db_pairing[0] == 'wol':
                alis[sample] = sam
                io_update(self, i_f=sam)
    return alis


def woltka_write_map(self, alis: dict) -> str:
    """Write the mapping file that servers as input to Woltka.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to the output folder
        .outputs: dict
            All outputs
    alis : dict
        Alignments

    Returns
    -------
    woltka_map : str
        Path to the output woltka samples file.
    """
    woltka_map = '%s/samples.map' % self.dir
    for idx, sam in enumerate(alis.keys()):
        echo = 'echo -e "%s\\t%s"' % (sam, alis[sam])
        if idx:
            self.outputs['cmds'].append('%s >> %s' % (echo, woltka_map))
        else:
            self.outputs['cmds'].append('%s > %s' % (echo, woltka_map))
    return woltka_map


def woltka_gotus(self, woltka_map: str) -> str:
    """Get the Woltka command for the gotu classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    woltka_map : str
        Path to the output woltka samples file.

    Returns
    -------
    genomes_out : str
        Path to the output genomes per sample table.
    """
    genomes_out = '%s/genomes.tsv' % self.dir
    if not isfile(genomes_out):
        cmd = 'woltka gotu'
        cmd += ' -i %s' % woltka_map
        cmd += ' -o %s' % genomes_out
        self.outputs['cmds'].append(cmd)
        self.outputs['outs'].append(genomes_out)
        io_update(self, o_f=genomes_out)
    return genomes_out


def woltka_tax_cmd(self, woltka_map: str, database: str) -> str:
    """Get the taxonomic classification outputs and prepare the
    Woltka commands for this classification.

    Parameters
    ----------
    self : Commands class instance
        .out_dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    woltka_map : str
        Path to the Woltka input file.
    database : str
        WOL database path

    Returns
    -------
    tax_outmap : str
        Path to the folder containing the taxonomic maps.
    """
    tax_out, tax_outmap = '%s/taxa' % self.dir, '%s/taxmap' % self.dir
    tax_outputs = ['phylum', 'genus', 'species', 'none']
    tax_todo = []
    for tdx, tax_output in enumerate(tax_outputs):
        cur_tax_output = '%s/%s.tsv' % (tax_out, tax_output)
        self.outputs['cmds'].append(cur_tax_output)
        if not isfile(cur_tax_output):
            tax_todo.append(tax_output)
            io_update(self, o_f=cur_tax_output)

    taxid = '%s/taxonomy/taxid.map' % database
    nodes = '%s/taxonomy/nodes.dmp' % database
    names = '%s/taxonomy/names.dmp' % database
    io_update(self, i_f=[taxid, nodes, names])
    if len(tax_todo):
        cur_cmd = '\n# taxonomic\n'
        cur_cmd += 'woltka classify'
        cur_cmd += ' -i %s' % woltka_map
        cur_cmd += ' --map %s' % taxid
        cur_cmd += ' --nodes %s' % nodes
        cur_cmd += ' --names %s' % names
        cur_cmd += ' --rank %s' % ','.join(tax_outputs)
        cur_cmd += ' --add-rank'
        cur_cmd += ' --add-lineage'
        cur_cmd += ' --name-as-id'
        cur_cmd += ' --to-tsv'
        cur_cmd += ' --outmap %s' % tax_outmap
        cur_cmd += ' -o %s' % tax_out
        self.outputs['cmds'].append(cur_cmd)
    return tax_outmap


def woltka_classif_go(self, woltka_map: str, tax_outmap: str, database: str):
    """Get the taxonomic classification outputs and prepare the Woltka
    commands for this classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    woltka_map : str
        Path to the Woltka input file.
    tax_outmap : str
        Path to the folder containing the taxonomic maps.
    database : str
        WOL database path
    """
    coords = '%s/proteins/coords.txt.xz' % database
    uniref_map = '%s/function/uniref/uniref.map.xz' % database
    uniref_names = '%s/function/uniref/uniref.names.xz' % database
    io_update(self, i_f=[coords, uniref_map, uniref_names])
    go_rt = '%s/function/go' % database
    gos = ['process', 'function', 'component']
    woltka_fun_out = '%s/go' % self.dir
    for go in gos:
        cmd = '\n# %s [no stratification]\n' % go
        cmd += 'woltka classify'
        cmd += ' -i %s' % woltka_map
        cmd += ' --coords %s' % coords
        cmd += ' --map-as-rank'
        cmd += ' --rank %s' % go
        cmd += ' --map %s' % uniref_map
        cmd += ' --to-tsv'
        cur_map = '%s/%s.map.xz' % (go_rt, go)
        cmd += ' --map %s' % cur_map
        io_update(self, i_f=cur_map)
        cur_out = '%s/%s.tsv' % (woltka_fun_out, go)
        cmd += ' -o %s' % cur_out
        if not isfile(cur_out):
            self.outputs['cmds'].append(cmd)
            io_update(self, o_f=cur_out)
        self.outputs['outs'].append(cur_out)

    stratifications = ['phylum', 'family', 'genus', 'species']
    for stratification in stratifications:
        woltka_fun_out = '%s/go_%s' % (self.dir, stratification)
        for go in gos:
            cmd = '\n# %s [%s]\n' % (go, stratification)
            cmd += 'woltka classify'
            cmd += ' -i %s' % woltka_map
            cmd += ' --coords %s' % coords
            cmd += ' --map-as-rank'
            cmd += ' --rank %s' % go
            cmd += ' --stratify %s/%s' % (tax_outmap, stratification)
            cmd += ' --map %s' % uniref_map
            cmd += ' --to-tsv'
            cur_map = '%s/%s.map.xz' % (go_rt, go)
            cmd += ' --map %s' % cur_map
            go_out = '%s/%s.tsv' % (woltka_fun_out, go)
            cmd += ' -o %s' % go_out
            if not isfile(go_out):
                self.outputs['cmds'].append(cmd)
                io_update(self, o_f=go_out)
            self.outputs['outs'].append(go_out)


def woltka_classif_genes(self, woltka_map: str, database: str):
    """Get the Woltka commands for the gene-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    woltka_map : str
        Path to the Woltka input file.
    database : str
        WOL database path

    Returns
    -------
    genes : str
        Path to the genes classification output.
    """
    coords = '%s/proteins/coords.txt.xz' % database
    genes = '%s/wol_genes.biom' % self.dir
    if not isfile(genes):
        cmd = '\n# per gene\n'
        cmd += 'woltka classify'
        cmd += ' -i %s' % woltka_map
        cmd += ' --coords %s' % coords
        cmd += ' -o %s' % genes
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=genes)
    else:
        io_update(self, i_f=genes)
    self.outputs['outs'].append(genes)
    return genes


def woltka_classif_uniref(self, genes: str, database: str):
    """Get the Woltka commands for the uniref-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    genes : str
        Path to the genes classification output.
    database : str
        WOL database path

    Returns
    -------
    uniref : str
        Path to the uniref classification.
    """
    uniref_map = '%s/function/uniref/uniref.map.xz' % database
    uniref_names = '%s/function/uniref/uniref.names.xz' % database
    uniref = '%s/wol_uniref.biom' % self.dir
    if not isfile(uniref):
        cmd = 'woltka tools collapse'
        cmd += ' --input %s' % genes
        cmd += ' --map %s' % uniref_map
        cmd += ' --names %s' % uniref_names
        cmd += ' --output %s' % uniref
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=uniref)
    else:
        io_update(self, i_f=uniref)
    self.outputs['outs'].append(uniref)
    return uniref


def woltka_classif_eggnog(self, uniref: str, database: str):
    """Get the Woltka commands for the eggnog-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    uniref : str
        Path to the uniref classification.
    database : str
        WOL database path
    """
    biom = '%s/eggnog/wol_eggnog.biom' % self.dir
    if not isfile(biom):
        cmd = 'woltka tools collapse'
        cmd += '--input %s' % uniref
        cmd += ' --map %s/function/eggnog/eggnog.map.xz' % database
        cmd += ' --output %s\n\n' % biom
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=biom)
    else:
        io_update(self, i_f=biom)
    tsv = '%s.tsv' % splitext(biom)[0]
    if not isfile(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=tsv)
    self.outputs['outs'].extend([tsv, biom])


def woltka_classif_cazy(self, genes: str, database: str):
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
    genes : str
        Path to the genes classification.
    database : str
        WOL database path
    """
    cazy_map = '%s/function/cazy/3tools.txt' % database
    biom = '%s/cazy/wol_cazy.biom' % self.dir
    if not isfile(biom):
        cmd = 'woltka tools collapse'
        cmd += '--input %s' % genes
        cmd += ' --map %s' % cazy_map
        cmd += ' --output %s\n\n' % biom
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=biom)
    else:
        io_update(self, i_f=biom)
    tsv = '%s.tsv' % splitext(biom)[0]
    if not isfile(tsv):
        cmd = 'biom convert -i %s -o %s.tmp --to-tsv\n' % (biom, tsv)
        cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
        cmd += 'rm %s.tmp\n' % tsv
        self.outputs['cmds'].append(cmd)
        io_update(self, o_f=tsv)
    self.outputs['outs'].extend([biom, tsv])


def woltka_classif_metacyc(self, genes: str, database: str):
    """Get the Woltka commands for the metacyc-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    genes : str
        Path to the genes classification.
    database : str
        WOL database path
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
    files = [genes]
    woltka_fun_out = '%s/metacyc' % self.dir
    cmd = ''
    for idx, (level, names, maps) in enumerate(metacyc):
        input_fp = files[idx]
        if not input_fp.endswith('.biom'):
            input_biom = '%s.biom' % splitext(input_fp)[0]
            if not isfile(input_biom):
                cmd += 'biom convert'
                cmd += ' -i %s' % input_fp
                cmd += ' -o %s' % input_biom
                cmd += ' --to-hdf5'
                cmd += ' --table-type="OTU table"\n'
                io_update(self, o_f=input_biom)
            else:
                io_update(self, i_f=input_biom)
        else:
            input_biom = input_fp
        cmd += '\n# %s [no stratification]\n' % level
        biom = '%s/metacyc-%s_pre.biom' % (woltka_fun_out, level)
        if not isfile(biom):
            cmd += 'woltka tools collapse'
            if names:
                cmd += ' --input %s' % input_biom
                cmd += ' --names %s/%s' % (metacyc_dir, names)
            else:
                if level == 'super_pathway':
                    cmd += ' --input %s' % files[4]
                elif level == 'regulation':
                    cmd += ' --input %s' % files[2]
                elif level == 'regulator':
                    cmd += ' --input %s' % files[6]
                elif level == 'ec':
                    cmd += ' --input %s' % files[3]
                else:
                    cmd += ' --input %s' % files[4]
            cmd += ' --map %s/%s' % (metacyc_dir, maps)
            cmd += ' --output %s\n' % biom
            io_update(self, o_f=biom)
        else:
            io_update(self, i_f=biom)
        tsv = '%s/metacyc-%s_pre.tsv' % (woltka_fun_out, level)
        if not isfile(tsv):
            cmd += 'biom convert'
            cmd += ' -i %s' % biom
            cmd += ' -o %s.tmp' % tsv
            cmd += ' --to-tsv\n'
            cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
            cmd += 'rm %s.tmp\n' % tsv
            io_update(self, o_f=tsv)
        else:
            io_update(self, i_f=tsv)
        files.extend([biom, tsv])
    if cmd:
        self.outputs['cmds'].append(cmd)
    self.outputs['outs'].extend(files)


def woltka_classif_kegg(self, uniref: str, database: str):
    """Get the Woltka commands for the kegg-level classification.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for SHOGUN.
        .outputs: dict
            All outputs
    uniref : str
        Path to the uniref classification.
    database : str
        WOL database path
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
        ('ko-pathway-compound', 'compound_name.txt',
         'pathway-to-compound.txt', 'ko-pathway'),
        ('ko-pathway-disease', 'disease_name.txt',
         'pathway-to-disease.txt', 'ko-pathway'),
        ('ko-pathway-ko', 'ko_name.txt', 'pathway-to-ko.txt', 'ko-pathway'),
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
    cmd = ''
    files = []
    kegg_maps = '%s/kegg_queried' % self.dir
    for (level, name, maps, prev) in ko_names_maps:
        if maps:
            biom = '%s/kegg/kegg-%s_pre.biom' % (self.dir, level)
            tsv = '%s/kegg/kegg-%s_pre.tsv' % (self.dir, level)
            if not prev:
                if not isfile(tsv):
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % uniref
                    cmd += ' --names %s/function/kegg/%s' % (database, name)
                    cmd += ' --map %s/function/kegg/%s' % (database, maps)
                    cmd += ' --output %s\n' % biom

                    cmd += ' biom convert -i %s' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += ' tail -n +2 %s.tmp\n' % tsv
                    cmd += ' > %s\n' % tsv
                    cmd += ' rm %s.tmp\n' % tsv
                    io_update(self, o_f=[biom, tsv])
                else:
                    io_update(self, i_f=biom)
            else:
                input_fp = '%s/kegg/kegg-%s.biom' % (self.dir, level)
                if not isfile(tsv):
                    cmd += 'woltka tools collapse'
                    cmd += ' --input %s' % input_fp
                    if name:
                        cmd += ' --names %s/%s' % (kegg_maps, name)
                    cmd += ' --map %s/%s' % (kegg_maps, maps)
                    cmd += ' --output %s\n' % biom
                    cmd += 'biom convert -i %s\n' % biom
                    cmd += ' -o %s.tmp --to-tsv\n' % tsv
                    cmd += 'tail -n +2 %s.tmp > %s\n' % (tsv, tsv)
                    cmd += 'rm %s.tmp\n\n' % tsv
                    io_update(self, o_f=[biom, tsv])
                else:
                    io_update(self, i_f=biom)
        else:
            if not isfile('%s/kegg_info.txt' % kegg_maps):
                cmd += 'cd %s\n' % kegg_maps
                cmd += 'cp %s %s/%s\n' % (tsv, kegg_maps, basename(tsv))
                kegg_query = '%s/wol/kegg_query.py' % RESOURCES
                cmd += 'python3 %s %s\n' % (kegg_query, basename(tsv))
                io_update(self, o_d=kegg_maps)
            else:
                io_update(self, i_d=kegg_maps)
    if cmd:
        self.outputs['cmds'].append(cmd)
    self.outputs['outs'].extend(files)


def woltka(self) -> None:
    """Prepare the command and collect outputs, io and dirs for Woltka.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .databases
            All databases class instance
    """
    database = self.databases.paths['wol']
    alis = woltka_aligments(self)
    woltka_map = woltka_write_map(self, alis)
    woltka_gotus(self, woltka_map)
    tax_outmap = woltka_tax_cmd(self, woltka_map, database)
    woltka_classif_go(self, woltka_map, tax_outmap, database)
    genes = woltka_classif_genes(self, woltka_map, database)
    uniref = woltka_classif_uniref(self, genes, database)
    woltka_classif_eggnog(self, uniref, database)
    woltka_classif_cazy(self, genes, database)
    woltka_classif_metacyc(self, genes, database)
    woltka_classif_kegg(self, uniref, database)


def get_midas_cmd(
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


def midas_species(self, focus_dir, db, tax) -> None:
    if db == self.databases.paths['midas']:
        species_out = '%s/species' % focus_dir
        species_profile = '%s/species_profile.txt' % species_out
        if not self.config.force and isfile(species_profile):
            io_update(self, i_d=species_out)
        else:
            self.outputs['cmds'].append(get_midas_cmd(self, focus_dir, db, tax))
            io_update(self, i_f=self.inputs[self.sam], o_d=species_out)
        self.outputs['outs'].append(species_out)
        self.outputs['dirs'].append(species_out)


def midas_genus(self, focus_dir, genes_out, db, tax, select):
    if self.config.force or not isfile('%s/readme.txt' % genes_out):
        self.outputs['cmds'].append(
            get_midas_cmd(self, focus_dir, db, tax, select))
        io_update(self, o_d=genes_out)
    self.outputs['outs'].append(genes_out)
    self.outputs['dirs'].append(genes_out)


def midas_snps(self, focus_dir, snps_out, db, tax, select):
    if self.config.force or not isfile('%s/readme.txt' % snps_out):
        self.outputs['cmds'].append(
            get_midas_cmd(self, focus_dir, db, tax, select))
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
        midas_species(self, focus_dir, db, 'species')

        select = set(get_species_select(self, db, species_list))

        genes_out = '%s/genes' % focus_dir
        midas_genus(self, focus_dir, genes_out, db, 'genus', select)

        snps_out = '%s/snps' % focus_dir
        midas_snps(self, focus_dir, snps_out, db, 'snps', select)


def get_kraken2_db(self, db):
    if db == 'default':
        db_path = self.databases.paths['kraken2']
    elif db in self.databases.paths:
        db_path = '%s/kraken2' % self.databases.paths[db]
        if not self.config.dev and not isdir(db_path):
            sys.exit('[kraken2] Not a database: %s' % db_path)
    else:
        sys.exit('[kraken2] Database not found: %s' % db)
    return db_path


def get_kraken2_cmd(self, out: str, db_path: str, report: str, result: str):
    cmd = 'kraken2 '
    cmd += ' -db %s' % db_path
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --report %s/report.tsv' % report
    cmd += ' --confidence %s' % self.soft.params['confidence']
    if len(self.inputs[self.sam]) > 1:
        unclass = ['%s/unclassified_%s.fastq' % (out, r) for r in [1, 2]]
        classif = ['%s/classified_%s.fastq' % (out, r) for r in [1, 2]]
        cmd += ' --unclassified-out %s/unclassified#.fastq' % out
        cmd += ' --classified-out %s/classified#.fastq' % out
        cmd += ' --paired'
    else:
        unclass = ['%s/unclassified.fastq' % out]
        classif = ['%s/classified.fastq' % out]
        cmd += ' --unclassified-out %s/unclassified.fastq' % out
        cmd += ' --classified-out %s/classified.fastq' % out
    if self.inputs[self.sam][0].endswith('.gz'):
        cmd += ' --gzip-compressed'
        io_update(self, o_f=([
            '%s.gz' % x for x in unclass] + ['%s.gz' % x for x in classif]))
    else:
        io_update(self, o_f=(unclass + classif))
    cmd += ' %s > %s' % (' '.join(self.inputs[self.sam]), result)
    return cmd


def kraken2(self) -> None:
    """Create command lines for kraken2.

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
    for db in self.soft.params['databases']:
        out = '%s/%s/%s' % (self.dir, self.sam, db)
        self.outputs['dirs'].append(out)
        report = '%s/report.tsv' % out
        result = '%s/result.tsv' % out
        if self.config.force or not isfile(result):
            db_path = get_kraken2_db(self, db)
            cmd = get_kraken2_cmd(self, out, db_path, report, result)
            self.outputs['outs'].append(result)
            self.outputs['cmds'].append(cmd)
            io_update(self, i_f=self.inputs[self.sam], o_f=[report, result])


def metaxa2(self) -> None:
    if self.soft.prev == 'kneaddata':
        fas = basename(self.inputs[self.sam][0]).split('_1.fastq')[0]
    else:
        fas = basename(self.inputs[self.sam][0]).split('_R1.fastq')[0]

    databases = {'greengenes': '/home/flejzerowicz/databases/metaxa2/gg',
                 'wol_ssu': '/home/flejzerowicz/databases/metaxa2/wl',
                 'wol_ssu_g': '/home/flejzerowicz/databases/metaxa2/wl_g'}
    io_update(self, i_f=self.inputs[self.sam])
    for db, db_path in databases.items():
        dir_db = self.dir + '/' + db
        if not isdir(dir_db):
            os.makedirs(dir_db)
        rad = dir_db + '/' + fas
        summary = '%s.summary.txt' % rad
        taxonomy = '%s.taxonomy.txt' % rad
        reltax = '%s.reltax.txt' % rad
        cmd = ''
        if self.config.force or not isfile(summary):
            cmd += 'metaxa2'
            for idx, fastq in enumerate(self.inputs[self.sam]):
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
            io_update(self, o_f=[summary, taxonomy, reltax])
        else:
            io_update(self, i_f=taxonomy)
        self.outputs['outs'].extend([summary, taxonomy, reltax])

        redist_rad = '%s.redist' % rad
        redist_taxonomy = '%s.taxonomy.summary.txt' % redist_rad
        if self.config.force or not isfile(redist_taxonomy):
            cmd += '\nmetaxa2_ttt'
            cmd += ' -i %s' % taxonomy
            cmd += ' -o %s' % redist_rad
            cmd += ' -r 0.8'
            cmd += ' -d 0.7'
            io_update(self, o_f=redist_taxonomy)

        if cmd:
            self.outputs['cmds'].append(cmd)
        self.outputs['outs'].append(redist_taxonomy)
