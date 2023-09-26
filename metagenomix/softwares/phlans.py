# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import dirname
from metagenomix.core.parameters import tech_params
from metagenomix._io_utils import (io_update, tech_specificity,
                                   to_do, status_update)


def get_read_count(
        self,
        tech: str,
        sam: str,
) -> str:
    """Get the number of reads in the current sample,
    or an empty string if the counting did not happen yet.

    Parameters
    ----------
    self : Commands class instance
        .softs : dict of software classes
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam: str
        Sample name

    Returns
    -------
    reads : str
        Number of reads
    """
    reads = ''
    if 'count' in self.softs:
        reads_fps = self.softs['count'].outputs.get((sam, tech), [])
        if reads_fps and not sum([to_do(x) for x in reads_fps]):
            reads_fp = reads_fps.replace('${SCRATCH_FOLDER}', '')
            with open(reads_fp) as f:
                for line in f:
                    reads = line[-1].strip().split(',')[-1]
                    break
        else:
            t = "marker_ab_table"
            if t in tech_params(self, tech)['t']:
                self.soft.messages.add('Run "count" to use read counts for' + t)
    return reads


def get_tax_levs(params):
    """Get taxa levels.

    Parameters
    ----------
    params : dict
        Parameters

    Returns
    -------
    tax_levs : list
        List of taxa levels
    """
    if 'a' in params:
        tax_levs = ['a']
    else:
        tax_levs = params['tax_lev']
    return tax_levs


def get_out_fp(self, tech, t, reads):
    """Get te output folder path.

    Parameters
    ----------
    self
    tech : str
        Technology
    t : str
        Value to the `t` argument
    reads : str
        Whether a reads file exists

    Returns
    -------
    out : str
        Path to the output folder
    """
    if t == 'marker_ab_table' and reads:
        out = '%s/%s/%s_reads' % (self.dir, tech, t)
    else:
        out = '%s/%s/%s' % (self.dir, tech, t)
    return out


def analyses(
        self,
        tech: str,
        sam: str,
        bowtie2out: str,
        cmds: list
) -> None:
    """Get the taxonomic analyses Metaphlan command.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam : str
        Name of the current sample
    bowtie2out : str
        Input type was a bowtie2 output of the same tool
    cmds : list
        Commands per tech
    """
    reads = get_read_count(self, tech, sam)
    params = tech_params(self, tech)
    tax_levs = get_tax_levs(params)

    for t in params['t']:
        out = get_out_fp(self, tech, t, reads)
        for tax_lev in tax_levs:
            rad = '%s/%s_t-%s' % (out, sam, tax_lev)
            cmd = 'metaphlan'
            cmd += ' %s' % bowtie2out
            cmd += ' -t %s' % t
            cmd += ' --tax_lev %s' % tax_lev
            cmd += ' --input_type bowtie2out'
            cmd += ' --sample_id %s' % sam
            cmd += ' --sample_id_key sample_name'
            cmd += ' --nproc %s' % params['cpus']
            for param in [
                'perc_nonzero', 'stat_q', 'min_cu_len', 'min_alignment_len',
            ]:
                cmd += ' --%s %s' % (param, params[param])
            for param in [
                'ignore_eukaryotes', 'ignore_bacteria', 'ignore_archaea',
                'add_viruses', 'avoid_disqm'
            ]:
                if params[param]:
                    cmd += ' --%s' % param
            if reads and t == 'marker_ab_table':
                cmd += ' --nreads %s' % reads
            elif t == 'marker_pres_table':
                cmd += ' --pres_th %s' % params['pres_th']
            elif t == 'clade_specific_strain_tracker':
                for strain in self.config.strains:
                    cms = cmd
                    cms += ' --clade %s' % strain.replace(' ', '_')
                    clade_out = '%s_%s.tsv' % (rad, strain.replace(' ', '_'))
                    cms += ' --output_file %s' % clade_out
                    if params['min_ab']:
                        cms += ' --min_ab %s' % params['min_ab']
                    if self.config.force or to_do(clade_out):
                        cmds.append(cmd)
                        io_update(self, o_f=clade_out, key=tech)
                        self.soft.add_status(
                            tech, self.sam_pool, 1, message=t, genome='strain')
                    else:
                        self.soft.add_status(
                            tech, self.sam_pool, 0, message=t, genome='strain')
            else:
                ab_out = '%s.tsv' % rad
                cmd += ' --output_file %s' % ab_out
                if self.config.force or to_do(ab_out):
                    cmds.append(cmd)
                    io_update(self, o_f=ab_out, key=tech)
                    self.soft.add_status(tech, self.sam_pool, 1, message=t)
                else:
                    self.soft.add_status(tech, self.sam_pool, 0, message=t)


def profiling(
        self,
        tech: str,
        sam: str,
        inputs: list,
        bowtie2out: str,
        tmpdir: str,
        cmds: list
) -> tuple:
    """Get the taxonomic profiling metaphlan command.

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
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam : str
        Name of the current sample
    inputs : list
        Paths for the input fastq files
    bowtie2out : str
        Input type was a bowtie2 output of the same tool
    tmpdir : str
        Path to a temporary directory
    cmds : list
        Commands per tech

    Returns
    -------
    outs : list
        sam_out : str
            Alignment output
        profile_out : str
            Metaphlan taxonomic profile output
    to_dos : list
    """
    sam_out = '%s/%s/sams/%s.sam.bz2' % (self.dir, tech, sam)
    profile_out = '%s/%s/profiles/%s.tsv' % (self.dir, tech, sam)

    params = tech_params(self, tech)
    if to_do(profile_out):
        io_update(self, o_f=profile_out, key=tech)
        cmd = 'metaphlan'
        if not to_do(bowtie2out):
            to_dos = status_update(self, tech, [bowtie2out])
            io_update(self, i_f=bowtie2out, key=tech)
            cmd += ' %s' % bowtie2out
            cmd += ' --input_type bowtie2out'
        else:
            to_dos = status_update(self, tech, inputs)
            io_update(self, i_f=inputs, key=tech)
            cmd += ' %s' % ','.join(inputs)
            cmd += ' --input_type fastq'
            cmd += ' --samout %s' % sam_out
            cmd += ' --bt2_ps %s' % params['bt2_ps']
            cmd += ' --min_mapq_val %s' % params['min_mapq_val']
            cmd += ' --read_min_len %s' % params['read_min_len']
            cmd += ' --bowtie2out %s' % bowtie2out
            cmd += ' --bowtie2db %s' % self.databases.paths['metaphlan']
        cmd += ' --tmp_dir %s' % tmpdir
        cmd += ' --output_file %s' % profile_out
        cmd += ' --nproc %s' % params['cpus']
        cmd += ' --sample_id_key %s' % sam
        cmds.append(cmd)
        self.soft.add_status(tech, self.sam_pool, 1)
    else:
        to_dos = []
        io_update(self, i_f=bowtie2out, key=tech)
        self.soft.add_status(tech, self.sam_pool, 0)
    outs = [sam_out, profile_out]
    return outs, to_dos


def metaphlan(self) -> None:
    """MetaPhlAn is a computational tool for profiling the composition of
    microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic
    shotgun sequencing data (i.e. not 16S) with species-level. With the newly
    added StrainPhlAn module, it is now possible to perform accurate
    strain-level microbial profiling.
    MetaPhlAn relies on ~1.1M unique clade-specific marker genes (the latest
    marker information file mpa_v296_CHOCOPhlAn_201901_marker_info.txt.bz2
    can be found here) identified from ~100,000 reference genomes (~99,
    500 bacterial and archaeal and ~500 eukaryotic).
    What's new in version 3:
    - New MetaPhlAn marker genes extracted with a newer version of ChocoPhlAn
      based on UniRef
    - Estimation of metagenome composed by unknown microbes with parameter
      --unknown_estimation
    - Automatic retrieval and installation of the latest MetaPhlAn database
      with parameter --index latest
    - Virus profiling with --add_viruses
    - Calculation of metagenome size for improved estimation of reads mapped
      to a given clade
    - Inclusion of NCBI taxonomy ID in the ouput file
    - CAMI (Taxonomic) Profiling Output Format included
    - Removal of reads with low MAPQ values

    References
    ----------
    Beghini, Francesco, et al. "Integrating taxonomic, functional,
    and strain-level profiling of diverse microbial communities with
    bioBakery 3." Elife 10 (2021): e65088.

    Notes
    -----
    GitHub  : https://github.com/biobakery/MetaPhlAn/tree/3.0/
    Docs    : http://segatalab.cibio.unitn.it/tools/metaphlan/index.html
    Paper   : https://doi.org/10.7554/eLife.65088

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .databases
            All databases
        .config
            Configurations
    """
    for (tech, sam), inputs in self.inputs[self.sam_pool].items():
        if tech_specificity(self, inputs, tech, sam, ['illumina']):
            continue

        tmpdir = '$TMPDIR/metaphlan_%s_%s' % (tech, self.sam_pool)
        bowtie2out = '%s/%s/bowtie2/%s.bowtie2.bz2' % (self.dir, tech, sam)

        cmds = []
        outs, to_dos = profiling(
            self, tech, sam, inputs, bowtie2out, tmpdir, cmds)
        self.outputs['outs'].setdefault(
            (tech, sam), []).extend([bowtie2out] + list(outs))
        self.outputs['dirs'] = [
            dirname(x) for x in self.outputs['outs'][(tech, sam)]]

        analyses(self, tech, sam, bowtie2out, cmds)
        if cmds:
            _cmd = ['mkdir -p %s\n' % tmpdir]
            if to_dos:
                self.outputs['cmds'][(tech,)] = [False]
            else:
                self.outputs['cmds'][(tech,)] = _cmd + cmds


def get_profile(self) -> list:
    """Get the taxonomic profiles if there are some, for humann to only run
    of the taxa from these profiles.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .databases
            All databases

    Returns
    -------
    profiles_dbs_outs : list
        Element are lists:
            profiles : list
                Paths to the taxonomic profile of the current sample.
            dbs : list
                Paths to the metaphlan database.
            outs : list
                Paths to the output folder.
    """
    profiles, dbs, outs = [], [], []
    db_path = self.databases.paths['humann']
    params = self.soft.params
    if params['profiles']:
        for profile_name, profile_fp in params['profiles'].items():
            # `profile_fp` must be "path/to/<sample>/profile_file.txt"
            profile_sams = {x.split('/')[-2]: x for x in glob.glob(profile_fp)}
            if self.sam_pool in profile_sams:
                profiles.append(profile_sams[self.sam_pool])
                io_update(self, i_f=profile_sams[self.sam_pool])
                dbs.append('')
                outs.append(
                    self.dir + '/profile_%s_%s_id%s_e%s_qC%s_sC%s_bypss%s' % (
                        profile_name,
                        params['uniref'],
                        params['id_thresh'],
                        params['evalue'],
                        params['query_coverage'],
                        params['subject_coverage'],
                        params['skip_translated']))
    else:
        profiles.append('')
        dbs.append(db_path)
        outs.append(self.dir + '/no_profile')
    profiles_dbs_outs = [profiles, dbs, outs]
    return profiles_dbs_outs


def get_cmd(
        self,
        profile: str,
        db: str,
        out: str,
        base: str,
        ali: str
) -> str:
    """Get the humann command line.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .soft.params
            Parameters for humann
    profile : str
        Path to the taxonomic profile of the current sample.
    db : str
        Path to the metaphlan database.
    out : str
        Path to the output folder.
    base : str
        Path to fastq file.
    ali : str
        Path to temporary, alignment sam file

    Returns
    -------
    cmd : str
        humann command line.
    """
    cmd = 'cat %s > %s\n' % (' '.join(self.inputs[self.sam_pool]), base)
    cmd += 'humann'
    cmd += ' --input %s' % base
    cmd += ' --output %s' % out
    cmd += ' --prescreen-threshold %s' % self.soft.params['prescreen_threshold']
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --search-mode %s' % self.soft.params['uniref']
    if len(db):
        cmd += ' --metaphlan-options '
        cmd += '" --mpa_pkl %s/mpa_v20_m200.pkl' % db
    cmd += ' --bowtie2db %s"' % db
    cmd += ' --identity-threshold %s' % self.soft.params['identity']
    cmd += ' --evalue %s' % self.soft.params['evalue']
    cmd += ' --translated-query-coverage-threshold'
    cmd += ' %s' % self.soft.params['subject_cov']
    cmd += ' --translated-query-coverage-threshold'
    cmd += ' %s' % self.soft.params['subject_cov']
    if profile:
        cmd += ' --taxonomic-profile %s' % profile
    if self.soft.params['skip_translated']:
        cmd += ' --bypass-translated-search'
    if not to_do(ali):
        cmd += ' -r'
    cmd += ' --nucleotide-database %s' % self.soft.params['nucleotide_db']
    cmd += ' --protein-database  %s\n' % self.soft.params['protein_db']
    return cmd


def renorm(
        input_fp: str,
        relab_fp: str
) -> str:
    """Normalize the reads counts into proportions.

    Parameters
    ----------
    input_fp : str
        Path to the input file.
    relab_fp : str
        Path to the output file.

    Returns
    -------
    cmd : str
        Normalization command line.
    """
    cmd = 'humann2_renorm_table'
    cmd += ' --input %s' % input_fp
    cmd += ' --output %s' % relab_fp
    cmd += ' -u relab\n'
    return cmd


def get_ali(
        self,
        o: str
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .outputs : dict
            All outputs
    o : str


    Returns
    -------
    ali : str
        Path to the alignment .sam file
    """
    ali = '%s/%s_tmp/%s_bowtie2_aligned.sam' % (o, self.sam_pool, self.sam_pool)
    if not self.config.force and not to_do(ali):
        io_update(self, i_f=ali)
    else:
        io_update(self, o_f=ali, o_d=o)
    return ali


def get_outputs(
        self,
        profile: str,
        db: str,
        out: str,
        ali: str,
        key: tuple
) -> tuple:
    """Get the profiling output files at genus, pathway and coverage levels.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .outputs : dict
            All outputs
        .config
            Configurations
    profile : str
        Path to the taxonomic profile of the current sample.
    db : str
        Path to the metaphlan database.
    out : str
        Path to the output folder.
    ali : str
        Path to temporary, alignment sam file
    key : tuple
        Technology and/or smaple name

    Returns
    -------
    gen : str
        Path to the genus abundance profile output file
    pwy : str
        Path to the pathway abundance profile output file
    cov : str
        Path to the pathway coverage profile output file
    """
    base = '%s/%s.fastq.gz' % (out, self.sam_pool)
    gen = '%s/%s_genefamilies.tsv' % (out, self.sam_pool)
    pwy = '%s/%s_pathabundance.tsv' % (out, self.sam_pool)
    cov = '%s/%s_pathcoverage.tsv' % (out, self.sam_pool)
    if not self.config.force and not sum([to_do(x) for x in [gen, pwy, cov]]):
        io_update(self, i_f=[gen, pwy])
    else:
        cmd = get_cmd(self, profile, db, out, base, ali)
        self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, o_f=[gen, pwy, cov])
    return gen, pwy, cov


def humann(self) -> None:
    """.

    References
    ----------


    Notes
    -----
    GitHub  :
    Docs    :
    Paper   :

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for humann
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters for humann
        .databases
            All databases
        .config
            Configurations
    """
    for profile, db, out in zip(*get_profile(self)):
        key = (tech, '_'.join(profile, db))
        ali = get_ali(self, out)
        gen, pwy, cov = get_outputs(self, ali, profile, db, out, key)
        gen_relab = '%s/%s_genefamilies_relab.tsv' % (out, self.sam_pool)
        pwy_relab = '%s/%s_pathabundance_relab.tsv' % (out, self.sam_pool)
        if self.config.force or to_do(gen_relab) or to_do(pwy_relab):
            cmd_gen = renorm(gen, gen_relab)
            cmd_pwy = renorm(pwy, pwy_relab)
            self.outputs['cmds'].setdefault(key, []).append(cmd_gen)
            self.outputs['cmds'].setdefault(key, []).append(cmd_pwy)
            io_update(self, o_f=[gen_relab, pwy_relab], key=key)
        self.outputs['dirs'].append(out)
        self.outputs['outs'].append([out, ali, gen, gen_relab, pwy,
                                     pwy_relab, cov])
