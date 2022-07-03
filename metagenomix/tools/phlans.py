# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import dirname, isfile
from metagenomix._io_utils import io_update


def metaphlan(self) -> None:
    """

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
    tmpdir = '$TMPDIR/metaphlan_%s' % self.sam
    tmp_cmd = ['mkdir -p %s\n' % tmpdir]

    bowtie2out = '%s/bowtie2/%s.bowtie2.bz2' % (self.dir, self.sam)

    outs = profiling(self, bowtie2out, tmpdir)
    self.outputs['outs'] = [bowtie2out] + list(outs)
    self.outputs['dir'] = [dirname(x) for x in self.outputs['outs']]

    reads = get_read_count(self)
    if reads:
        analyses(self, bowtie2out, reads)
    if self.outputs['cmds']:
        self.outputs['cmds'] = tmp_cmd + self.outputs['cmds']


def get_read_count(self) -> str:
    """Get the number of reads in the current sample,
    or an empty string if the counting did not happen yet.

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name

    Returns
    -------
    reads : str
        Number of reads.
    """
    reads = ''
    reads_fps = self.softs['count'].outputs
    if self.sam in reads_fps and isfile(reads_fps[self.sam][0]):
        with open(reads_fps[self.sam][0]) as f:
            for line in f:
                reads = line[-1].strip().split(',')[-1]
                break
    return reads


def analyses(self, bowtie2out: str, reads: str) -> None:
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
    bowtie2out : str
        Input type was a bowtie2 output of the same tool
    reads : str
        Number of reads in the current sample
    """
    tax_levs = ['a', 'k', 'p', 'c', 'o', 'f', 'g', 's']
    for analysis in ['rel_ab',
                     'rel_ab_w_read_stats',
                     'clade_profiles',
                     'marker_ab_table',
                     'marker_pres_table',
                     'clade_specific_strain_tracker']:
        out = '%s/%s' % (self.dir, analysis)
        for tax_lev in tax_levs:
            rad = '%s/%s_t-%s' % (out, self.sam, tax_lev)
            cmd = 'metaphlan'
            cmd += ' %s' % bowtie2out
            cmd += ' -t %s' % analysis
            cmd += ' --tax_lev %s' % tax_lev
            cmd += ' --input_type bowtie2out'
            cmd += ' --sample_id %s' % self.sam
            cmd += ' --sample_id_key sample_name'
            cmd += ' --nproc %s' % self.soft.params['cpus']
            if reads and analysis == 'marker_ab_table':
                cmd += ' --nreads %s' % reads
            if analysis == 'clade_specific_strain_tracker':
                for strain in self.config.strains:
                    strain_cmd = cmd
                    strain_cmd += ' --clade %s' % strain.replace(' ', '_')
                    clade_out = '%s_%s.tsv' % (rad, strain.replace(' ', '_'))
                    strain_cmd += ' --output_file %s' % clade_out
                    if self.config.force or not isfile(clade_out):
                        self.outputs['cmds'].append(strain_cmd)
                        io_update(self, o_f=clade_out)
            else:
                ab_out = '%s.tsv' % rad
                cmd += ' --output_file %s' % ab_out
                if self.config.force or not isfile(ab_out):
                    self.outputs['cmds'].append(cmd)
                    io_update(self, o_f=ab_out)


def profiling(self, bowtie2out: str, tmpdir: str) -> tuple:
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
    bowtie2out : str
        Input type was a bowtie2 output of the same tool
    tmpdir : str
        Path to a temporary directory

    Returns
    -------
    sam_out : str
        Alignment output
    profile_out : str
        Metaphlan taxonomic profile output
    """
    sam_out = '%s/sams/%s.sam.bz2' % (self.dir, self.sam)
    profile_out = '%s/profiles/%s.tsv' % (self.dir, self.sam)

    if isfile(profile_out):
        io_update(self, i_f=bowtie2out)
    else:
        io_update(self, o_f=profile_out)
        cmd = 'metaphlan'
        if isfile(bowtie2out):
            io_update(self, i_f=bowtie2out)
            cmd += ' %s' % bowtie2out
            cmd += ' --input_type bowtie2out'
        else:
            io_update(self, i_f=self.inputs[self.sam])
            cmd += ' %s' % ','.join(self.inputs[self.sam])
            cmd += ' --input_type fastq'
            cmd += ' --samout %s' % sam_out
            cmd += ' --bowtie2out %s' % bowtie2out
            cmd += ' --bowtie2db %s' % self.databases.paths['metaphlan']
        cmd += ' --tmp_dir %s' % tmpdir
        cmd += ' --output_file %s' % profile_out
        cmd += ' --nproc %s' % self.soft.params['cpus']
        cmd += ' --sample_id_key %s' % self.sam
        self.outputs['cmds'].append(cmd)

    return sam_out, profile_out


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
            if self.sam in profile_sams:
                profiles.append(profile_sams[self.sam])
                io_update(self, i_f=profile_sams[self.sam])
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


def get_cmd(self, profile: str, db: str, out: str, base: str, ali: str):
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
    cmd = 'cat %s > %s\n' % (' '.join(self.inputs[self.sam]), base)
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
    if isfile(ali):
        cmd += ' -r'
    cmd += ' --nucleotide-database %s' % self.soft.params['nucleotide_db']
    cmd += ' --protein-database  %s\n' % self.soft.params['protein_db']
    return cmd


def renorm(input_fp: str, relab_fp: str) -> str:
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


def get_ali(self, out: str):
    """

    Parameters
    ----------
    self : Commands class instance
        .sam : str
            Sample name
        .outputs : dict
            All outputs
    out : str


    Returns
    -------
    ali : str
        Path to the alignment .sam file
    """
    ali = '%s/%s_tmp/%s_bowtie2_aligned.sam' % (out, self.sam, self.sam)
    if not self.config.force and isfile(ali):
        io_update(self, i_f=ali)
    else:
        io_update(self, o_f=ali, o_d=out)
    return ali


def get_outputs(self, profile: str, db: str, out: str, ali: str) -> tuple:
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

    Returns
    -------
    gen : str
        Path to the genus abundance profile output file
    pwy : str
        Path to the pathway abundance profile output file
    cov : str
        Path to the pathway coverage profile output file
    """
    base = '%s/%s.fastq.gz' % (out, self.sam)
    gen = '%s/%s_genefamilies.tsv' % (out, self.sam)
    pwy = '%s/%s_pathabundance.tsv' % (out, self.sam)
    cov = '%s/%s_pathcoverage.tsv' % (out, self.sam)
    if not self.config.force and isfile(gen) and isfile(pwy) and isfile(cov):
        io_update(self, i_f=[gen, pwy])
    else:
        self.outputs['cmds'].append(get_cmd(self, profile, db, out, base, ali))
        io_update(self, o_f=[gen, pwy, cov])
    return gen, pwy, cov


def humann(self) -> None:
    """Create command lines for humann for the current database's species focus.

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
        ali = get_ali(self, out)
        gen, pwy, cov = get_outputs(self, ali, profile, db, out)
        gen_relab = '%s/%s_genefamilies_relab.tsv' % (out, self.sam)
        pwy_relab = '%s/%s_pathabundance_relab.tsv' % (out, self.sam)
        if self.config.force or not isfile(gen_relab) or not isfile(pwy_relab):
            self.outputs['cmds'].append(renorm(gen, gen_relab))
            self.outputs['cmds'].append(renorm(pwy, pwy_relab))
            io_update(self, o_f=[gen_relab, pwy_relab])
        self.outputs['dirs'].append(out)
        self.outputs['outs'].append([out, ali, gen, gen_relab,
                                     pwy, pwy_relab, cov])


def get_sample_to_marker(self, sam_dir, markers_dir):
    if self.config.force or len(glob.glob('%s/*.sam.bz2' % sam_dir)):
        cmd = 'sample2markers.py'
        cmd += ' -i %s/*.sam.bz2' % sam_dir
        cmd += ' -o %s' % markers_dir
        cmd += ' -n %s' % self.soft.params['cpus']
        self.outputs['cmds'].append(cmd)


def extract_markers(self, st, strain_dir, db_dir):
    if self.config.force or not isfile('%s/%s.fna' % (strain_dir, st)):
        cmd = 'extract_markers.py'
        cmd += ' -c %s' % st
        cmd += ' -o %s' % db_dir
        self.outputs['cmds'].append(cmd)


def strainphlan(self) -> None:
    """Create command lines for strainphlan.

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
    bt2_fp, sam_fp = self.inputs[list(self.inputs.keys())[0]][:2]
    sam_dir = dirname(sam_fp)
    tmp_dir = '$TMPDIR/strainphlan'
    db_dir = '%s/db_markers' % self.dir
    meta_dir = '%s/metadata' % self.dir
    markers_dir = '%s/consensus_markers' % self.dir
    wol_pd = self.databases.paths['wol']['taxonomy']

    self.outputs['outs'] = dict({})
    self.outputs['dirs'].extend([db_dir, meta_dir])
    self.outputs['cmds'].append('mkdir -p %s' % tmp_dir)
    io_update(self, i_d=sam_dir, o_d=[markers_dir, db_dir, meta_dir])

    get_sample_to_marker(self, sam_dir, markers_dir)

    for strain_group, strains in self.config.strains.items():
        strain_dir = '%s/%s' % (db_dir, strain_group)
        self.outputs['dirs'].append(strain_dir)
        for strain in [x.replace(' ', '_') for x in strains]:
            if strain[0] != 's':
                continue
            extract_markers(self, strain, strain_dir, db_dir)
            odir = '%s/output/%s' % (self.dir, strain_group)
            self.outputs['dirs'].append(odir)
            self.outputs['outs'][(strain_group, strain)] = odir
            io_update(self, o_d=odir)
            tree = '%s/RAxML_bestTree.%s.StrainPhlAn3.tre' % (odir, strain)
            if self.config.force or not isfile(tree):
                cmd = 'strainphlan'
                cmd += ' -s %s/*.pkl' % markers_dir
                cmd += ' -m %s/%s.fna' % (db_dir, strain)
                if strain in wol_pd['species']:
                    fna = '%s/%s.fna.bz2' % (
                        self.databases.paths['wol']['fna'],
                        wol_pd.loc[wol_pd['species'] == strain, 0])
                    io_update(self, i_f=fna)
                    cmd += ' -r %s' % fna
                cmd += ' -o %s' % odir
                cmd += ' -n %s' % self.soft.params['cpus']
                cmd += ' -c %s' % strain
                cmd += ' --mutation_rates'
                self.outputs['cmds'].append(cmd)
