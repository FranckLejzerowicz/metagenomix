# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import dirname, isfile


def metaphlan(out_dir: str, sam: str, inputs: dict, path: str, params: dict,
              strains: list, reads_fps: dict) -> tuple:
    """

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    path : str
        Path to the metaphlan database.
    params : dict
        Run parameters.
    strains : list
        Names of the strains to focus on.
    reads_fps : dict
        Paths to the read counts files.

    Returns
    -------
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.
    outputs : list
        All outputs paths.
    """
    io, cmds = {'I': [], 'O': []}, []
    tmpdir = '$TMPDIR/metaphlan_%s' % sam
    tmp_cmd = 'mkdir -p %s\n' % tmpdir

    bowtie2out = '%s/bowtie2/%s.bowtie2.bz2' % (out_dir, sam)
    outputs = list([bowtie2out])
    outputs.extend(list(profiling(
        out_dir, sam, inputs, path, params, bowtie2out, tmpdir, io, cmds)))

    reads = get_read_count(sam, reads_fps)
    if reads:
        analyses(out_dir, sam, params, strains, bowtie2out, reads, io, cmds)
    if cmds:
        cmds = [tmp_cmd] + cmds

    return io, cmds, outputs


def get_read_count(sam, reads_fps) -> str:
    """Get the number of reads in the current sample,
    or an empty string if the counting did not happen yet.

    Parameters
    ----------
    sam : str
        Sample name.
    reads_fps : dict
        Paths to the read counts files.

    Returns
    -------
    reads : str
        Number of reads.
    """
    reads = ''
    if sam in reads_fps and isfile(reads_fps[sam][0]):
        with open(reads_fps[sam][0]) as f:
            for line in f:
                reads = line[-1].strip().split(',')[-1]
                break
    return reads


def analyses(out_dir: str, sam: str, params: dict, strains: list,
             bowtie2out: str, reads: str, io: dict, cmds: list) -> None:
    """Get the taxonomic analyses Metaphlan command.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    params : dict
        Run parameters.
    strains : list
        Names of the strains to focus on.
    bowtie2out : str
        Input type was a bowtie2 output of the same tool.
    reads : str
        Number of reads in the current sample.
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.
    """
    tax_levs = ['a', 'k', 'p', 'c', 'o', 'f', 'g', 's']
    for analysis in ['rel_ab',
                     'rel_ab_w_read_stats',
                     'clade_profiles',
                     'marker_ab_table',
                     'marker_pres_table',
                     'clade_specific_strain_tracker']:
        dir_out = '%s/%s' % (out_dir, analysis)
        for tax_lev in tax_levs:
            rad = '%s/%s_t-%s' % (dir_out, sam, tax_lev)
            cmd = 'metaphlan'
            cmd += ' %s' % bowtie2out
            cmd += ' --input_type bowtie2out'
            cmd += ' -t %s' % analysis
            cmd += ' --nproc %s' % params['cpus']
            cmd += ' --sample_id_key sample_name'
            cmd += ' --sample_id %s' % sam
            cmd += ' --tax_lev %s' % tax_lev
            if reads and analysis == 'marker_ab_table':
                cmd += ' --nreads %s' % reads
            if analysis == 'clade_specific_strain_tracker':
                for strain in strains:
                    strain_cmd = cmd
                    strain_cmd += ' --clade %s' % strain.replace(' ', '_')
                    clade_out = '%s_%s.tsv' % (rad, strain.replace(' ', '_'))
                    strain_cmd += ' -o %s' % clade_out
                    if not isfile(clade_out):
                        io['O'].append(clade_out)
                        cmds.append(strain_cmd)
            else:
                ab_out = '%s.tsv' % rad
                cmd += ' -o %s' % ab_out
                if not isfile(ab_out):
                    io['O'].append(ab_out)
                    cmds.append(cmd)


def profiling(out_dir: str, sam: str, inputs: dict, path: str, params: dict,
              bowtie2out: str, tmpdir: str, io: dict, cmds: list) -> tuple:
    """Get the taxonomic profiling metaphlan command.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    path : str
        Path to the metaphlan database.
    params : dict
        Run parameters.
    bowtie2out : str
        Input type was a bowtie2 output of the same tool.
    tmpdir : str
        Path to a temporary directory.
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.

    Returns
    -------
    sam_out : str
        Alignment output.
    profile_out : str
        Metaphlan taxonomic profile output.
    """
    sam_out = '%s/sams/%s.sam.bz2' % (out_dir, sam)
    profile_out = '%s/profiles/%s.tsv' % (out_dir, sam)

    if isfile(profile_out):
        io['I'].append(bowtie2out)
    else:
        io['O'].append(profile_out)
        cmd = 'metaphlan'
        if isfile(bowtie2out):
            io['I'].append(bowtie2out)
            cmd += ' %s' % bowtie2out
            cmd += ' --input_type bowtie2out'
        else:
            io['I'].extend(inputs[sam])
            cmd += ' %s' % ','.join(inputs[sam])
            cmd += ' --input_type fastq'
            cmd += ' --samout %s' % sam_out
            cmd += ' --bowtie2out %s' % bowtie2out
            cmd += ' --bowtie2db %s' % path
        cmd += ' --tmp_dir %s' % tmpdir
        cmd += ' -o %s' % profile_out
        cmd += ' --nproc %s' % params['cpus']
        cmd += ' --sample_id_key %s' % sam
        cmds.append(cmd)

    return sam_out, profile_out


def get_profile(out_dir: str, sam: str, path: str, params: dict,
                profile: dict) -> tuple:
    """Get the taxonomic profiles if there are some, for humann to only run
    of the taxa from these profiles.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    path : str
        Path to the humann database.
    params : dict
        humann run paramters
    profile : dict
        Taxonomic profiles for species to focus humann on.

    Returns
    -------
    profiles : list
        Paths to the taxonomic profile of the current sample.
    dbs : list
        Paths to the metaphlan database.
    outs : list
        Paths to the output folder.
    """
    profiles, dbs, outs = [], [], []
    if profile:
        for profile_name, profile_fp in profile.items():
            # `profile_fp` must be "path/to/<sample>/profile_file.txt"
            profile_sams = {x.split('/')[-2]: x for x in glob.glob(profile_fp)}
            if sam in profile_sams:
                profiles.append(profile_sams[sam])
                dbs.append('')
                outs.append(
                    out_dir + '/profile_%s_%s_id%s_e%s_qC%s_sC%s_bypss%s' % (
                        profile_name, params['uniref'], params['id_thresh'],
                        params['evalue'], params['query_coverage'],
                        params['subject_coverage'], params['skip_translated']))
    else:
        profiles.append('')
        dbs.append(path)
        outs.append(out_dir + '/no_profile')
    return profiles, dbs, outs


def get_cmd(sam: str, inputs: dict, params: dict, profile: str, db: str,
            out: str, base_file: str, b2_ali_fp: str):
    """Get the humann command line.

    Parameters
    ----------
    sam : str
        Sample name.
    inputs : dict
        Input files.
    params : dict
        Run parameters.
    profile : str
        Path to the taxonomic profile of the current sample.
    db : str
        Path to the metaphlan database.
    out : str
        Path to the output folder.
    base_file : str
        Path to fastq file.
    b2_ali_fp : str
        Path to temporary, alignment sam file

    Returns
    -------
    cmd : str
        humann command line.
    """
    cmd = 'cat %s %s > %s\n' % (inputs[sam][0], inputs[sam][1], base_file)
    cmd += 'humann'
    cmd += ' --input %s' % base_file
    cmd += ' --output %s' % out
    cmd += ' --prescreen-threshold 0.001'
    cmd += ' --threads %s' % params['cpus']
    cmd += ' --search-mode %s' % params['uniref']
    if len(db):
        cmd += ' --metaphlan-options '
        cmd += '" --mpa_pkl %s/mpa_v20_m200.pkl --bowtie2db %s "' % (db, db)
    cmd += ' --identity-threshold %s' % params['id_thresh']
    cmd += ' --evalue %s' % params['evalue']
    cmd += ' --translated-query-coverage-threshold %s' % params['query_cov']
    cmd += ' --translated-subject-coverage-threshold %s' % params['subject_cov']
    if profile:
        cmd += ' --taxonomic-profile %s' % profile
    if params['skip_translated']:
        cmd += ' --bypass-translated-search'
    if isfile(b2_ali_fp):
        cmd += ' -r'
    cmd += ' --nucleotide-database %s' % params['nucleotide_db']
    cmd += ' --protein-database  %s\n' % params['protein_db']
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


def humann(out_dir: str, sam: str, inputs: dict, path: str,
           params: dict, profile: dict) -> tuple:
    """Create command lines for humann for the current database's species focus.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    sam : str
        Sample name.
    inputs : dict
        Input files.
    path : str
        Path to the humann database.
    params : dict
        Run parameters.
    profile : dict
        Taxonomic profiles for species to focus humann on.

    Returns
    -------
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.
    outputs : list
        All outputs paths.
    """
    io = {'I': [], 'O': {'f': [], 'd': []}}
    profiles, dbs, outs = get_profile(out_dir, sam, path, params, profile)
    if profiles:
        io['I'].extend(profiles)

    cmds = []
    outputs = []
    for prof, db, out in zip(*[profiles, dbs, outs]):

        b2_ali_fp = '%s/%s_humann_tmp/%s_bowtie2_aligned.sam' % (out, sam, sam)
        if isfile(b2_ali_fp):
            io['I'].append(b2_ali_fp)
        else:
            io['O']['f'].append(b2_ali_fp)
            io['O']['d'].append(out)

        base_file = '%s/%s.fastq.gz' % (out, sam)
        gene_fp = '%s/%s_genefamilies.tsv' % (out, sam)
        path_fp = '%s/%s_pathabundance.tsv' % (out, sam)
        cov_fp = '%s/%s_pathcoverage.tsv' % (out, sam)
        if not isfile(gene_fp) or not isfile(path_fp) or not isfile(cov_fp):
            cmds.append(get_cmd(
                sam, inputs, params, prof, db, out, base_file, b2_ali_fp))
            io['O']['f'].extend([gene_fp, path_fp, cov_fp])
        else:
            io['I'].extend([gene_fp, path_fp])

        gene_relab_fp = '%s/%s_genefamilies_relab.tsv' % (out, sam)
        path_relab_fp = '%s/%s_pathabundance_relab.tsv' % (out, sam)
        if not isfile(gene_relab_fp) or not isfile(path_relab_fp):
            io['O']['f'].extend([gene_relab_fp, path_relab_fp])
            cmds.append(renorm(gene_fp, gene_relab_fp))
            cmds.append(renorm(path_fp, path_relab_fp))

        outputs.append([out, b2_ali_fp, gene_fp, gene_relab_fp,
                        path_fp, path_relab_fp, cov_fp])

    return io, cmds, outputs


def strainphlan(out_dir: str, inputs: dict, params: dict,
                wol: dict, strains: dict) -> tuple:
    """Create command lines for strainphlan.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS.
    inputs : dict
        Input files.
    params : dict
        Run parameters.
    wol : dict
        Web Of Life database files
    strains : dict
        Species to focus on for strain analyses.

    Returns
    -------
    io : dict
        Inputs and outputs to potentially move to scratch and back.
    cmds : list
        All MIDAS command lines.
    outputs : list
        All outputs paths.
    dirs : list
        List of output folders to create.
    """
    io, cmds, outputs, dirs = {'I': {'f': [], 'd': []}, 'O': []}, [], {}, []
    bt2_fp, sam_fp = inputs[list(inputs.keys())[0]][:2]
    sam_dir = dirname(sam_fp)
    io['I']['d'].append(sam_dir)

    tmp_dir = '$TMPDIR/strainphlan'
    cmds.append('mkdir -p %s' % tmp_dir)

    db_dir = '%s/db_markers' % out_dir
    meta_dir = '%s/metadata' % out_dir
    markers_dir = '%s/consensus_markers' % out_dir
    dirs.extend([db_dir, meta_dir])
    io['O'].extend([markers_dir, db_dir, meta_dir])

    if len(glob.glob('%s/*.sam.bz2' % sam_dir)):
        cmd = 'sample2markers.py'
        cmd += ' -i %s/*.sam.bz2' % sam_dir
        cmd += ' -o %s' % markers_dir
        cmd += ' -n %s' % params['cpus']
        cmds.append(cmd)

    wol_genomes_pd = wol['taxonomy']
    for strain_group, strains in strains.items():
        strain_dir = '%s/%s' % (db_dir, strain_group)
        dirs.append(strain_dir)
        for st in [x.replace(' ', '_') for x in strains]:
            if st[0] != 's':
                continue
            if not isfile('%s/%s.fna' % (strain_dir, st)):
                cmd = 'extract_markers.py'
                cmd += ' -c %s' % st
                cmd += ' -o %s' % db_dir
                cmds.append(cmd)
            odir = '%s/output/%s' % (out_dir, strain_group)
            dirs.append(odir)
            io['O'].append(odir)
            outputs[(strain_group, st)] = odir
            if not isfile('%s/RAxML_bestTree.%s.StrainPhlAn3.tre' % (odir, st)):
                cmd = 'strainphlan'
                cmd += ' -s %s/*.pkl' % markers_dir
                cmd += ' -m %s/%s.fna' % (db_dir, st)
                if st in wol_genomes_pd['species']:
                    fna = '%s/%s.fna.bz2' % (wol['fna'], wol_genomes_pd.loc[
                        wol_genomes_pd['species'] == st, 0])
                    io['I']['f'].append(fna)
                    cmd += ' -r %s' % fna
                cmd += ' -o %s' % odir
                cmd += ' -n %s' % params['cpus']
                cmd += ' -c %s' % st
                cmd += ' --mutation_rates'
                cmds.append(cmd)

    return io, cmds, outputs, dirs
