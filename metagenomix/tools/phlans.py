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
              strains: list, reads_fps: dict, config) -> dict:
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
    config
        Configurations

    Returns
    -------
    outputs : dict
        All outputs
    """
    outputs = {'io': {'I': {'f': set()}, 'O': {'f': set()}},
               'cmds': [], 'dirs': [], 'outs': []}
    tmpdir = '$TMPDIR/metaphlan_%s' % sam
    tmp_cmd = ['mkdir -p %s\n' % tmpdir]

    b2o = '%s/bowtie2/%s.bowtie2.bz2' % (out_dir, sam)
    outs = profiling(out_dir, sam, inputs, path, params, b2o, tmpdir, outputs)
    outputs['outs'] = [b2o] + list(outs)
    outputs['dir'] = [dirname(x) for x in outputs['outs']]

    reads = get_read_count(sam, reads_fps)
    if reads:
        analyses(out_dir, sam, params, strains, b2o, reads, outputs, config)
    if outputs['cmds']:
        outputs['cmds'] = tmp_cmd + outputs['cmds']

    return outputs


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
             b2o: str, reads: str, outputs: dict, config) -> None:
    """Get the taxonomic analyses Metaphlan command.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS
    sam : str
        Sample name.
    params : dict
        Run parameters.
    strains : list
        Names of the strains to focus on.
    b2o : str
        Input type was a bowtie2 output of the same tool
    reads : str
        Number of reads in the current sample
    outputs : dict
        All outputs
    config
        Configurations
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
            cmd += ' %s' % b2o
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
                    if config.force or not isfile(clade_out):
                        outputs['io']['O']['f'].append(clade_out)
                        outputs['cmds'].append(strain_cmd)
            else:
                ab_out = '%s.tsv' % rad
                cmd += ' -o %s' % ab_out
                if config.force or not isfile(ab_out):
                    outputs['io']['O']['f'].append(ab_out)
                    outputs['cmds'].append(cmd)


def profiling(out_dir: str, sam: str, inputs: dict, path: str, params: dict,
              bowtie2out: str, tmpdir: str, outputs: dict) -> tuple:
    """Get the taxonomic profiling metaphlan command.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS
    sam : str
        Sample name
    inputs : dict
        Input files
    path : str
        Path to the metaphlan database
    params : dict
        Run parameters
    bowtie2out : str
        Input type was a bowtie2 output of the same tool
    tmpdir : str
        Path to a temporary directory
    outputs : dict
        All Outputs

    Returns
    -------
    sam_out : str
        Alignment output
    profile_out : str
        Metaphlan taxonomic profile output
    """
    sam_out = '%s/sams/%s.sam.bz2' % (out_dir, sam)
    profile_out = '%s/profiles/%s.tsv' % (out_dir, sam)

    if isfile(profile_out):
        outputs['io']['I']['f'].append(bowtie2out)
    else:
        outputs['io']['O']['f'].append(profile_out)
        cmd = 'metaphlan'
        if isfile(bowtie2out):
            outputs['io']['I']['f'].append(bowtie2out)
            cmd += ' %s' % bowtie2out
            cmd += ' --input_type bowtie2out'
        else:
            outputs['io']['I']['f'].extend(inputs[sam])
            cmd += ' %s' % ','.join(inputs[sam])
            cmd += ' --input_type fastq'
            cmd += ' --samout %s' % sam_out
            cmd += ' --bowtie2out %s' % bowtie2out
            cmd += ' --bowtie2db %s' % path
        cmd += ' --tmp_dir %s' % tmpdir
        cmd += ' -o %s' % profile_out
        cmd += ' --nproc %s' % params['cpus']
        cmd += ' --sample_id_key %s' % sam
        outputs['cmds'].append(cmd)

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
            out: str, base_file: str, ali: str):
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
    ali : str
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
    if isfile(ali):
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
           params: dict, profile: dict, conf) -> dict:
    """Create command lines for humann for the current database's species focus.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder for MIDAS
    sam : str
        Sample name
    inputs : dict
        Input files
    path : str
        Path to the humann database
    params : dict
        Run parameters
    profile : dict
        Taxonomic profiles for species to focus humann on
    conf
        Configurations

    Returns
    -------
    outputs : dict
        All outputs
    """
    outputs = {'io': {'I': {'f': set()}, 'O': {'d': set(), 'f': set()}},
               'cmds': [], 'dirs': [], 'outs': []}

    profiles, dbs, outs = get_profile(out_dir, sam, path, params, profile)
    if profiles:
        outputs['io']['I']['f'].update(profiles)

    for prof, db, out in zip(*[profiles, dbs, outs]):

        ali = '%s/%s_humann_tmp/%s_bowtie2_aligned.sam' % (out, sam, sam)
        if not conf.force and isfile(ali):
            outputs['io']['I']['f'].add(ali)
        else:
            outputs['io']['O']['f'].add(ali)
            outputs['io']['O']['d'].add(out)

        base_file = '%s/%s.fastq.gz' % (out, sam)
        gen = '%s/%s_genefamilies.tsv' % (out, sam)
        pwy = '%s/%s_pathabundance.tsv' % (out, sam)
        cov = '%s/%s_pathcoverage.tsv' % (out, sam)
        if conf.force or not isfile(gen) or not isfile(pwy) or not isfile(cov):
            outputs['cmds'].append(get_cmd(
                sam, inputs, params, prof, db, out, base_file, ali))
            outputs['io']['O']['f'].update([gen, pwy, cov])
        else:
            outputs['io']['I']['f'].update([gen, pwy])

        gen_relab = '%s/%s_genefamilies_relab.tsv' % (out, sam)
        pwy_relab = '%s/%s_pathabundance_relab.tsv' % (out, sam)
        if conf.force or not isfile(gen_relab) or not isfile(pwy_relab):
            outputs['io']['O']['f'].update([gen_relab, pwy_relab])
            outputs['cmds'].append(renorm(gen, gen_relab))
            outputs['cmds'].append(renorm(pwy, pwy_relab))

        outputs['dirs'].append(out)
        outputs['outs'].append([out, ali, gen, gen_relab, pwy, pwy_relab, cov])

    return outputs


def strainphlan(out_dir: str, inputs: dict, params: dict,
                wol: dict, strains: dict, config) -> dict:
    """Create command lines for strainphlan.

    Parameters
    ----------
    out_dir : str
        Path to pipeline output folder
    inputs : dict
        Input files
    params : dict
        Run parameters
    wol : dict
        Web Of Life database files
    strains : dict
        Species to focus on for strain analyses
    config
        Configurations

    Returns
    -------
    outputs : dict
        All outputs
    """
    outputs = {'io': {'I': {'f': set(), 'd': set()}, 'O': {'d': set()}},
               'cmds': [], 'dirs': [], 'outs': dict({})}
    bt2_fp, sam_fp = inputs[list(inputs.keys())[0]][:2]
    sam_dir = dirname(sam_fp)
    outputs['io']['I']['d'].add(sam_dir)

    tmp_dir = '$TMPDIR/strainphlan'
    outputs['cmds'].append('mkdir -p %s' % tmp_dir)

    db_dir = '%s/db_markers' % out_dir
    meta_dir = '%s/metadata' % out_dir
    markers_dir = '%s/consensus_markers' % out_dir
    outputs['dirs'].extend([db_dir, meta_dir])
    outputs['io']['O']['d'].update([markers_dir, db_dir, meta_dir])

    if config.force or len(glob.glob('%s/*.sam.bz2' % sam_dir)):
        cmd = 'sample2markers.py'
        cmd += ' -i %s/*.sam.bz2' % sam_dir
        cmd += ' -o %s' % markers_dir
        cmd += ' -n %s' % params['cpus']
        outputs['cmds'].append(cmd)

    wol_genomes_pd = wol['taxonomy']
    for strain_group, strains in strains.items():
        strain_dir = '%s/%s' % (db_dir, strain_group)
        outputs['dirs'].append(strain_dir)
        for st in [x.replace(' ', '_') for x in strains]:
            if st[0] != 's':
                continue
            if config.force or not isfile('%s/%s.fna' % (strain_dir, st)):
                cmd = 'extract_markers.py'
                cmd += ' -c %s' % st
                cmd += ' -o %s' % db_dir
                outputs['cmds'].append(cmd)
            odir = '%s/output/%s' % (out_dir, strain_group)
            outputs['dirs'].append(odir)
            outputs['io']['O']['d'].add(odir)
            outputs['outs'][(strain_group, st)] = odir
            tree = '%s/RAxML_bestTree.%s.StrainPhlAn3.tre' % (odir, st)
            if config.force or not isfile(tree):
                cmd = 'strainphlan'
                cmd += ' -s %s/*.pkl' % markers_dir
                cmd += ' -m %s/%s.fna' % (db_dir, st)
                if st in wol_genomes_pd['species']:
                    fna = '%s/%s.fna.bz2' % (wol['fna'], wol_genomes_pd.loc[
                        wol_genomes_pd['species'] == st, 0])
                    outputs['io']['I']['f'].add(fna)
                    cmd += ' -r %s' % fna
                cmd += ' -o %s' % odir
                cmd += ' -n %s' % params['cpus']
                cmd += ' -c %s' % st
                cmd += ' --mutation_rates'
                outputs['cmds'].append(cmd)

    return outputs
