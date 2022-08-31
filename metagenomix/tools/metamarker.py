# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import basename, splitext

from metagenomix._io_utils import io_update


def get_meta_groupings(self) -> dict:
    """

    Returns
    -------
    meta_groups : dict
        Samples per group's unique value per group variable
    """
    meta_groups = {}
    meta = self.config.meta.copy()
    for group in self.soft.params['groups']:  # for each group from user params
        meta_groups[group] = {}
        group_type = str(meta[group].dtype)
        # can not use group/variables that are numerical
        if group_type != 'object':
            continue
        for factor in set(meta[group]):
            sams = meta.loc[meta[group] == factor, :].index.tolist()
            meta_groups[group][factor] = sams
        sample_folder = '%s/%s' % (self.dir, group)
        self.outputs['cmds'].append('mkdir -p %s' % sample_folder)
    return meta_groups


def get_input_files(self):
    input_files = {}
    all_input_files = []
    # fdx = 0

    #
    # GET THE RIGHT OUTPUT
    #
    for root, dirs, files in os.walk('%s/%s' % (self.dir, self.soft.prev)):
        for fil in files:
            if 'bin' in fil and 'metawrap_' in root.split('/')[-1] and '_bins' in root.split('/')[-1]:
                # fdx += 1
                # if fdx > 5:
                #     continue
                path = root + '/' + fil
                sam = root.split(self.soft.prev)[-1].split('/')[3]
                input_files.setdefault(sam, []).append(path)
                all_input_files.append(path)
                io_update(self, i_f=all_input_files)
    return input_files


def metamarker(self) -> None:
    """Prepare the command and collect outputs, io and dirs for Metamarker.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .outputs : dict
            All outputs
        .soft
            Software class instances and attributes
        .config
            Configuration class instance
    """
    # note that metamarker can run on genomes data after e.g.:
    #   - metawrap_ref
    #   - other
    #  which should be contained in some user_params

    meta_groups = get_meta_groupings(self)
    input_files = get_input_files(self)
    data = metamarker_step1(self, meta_groups, input_files)
    metamarker_step2(self, data)
    metamarker_step3(self, data, meta_groups)
    metamarker_step4(self, data)
    metamarker_step5(self, data, meta_groups)
    metamarker_step6(self, data)


def metamarker_step1(self, meta_groups: dict, input_files: dict) -> dict:

    data = {'npzs': [], 'nsamples': [], 'sample_folders': []}
    for meta_group, samples in meta_groups.items():
        self.outputs['cmds'].append('\n')
        sample_folder = '%s/%s' % (self.dir, meta_group)
        n_samples = 0
        meta_group_fas = []
        for sdx, sam in enumerate(samples):
            if sam not in input_files:
                continue
            for idx, input_file in enumerate(input_files[sam]):
                cur_fas = '%s/%s.%s.fasta' % (
                    sample_folder, sam, splitext(basename(input_file))[0])
                cur_cmd = '/bin/cp %s %s' % (input_file, cur_fas)
                n_samples += 1
                self.outputs['cmds'].append(cur_cmd)
                meta_group_fas.append(cur_fas)

        out_npz = '%s/outStep1/%s/%s.npz' % (self.dir, meta_group, meta_group)
        out_folder = '%s/outStep1/%s' % (self.dir, meta_group)
        cur_cmd = 'Step1_MakeHashTable.py'
        cur_cmd += ' --sample %s' % sample_folder
        cur_cmd += ' --output %s' % out_folder
        cur_cmd += ' --npz %s' % out_npz
        mkdir_cmd = 'mkdir -p %s' % out_folder

        cmd = '\n\necho "%s"\n' % mkdir_cmd
        cmd += mkdir_cmd
        cmd += '\necho "%s"\n' % cur_cmd
        cmd += cur_cmd
        self.outputs['cmds'].append(cmd)

        data['npzs'].append(out_npz)
        data['nsamples'].append(n_samples)
        data['sample_folders'].append(sample_folder)

    return data


def metamarker_step2(self, data: dict):
    out_folder = '%s/out_Step2' % self.dir
    cur_cmd = 'Step2_MakeShortMarker.py'
    cur_cmd += ' --caseHash %s' % data['npzs'][0]
    cur_cmd += ' --controlHash %s' % data['npzs'][1]
    cur_cmd += ' --spade /home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin'
    cur_cmd += ' --output %s' % out_folder
    cur_cmd += ' --threads %s' % self.soft.params['cpus']
    mkdir_cmd = 'mkdir -p %s' % out_folder
    cmd = '\n\necho "%s"\n' % mkdir_cmd
    cmd += mkdir_cmd
    cmd += '\necho "%s"\n' % cur_cmd
    cmd += cur_cmd
    self.outputs['cmds'].append(cmd)
    data['step2_outs'] = ['%s/Short_Markers_case.fasta' % out_folder,
                          '%s/Short_Markers_control.fasta' % out_folder]


def metamarker_step3(self, data: dict, meta_groups: dict):
    data['step3_outs'] = []
    self.outputs['cmds'].append('\n')
    for mdx, meta_group in enumerate(meta_groups.keys()):
        cur_out = '%s/out_Step3/%s' % (self.dir, meta_group)
        data['step3_outs'].append(cur_out)
        tmpdir = '$TMPDIR/MetaMarker/%s' % meta_group
        cur_cmd = 'Step3_MakeLongMarker.py'
        cur_cmd += ' --short_marker %s' % data['step2_outs'][mdx]
        cur_cmd += ' --sample %s' % data['sample_folders'][mdx]
        cur_cmd += ' --tmp %s' % tmpdir
        cur_cmd += ' --output %s' % cur_out
        cur_cmd += ' --spade /home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin'
        cur_cmd += ' --usearch /home/flejzerowicz/softs'
        cur_cmd += ' --gt /home/flejzerowicz/usr/local/genometools/bin'
        cur_cmd += ' --threads %s' % self.soft.params['cpus']
        cur_cmd += '  --id %s' % self.soft.params['identity']
        mkdir_cmd = 'mkdir -p %s' % cur_out
        cmd = '\n\necho "%s"\n' % mkdir_cmd
        cmd += mkdir_cmd
        cmd += '\necho "%s"\n' % cur_cmd
        cmd += cur_cmd
        self.outputs['cmds'].append(cmd)


def metamarker_step4(self, data: dict):
    cur_out = '%s/out_Step4' % self.dir
    cur_cmd = 'Step4_MarkerCleaning.py'
    cur_cmd += ' --caseMarker %s' % data['step3_outs'][0]
    cur_cmd += ' --caseSize %s' % data['nsamples'][0]
    cur_cmd += ' --controlMarker %s' % data['step3_outs'][1]
    cur_cmd += ' --controlSize %s' % data['nsamples'][1]
    cur_cmd += ' --output %s' % cur_out
    cur_cmd += ' --cdhit /home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208'
    mkdir_cmd = 'mkdir -p %s' % cur_out
    cmd = '\n\necho "%s"\n' % mkdir_cmd
    cmd += mkdir_cmd
    cmd += '\necho "%s"\n' % cur_cmd
    cmd += cur_cmd
    self.outputs['cmds'].append(cmd)


def metamarker_step5(self, data: dict, meta_groups: dict):
    data['step5_outs'] = []
    self.outputs['cmds'].append('\n')
    for mdx, meta_group in enumerate(meta_groups.keys()):
        cur_out = '%s/out_Step5/%s' % (self.dir, meta_group)
        tmpdir = '$TMPDIR/MetaMarker/Step5/%s' % meta_group
        data['step5_outs'].append(cur_out)
        cur_cmd = 'Step5_MarkerAbundace.py'
        cur_cmd += ' --marker %s/out_Step4/Selected_Markers.fasta' % self.dir
        cur_cmd += ' --sample %s' % data['sample_folders'][mdx]
        cur_cmd += ' --output %s' % cur_out
        cur_cmd += ' --tmp %s' % tmpdir
        cur_cmd += ' --alnlength 50'
        cur_cmd += ' --id %s' % self.soft.params['identity']
        cur_cmd += ' --bowtie2 /home/flejzerowicz/usr/miniconda3/bin'
        cur_cmd += ' --threads %s' % self.soft.params['cpus']
        cur_cmd += ' --removeTemp'
        mkdir_cmd = 'mkdir -p %s' % cur_out
        cmd = '\n\necho "%s"\n' % mkdir_cmd
        cmd += mkdir_cmd
        cmd += '\necho "%s"\n' % cur_cmd
        cmd += cur_cmd
        self.outputs['cmds'].append(cmd)


def metamarker_step6(self, data: dict):
    cur_out = '%s/out_Step6' % self.dir
    cur_cmd = 'Step6_MarkerRank.py'
    cur_cmd += ' --case %s' % data['step5_outs'][0]
    cur_cmd += ' --control %s' % data['step5_outs'][1]
    cur_cmd += ' --marker %s/out_Step4/Selected_Markers.fasta' % self.dir
    cur_cmd += ' --cdhit /home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208'
    cur_cmd += ' --output %s' % cur_out
    mkdir_cmd = 'mkdir -p %s' % cur_out
    cmd = '\n\necho "%s"\n' % mkdir_cmd
    cmd += mkdir_cmd
    cmd += '\necho "%s"\n' % cur_cmd
    cmd += cur_cmd
    self.outputs['cmds'].append(cmd)
    step6_outs = ['%s/%s' % (cur_out, x) for x in [
        'Final_Marker_Case.csv', 'Final_Marker_Case.fasta',
        'Final_Marker_Control.csv', 'Final_Marker_Control.fasta',
        'Marker_Stats.csv']]
    io_update(self, o_f=step6_outs)
    self.outputs['outs'].extend(step6_outs)
