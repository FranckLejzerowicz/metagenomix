# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

#
# def metamarker(soft, out_dir, inputs):
#
#     ald_soft_prev = soft_prev
#     soft_prev, grouping = stand_alone_groupings[soft]
#     outp_dir = outp_dir.replace('/%s/after_%s' % (soft, ald_soft_prev),
#                                 '/%s/after_%s' % (soft, soft_prev))
#     mkdr(outp_dir)
#
#     metamarker_cmds = []
#
#     meta_groups = dict \
#         ([x, meta_pd.loc[meta_pd[grouping ]= =x, :].index.tolist()] for x in meta_pd[grouping].unique())
#     for meta_group, samples in meta_groups.items():
#         sample_folder = '%s/%s' % (outp_dir, meta_group)
#         metamarker_cmds.append('mkdir -p %s' % sample_folder)
#
#     input_files = {}
#     all_input_files = []
#     # fdx = 0
#     for root, dirs, files in os.walk('%s/pipeline/%s' % (rad, soft_prev)):
#         for fil in files:
#             if 'bin' in fil and 'metawrap_' in root.split('/')
#                 [-1] and '_bins' in root.split('/')[-1]:
#                 # fdx += 1
#                 # if fdx > 5:
#                 #     continue
#                 path = root + '/' + fil
#                 sam = root.split(soft_prev)[-1].split('/')[3]
#                 input_files.setdefault(sam, []).append(path)
#                 all_input_files.append(path)
#
#     npzs = []
#     nSamples = []
#     sample_folders = []
#     metamarker_outs = []
#     for meta_group, samples in meta_groups.items():
#         metamarker_cmds.append('\n')
#         sample_folder = '%s/%s' % (outp_dir, meta_group)
#         cur_nSamples = 0
#         meta_group_fas = []
#         for sdx, sam in enumerate(samples):
#             if sam not in input_files:
#                 continue
#             for idx, input_file in enumerate(input_files[sam]):
#                 cur_fas = '%s/%s.%s.fasta' % \
#                 (sample_folder, sam, splitext(basename(input_file))[0])
#                 cur_cmd = '/bin/cp %s %s' % (input_file, cur_fas)
#                 cur_nSamples += 1
#                 metamarker_cmds.append(cur_cmd)
#                 meta_group_fas.append(cur_fas)
#         nSamples.append(cur_nSamples)
#
#         sample_folders.append(sample_folder)
#         out_npz = '%s/outStep1/%s/%s.npz' % (outp_dir, meta_group, meta_group)
#         out_folder = '%s/outStep1/%s' % (outp_dir, meta_group)
#         cur_cmd = 'Step1_MakeHashTable.py'
#         cur_cmd += ' --sample %s' % sample_folder
#         cur_cmd += ' --output %s' % out_folder
#         cur_cmd += ' --npz %s' % out_npz
#         mkdir_cmd = 'mkdir -p %s' % out_folder
#         metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#         metamarker_cmds.append(mkdir_cmd)
#         metamarker_cmds.append('\necho "%s"' % cur_cmd)
#         metamarker_cmds.append(cur_cmd)
#         npzs.append(out_npz)
#
#     out_folder = '%s/out_Step2' % outp_dir
#     cur_cmd = 'Step2_MakeShortMarker.py'
#     cur_cmd += ' --caseHash %s' % npzs[0]
#     cur_cmd += ' --controlHash %s' % npzs[1]
#     cur_cmd += ' --spade /home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin'
#     cur_cmd += ' --output %s' % out_folder
#     cur_cmd += ' --threads %s' % procs
#     metamarker_cmds.append('\n')
#     mkdir_cmd = 'mkdir -p %s' % out_folder
#     metamarker_cmds.append('\n')
#     metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#     metamarker_cmds.append(mkdir_cmd)
#     metamarker_cmds.append('\necho "%s"' % cur_cmd)
#     metamarker_cmds.append(cur_cmd)
#
#     step2_outs = ['%s/Short_Markers_case.fasta' % out_folder,
#                   '%s/Short_Markers_control.fasta' % out_folder]
#
#     step3_outs = []
#     metamarker_cmds.append('\n')
#     for mdx, meta_group in enumerate(meta_groups.keys()):
#         cur_out = '%s/out_Step3/%s' % (outp_dir, meta_group)
#         step3_outs.append(cur_out)
#         tmpdir = '$TMPDIR/MetaMarker/%s' % meta_group
#         cur_cmd = 'Step3_MakeLongMarker.py'
#         cur_cmd += ' --short_marker %s' % step2_outs[mdx]
#         cur_cmd += ' --sample %s' % sample_folders[mdx]
#         cur_cmd += ' --tmp %s' % tmpdir
#         cur_cmd += ' --output %s' % cur_out
#         cur_cmd += ' --spade /home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin'
#         cur_cmd += ' --usearch /home/flejzerowicz/softs'
#         cur_cmd += ' --gt /home/flejzerowicz/usr/local/genometools/bin'
#         cur_cmd += ' --threads %s' % procs
#         cur_cmd += '  --id 0.9'
#         mkdir_cmd = 'mkdir -p %s' % cur_out
#         metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#         metamarker_cmds.append(mkdir_cmd)
#         metamarker_cmds.append('\necho "%s"' % cur_cmd)
#         metamarker_cmds.append(cur_cmd)
#
#     cur_out = '%s/out_Step4' % outp_dir
#     cur_cmd = 'Step4_MarkerCleaning.py'
#     cur_cmd += ' --caseMarker %s' % step3_outs[0]
#     cur_cmd += ' --caseSize %s' % nSamples[0]
#     cur_cmd += ' --controlMarker %s' % step3_outs[1]
#     cur_cmd += ' --controlSize %s' % nSamples[1]
#     cur_cmd += ' --output %s' % cur_out
#     cur_cmd += ' --cdhit /home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208'
#     mkdir_cmd = 'mkdir -p %s' % cur_out
#     metamarker_cmds.append('\n')
#     metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#     metamarker_cmds.append(mkdir_cmd)
#     metamarker_cmds.append('\necho "%s"' % cur_cmd)
#     metamarker_cmds.append(cur_cmd)
#
#     step5_outs = []
#     metamarker_cmds.append('\n')
#     for mdx, meta_group in enumerate(meta_groups.keys()):
#         cur_out = '%s/out_Step5/%s' % (outp_dir, meta_group)
#         tmpdir = '$TMPDIR/MetaMarker/Step5/%s' % meta_group
#         step5_outs.append(cur_out)
#         cur_cmd = 'Step5_MarkerAbundace.py'
#         cur_cmd += ' --marker %s/out_Step4/Selected_Markers.fasta' % outp_dir
#         cur_cmd += ' --sample %s' % sample_folders[mdx]
#         cur_cmd += ' --output %s' % cur_out
#         cur_cmd += ' --tmp %s' % tmpdir
#         cur_cmd += ' --alnlength 50'
#         cur_cmd += ' --id 0.9'
#         cur_cmd += ' --bowtie2 /home/flejzerowicz/usr/miniconda3/bin'
#         cur_cmd += ' --threads %s' % procs
#         cur_cmd += ' --removeTemp'
#         mkdir_cmd = 'mkdir -p %s' % cur_out
#         metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#         metamarker_cmds.append(mkdir_cmd)
#         metamarker_cmds.append('\necho "%s"' % cur_cmd)
#         metamarker_cmds.append(cur_cmd)
#
#     cur_out = '%s/out_Step6' % outp_dir
#     cur_cmd = 'Step6_MarkerRank.py'
#     cur_cmd += ' --case %s' % step5_outs[0]
#     cur_cmd += ' --control %s' % step5_outs[1]
#     cur_cmd += ' --marker %s/out_Step4/Selected_Markers.fasta' % outp_dir
#     cur_cmd += ' --cdhit /home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208'
#     cur_cmd += ' --output %s' % cur_out
#     mkdir_cmd = 'mkdir -p %s' % cur_out
#     metamarker_cmds.append('\n')
#     metamarker_cmds.append('\necho "%s"' % mkdir_cmd)
#     metamarker_cmds.append(mkdir_cmd)
#     metamarker_cmds.append('\necho "%s"' % cur_cmd)
#     metamarker_cmds.append(cur_cmd)
#
#     step6_outs = ['%s/%s' % (cur_out, x) for x in ['Final_Marker_Case.csv',
#                                                    'Final_Marker_Case.fasta',
#                                                    'Final_Marker_Control.csv',
#                                                    'Final_Marker_Control.fasta',
#                                                    'Marker_Stats.csv']]
#
#     IO['I']['f'].extend(all_input_files)
#     IO['O']['f'].extend(step6_outs)
#     metamarker_outs.extend(step6_outs)
#
#     return metamarker_cmds, metamarker_outs, procs, mem_n_u, IO

