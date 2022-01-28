# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def prep_midas():

    procs = 24
    mem_n_u = ('50', 'gb')

    midas_cmds = []
    midas_outs = []

    done = []

    self.tool.io['I']['f'].extend(input_paths)
    for midas_case, (midas_db, midas_species_list) in midas_focus.items():

        if midas_case:
            sam_outp_dir = '%s/%s/%s' % (outp_dir, midas_case, sam)
        else:
            sam_outp_dir = '%s/%s' % (outp_dir, sam)
        if midas_db:
            cur_db = midas_db
        else:
            cur_db = '/home/flejzerowicz/databases/midas_db_v1.2'

            cmd = 'run_midas.py species'
            cmd += ' %s' % sam_outp_dir
            cmd += ' -1 %s' % input_paths[0]
            cmd += ' -2 %s' % input_paths[1]
            cmd += ' -d %s' % cur_db
            cmd += ' -t %s' % procs
            cmd += ' --remove_temp'
            species_out = '%s/species' % sam_outp_dir
            # species_out = '%s/%s/species' % (sam_outp_dir, sam)
            if isfile('%s/species_profile.txt' % species_out):
                done.append('species')
            else:
                midas_cmds.append(cmd)
                IO['O']['d'].append(species_out)
            midas_outs.append(species_out)

        midas_species_select = set()
        if midas_species_list:
            midas_species = [x.strip() for x in
                             open(midas_species_list).readlines()]
            with open('%s/species_info.txt' % cur_db) as f:
                for line in f:
                    for genus_species in midas_species:
                        if genus_species.replace(' ', '_') in line:
                            midas_species_select.add(line.split('\t')[0])

        cmd = 'run_midas.py genes'
        # if 'species' in done:
        #     cmd += ' %s/%s' % (sam_outp_dir, sam)
        # else:
        cmd += ' %s' % sam_outp_dir
        cmd += ' -1 %s' % input_paths[0]
        cmd += ' -2 %s' % input_paths[1]
        cmd += ' -d %s' % cur_db
        cmd += ' -t %s' % procs
        cmd += ' --remove_temp'
        cmd += ' --species_cov 1'
        if midas_species_select:
            cmd += ' --species_id %s' % ','.join(list(midas_species_select))
        genes_out = '%s/genes' % sam_outp_dir
        # genes_out = '%s/%s/genes' % (sam_outp_dir, sam)
        if isfile('%s/readme.txt' % genes_out):
            done.append('genes')
        else:
            midas_cmds.append(cmd)
            # cmd = 'mv %s/genes/* %s/.' % (genes_out, genes_out)
            # midas_cmds.append(cmd)
            # cmd = 'rm -rf %s/genes' % genes_out
            # midas_cmds.append(cmd)
            if 'species' in done:
                IO['I']['d'].append(species_out)
            IO['O']['d'].append(genes_out)
        midas_outs.append(genes_out)

        cmd = 'run_midas.py snps'
        # if 'species' in done:
        #     cmd += ' %s/%s' % (sam_outp_dir, sam)
        # else:
        cmd += ' %s' % sam_outp_dir
        cmd += ' -1 %s' % input_paths[0]
        cmd += ' -2 %s' % input_paths[1]
        cmd += ' -d %s' % cur_db
        cmd += ' -t %s' % procs
        cmd += ' --remove_temp'
        cmd += ' --species_cov 1'
        if midas_species_select:
            cmd += ' --species_id %s' % ','.join(list(midas_species_select))
        snps_out = '%s/snps' % sam_outp_dir
        # snps_out = '%s/%s/snps' % (sam_outp_dir, sam)
        if isfile('%s/readme.txt' % snps_out):
            done.append('snps')
        else:
            midas_cmds.append(cmd)
            # cmd = 'mv %s/snps/* %s/.' % (snps_out, snps_out)
            # midas_cmds.append(cmd)
            # cmd = 'rm -rf %s/snps' % snps_out
            # midas_cmds.append(cmd)
            if 'species' in done:
                IO['I']['d'].append(species_out)
            IO['O']['d'].append(snps_out)
        midas_outs.append(snps_out)
    return midas_cmds, midas_outs, procs, mem_n_u, IO
