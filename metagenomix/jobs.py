# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import subprocess
import numpy as np
from os.path import dirname, splitext

from metagenomix._io_utils import mkdr, get_roundtrip


class CreateScripts(object):

    def __init__(self, config):
        self.config = config
        self.cmd = []
        self.cmds = {}
        self.chunks = {}
        self.soft = ''
        self.sh = ''
        self.run = {'database': {}, 'software': {}}
        self.job_fps = []
        self.job_name = ''
        self.modules = {}
        self.pjct = self.get_prjct()
        self.scheduler = self.get_scheduler()

    def get_scheduler(self):
        if self.config.jobs:
            if self.config.torque:
                return 'qsub'
            else:
                return 'sbatch'
        else:
            return 'sh'

    def get_prjct(self):
        prjct = [x for x in self.config.project if x.lower() not in 'aeiouy']
        if prjct:
            return ''.join(prjct)
        else:
            return self.config.project

    def get_chunks(self, n_chunks):
        self.chunks = {}
        if n_chunks and len(self.cmds) >= n_chunks:
            cmds_d = dict(enumerate(self.cmds))
            chunks = [list(x) for x in np.array_split(list(cmds_d), n_chunks)]
            for cdx, chunk in enumerate(chunks):
                if len(chunk):
                    self.chunks[str(cdx)] = [cmds_d[x] for x in chunk]
        else:
            for key in self.cmds:
                if isinstance(key, str):
                    self.chunks[key] = [key]
                else:
                    self.chunks['_'.join([k for k in key if k])] = [key]

    def get_main_sh(self, name, soft=None) -> str:
        main = '%s/run_%s' % (self.sh.rsplit('/', 2)[0], name)
        if soft:
            main += '_after_%s.sh' % soft.prev
            self.run['software'][name] = main
        else:
            main += '.sh'
            self.run['database'][name] = main
        return main

    def write_main(self, name, soft=None) -> None:
        main = self.get_main_sh(name, soft)
        with open(main, 'w') as o:
            if len(self.job_fps):
                job_dir = dirname(self.job_fps[0])
                o.write('mkdir -p %s/output\n' % job_dir)
                o.write('cd %s/output\n' % job_dir)
                for job_fp in self.job_fps:
                    o.write('%s %s\n' % (self.scheduler, job_fp))
                self.job_fps = []

    def database_cmds(self, databases):
        for db, cmds in databases.commands.items():
            self.cmds = cmds
            self.get_chunks(self.config.params['chunks'])
            self.write_jobs(db)
            self.write_main(db)

    def get_modules(self, name: str):
        self.modules = self.config.modules.get(name, [])

    def scratch(self, soft, key, cmds):
        if soft.params['scratch'] and self.config.jobs and key in soft.io:
            roundtrip = get_roundtrip(soft.io[key])
            scratch_cmds = ['\n# Move to SCRATCH_FOLDER']
            scratch_cmds += roundtrip['to']
            scratch_cmds += ['\n# %s commands (%s)' % (soft.name, key)]
            scratch_cmds += cmds
            if self.config.move_back:
                scratch_cmds += ['\n# Move from SCRATCH_FOLDER']
                scratch_cmds += roundtrip['from']
            self.cmds[key] = scratch_cmds
        else:
            self.cmds[key] = cmds

    def get_cmds(self, soft, commands):
        self.cmds = {}
        for sam_or_pool, cmds in soft.cmds.items():
            if isinstance(cmds, list):
                self.scratch(soft, sam_or_pool, cmds)
            elif isinstance(cmds, dict):
                if sam_or_pool in commands.pools:
                    # for group in commands.pools[sam_or_pool]:
                    for group, group_cmds in cmds.items():
                        self.scratch(soft, (sam_or_pool, group), group_cmds)
                else:
                    self.scratch(soft, sam_or_pool, cmds)
            else:
                sys.exit('The collected commands are neither list of dict!')

    @staticmethod
    def print_status(soft):
        if soft.status == {'Done'}:
            print('-> Done')
        elif len(soft.status):
            print()
            for stat in soft.status:
                if stat != 'Done':
                    print('\t-> %s' % stat)
        else:
            print()

    def software_cmds(self, commands):
        for sdx, (name, soft) in enumerate(commands.softs.items()):
            print('[Writing commands] #%s: %s' % (sdx, name), end=' ')
            self.print_status(soft)
            if not len(soft.cmds):
                continue
            self.get_modules(name)
            self.get_cmds(soft, commands)
            self.get_chunks(soft.params['chunks'])
            self.write_jobs(name, soft)
            self.write_main(name, soft)

    def get_sh(self, name: str, chunk_name: str, soft=None) -> None:
        """

        Parameters
        ----------
        name : str
            Name of the current software of the pipeline workflow
        chunk_name : str
            Name of the current chunk of commands
        soft
        """
        if soft:
            self.sh = '%s/%s/jobs/run_%s_after_%s_%s.sh' % (
                self.config.dir, name, name, soft.prev, chunk_name)
        else:
            self.sh = '%s/databases/jobs/build_%s_%s.sh' % (
                self.config.dir, name, chunk_name)
        mkdr(dirname(self.sh))

    def prep_script(self, params: dict) -> None:
        # mandatory options
        self.cmd = [
            'Xhpc',
            '-j', self.job_name,
            '-t', str(params['time']),
            '-c', str(params['cpus']),
            '-M', str(params['mem_num']), params['mem_dim'],
            '--no-stat',
            '-i', self.sh]
        # whether the cpus requests is per node
        if params['nodes']:
            self.cmd.extend(['-n', str(params['nodes'])])
        # always provide an account
        if self.config.account:
            self.cmd.extend(['-a', self.config.account])
        # get the job script file path and use it
        job_script = '%s.slm' % splitext(self.sh)[0]
        if self.config.torque:
            self.cmd.append('--torque')
            job_script = '%s.pbs' % splitext(self.sh)[0]
        self.job_fps.append(job_script)
        self.cmd.extend(['-o', job_script])
        # whether an environment must be used for the current software
        if not self.modules and params['env']:
            self.cmd.extend(['-e', params['env']])
        # setup the scratch location to be used for the current software
        if isinstance(params['scratch'], int):
            self.cmd.extend(['--localscratch', str(params['scratch'])])
        elif params['scratch'] == 'scratch':
            self.cmd.append('--scratch')
        elif params['scratch'] == 'userscratch':
            self.cmd.append('--userscratch')
        if params['machine']:
            self.cmd.append('--%s' % params['machine'])
        # quiet Xhpc if metagenomix is not supposed to be verbose
        if not self.config.verbose:
            self.cmd.append('--quiet')
        # print(' '.join(self.cmd))

    def call_cmd(self):
        cmd = ' '.join(self.cmd)
        if self.config.verbose:
            print('[Running]', cmd)
        subprocess.call(cmd.split())
        os.remove(self.sh)

    def write_chunks(self, chunk_keys: list, soft):
        with open(self.sh, 'w') as sh:
            if self.modules:
                sh.write('module purge\n')
            for module in self.modules:
                sh.write('module load %s\n' % module)
            for chunk_key in chunk_keys:
                for cmd in self.cmds[chunk_key]:
                    if soft.params['scratch'] and self.config.jobs:
                        sh.write('%s\n' % cmd)
                    else:
                        sh.write('%s\n' % cmd.replace('${SCRATCH_FOLDER}', ''))

    def get_job_name(self, name: str, chunk_name: str):
        self.job_name = name + '.' + self.pjct + '.' + chunk_name

    def write_script(self, soft=None):
        if self.config.jobs:
            if soft:
                self.prep_script(soft.params)
            else:
                self.prep_script(self.config.params)
            self.call_cmd()
        else:
            self.job_fps.append('%s.sh' % splitext(self.sh)[0])

    def write_jobs(self, name: str, soft=None):
        for chunk_name, chunk_keys in self.chunks.items():
            self.get_sh(name, chunk_name, soft)
            self.write_chunks(chunk_keys, soft)
            self.get_job_name(name, chunk_name)
            self.write_script(soft)

    def display(self):
        if len(self.run['database']) or len(self.run['software']):
            print()
            print('< PLEASE CONSIDER CHECKING THE SCRIPTS MANUALLY >')
        for database_software, name_main in self.run.items():
            print()
            print('# ========== #')
            print('  ', database_software)
            print('# ========== #')
            for name, main in name_main.items():
                print()
                print('>', name)
                print('sh', main)
