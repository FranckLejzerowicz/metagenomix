# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess
import numpy as np
from os.path import dirname, splitext

from metagenomix._io_utils import mkdr


class CreateScripts(object):

    def __init__(self, config):
        self.config = config
        self.cmd = []
        self.cmds = {}
        self.cmds_chunks = []
        self.soft = ''
        self.sh = ''
        self.run = []
        self.main_fps = []
        self.job_fps = []
        self.job_name = ''
        self.module = ''
        self.pjct = self.get_prjct()
        self.scheduler = self.get_scheduler()

    def get_scheduler(self):
        if self.config.jobs:
            if self.config.torque:
                return 'qsub'
            else:
                return 'sbatch'
        return ''

    def get_prjct(self):
        prjct = [x for x in self.config.project if x.lower() not in 'aeiouy']
        if prjct:
            return ''.join(prjct)
        else:
            return self.config.project

    def get_cmds_chunks(self, n_chunks):
        if len(self.cmds) > int(n_chunks):
            cmds_split = np.array_split(list(self.cmds), n_chunks)
            self.cmds_chunks = [list(x) for x in cmds_split if len(x)]
        else:
            self.cmds_chunks = [[x] for x in self.cmds]

    def write_main(self, name, soft=None) -> None:
        main = '%s/run_%s' % (self.sh.rsplit('/', 2)[0], name)
        if soft:
            main += '_after_%s' % soft.prev
        main += '.sh'
        with open(main, 'w') as o:
            if self.scheduler:
                for job_fp in self.job_fps:
                    o.write('%s %s\n' % (self.scheduler, job_fp))
                self.job_fps = []
        self.run[name] = main

    def database_cmds(self, databases):
        for db, cmds in databases.commands.items():
            self.cmds = cmds
            self.get_cmds_chunks(self.config.params['chunks'])
            self.write_jobs(db)
            self.write_main(db)

    def get_module(self, name: str):
        self.module = self.config.modules.get(name)

    def software_cmds(self, commands):
        for name, soft in commands.softs.items():
            self.get_module(name)
            self.cmds = soft.cmds
            self.get_cmds_chunks(soft.params['chunks'])
            self.write_jobs(name, soft)
            self.write_main(name, soft)

    def get_sh(self, name: str, cdx: int, soft=None) -> None:
        if soft:
            self.sh = '%s/%s/jobs/run_%s_after_%s_%s.sh' % (
                self.config.dir, name, name, soft.prev, cdx)
        else:
            self.sh = '%s/databases/jobs/build_%s_%s.sh' % (
                self.config.dir, name, cdx)
        mkdr(dirname(self.sh))

    def prep_script(self, params: dict) -> None:
        self.cmd = [
            'Xhpc',
            '-j', self.job_name,
            '-t', params['time'],
            '-c', str(params['cpus']),
            '-n', str(params['nodes']),
            '-M', str(params['mem_num']), params['mem_dim'],
            '--no-stat',
            '-i', self.sh]

        job_script = '%s.slm' % splitext(self.sh)[0]
        if self.config.torque:
            self.cmd.append('--torque')
            job_script = '%s.pbs' % splitext(self.sh)[0]
        self.job_fps.append(job_script)
        self.cmd.extend(['-o', job_script])

        if not self.module and params['env']:
            self.cmd.extend(['-e', params['env']])

        if self.config.userscratch:
            self.cmd.append('--userscratch')
        if self.config.scratch:
            self.cmd.append('--scratch')
        if self.config.localscratch:
            self.cmd.extend(['--localscratch', self.config.localscratch])

    def call_cmd(self):
        cmd = ' '.join(map(str, self.cmd))
        if self.config.verbose:
            print('[Running]', cmd)
        subprocess.call(cmd.split())

    def write_chunks(self, chunks: list):
        with open(self.sh, 'w') as sh:
            if self.module:
                sh.write('module load %s\n' % self.module)
            for chunk in chunks:
                for cmd in self.cmds[chunk]:
                    sh.write('%s\n' % cmd)

    def get_job_name(self, name: str, adx: int):
        self.job_name = name + '.' + self.pjct + '.' + str(adx)

    def write_script(self, soft=None):
        if self.config.jobs:
            if soft:
                self.prep_script(soft.params)
            else:
                self.prep_script(self.config.params)
            self.call_cmd()

    def write_jobs(self, name: str, soft=None):
        for cdx, chunks in enumerate(self.cmds_chunks):
            self.get_sh(name, cdx, soft)
            self.write_chunks(chunks)
            self.get_job_name(name, cdx)
            self.write_script(soft)
