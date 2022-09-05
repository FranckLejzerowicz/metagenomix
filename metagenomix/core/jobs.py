# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import json
import yaml
import random
import hashlib
import subprocess
import numpy as np
import datetime as dt
from os.path import dirname, isdir, splitext

from metagenomix._io_utils import (
    mkdr, get_roundtrip, print_status_table, get_runs)


class CreateScripts(object):

    def __init__(self, config, databases, workflow, commands):
        self.config = config
        self.databases = databases
        self.commands = commands
        self.graph = workflow.graph
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
        self.time = dt.datetime.now().strftime("%d/%m/%Y-%H:%M")
        self.scripts = []
        self.hash = ''
        self.hashed = ''

    def make_dirs(self):
        for name, soft in self.commands.softs.items():
            for directory in sorted(soft.dirs):
                if not isdir(directory):
                    os.makedirs(directory)

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
        if soft:
            main = '%s/run.sh' % soft.dir
            self.run['software'][name] = main
        else:
            main = '%s/run_%s.sh' % (self.sh.rsplit('/', 2)[0], name)
            self.run['database'][name] = main
        return main

    def write_main(self, name, soft=None) -> None:
        main = self.get_main_sh(name, soft)
        with open(main, 'w') as o:
            if len(self.job_fps):
                job_dir = dirname(self.job_fps[0])
                o.write('mkdir -p %s/output\n' % job_dir)
                o.write('cd %s/output\n' % job_dir)
                for jdx, job_fp in enumerate(self.job_fps):
                    o.write('%s %s\n' % (self.scheduler, job_fp))
                    if self.config.dev:
                        if not isdir('%s/output' % job_dir):
                            os.makedirs('%s/output' % job_dir)
                        f = 'slurm-%s_000000%s' % (self.job_name, jdx)
                        for e in ['.o', '.e']:
                            job_std = '%s/output/%s%s' % (job_dir, f, e)
                            with open(job_std, 'w') as o_dev:
                                for r in range(10, random.randrange(11, 100)):
                                    o_dev.write('random %s\n' % r)
                                o_dev.write('Done!\n')
                self.job_fps = []

    def get_provenance_fp(self, name, soft=None) -> str:
        if soft:
            fp = '%s/provenance_%s.txt' % (soft.dir, self.hashed)
        else:
            fp = '%s/provenance_%s_%s.txt' % (self.sh.rsplit('/', 2)[0],
                                              name, self.hashed)
        return fp

    def write_provenance(self, name, soft=None) -> None:
        steps = self.graph.paths[soft.name][0]
        step_max = max([len(x) for x in steps])
        provenance_fp = self.get_provenance_fp(name, soft)
        with open(provenance_fp, 'w') as o:
            o.write("Pipeline steps to this output (and analysis type):\n")
            for sdx, step in enumerate(steps):
                o.write('%s\t%s%s\t: %s\n' % (
                    sdx, step, (' ' * (step_max-len(step))),
                    self.config.tools.get(step, '')))
            o.write("\nParameters for these steps' outputs:\n")
            for sdx, step in enumerate(steps):
                if step not in self.commands.softs:
                    o.write('\n ===== %s: %s (no param) =====\n' % (sdx, step))
                    continue
                o.write('\n========== %s: %s ==========\n' % (sdx, step))
                params_show = self.commands.softs[step].params.copy()
                if 'databases' in params_show:
                    databases = params_show['databases']
                    del params_show['databases']
                    params = {'search': {soft.name.split('_')[-1]: params_show,
                                         'databases': databases}}
                    o.write('%s\n' % yaml.dump(params))
                else:
                    o.write('%s\n' % yaml.dump(params_show))

    def database_cmds(self):
        for db, cmds in self.databases.commands.items():
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

    def get_cmds(self, soft):
        self.cmds = {}
        for sam_or_pool, cmds in soft.cmds.items():
            if isinstance(cmds, list):
                self.scratch(soft, sam_or_pool, cmds)
            elif isinstance(cmds, dict):
                if sam_or_pool in self.commands.pools:
                    for group, group_cmds in cmds.items():
                        self.scratch(soft, (sam_or_pool, group), group_cmds)
                else:
                    self.scratch(soft, sam_or_pool, cmds)
            else:
                sys.exit('The collected commands are neither list of dict!')

    def print_status(self, m, sdx, name, soft):
        gap = (m - len(name) - len(str(sdx))) + 1
        print('\t%s [%s] %s%s' % (sdx, name, ('.' * gap), ('.' * 8)), end=' ')
        print_status_table(soft, self.config.show_status)

    def get_hash(self):
        h = hashlib.blake2b(digest_size=10)
        h.update(self.hash.encode('utf-8'))
        hashed = str(h.hexdigest())
        self.hashed = hashed

    def get_soft_hash(self, soft=None):
        dflts = ['time', 'nodes', 'mem', 'mem_dim', 'env', 'chunks',
                 'scratch', 'machine', 'partition']
        if soft:
            steps = self.graph.paths[soft.name]
            params = dict(x for x in soft.params.items() if x[0] not in dflts)
            self.hash += ''.join([
                json.dumps(step) + json.dumps(params) for step in steps])
            self.get_hash()

    def versioning(self):
        runs_dir = '%s/_created' % self.config.dir
        if not isdir(runs_dir):
            os.makedirs(runs_dir)
        self.write_runs(runs_dir)

    def write_runs(self, runs_dir):
        runs_fp = runs_dir + '/' + self.hashed + '.txt'
        runs = get_runs(runs_fp)
        with open(runs_fp, 'w') as o:
            for script in self.scripts:
                o.write('%s\n' % script)
            o.write('\n--------- runs ---------\n')
            for run in runs:
                o.write('Date: %s\n' % run)
            o.write('Date: %s\n' % self.time)
            o.write('\n----------------------------\n')
        print('\n[config run] Written: %s' % runs_fp)

    def software_cmds(self):
        m = max(len(x) for x in self.commands.softs) + 1
        for sdx, (name, soft) in enumerate(self.commands.softs.items()):
            self.print_status(m, sdx, name, soft)
            if not len(soft.cmds):
                continue
            self.get_modules(name)
            self.get_cmds(soft)
            self.get_chunks(soft.params['chunks'])
            self.write_jobs(name, soft)
            self.get_soft_hash(soft)
            self.write_main(name, soft)
            self.write_provenance(name, soft)

    def get_sh(self, name: str, chunk: str, soft=None) -> None:
        """

        Parameters
        ----------
        name : str
            Name of the current software of the pipeline workflow
        chunk : str
            Name of the current chunk of commands
        soft
        """
        if soft:
            self.sh = '%s/jobs/run_%s.sh' % (
                soft.dir, chunk.replace('/', '_'))
        else:
            self.sh = '%s/databases/jobs/build_%s_%s.sh' % (
                self.config.dir, name, chunk)
        mkdr(dirname(self.sh))

    def prep_script(self, params: dict) -> None:
        # mandatory options
        self.cmd = [
            'Xhpc',
            '-i', self.sh,
            '-j', self.job_name,
            '-t', str(params['time']),
            '-c', str(params['cpus']),
            '-M', str(params['mem']), params['mem_dim'],
            '--no-stat']
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
        # machine-specific setup: env activating and slurm partition
        if params['machine']:
            self.cmd.append('--%s' % params['machine'])
        if params['partition']:
            self.cmd.append('--p-partition %s' % params['partition'])
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
        if not self.config.verbose:
            self.cmd.append('--quiet')

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

    def get_job_name(self, name: str, chunk: str, soft=None):
        self.job_name = name + '.' + self.pjct
        if soft:
            self.job_name += '.' + ''.join([
                x for x in soft.dir.split('after_')[-1] if x not in 'aeiuoy'])
        self.job_name += '.' + chunk.replace('/', '_')

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
        for chunk, chunk_keys in self.chunks.items():
            self.get_sh(name, chunk, soft)
            self.write_chunks(chunk_keys, soft)
            self.get_job_name(name, chunk, soft)
            self.write_script(soft)

    def display(self):
        for db_soft, name_main in self.run.items():
            if len(self.run[db_soft]):
                print()
                soft_print = '========== %s scripts ========== ' % db_soft
                print(soft_print)
                self.scripts.append(soft_print)
            for name, main in name_main.items():
                main_print = '>%s\nsh %s' % (name, main)
                print(main_print)
                self.scripts.append(main_print)
