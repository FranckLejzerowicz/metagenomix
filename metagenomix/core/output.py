# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
import glob
import pandas as pd
import pkg_resources


class Output(object):

    def __init__(self, folder, soft):
        self.soft = soft
        self.path = folder + '/' + soft
        self.afters = self.get_afters()
        self.table = self.init_table()
        # summary information about
        self.outputs = {}  # software outputs
        self.jobs = {}  # software jobs
        self.sbatch = {'account', 'partition', 'gres', 'time', 'mem', 'ntasks'}
        self.res = {}
        self.stdoes = {}
        self.manage()

    def get_afters(self) -> dict:
        afters = {}
        for folder in os.listdir(self.path):
            if folder.startswith('after_'):
                afters[folder[6:]] = {}
        return afters

    def manage(self):
        """Collect info to submit to user to obtain management decisions."""
        self.manage_jobs()  # for the job files (jobs and stdout/err files)
        # self.manage_outputs()  # for the output files (results)

    def init_table(self) -> pd.DataFrame:
        tab_d = dict((x, []) for x in ['after', 'content'])
        for attr in ['jobs', 'remove', 'rename']:
            if hasattr(self, attr):
                tab_d[attr] = []
        return pd.DataFrame(tab_d)

    def job_outputs(self, job_name):
        stdouts = glob.glob('%s/jobs/output/*/*%s*.o' % (self.path, job_name))
        for sdx, stdout in enumerate(stdouts):
            self.stdoes.setdefault(sdx, []).append(stdout[:-1])

    def update_jobs(self, job, ext):
        self.jobs.setdefault('file', []).append(job)
        self.jobs.setdefault('job', []).append(job.split('run_')[-1].strip(ext))
        for key in self.sbatch:
            self.jobs.setdefault(key, []).append(self.res.get(key, ''))
        for stdoe in self.stdoes:
            self.jobs.setdefault(stdoe, []).append(self.res.get(stdoe, ''))

    def manage_jobs(self):
        rad = self.path + '/jobs/*'
        for ext in ['.sh', '.slm']:
            print(rad + ext)
            jobs = glob.glob(rad + ext)
            print(jobs)
            for job in jobs:
                if ext == '.slm':
                    self.parse_slm(job)
                self.update_jobs(job, ext)
                self.stdoes = {}
                self.res = {}

    def parse_slm(self, job):
        job_name = None
        print(job)
        with open(job) as f:
            for line in f:
                if '--job-name' in line:
                    job_name = line.strip().split('=')[-1]
                elif line.startswith('#SBATCH') and '=' in line:
                    key, value = line.strip().split('SBATCH ')[-1].split('=')
                    if key.strip('-') in self.sbatch:
                        self.res[key.strip('-')] = value
        self.job_outputs(job_name)
        print(self.stdoes)
        print(self.stdoesfd)

    def parse_after(self, folder):
        for root, folders, files in os.walk(folder):
            if len(files):
                print(root)
                print(folders)
                print(files)
                print(filesfdsa)

    def manage_outputs(self):
        folders = [x for x in os.listdir(self.path) if x != 'jobs']
        for folder in folders:
            self.parse_after(folder)


class Softwares(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.inputs = {'res': {}, 'dir': set(), 'pip': set(), 'usr': set()}
        self.names = set()
        self.roles = {}
        self.softs = {}
        self.get_softs()

    def _from_dir(self):
        self.inputs['dir'].update([soft for soft in os.listdir(self.folder)])

    def _from_usr(self):
        for soft in self.softwares:
            if soft in self.inputs['dir']:
                self.inputs['usr'].add(soft)
            else:
                r = re.compile(soft)
                self.inputs['usr'].update(
                    set(filter(r.match, list(self.inputs['dir']))))

    def _from_res(self):
        res = pkg_resources.resource_filename("metagenomix", "resources")
        with open('%s/softwares.txt' % res) as f:
            for line in f:
                soft, role = line.strip().split('\t')
                self.roles.setdefault(role.split(' (')[0],set()).add(soft)
                self.inputs['res'][soft] = role

    def _from_pipeline(self):
        if self.pipeline:
            with open(self.pipeline) as f:
                self.inputs['pip'].update([x.strip().split()[-1] for x in
                                           f if x.strip() and x[0] != '#'])

    def _intersection(self) -> set:
        softs = self.inputs['dir'] & set(self.inputs['res'])
        if self.inputs['pip'] | self.inputs['usr']:
            softs = softs & (self.inputs['pip'] | self.inputs['usr'])
        return softs

    def _to_manage(self) -> dict:
        softs_intersection = self._intersection()
        for role, softs in self.roles.items():
            cur_softs = softs & softs_intersection
            if cur_softs:
                self.softs[role] = {t: self.inputs['res'][t] for t in cur_softs}
                self.names.update(cur_softs)

    def show(self):
        print('Tools to %s per role:' % self.command)
        for role, softs in self.softs.items():
            print('  [Role: %s]' % role)
            n = max([len(soft) for soft in softs])
            roles = set(softs.values())
            for soft, role_ in softs.items():
                print('\t- %s' % soft, end='')
                if role_ == role or len(roles) == 1:
                    print()
                else:
                    print(' ' * (n-len(soft)), '\t:', role_.split('(')[-1][:-1])

    def get_softs(self):
        print('* Getting softwares to %s' % self.command)
        self._from_dir()
        self._from_usr()
        self._from_res()
        self._from_pipeline()
        self._to_manage()
        if self.verbose:
            self.show()
