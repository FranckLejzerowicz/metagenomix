# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import pandas as pd


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
