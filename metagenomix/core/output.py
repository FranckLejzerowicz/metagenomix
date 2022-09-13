# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import os
from glob import glob
import pandas as pd
import pkg_resources
from os.path import getsize, isdir, isfile, splitext


class Output(object):

    def __init__(self, folder, soft):
        self.soft = soft
        self.soft_dir = folder + '/' + soft
        # summary information about
        self.name = None
        self.after = None
        self.script = None
        self.after_dir = None
        self.directives = {'account', 'partition', 'gres',
                           'time', 'mem', 'ntasks', 'job-name'}
        self.oe = []
        self.hpc = []
        self.input = []
        self.outputs = {}
        self.sample_pool = {}
        self.run_scripts = {}
        self.dir_scripts = {}

    def get_afters(self):
        for folder in os.listdir(self.soft_dir):
            if folder.startswith('after_'):
                key = tuple(folder.lstrip('after_').rsplit('_', 1))
                self.outputs[key] = {'stdouts': {}, 'results': {}, 'sizes': {}}

    def get_outputs(self):
        for after in list(self.outputs):
            self.after = after
            self.get_contents()
            self.get_results()
            self.parse_provenance()
            self.parse_run_sh()
            self.make_table(after)

    def get_contents(self):
        self.after_dir = '%s/after_%s' % (self.soft_dir, '_'.join(self.after))
        all_scripts = set(glob('%s/jobs/*.*' % self.after_dir))
        self.outputs[self.after]['all_scripts'] = all_scripts

    def get_results(self):
        for data_dir in [x for x in os.listdir(self.after_dir) if x != 'jobs']:
            cur_dir = self.after_dir + '/' + data_dir
            if not isdir(cur_dir):
                continue
            self.outputs[self.after]['results'][data_dir] = {}
            self.outputs[self.after]['sizes'][data_dir] = 0
            for root_, dirs, fs in os.walk(cur_dir):
                if root_ == cur_dir:
                    continue
                root = root_.split('%s/' % cur_dir)[-1]
                sizes = sum([getsize('%s/%s' % (root_, x)) for x in fs])
                self.outputs[self.after]['sizes'][data_dir] += sizes
                self.outputs[self.after]['results'][data_dir].update({root: fs})

    def parse_provenance(self):
        self.sample_pool[self.after] = 'sample'
        provenance = '%s/provenance.txt' % self.after_dir
        if isfile(provenance):
            with open(provenance) as provenance_handle:
                for line in provenance_handle:
                    if not line.strip():
                        break
                    if line.split()[1] == 'pooling':
                        self.sample_pool[self.after] = 'pool'
                        break

    def parse_run_sh(self):
        run_sh = '%s/run.sh' % self.after_dir
        if isfile(run_sh):
            scripts = set()
            with open(run_sh) as run_handle:
                for ldx, line in enumerate(run_handle):
                    if ldx < 2:
                        continue
                    self.script = line.strip().split()[-1]
                    scripts.add(self.script)
                    if self.script.endswith('.slm'):
                        self.parse_slm()
                        self.job_outputs()
                        self.outputs[self.after][self.name] = self.script
                    self.get_inputs()
            self.outputs[self.after]['run_scripts'] = scripts

    def make_table(self, after):
        for (tab, attr) in [
            ('oe', self.oe),
            ('hpc', self.hpc),
            ('input', self.input),
        ]:
            self.outputs[after][tab] = pd.DataFrame(attr)
        self.oe = []
        self.hpc = []
        self.input = []

    def parse_slm(self):
        directives = {}
        with open(self.script) as f:
            for line in f:
                if line.startswith('#SBATCH') and '=' in line:
                    key, value = line.strip().split('SBATCH ')[-1].split('=')
                    if key.strip('-') in self.directives:
                        directives[key.split('-')[-1]] = value
        self.name = directives.get('name')
        self.hpc.append(directives)

    def job_outputs(self):
        o_fps = '%s/jobs/output/*%s*.o' % (self.after_dir, self.name)
        for sdx, stdout in enumerate(sorted(glob(o_fps))):
            job_id = stdout[:-2].split('_')[-1]
            self.outputs[self.after]['stdouts'][(self.name, job_id)] = stdout

            status = {'name': self.name, 'done': 'N', 'id': job_id,
                      'cpu_usage': None, 'memory_usage': None, 'error': None}

            ls = open(stdout).readlines()[-25:]
            if 'Done!' in ls:
                status['done'] = 'Y'
            if 'Task and CPU usage stats:' in ls:
                line = ls[ls.index('Task and CPU usage stats:') + 4].split()
                status['cpu_usage'] = ' '.join([line[2], line[6]])
            if 'Memory usage stats:' in ls:
                line = ls[ls.index('Memory usage stats:') + 4].split()
                status['memory_usage'] = line[3]

            stderr = '%s.e' % splitext(stdout)[0]
            ls = open(stderr).readlines()[-1]
            if 'error' in ls:
                status['error'] = ls
            self.oe.append(status)

    def get_inputs(self):
        cols = [self.sample_pool[self.after], 'tech', 'group', 'genome_taxon']
        if isfile(self.script):
            with open(self.script) as script_handle:
                for line in script_handle:
                    if line.startswith('# data: '):
                        data = line[8:].strip().split()
                        if len(data) < 4:
                            data += [None] * (4 - len(data))
                        data = {col: data[x] for x, col in enumerate(cols)}
                        data['name'] = self.name
                        self.input.append(data)


class Softwares(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.inputs = {'res': {}, 'dir': set(), 'pip': set(), 'usr': set()}
        self.names = set()
        self.roles = {}
        self.softs = {}
        self.get_softs()

    def _from_dir(self):
        self.inputs['dir'].update([soft for soft in os.listdir(self.dir)])

    def _from_usr(self):
        for soft in self.softwares:
            if soft in self.inputs['dir']:
                self.inputs['usr'].add(soft)
            else:
                r = re.compile(soft.replace('*', '.*').replace('..*', '.*'))
                self.inputs['usr'].update(
                    set(filter(r.match, list(self.inputs['dir']))))

    def _from_res(self):
        res = pkg_resources.resource_filename("metagenomix", "resources")
        with open('%s/softwares.txt' % res) as f:
            for line in f:
                soft, role = line.strip().split('\t')
                self.roles.setdefault(role.split(' (')[0], set()).add(soft)
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
        self.show()
