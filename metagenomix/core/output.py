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
from os.path import basename, getsize, isdir, isfile, splitext


class Output(object):
    """

    Attributes
    ----------
    soft : str

    soft_dir : str
        Path to the software output folder (without after folder)
    name : str
        Name of the job
    after : tuple
        (previous software, hash value)
    script : str
        Path to the job script file (.slm or .pbs)
    after_dir : str
        Path of a current after folder
    oe : list
        All job standard output files infos
        [{'name': name, 'id': job_id, 'stdout': stdout, 'done': 'N',
          'AllocCPUS': 1, 'MinCPU': None, 'AveCPU': None,
          'MaxRSS': None, 'AveRSS': None, 'error': None}, ...]
    hpc : list

    input : list

    outputs : dict

    """
    def __init__(self, folder, soft):
        self.soft = soft
        self.soft_dir = folder + '/' + soft
        self.name = None
        self.after = None
        self.script = None
        self.after_dir = None
        self.oe = []
        self.hpc = []
        self.input = []
        self.outputs = {}

    def get_outputs(self):
        for folder in glob('%s/after_*' % self.soft_dir):
            self.after = tuple(folder.split('after_')[-1].rsplit('_', 1))
            self.after_dir = folder
            self.outputs[self.after] = {
                'sample_pool': 'sample',
                'provenance': {},
                'results': {},
                'sizes': {}
            }
            self.get_results()
            self.parse_provenance()
            self.parse_run_sh()
            self.job_outputs()
            self.make_table()

    def get_results(self):
        """Collect the output files (values) of each folder (keys) located
        inside each software's non-jobs output folder (and their sizes)."""
        for tech in os.listdir(self.after_dir):
            cur_dir = self.after_dir + '/' + tech
            if not isdir(cur_dir) or tech == 'jobs':
                continue
            self.outputs[self.after]['sizes'][tech] = 0
            for root_, dirs, fs in os.walk(cur_dir):
                if root_ == cur_dir:
                    continue
                sizes = sum([getsize('%s/%s' % (root_, x)) for x in fs])
                self.outputs[self.after]['sizes'][tech] += sizes
                if tech not in self.outputs[self.after]['results']:
                    self.outputs[self.after]['results'][tech] = {}
                folder = root_.split('%s/%s/' % (self.after_dir, tech))[-1]
                self.outputs[self.after]['results'][tech][folder] = fs

    def parse_provenance(self):
        """Collect the provenance as step (key) - software (value) dict and set
        'sample_pool' to 'pool' for post-pooling steps (or set to 'sample')."""
        provenance = '%s/provenance.txt' % self.after_dir
        if isfile(provenance):
            with open(provenance) as provenance_handle:
                for ldx, line in enumerate(provenance_handle):
                    if ldx:
                        if not line.strip():
                            break
                        if line.split()[1] == 'pooling':
                            self.outputs[self.after]['sample_pool'] = 'pool'
                        n, soft = line.split()[:2]
                        self.outputs[self.after]['provenance'][n] = soft

    def parse_run_sh(self):
        """Parse the latest run.sh script to access its run_*.slm job scripts
        and collect their directives, status and outputs."""
        run_sh = '%s/run.sh' % self.after_dir
        if isfile(run_sh):
            with open(run_sh) as run_handle:
                for ldx, line in enumerate(run_handle):
                    if ldx < 2:
                        continue
                    self.script = line.strip().split()[-1]
                    if self.script.endswith('.slm') and isfile(self.script):
                        self.parse_scripts()
                    self.get_inputs()

    def make_table(self):
        for (tab, attr) in [
            ('oe', self.oe),
            ('hpc', self.hpc),
            ('input', self.input),
        ]:
            self.outputs[self.after][tab] = pd.DataFrame(attr)
        self.oe = []
        self.hpc = []
        self.input = []

    def parse_scripts(self):
        """Set to the value of class attribute 'name' to the job name and add
        the directives as a dict to the 'hpc' attribute (to cast a table)."""
        directives = {'script': self.script}
        with open(self.script) as f:
            for line in f:
                if not line.strip():
                    break
                if line.startswith('#SBATCH') and '=' in line:
                    k, v = line.strip().split('SBATCH ')[-1].split('=')
                    if k.strip('-') in {
                        'account',
                        'partition',
                        'gres',
                        'time',
                        'mem',
                        'ntasks',
                        'job-name'
                    }:
                        directives[k.split('-')[-1]] = v
        self.name = directives.get('name')
        self.hpc.append(directives)

    def job_outputs(self):
        o_fps = '%s/jobs/output/*.o' % self.after_dir
        for sdx, stdout in enumerate(sorted(glob(o_fps))):
            name, job_id = basename(stdout[:-2]).split('-', 1)[1].rsplit('_', 1)
            status = {
                'name': name, 'id': job_id,
                'stdout': stdout, 'done': 'N',
                'AllocCPUS': 1, 'MinCPU': None,
                'AveCPU': None, 'MaxRSS': None,
                'AveRSS': None, 'error': None
            }
            ls = open(stdout).readlines()[-25:]
            uses = {
                'Memory usage stats:\n': ['MaxRSS', 'AveRSS'],
                'Task and CPU usage stats:\n': ['AllocCPUS', 'MinCPU', 'AveCPU']
            }
            for use, vals in uses.items():
                if use in ls:
                    head, stat = ls[ls.index(use)+1], ls[ls.index(use)+4]
                    for v in vals:
                        if stat[head.index(v)+len(v)-1].strip():
                            status[v] = stat[:head.index(v)+len(v)].split()[-1]
            if 'Done!\n' in ls:
                status['done'] = 'Y'

            stderr = '%s.e' % splitext(stdout)[0]
            ls = open(stderr).readlines()
            if ls and 'error' in ls[-1]:
                status['error'] = ls[-1]
            self.oe.append(status)

    def get_inputs(self):
        cols = [self.outputs[self.after]['sample_pool'],
                'tech', 'group', 'genome_taxon']
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
    """This class gets the softwares to export, manage or monitor.
    For these softwares, the class will have two attribues:

    This class is called by:
        class Exported(object)
        class Manage(object)
        class Monitored(object)

    Attributes
    ----------
    names : set
        All softwares to export, manage or monitor
    roles : dict
        Set of softwares (values) per role (keys)
    softs : dict
        Dictionaries of {software: full_role} (values) per short role (keys)
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.inputs = {'res': {}, 'dir': set(), 'pip': set(), 'usr': set()}
        self.names = set()
        self.roles = {}
        self.softs = {}
        self.get_softs()

    def _from_dir(self):
        """Get the softwares that have a folder in the output directory."""
        self.inputs['dir'].update([soft for soft in os.listdir(self.dir)])

    def _from_usr(self):
        """Get user-defined softwares but only those that have a folder in
        the output directory."""
        # for every command-line software(s) (`-s, --software`)
        for soft in self.softwares:
            if soft in self.inputs['dir']:
                # add it the user softwares if it has as folder
                self.inputs['usr'].add(soft)
            else:
                # add all its regex matches it this regex match to folder(s)
                r = re.compile(soft.replace('*', '.*').replace('..*', '.*'))
                self.inputs['usr'].update(
                    set(filter(r.match, list(self.inputs['dir']))))

    def _from_res(self):
        """Get the softwares in the `softwares.txt` resource file."""
        # get the paths resources folder
        res = pkg_resources.resource_filename("metagenomix", "resources")
        with open('%s/softwares.txt' % res) as f:
            # parse the `softwares.txt` file to collect:
            for line in f:
                soft, role = line.strip().split('\t')
                # set of software names per role (category)
                self.roles.setdefault(role.split(' (')[0], set()).add(soft)
                # role (in full) or each software name
                self.inputs['res'][soft] = role

    def _from_pip(self):
        """Get the softwares scheduled to run the pipeline."""
        # it could be that for exporting, there is no pipeline configuration
        if self.pipeline_tsv:
            # but if there is a pipeline configuration
            with open(self.pipeline_tsv) as f:
                # collect then softwares scheduled to run the pipeline
                self.inputs['pip'].update([x.strip().split()[-1] for x in
                                           f if x.strip() and x[0] != '#'])

    def _intersection(self) -> set:
        """Get the set of softwares in the intersection of softwares from the
        various softwares found in the output folder, the resources file,
        and either the pipeline and user-defined softwares (if monitoring:
        restrictive scope) ot the pipeline or user-defined softwares (if
        managing: broader scope).
        """
        # softwares have to be present in both the output folder and resource
        softs = self.inputs['dir'] & set(self.inputs['res'])
        # get the union of softwares from the pipeline and passed by the user
        pip_usr = self.inputs['pip'] | self.inputs['usr']
        # if at least one software was given by user (command-line or pipeline)
        if pip_usr:
            if self.command == 'monitor':
                # restrict monitoring to softwares given by all these inputs
                if self.inputs['pip']:
                    softs = softs & self.inputs['pip']
                if self.inputs['usr']:
                    softs = softs & self.inputs['usr']
            else:
                # restrict exporting or management to softwares across inputs
                softs = softs & pip_usr
        return softs

    def _in_scope(self) -> dict:
        """Collect softwares per role and the full list of softwares in scope"""
        # get softwares in the intersection of softwares
        softs_intersection = self._intersection()
        for role, softs in self.roles.items():
            cur_softs = softs & softs_intersection
            if cur_softs:
                # e.g., {'preprocessing': {'fastp': 'preprocessing'}}
                self.softs[role] = {t: self.inputs['res'][t] for t in cur_softs}
                # add to the full set of softwares in scope
                self.names.update(cur_softs)

    def show(self):
        """Print the softwares in the scope of softwares to export, manage
        or monitor. For example:
        ```
        Tools to monitor per role:
          [Role: preprocessing]
                - fastp
          [Role: paired-read merging]
                - bbmerge
                - flash
          [Role: coassembly setup]
                - pooling
          [Role: MAG]
                - drep
        ```
        """
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
        """Get the softwares to export, manage or monitor."""
        print('* Getting softwares to %s' % self.command)
        self._from_dir()  # fills self.inputs['dir'] with folder'd software
        self._from_usr()  # fills self.inputs['usr'] with user-defined softwares
        self._from_res()  # fills self.inputs['res'] with resource softwares
        self._from_pip()  # fills self.inputs['pip'] with pipeline softwares
        self._in_scope()  # restrict scope of softwares
        self.show()
