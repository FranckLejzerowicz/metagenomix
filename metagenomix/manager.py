# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import pkg_resources
import pandas as pd
import datetime as dt
from os.path import basename
from metagenomix.core.output import Output


def manager(**kwargs):
    """Run metagenomix to assist user renaming or discarding outputs.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix manager ===\n')
    managing = Managed(**kwargs)

    print('* setting up the output file name')
    managing.get_out()
    print(' -> %s' % managing.out)

    print('* getting softwares to manage')
    managing.get_softs()

    print('* parsing folders to collect info about jobs and outputs')
    managing.parse_softs()

    for name, soft in managing.data.items():
        print()
        print('software:', name)
        print(pd.DataFrame(soft.jobs))


class Managed(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.softs = {'res': {}, 'dir': set(), 'pip': set(), 'usr': set()}
        self.time = dt.datetime.now().strftime("%d/%m/%Y, %H") + 'h'
        self.manage_dir = '%s/_management' % self.folder
        self.roles = {}
        self.data = {}

    def _from_dir(self):
        self.softs['dir'].update([soft for soft in os.listdir(self.folder)])

    def _from_usr(self):
        for soft in self.softwares:
            if soft in self.softs['dir']:
                self.softs['usr'].add(soft)
            else:
                r = re.compile(soft)
                self.softs['usr'].update(
                    set(filter(r.match, list(self.softs['dir']))))

    def _from_res(self):
        res = pkg_resources.resource_filename("metagenomix", "resources")
        with open('%s/softwares.txt' % res) as f:
            for line in f:
                soft, role = line.strip().split('\t')
                self.roles.setdefault(role.split(' (')[0],set()).add(soft)
                self.softs['res'][soft] = role

    def _from_pipeline(self):
        if self.pipeline:
            with open(self.pipeline) as f:
                self.softs['pip'].update([line.strip().split()[-1] for line in
                                          f if line.strip() and line[0] != '#'])

    def _intersection(self) -> set:
        softs = self.softs['dir'] & set(self.softs['res'])
        if self.softs['pip'] | self.softs['usr']:
            softs = softs & (self.softs['pip'] | self.softs['usr'])
        return softs

    def _to_manage(self) -> dict:
        softs_intersection = self._intersection()
        role_softs = {}
        for role, softs in self.roles.items():
            cur_softs = softs & softs_intersection
            if cur_softs:
                role_softs[role] = {t: self.softs['res'][t] for t in cur_softs}
        return role_softs

    def show(self):
        print('Tools to manage per role:')
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
        self._from_dir()
        self._from_usr()
        self._from_res()
        self._from_pipeline()
        self.softs = self._to_manage()
        if self.verbose:
            self.show()

    def get_out(self):
        if self.out is None:
            path = self.time + '.txt'
        else:
            path = basename(self.out)
            if '/' in self.out:
                print('Using "%s" to write in "%s"' % (path, self.manage_dir))
        self.out = self.manage_dir + '/' + path

    def parse_softs(self):
        """An Output class instance is created for each software to manage,
        and placed as value to the dict with the software name of key."""
        for role, softs in self.softs.items():
            for soft in softs:
                self.data[soft] = Output(self.folder, soft)
