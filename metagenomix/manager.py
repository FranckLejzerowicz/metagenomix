# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import datetime as dt
from os.path import basename
from metagenomix.core.output import Softwares, Output


def manager(**kwargs):
    """Run metagenomix to assist user renaming or discarding outputs.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix manager ===\n')
    kwargs['command'] = 'manage'
    managing = Manage(**kwargs)

    print('* Setting up the output file name')
    managing.get_out()
    print(' -> %s' % managing.out)

    print('* Parsing folders to collect info about jobs and outputs')
    managing.parse_softs()

    for name, soft in managing.data.items():
        print()
        print('Software:', name)
        print(pd.DataFrame(soft.jobs))


class Manage(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.softs = Softwares(**kwargs)
        self.time = dt.datetime.now().strftime("%d/%m/%Y, %H") + 'h'
        self.manage_dir = '%s/_menagament' % self.folder
        self.data = {}

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
