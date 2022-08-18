# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import pandas as pd
import datetime as dt
from os.path import basename


def manager(**kwargs):
    """Run metagenomix to assist user renaming or discarding outputs.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    managing = Managed(**kwargs)
    print(managing.__dict__)


class Managed(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.get_out()
        self.table = self.init_table()

    def init_table(self):
        tab_d = dict((x, []) for x in ['software', 'after', 'content'])
        for attr in ['jobs', 'remove', 'rename']:
            if hasattr(self, attr):
                tab_d[attr] = []
        return pd.DataFrame(tab_d)

    def get_out(self):
        out_dir = '%s/_management' % self.dir
        if self.out is None:
            time = re.sub('[. :]', '-', str(dt.datetime.now()).split('.')[0])
            self.out = out_dir + '/' + time + '.txt'
        else:
            base = basename(self.out)
            if '/' in self.out:
                print('Using basename "%s" to write in "%s"' % (out_dir, base))
            self.out = out_dir + '/' + base


