# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import datetime as dt
from os.path import basename, isdir

from metagenomix.metagenomix import metagenomix
from metagenomix._io_utils import print_status_table
from metagenomix.core.output import Output


def monitoring(**kwargs):
    """Show the status of the planned outputs."""
    print('\n === metagenomix checker ===\n')
    # Collect all command and init the script creating instance
    monitored = Monitored(**kwargs)

    print('* setting up the output file name')
    monitored.make_status_dir()
    monitored.get_out()

    print('* collecting and showing the current status of the analyses')
    monitored.monitor_status()
    monitored.write_status()

    monitored.parse_softs()
    monitored.monitor_softs()

    for name, soft in monitored.data.items():
        print()
        print('software:', name)
        print(pd.DataFrame(soft.jobs))


class Monitored(object):

    def __init__(self, **kwargs):
        kwargs['localscratch'] = None
        kwargs['userscratch'] = False
        kwargs['purge_pfams'] = None
        kwargs['show_params'] = None
        kwargs['show_pfams'] = None
        kwargs['scratch'] = False
        kwargs['chunks'] = None
        kwargs['jobs'] = None
        config, databases, workflow, commands = metagenomix(**kwargs)
        self.__dict__.update(kwargs)
        self.config = config
        self.databases = databases
        self.graph = workflow.graph
        self.commands = commands
        # self.softs = {'res': {}, 'dir': set(), 'pip': set(), 'usr': set()}
        # USE THE SOFTWARES OF THE PARSED COMMANDS----
        self.time = dt.datetime.now().strftime("%d/%m/%Y, %H") + 'h'
        self.status_dir = '%s/_status' % self.output_dir
        self.roles = {}
        self.data = {}

    def monitor_status(self):
        m = max(len(x) for x in self.commands.softs) + 1
        for sdx, (name, soft) in enumerate(self.commands.softs.items()):
            n = (m - len(name) - len(str(sdx))) + 1
            cur_soft = '%s [%s]' % (sdx, name)
            soft.tables.append(cur_soft)
            print('\t%s %s%s' % (cur_soft, ('.' * n), ('.' * 8)), end=' ')
            print_status_table(soft, True)

    def make_status_dir(self):
        if not isdir(self.status_dir):
            os.makedirs(self.status_dir)

    def get_out(self):
        if self.summary_fp is None:
            path = self.time + '.txt'
        else:
            path = basename(self.summary_fp)
            if '/' in self.summary_fp:
                print('Using "%s" to write in "%s"' % (path, self.status_dir))
        self.summary_fp = self.status_dir + '/' + path

    def get_hash(self, soft) -> str:
        steps = self.graph.paths[soft.name]

        hash_string = ''.join([str(step) + str(soft.params) for step in
                               steps])
        hashed = abs(hash(hash_string)) % (10 ** 8)
        return hashed

    def write_status(self):
        with open(self.summary_fp, 'w') as o:
            o.write('# Summary of the data that is currently needed as input\n')
            o.write('# or not yet produced as output.\n')
            o.write('# This file format will evolve...\n')
            o.write('# Date of status summary: %s\n' % self.time)
            for sdx, (name, soft) in enumerate(self.commands.softs.items()):
                # hashed = self.get_hash(soft)
                o.write('\n%s\n' % '\t'.join(soft.tables[:2]))
                for table in soft.tables[2:]:
                    if table is None:
                        o.write(' -> All necessary data available\n')
                    else:
                        o.write('%s\n' % table)
        print('Witten: %s' % self.summary_fp)

    def parse_softs(self):
        """An Output class instance is created for each software to manage,
        and placed as value to the dict with the software name of key."""
        for name, soft in self.commands.softs.items():
            if isdir(self.output_dir + '/' + name):
                self.data[name] = Output(self.output_dir, name)
