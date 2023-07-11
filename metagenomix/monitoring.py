# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import datetime as dt
from os.path import abspath, isdir

from metagenomix.metagenomix import metagenomix
from metagenomix._io_utils import print_status_table
from metagenomix.core.output import Output, Softwares


def monitoring(**kwargs):
    """Show the status of the planned outputs."""

    print('\n>>> `metagenomix monitor` started >>>\n')

    kwargs['command'] = 'monitor'
    # Collect all command and init the script creating instance
    monitor = Monitored(**kwargs)
    print("monitor.dir:", monitor.dir)
    print("monitor.time:", monitor.time)
    print("monitor.log_dir:", monitor.log_dir)
    print("monitor.roles:", monitor.roles)
    print("monitor.monitored:", monitor.monitored)
    print("monitor.softwares.names:", monitor.softwares.names)
    print("monitor.softwares.softs:", monitor.softwares.softs)
    print('* setting up the output file name')
    monitor.make_status_dir()
    monitor.get_out()
    print("monitor.result_fp:", monitor.result_fp)

    print('* collecting and showing the current status of the analyses')
    print('monitor.monitor_status()')
    monitor.monitor_status()
    print('monitor.get_status()')
    monitor.get_status()
    print('monitor.write_status()')
    monitor.write_status()
    print('monitor.parse_softs()')
    monitor.parse_softs()
    print('monitor.monitor_softs()')
    monitor.monitor_softs()

    print('\n<<< `metagenomix monitor` completed <<<\n')


class Monitored(object):

    def __init__(self, **kwargs):
        kwargs['chunks'] = None
        kwargs['force'] = False
        kwargs['verbose'] = False
        kwargs['scratch'] = False
        kwargs['show_pfams'] = None
        kwargs['purge_pfams'] = None
        kwargs['userscratch'] = False
        kwargs['localscratch'] = None
        kwargs['dir'] = abspath(kwargs['output_dir'])
        config, databases, workflow, commands = metagenomix(**kwargs)
        self.__dict__.update(kwargs)
        self.softwares = Softwares(**kwargs)
        self.config = config
        self.databases = databases
        self.graph = workflow.graph
        self.commands = commands
        self.time = dt.datetime.now().strftime("%d-%m-%Y_%H-%M")
        self.log_dir = '%s/_monitored' % config.dir
        self.roles = {}
        self.monitored = {}

    def monitor_status(self):
        """Print enumerated pipeline softwares and
        """
        # get a text padding length
        m = max(len(x) for x in self.commands.softs) + 1
        # for each software in enumeration
        for sdx, (name, soft) in enumerate(self.commands.softs.items()):
            # skip softwares that not registered in the resources/softwares.txt
            print(self.softwares.names)
            if name not in self.softwares.names:
                continue
            # get the current software name length per padding
            n = (m - len(name) - len(str(sdx))) + 1
            cur_soft = '%s [%s]' % (sdx, name)
            # print and collect the print for pretty´-table printing
            print('\t%s %s%s' % (cur_soft, ('.' * n), ('.' * 8)), end=' ')
            soft.tables.append(cur_soft)
            print('--------------1')
            print_status_table(soft, True)
            print('--------------2')
            print(gfdsa)

    def make_status_dir(self):
        if not isdir(self.log_dir):
            os.makedirs(self.log_dir)

    def get_out(self):
        if self.result_fp is None:
            self.result_fp = self.log_dir + '/' + self.time + '.txt'
        elif '/' in self.result_fp:
            self.result_fp = abspath(self.result_fp)
        else:
            self.result_fp = self.log_dir + '/' + self.result_fp

    def get_status(self):
        for name, soft in self.commands.softs.items():
            if name not in self.softwares.names:
                continue
            print()
            print("name:", name)
            print("soft.io")
            print(soft.io)
            print("soft.status")
            print(soft.status)
            print("soft.links")
            print(soft.links)

    def write_status(self):
        with open(self.result_fp, 'w') as o:
            o.write('# Summary of the data that is currently needed as input\n')
            o.write('# or not yet produced as output.\n')
            o.write('# This file format will evolve...\n')
            o.write('# Date of status results: %s\n' % self.time)
            for sdx, (name, soft) in enumerate(self.commands.softs.items()):
                if name not in self.softwares.names:
                    continue
                o.write('\n%s\n' % '\t'.join(soft.tables[:2]))
                for table in soft.tables[2:]:
                    if table is None:
                        o.write(' -> All necessary data available\n')
                    else:
                        o.write('%s\n' % table)
        print('Written: %s' % self.result_fp)

    def parse_softs(self):
        """An Output class instance is created for each software to manage,
        and placed as value to the dict with the software name of key."""
        for name, soft in self.commands.softs.items():
            if name not in self.softwares.names:
                continue
            if isdir(self.output_dir + '/' + name):
                output = Output(self.output_dir, name)
                output.get_outputs()
                self.monitored[name] = output.outputs

    def monitor_softs(self):
        for name, outputs in self.monitored.items():
            print()
            print()
            print('#' * 40)
            print('name:', name)
            print('#' * 40)
            for (after, hash_value), data in outputs.items():
                print()
                print('name:', name)
                print("after:", after)
                print("hash_value:", hash_value)
                for k, d in data.items():
                    print()
                    print('-----------------------')
                    print("        ", k)
                    print('-----------------------')
                    print(d)
                    if k in ['oe', 'hpc', 'input']:
                        print()
                        print(d.values)
                        print()
