# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import datetime as dt

import pandas as pd
from tabulate import tabulate
from os.path import abspath, basename, dirname
from metagenomix.core.output import Softwares, Output
from metagenomix._io_utils import get_input_info, get_folder_info, get_size_info


def manager(**kwargs):
    """Run metagenomix to assist user renaming or discarding outputs.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix manager ===\n')

    tasks = ['jobs', 'remove', 'rename', 'store']
    if not sum(kwargs[task] for task in tasks):
        sys.exit('At least 1 management task needed (--%s)' % ',--'.join(tasks))

    kwargs['command'] = 'manage'
    managing = Manage(**kwargs)

    print('* Setting up the output file name')
    managing.get_output_path()
    print(' -> %s' % managing.out)

    print('* Parsing folders to collect info about jobs and outputs')
    managing.get_softs_data()

    print('* Managing')
    managing.manage()

    print('* Applying management decisions:')
    if managing.remove or managing.jobs:
        print(' - Removing')
        managing.removing()

    if managing.rename:
        print(' - Renaming')
        managing.renaming()

    if managing.store:
        print(' - Storing')
        managing.storing()


class Manage(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.dir = abspath(self.dir)
        self.softwares = Softwares(**kwargs)
        self.time = dt.datetime.now().strftime("%d/%m/%Y, %H") + 'h'
        self.manage_dir = '%s/_managament' % self.dir
        self.managed = {}
        self.done = {}
        self.h = None
        self.tab = None
        self.data = None
        self.soft = None
        self.role = None
        self.name = None
        self.after = None
        self.removes = []
        self.renames = {}
        self.rename_pd = pd.DataFrame()

    def get_output_path(self):
        if self.out is None:
            path = self.time + '.txt'
        else:
            path = basename(self.out)
            if '/' in self.out:
                print('Using "%s" to write in "%s"' % (path, self.manage_dir))
        self.out = self.manage_dir + '/' + path

    def get_softs_data(self):
        """An Output class instance is created for each software to manage,
        and placed as value to the dict with the software name of key."""
        for role, softs in self.softwares.softs.items():
            self.role = role
            for soft in softs:
                output = Output(self.dir, soft)
                output.get_afters()
                output.get_outputs()
                self.managed[soft] = output.outputs

    def manage(self):
        for soft, afters in self.managed.items():
            self.soft = soft
            term = 'after %s software' % len(afters)
            if len(afters) > 1:
                term += 's'
            m = '[%s] %s (%s)' % (soft, term, '; '.join([x[0] for x in afters]))
            print('\n%s\n%s\n%s' % (('=' * len(m)), m, ('=' * len(m))))
            # if self.jobs:
            #     self.manage_jobs(afters)
            # if self.rename:
            #     self.rename_folders(afters)
            if self.store:
                self.store_data(afters)

    def store_data(self, afters):
        for (after, h), data in afters.items():
            self.data, self.after, self.h = data, after, h
            self.print_after('sizes')

    def manage_jobs(self, afters):
        for (after, h), data in afters.items():
            self.data, self.after, self.h = data, after, h
            self.print_after('input')
            # self.remove_scripts()
            self.get_oes()
            self.print_oes()
            self.manage_oes()

    def get_renaming_table(self):
        rename_pd = []
        for ix, (t_dir, ffs) in enumerate(self.data['results'].items()):
            idx = str(ix)
            rename_pd.append([idx, 'folder', '.', t_dir])
            for jx, (ff, files) in enumerate(ffs.items()):
                idx = '%s.%s' % (ix, jx)
                if '/' in ff:
                    rename_pd.append([
                        idx, 'folder', t_dir + '/' + dirname(ff), basename(ff)])
                else:
                    rename_pd.append([idx, 'folder', t_dir, ff])
                for fx, fil in enumerate(files):
                    idx = '%s.%s.%s' % (ix, jx, fx)
                    rename_pd.append([idx, 'file', t_dir + '/' + ff, fil])
        path = '%s/after_%s_%s/' % (self.soft, self.after, self.h)
        self.rename_pd = pd.DataFrame(rename_pd, columns=[
            'index', 'dir/file', 'path (in:)\n%s' % path, 'name (to rename)'])

    def rename_folders(self, afters):
        for (after, h), data in afters.items():
            self.data, self.after, self.h = data, after, h
            if not self.data['results']:
                continue
            self.print_after('results')
            self.get_renaming_table()
            self.pretty_print(self.rename_pd)
            self.get_renaming()

    def get_renaming(self):
        index = ''
        text = '\tindex  : '
        indices = set(self.rename_pd['index'])
        print('\n\tEnter index and text to rename these paths (Enter to exit):')
        while 1:
            inp = str(input(text))
            if text == '\tindex  : ':
                if inp and inp not in indices:
                    print('\t -> "%s" is not an index above' % inp)
                    continue
                if index:
                    self.collect_rename(index, rename)
                text = '\trename : '
                index = inp
            else:
                text = '\tindex  : '
                rename = inp
            if not inp:
                break

    def collect_rename(self, index, rename):
        path = '%s/%s/after_%s_%s' % (self.dir, self.soft,
                                      self.after, self.h)
        if path not in self.renames:
            self.renames[path] = {}
        row = self.rename_pd.loc[self.rename_pd['index'] == index].values[0]
        self.renames[path].setdefault(len(index), []).append(
            (row[-2], row[-1], rename))

    def get_tables(self):
        for tab in ['oe', 'hpc', 'input']:
            d = {tab: self.data[tab].loc[self.data[tab]['name'] == self.name]}
            self.__dict__.update(d)

    def print_name(self, idx):
        self.tab = self.input.iloc[:, :-1]
        tab = self.tab.dropna(axis=1)
        unit = 'input unit'
        if tab.shape[0] > 1:
            unit += 's'
        print('\n\t[%s] %s %s for "%s"' % (idx, tab.shape[0], unit, self.name))
        done = self.done.get(self.name)
        if done['N'] or done['Y'] and self.remove:
            self.pretty_print(tab, '\t    ', 'presto')

    def manage_oe(self, idx):
        done = self.done.get(self.name)
        if done:
            self.print_name(idx)
            self.remove_jobs('N')
            if self.remove:
                self.remove_jobs('Y')

    def remove_jobs(self, y_n):
        ids = self.done[self.name][y_n]
        if ids:
            fps = [self.data['stdouts'][(self.name, x)] for x in ids]
            n = len(fps)
            fps = self.get_paths(fps, True)
            if y_n == 'Y':
                status = 'ERRONEOUS'
            else:
                status = 'COMPLETED'
            if input('\n\t     Remove %s %s jobs?%s' % (n, status, fps)) == 'y':
                self.removes.extend(fps)

    def manage_oes(self):
        for idx, name in enumerate(self.data['input']['name'].unique()):
            self.name = name
            self.get_tables()
            self.manage_oe(idx)

    def get_oes(self):
        self.done = {}
        if len(self.data['oe']):
            oe = self.data['oe'].iloc[:, :3].pivot_table(
                index=['name'], columns=['done'], aggfunc=list)
            oe.columns = oe.columns.droplevel()
            self.done = oe.fillna('').T.to_dict()

    def print_oes(self):
        data = self.data['input'].dropna(axis=1).copy()
        data['done'] = ['Y' if self.done and self.done[x]['Y']
                        else 'N' for x in data['name']]
        data['job name'] = data['name']
        data['job IDs'] = [';'.join(self.done[x]['Y']) if self.done[x]['Y'] else
                           ';'.join(self.done[x]['N']) for x in data['name']]
        data = data.drop(columns=['name'])
        self.pretty_print(data)

    @staticmethod
    def pretty_print(data, tab='\t', fmt='pretty', showindex=False):
        t = tabulate(data, headers='keys', tablefmt=fmt,
                     showindex=showindex, stralign="left")
        print('\n%s%s' % (tab, t.replace('\n', '\n%s' % tab)))

    def print_after(self, data):
        print_ = '\n    * after "%s"' % self.after
        dat = self.data[data]
        print(dat)
        if data == 'input':
            print_ = '\n  [jobs] * %s' % get_input_info(dat)
        elif data == 'results':
            print_ = '\n  [results] * %s' % get_folder_info(dat)
        elif data == 'sizes':
            print_ = '\n  [sizes] * %s' % get_size_info(dat)
        print(print_)

    def get_paths(self, paths, jobs=False):
        if jobs:
            paths = ['\t\t- %s{o,e}' % x.split(self.h)[-1][1:-1] for x in paths]
        else:
            paths = ['\t\t- %s' % x.split(self.h)[-1][1:] for x in paths]
        paths = '\n  %s\n\ty/[n] ' % '\n'.join(paths)
        return paths

    def remove_scripts(self):
        if self.remove:
            all_scripts = self.data['all_scripts']
            run_scripts = self.data['run_scripts']
            scripts = all_scripts.difference(run_scripts)
            if scripts:
                fps = self.get_paths(scripts)
                if input('\n\tRemove previous scripts?%s\n' % fps) == 'y':
                    self.removes.extend(scripts)

    def removing(self):
        if self.removes:
            remove = 'y'
            if self.confirm:
                remove = input('\tRemove?%s' % self.get_paths(self.removes))
            if remove == 'y':
                for path in self.removes:
                    os.remove(path)
            else:
                print('\t-> Removal aborted.')

    def renaming(self):
        renames = []
        for after_dir, index_renames in self.renames.items():
            soft, after = after_dir.rsplit('/', 2)[1:]
            if self.confirm:
                print('\t[%s]' % soft)
                shows = []
            for index in sorted(index_renames.keys())[::-1]:
                for (folder, name, rename) in index_renames[index]:
                    renames.append([after_dir + '/' + folder + '/' + name,
                                    after_dir + '/' + folder + '/' + rename])
                    if self.confirm:
                        if folder == '.':
                            cur = '%s/%s/%s' % (soft, after, name)
                        else:
                            cur = '%s/%s/%s/%s' % (soft, after, folder, name)
                        shows.append([cur, rename])

            if self.confirm:
                t = tabulate(shows, headers='', tablefmt='presto',
                             showindex=False).replace('|', '->')
                print('\t - %s' % t.replace('\n', '\n\t - '))

        if renames and self.confirm:
            if input('\ty/[n] ') == 'y':
                for (name, rename) in renames:
                    os.rename(name, rename)

    def storing(self):
        pass
