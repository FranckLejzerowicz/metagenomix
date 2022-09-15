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

import numpy as np
import pandas as pd
from tabulate import tabulate
from os.path import abspath, basename, dirname, isdir, splitext
from metagenomix.core.output import Softwares, Output
from metagenomix._io_utils import get_input_info, get_res_info, get_size_info


def manager(**kwargs):
    """Deal with the contents of your pipeline output folder.

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

    print('* Parsing folders to collect info about jobs and outputs')
    managing.get_softs_data()

    print('* Managing')
    managing.manage()

    print('\n\n%s\n* Applying management decisions:' % ('=' * 50))
    if managing.remove or managing.jobs:
        managing.removing()
    if managing.rename:
        managing.renaming()
    if managing.store:
        managing.storing()


class Manage(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.check_disk()
        self.disk = abspath(self.disk)
        self.dir = abspath(self.dir)
        self.softwares = Softwares(**kwargs)
        self.time = dt.datetime.now().strftime("%d/%m/%Y, %H") + 'h'
        self.manage_dir = '%s/_managed' % self.dir
        self.managed = {}
        self.done = {}
        self.h = None
        self.tab = None
        self.data = None
        self.soft = None
        self.role = None
        self.name = None
        self.task = None
        self.sizes = None
        self.after = None
        self.removes = []
        self.stores = []
        self.mkdir = {}
        self.rsync = {}
        self.ln = {}
        self.rm = {}
        self.folders = {}
        self.renames = {}
        self.rename_pd = pd.DataFrame()

    def check_disk(self):
        disk_root = (self.disk, None)
        if '/' in self.disk:
            disk_root = self.disk.rsplit('/', 1)

        if self.store and not isdir(disk_root[0]):
            error = 'Error:`--store` needs a valid path'
            if disk_root[0]:
                error += ': "%s" not found' % self.disk
                if 'SLURM_JOB_ID' in os.environ or 'PBS_JOBID' in os.environ:
                    error += '\n(Check that path is accessible from a job/node)'
            sys.exit(error)

    def get_output_path(self):
        if self.out is None:
            name = self.time
        else:
            name = splitext(basename(self.out))[0]
        self.out = self.manage_dir + '/' + name

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
            sep = ('*' * len(m))
            print('\n\n\n%s\n%s\n%s' % (sep, m, sep))
            if self.jobs:
                self.task = 'jobs'
                self.manage_jobs(afters)
            if self.rename:
                self.task = 'rename'
                self.rename_folders(afters)
            if self.store:
                self.task = 'store'
                self.store_data(afters)

    def store_data(self, afters):
        for (after, h), data in afters.items():
            self.data, self.after, self.h, self.stores = data, after, h, []
            self.folders = self.data['results']
            self.print_after()
            self.manage_storage()
            if self.stores:
                self.set_stores()

    def set_stores(self):
        stores_pd = pd.DataFrame(self.stores, columns=[
            'local', 'output', 'root', 'tech_dir', 'folder', 'fil'])
        for row in stores_pd.values:
            fil = '/'.join(row)
            folder = dirname(fil)
            dest = '%s/%s' % (self.disk, '/'.join(row[1:-1]))
            out = dest + '/' + row[-1]
            if not isdir(dest):
                self.mkdir.setdefault(folder, set()).add((dest, ))
            self.rsync.setdefault(folder, set()).add((folder, dest))
            self.rm.setdefault(folder, set()).add((fil,))
            self.ln.setdefault(folder, set()).add((out, fil))

    def get_store_input(self, details=True):
        if details and len(self.sizes) > 1:
            inp = input('y/n/[d]etails: ')
        else:
            inp = input('y/[n]: ')
        return inp

    def manage_storage(self, details=True):
        inp = self.get_store_input(details)
        self.store_level(inp, details)

    def store_level(self, inp, details=True):
        if inp == 'y':
            self.store_all()
        elif inp in ['', 'd']:
            if details and len(self.sizes) > 1:
                self.store_details()
        elif inp != 'n':
            sys.exit('Error: "%s" is unknown (not in ["y", "n", "d"])' % inp)

    def store_all(self):
        local, out = self.dir.rsplit('/', 1)
        root = '%s/after_%s_%s' % (self.soft, self.after, self.h)
        for tech_dir, folders_files in self.folders.items():
            for folder, files in folders_files.items():
                for fp in files:
                    self.stores.append([local, out, root, tech_dir, folder, fp])

    def store_details(self):
        for folder in self.data['sizes']:
            print('\t* %s (%s): Store?' % (folder, self.sizes[folder]), end=' ')
            self.folders = {folder: self.data['results'][folder]}
            self.manage_storage(False)

    def manage_jobs(self, afters):
        for (after, h), data in afters.items():
            self.data, self.after, self.h = data, after, h
            self.print_after()
            self.get_oes()
            self.print_oes()
            self.manage_oes()
            self.remove_scripts()

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
            self.print_after()
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
        print('\n\n%s%s' % (tab, t.replace('\n', '\n%s' % tab)))

    def print_after(self):
        print_ = '\n    * after "%s"' % self.after
        if self.task == 'jobs':
            task = ' Management task: Jobs '
            print_ = '%s' % get_input_info(self.data['input'])
        elif self.task == 'rename':
            task = ' Management task: Rename '
            print_ = '%s' % get_res_info(self.data['results'])
        elif self.task == 'store':
            task = ' Management task: Storage '
            print_, self.sizes = get_size_info(self.data['sizes'])
            print_ = '%s: Store' % print_
            if len(self.sizes) > 1:
                print_ += ' all'
            print_ += '?'
        sep = '  %s' % ('-' * len(task))
        print('\n\n%s\n  %s\n%s\n  > %s' % (sep, task, sep, print_), end=' ')

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
        print('  - Removing')
        if self.removes:
            remove = 'y'
            if self.confirm:
                remove = input('\tRemove?%s' % self.get_paths(self.removes))
            if remove == 'y':
                for path in self.removes:
                    os.remove(path)
                print('\t-> done')
            else:
                print('\t-> aborted')
        else:
            print('\t-> nothing to remove')

    def renaming(self):
        print('  - Renaming')
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
            print('\t-> done')
        else:
            print('\t-> nothing to rename')

    def get_chunks(self):
        cmds = {}
        if self.chunks and 1 < self.chunks <= len(self.rsync):
            cmds_ = dict(enumerate(self.rsync))
            chunks = [list(x) for x in np.array_split(list(cmds_), self.chunks)]
            for cdx, c in enumerate(chunks):
                if len(c):
                    cmds['-%s' % str(cdx + 1)] = [cmds_[x] for x in c]
        elif self.rsync:
            cmds[''] = sorted(self.rsync)
        return cmds

    def storing(self):
        self.make_disk_dir()
        self.make_scripts_dir()
        scripts = self.write_scripts()
        if scripts:
            sh = self.write_screen_jobs(scripts)
            print('  - Storing')
            message = 'please run the following script to spawn screen session'
            if len(scripts) > 1:
                message += 's'
            print('\t-> %s' % message)
            print('\t   sh %s' % sh)
        else:
            print('\t-> nothing to store')

    def write_screen_jobs(self, scripts):
        sh = '%s/store.sh' % self.out
        with open(sh, 'w') as o:
            for part, script in sorted(scripts.items()):
                if part:
                    name = 'store_%s_of_%s' % (part.split()[1], part.split()[3])
                else:
                    name = 'store'
                echo = 'To monitor storage, please run: `screen -r %s`' % name
                screen = 'screen -dm -S %s /bin/bash "%s"' % (name, script)
                o.write('%s\n' % screen)
                o.write('echo "%s"\n' % echo)
            o.write('echo "`screen -ls` to list running screen session(s)"\n')
            o.write('echo "<ctrl-d> to detach when within screen session"\n')
            o.write('echo "<ctrl-k> to kill a screen session from within"\n')
        return sh

    def write_scripts(self):
        scripts = {}
        chunks = self.get_chunks()
        for chunk, keys in chunks.items():
            part = ''
            if chunk:
                part += ' [ %s / %s ]' % (chunk[1:], len(chunks))
            sh = '%s/scripts/store%s.sh' % (self.out, chunk)
            scripts[part] = sh
            with open(sh, 'w') as o:
                message = 'Storing pipeline data to: %s%s' % (self.disk, part)
                o.write('echo "%s"\n' % message)
                for key in keys:
                    o.write('# Storing folder: "%s"\n' % key)
                    for (step, args, j) in [
                        ('mkdir', 'p', ' '),
                        ('rsync', 'aurqv', '/ '),
                        ('rm', 'rf', ' '),
                        ('ln', 's', ' '),
                    ]:
                        for paths in self.__dict__[step].get(key, []):
                            o.write('%s -%s %s\n' % (step, args, j.join(paths)))
                o.write('echo "done"\n')
        return scripts

    def make_disk_dir(self):
        disk = abspath(self.disk)
        if not isdir(disk):
            os.makedirs(disk)

    def make_scripts_dir(self):
        folder = '%s/scripts' % self.out
        if not isdir(folder):
            os.makedirs(folder)

