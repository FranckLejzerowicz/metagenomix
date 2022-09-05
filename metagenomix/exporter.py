# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import socket
import subprocess
import datetime as dt
from os.path import (
    abspath, basename, dirname, expanduser, expandvars, isdir, isfile, splitext)

from metagenomix.core.output import Softwares


def exporter(**kwargs):
    """Prepare an archive for specific pipeline outputs.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix exporter ===\n')
    kwargs['command'] = 'export'
    exporting = Exported(**kwargs)
    exporting.get_folder_exp()
    exporting.get_extensions()
    exporting.get_inputs()
    exporting.get_output()
    print('\t>', exporting.out)
    exporting.get_commands()
    exporting.write_job()
    exporting.showdown()


class Exported(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.softs = Softwares(**kwargs)
        self.time = dt.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        self.export_dir = abspath('%s/_exports' % self.folder)
        self.dir = ''
        self.extensions = []
        self.to_exports = []
        self.commands = []

    def get_folder_exp(self):
        self.dir = '%s/%s' % (self.export_dir, self.time)
        if not self.local and self.location in os.environ:
            self.dir = '%s/exports_%s' % (os.environ[self.location], self.time)

    def get_extensions(self):
        if self.exts:
            print('* Getting extensions to target specific files')
            self.exts = ['.%s' % x if x[0] != '.' else x for x in self.exts]

    def get_inputs(self):
        print('\n* Getting files to export:')
        if self.regex:
            regex = re.compile(
                r'%s' % '|'.join(list(self.regex)), flags=re.IGNORECASE)

        m = ''
        for root, dirs, files in os.walk(self.folder):
            if root == self.folder:
                continue
            soft = root.split('%s/' % self.folder)[-1].split('/')[0]
            if self.softs.names and soft not in self.softs.names:
                continue

            for fil in files:
                if self.regex:
                    if re.search(regex, fil):
                        self.to_exports.append('%s/%s' % (root, fil))
                        continue
                ext = splitext(fil)[1]
                if self.exts and ext in self.exts:
                    m = '(with extensions "%s" )' % '", "'.join(self.exts)
                    self.to_exports.append('%s/%s' % (root, fil))
                else:
                    self.to_exports.append('%s/%s' % (root, fil))

        print('\t%s files %swill be exported' % (len(self.to_exports), m))

    def get_output(self):
        print('* Getting output archive path')
        self.out = abspath(self.out)
        if not self.local and self.location in os.environ:
            self.out = '%s/%s' % (os.environ[self.location], basename(self.out))
        self.out = '%s.tar.gz' % splitext(expandvars(self.out))[0]

    def get_commands(self):
        print('* Getting copying commmands')
        self.get_copy_commands()
        print('* Getting archiving commmands')
        self.get_archiving_commands()

    def get_copy_commands(self):
        for to_export in self.to_exports:
            exp = to_export.replace(self.folder.rstrip('/'), self.dir)
            if isfile(exp):
                continue
            if not isdir(dirname(exp)):
                self.commands.append('mkdir -p %s' % dirname(exp))
            self.commands.append('cp %s %s' % (to_export, exp))

    def get_archiving_commands(self):
        tar_cmd = 'tar czf %s -C %s .' % (self.out, self.dir)
        rm_cmd = 'rm -rf %s' % self.dir
        self.commands.append(tar_cmd)
        self.commands.append(rm_cmd)

    def write_sh(self):
        job_sh = '%s/jobs/export_%s.sh' % (self.export_dir, self.time)
        if not isdir(dirname(job_sh)):
            os.makedirs(dirname(job_sh))
        with open(job_sh, 'w') as o:
            for command in self.commands:
                o.write('%s\n' % command)
        return job_sh

    def run_xhpc(self, sh):
        if self.torque:
            hpc = '%s.pbs' % splitext(sh)[0]
        else:
            hpc = '%s.slm' % splitext(sh)[0]
        cmd = 'Xhpc -i %s -o %s -j xpt_%s -t 2 --no-stat' % (sh, hpc, self.time)
        if self.account:
            cmd += ' -a %s' % self.account
        cmd += ' --quiet'
        subprocess.call(cmd.split())
        os.remove(sh)
        return hpc

    def write_job(self):
        sh = self.write_sh()
        if self.jobs:
            hpc = self.run_xhpc(sh)
            with open(sh, 'w') as o:
                o.write('mkdir %s/jobs/output\n' % self.export_dir)
                o.write('cd %s/jobs/output\n' % self.export_dir)
                if self.torque:
                    o.write('squeue %s\n' % hpc)
                else:
                    o.write('sbatch %s\n' % hpc)
        print('\n* To export, run the following:\nsh %s\n' % sh)

    def showdown(self):
        user = expanduser('~').split('/')[-1]
        hostname = socket.gethostname()
        print('* Then, copy(-edit)-paste to download from server to local:')
        print('scp %s@%s:%s .' % (user, hostname, self.out))
