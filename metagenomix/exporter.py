# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import timeit
import socket
import subprocess
import datetime as dt
import sys
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
    print('\n>>> `metagenomix export` started >>>\n')
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
    print('\n<<< `metagenomix export` completed <<<\n')


class Exported(object):

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.softs = Softwares(**kwargs)
        self.time = dt.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        self.dir = abspath(self.dir)
        self.export_dir = abspath('%s/_exports' % self.dir)
        self.out = ''
        self.loc = ''
        self.provenances = []
        self.extensions = []
        self.to_exports = []
        self.commands = []

    def get_folder_exp(self):
        self.out = '%s/%s' % (self.export_dir, self.time)
        if self.location in os.environ:
            self.loc = os.environ[self.location]
        elif self.location and self.location[0] == '/':
            self.loc = self.location
        self.loc = self.loc.rstrip('/')
        if self.loc and self.loc[0] != '/':
            sys.exit('Location "%s" not an absolute path: %s' % (
                self.location, self.loc))
        if not self.local and self.loc:
            self.out = '%s/exports_%s' % (self.loc, self.time)

    def get_extensions(self):
        if self.exts:
            print('* Getting extensions to target specific files')
            self.exts = ['.%s' % x if x[0] != '.' else x for x in self.exts]
            self.exts.extend(['%s.gz' % x for x in self.exts])

    def get_inputs(self):
        print('\n* Getting files to export:')
        m = ''
        for soft in self.softs.names:
            soft_dir = self.dir + '/' + soft
            added, prov = 0, []
            for idx, (root, dirs, files) in enumerate(os.walk(soft_dir)):
                if not idx:
                    prov = ['%s/%s/provenance.txt' % (root, d) for d in dirs]
                    continue
                for fil in files:
                    pref, ext = splitext(fil)
                    if fil.endswith('gz'):
                        ext = '%s.gz' % splitext(pref)[1]
                    if self.exts:
                        if ext in self.exts:
                            m = '(with extension "%s")' % '", "'.join(self.exts)
                            self.to_exports.append('%s/%s' % (root, fil))
                            added = 1
                    else:
                        self.to_exports.append('%s/%s' % (root, fil))
                        added = 1
            if added and prov:
                self.to_exports.extend(prov)

        if self.regex:
            regexes = r'%s' % '|'.join(list(self.regex))
            rgx = re.compile(regexes, flags=re.IGNORECASE)
            self.to_exports = [x for x in self.to_exports if re.search(rgx, x)
                               or 'provenance.txt' in x]

        print('\t%s files %s will be exported' % (len(self.to_exports), m))

    # def get_inputs(self):
    #     print('\n* Getting files to export:')
    #     m = ''
    #     for root, dirs, files in os.walk(self.dir):
    #         if root == self.dir:
    #             continue
    #         soft = root.split('%s/' % self.dir)[-1].split('/')[0]
    #         if self.softs.names and soft not in self.softs.names:
    #             continue
    #         for fil in files:
    #             ext = splitext(fil)[1]
    #             if fil.endswith('gz'):
    #                 ext = '%s.gz' % splitext(splitext(fil)[0])[1]
    #             if self.exts:
    #                 if ext in self.exts:
    #                     m = '(with extensions "%s" )' % '", "'.join(self.exts)
    #                     self.to_exports.append('%s/%s' % (root, fil))
    #             else:
    #                 self.to_exports.append('%s/%s' % (root, fil))
    #
    #     if self.regex:
    #         regexes = r'%s' % '|'.join(list(self.regex))
    #         rgx = re.compile(regexes, flags=re.IGNORECASE)
    #         self.to_exports = [x for x in self.to_exports if re.search(rgx, x)]
    #
    #     print('\t%s files %swill be exported' % (len(self.to_exports), m))

    def get_output_name(self):
        if not self.output.endswith('.tar.gz'):
            if self.output.endswith('.gz'):
                self.output = '%s.tar.gz' % splitext(expandvars(self.output))[0]
            elif self.output.endswith('.tar'):
                self.output = '%s.gz' % self.output

    def get_output(self):
        print('* Getting output archive path')
        if self.output:
            self.output = abspath(self.output)
        else:
            self.output = '%s/%s.tar' % (self.dir, '_'.join(self.softs.names))
        if not self.local and self.loc:
            self.output = '%s/%s' % (self.loc, basename(self.output))
        self.get_output_name()

    def get_commands(self):
        print('* Getting copying commmands')
        self.get_copy_commands()
        print('* Getting archiving commmands')
        self.get_archiving_commands()

    def get_copy_commands(self):
        for to_export in self.to_exports:
            exp = to_export.replace(self.dir, self.out)
            if isfile(exp):
                continue
            if not isdir(dirname(exp)):
                self.commands.append('mkdir -p %s' % dirname(exp))
            self.commands.append('cp %s %s' % (to_export, exp))

    def get_archiving_commands(self):
        tar_cmd = 'tar czf %s -C %s .' % (self.output, self.out)
        rm_cmd = 'rm -rf %s' % self.out
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
            pbs = ' --torque'
        else:
            hpc = '%s.slm' % splitext(sh)[0]
            pbs = ''
        cmd = 'Xhpc -i %s -o %s -j xpt_%s -t 2%s' % (sh, hpc, self.time, pbs)
        if self.account:
            cmd += ' -a %s' % self.account
        cmd += ' --quiet'
        cmd += ' --no-abspath'
        subprocess.call(cmd.split())
        # os.remove(sh)
        return hpc

    def write_job(self):
        sh = self.write_sh()
        print('\n* To export, run one of the following:')
        if self.jobs:
            hpc = self.run_xhpc(sh)
            run = '%s_job.sh' % splitext(sh)[0]
            with open(run, 'w') as o:
                o.write('mkdir -p %s/jobs/output\n' % self.export_dir)
                o.write('cd %s/jobs/output\n' % self.export_dir)
                if self.torque:
                    o.write('squeue %s\n' % hpc)
                else:
                    o.write('sbatch %s\n' % hpc)
        print('[jobs: NIRD files unreachable] sh %s' % run)
        print('[no job: NIRD files reachable] sh %s\n' % sh)

    def showdown(self):
        user = expanduser('~').split('/')[-1]
        hostname = socket.gethostname()
        ip_address = socket.gethostbyname(hostname)
        if hostname.startswith('login') and hostname.endswith('saga'):
            hostname = 'saga.sigma2.no'
        print('* Then, copy(-edit)-paste on your computer to download from '
              'server:')
        print('scp %s@%s:%s .' % (user, ip_address, self.output))
        print(' or')
        print('scp %s@%s:%s .' % (user, hostname, self.output))

