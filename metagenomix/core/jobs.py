# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import yaml
import glob
import shutil
import random
import subprocess
import pkg_resources
import numpy as np
import datetime as dt
from tabulate import tabulate
from os.path import dirname, exists, isdir, islink, splitext

from metagenomix._io_utils import (
    mkdr, get_roundtrip, print_status_table, compute_hash, get_md5, get_dates)
from metagenomix.core.slurm import (
    set_directives, set_preamble, set_scratching, set_tmpdir)

STORAGE = pkg_resources.resource_filename("metagenomix", "storage")
TESTS = pkg_resources.resource_filename("metagenomix", "tests")


class Created(object):

    def __init__(self, config, databases, workflow, commands):
        self.config = config
        self.databases = databases
        self.commands = commands
        self.graph = workflow.graph
        self.cmd = []
        self.cmds = {}
        self.links = {}
        self.soft_links = {}
        self.links_stats = {}
        self.avail = {}
        self.chunks = {}
        self.soft = ''
        self.hash = ''
        self.sh = ''
        self.main_sh = ''
        self.move_sh = ''
        self.move_shs = {}
        self.run = {'database': {}, 'software': {}}
        self.job_fps = []
        self.job_name = ''
        self.modules = {}
        self.pjct = self.get_prjct()
        self.scheduler = self.get_scheduler()
        self.time = dt.datetime.now().strftime("%d/%m/%Y-%H:%M")
        self.scripts = []
        self.scripts_parts = []
        self.log_dir = '%s/_created' % config.dir
        self.link_script = None

    def make_dirs(self):
        for name, hashes_softs in self.commands.softs.items():
            for hash, soft in hashes_softs.items():
                for directory in sorted(soft.dirs):
                    if not islink(directory) and not isdir(directory):
                        os.makedirs(directory)

    def get_scheduler(self):
        if self.config.jobs:
            if self.config.torque:
                return 'qsub'
            else:
                return 'sbatch'
        else:
            return 'sh'

    def get_prjct(self):
        prjct = [x for x in self.config.project if x.lower() not in 'aeiouy']
        if prjct:
            return ''.join(prjct)
        else:
            return self.config.project

    def get_chunks(self, n_chunks):
        self.chunks = {}
        if n_chunks and len(self.cmds) >= n_chunks:
            cmds_d = dict(enumerate(self.cmds))
            chunks = [list(x) for x in np.array_split(list(cmds_d), n_chunks)]
            for cdx, chunk in enumerate(chunks):
                if len(chunk):
                    self.chunks[str(cdx)] = [cmds_d[x] for x in chunk]
        else:
            for key in self.cmds:
                chunkey = '%s_%s' % (key[0], '-'.join(key[1]).replace('/', '_'))
                self.chunks[chunkey] = [key]

    def get_main_sh(self, name, hashed='', soft=None, local='') -> None:
        if soft:
            n = 'software'
            key = '%s: %s%s' % (name, hashed, local)
            if self.config.show_paths:
                key += '\n%s' % " -> ".join(soft.path)
            if self.cmds or local:
                self.main_sh = '%s/run%s.sh' % (soft.dir, local)
                self.run[n][key] = self.main_sh
            else:
                self.main_sh = ''
                self.run[n][key] = 'Input needed (output from %s)' % soft.prev
        else:
            self.main_sh = '%s/run_%s.sh' % (self.sh.rsplit('/', 2)[0], name)
            self.run['database'][name] = self.main_sh

    def write_main(self) -> None:
        with open(self.main_sh, 'w') as o:
            if len(self.job_fps):
                job_dir = dirname(self.job_fps[0])
                o.write('mkdir -p %s/output\n' % job_dir)
                o.write('cd %s/output\n' % job_dir)
                for jdx, job_fp in enumerate(self.job_fps):
                    o.write('%s %s\n' % (self.scheduler, job_fp))
                self.job_fps = []

    def write_bash(self, name, hashed, soft) -> None:
        if soft.bash:
            self.get_main_sh(name, hashed, soft, '_locally')
            if self.main_sh:
                with open(self.main_sh, 'w') as o:
                    for cmd in soft.bash:
                        o.write('%s\n' % cmd)

    def get_provenance_fp(self, name, soft=None) -> str:
        if soft:
            fp = '%s/provenance.txt' % soft.dir
        else:
            fp = '%s/provenance_%s.txt' % (self.sh.rsplit('/', 2)[0], name)
        return fp

    def write_provenance(self, name, soft) -> None:
        step_max = max([len(x) for x in soft.path])
        provenance_fp = self.get_provenance_fp(name, soft)
        with open(provenance_fp, 'w') as o:
            o.write("Pipeline steps to this output (and analysis type):\n")
            for sdx, step in enumerate(soft.path):
                o.write('%s\t%s%s\t: %s\n' % (
                    sdx, step, (' ' * (step_max-len(step))),
                    self.config.tools.get(step, '')))
            o.write("\nParameters for these steps' outputs:\n")
            for sdx, step in enumerate(soft.path):
                if not sdx:
                    continue
                o.write('\n========== %s: %s ==========\n' % (sdx, step))
                params_show = self.commands.softs[step][self.commands.hashes[
                    tuple(soft.path[:(sdx + 1)])]].params.copy()
                if step.startswith('search'):
                    databases = params_show['databases']
                    del params_show['databases']
                    params = {'search': {soft.name.split('_')[-1]: params_show,
                                         'databases': databases}}
                    o.write('%s\n' % yaml.dump(params))
                else:
                    o.write('%s\n' % yaml.dump(params_show))

    def database_cmds(self):
        for db, cmds in self.databases.commands.items():
            self.cmds = cmds
            self.get_chunks(self.config.params['chunks'])
            self.write_jobs(db)
            self.write_main(db)

    def get_modules(self, name: str):
        self.modules = set(self.config.modules.get(name.split('_')[0], set()))
        self.modules.update(set(self.config.modules.get(name, set())))

    def scratch(self, soft, key, cmds):
        if soft.params['scratch'] and self.config.jobs and key in soft.io:
            roundtrip = get_roundtrip(soft.io[key])
            extended_cmds = ['\n# Move to SCRATCH_FOLDER']
            extended_cmds += roundtrip['to']
            extended_cmds += ['\n# data: %s %s' % (key[0], ' '.join(key[1]))]
            extended_cmds += cmds
            if self.config.move_back:
                extended_cmds += ['\n# Move from SCRATCH_FOLDER']
                extended_cmds += roundtrip['from']
        else:
            extended_cmds = ['\n# data: %s %s' % (key[0], ' '.join(key[1]))]
            extended_cmds += cmds
        self.cmds[key] = extended_cmds
        self.dummy_outputs(soft, key)

    def dummy_outputs(self, soft, key):
        if self.config.dev:

            for fil_ in soft.io[key].get(('O', 'f'), []):
                if fil_.endswith('*'):
                    continue
                fil_loc = fil_.replace('${SCRATCH_FOLDER}', '')
                fil_sto = fil_loc.replace(dirname(self.config.dir), STORAGE)
                if exists(fil_loc):
                    os.remove(fil_loc)
                if exists(fil_sto):
                    os.remove(fil_sto)

            for fol_ in soft.io[key].get(('O', 'd'), []):
                if fol_.endswith('*'):
                    continue
                fol_loc = fol_.replace('${SCRATCH_FOLDER}', '')
                fol_sto = fol_loc.replace(dirname(self.config.dir), STORAGE)
                if islink(fol_loc):
                    os.unlink(fol_loc)
                elif isdir(fol_loc):
                    shutil.rmtree(fol_loc)
                if isdir(fol_sto):
                    shutil.rmtree(fol_sto)

            if random.choice([0, 1]):

                for fil_ in soft.io[key].get(('O', 'f'), []):
                    if fil_.endswith('*'):
                        continue
                    fil_loc = fil_.replace('${SCRATCH_FOLDER}', '')
                    fil_sto = fil_loc.replace(dirname(self.config.dir), STORAGE)
                    if not isdir(dirname(fil_loc)):
                        os.makedirs(dirname(fil_loc))
                    if not isdir(dirname(fil_sto)):
                        os.makedirs(dirname(fil_sto))
                    with open(fil_sto, 'w') as o:
                        for r in range(10, random.randrange(11, 100)):
                            o.write('test\n')
                    if not exists(fil_loc):
                        os.symlink(fil_sto, fil_loc)

                for fol_ in soft.io[key].get(('O', 'd'), []):
                    if fol_.endswith('*'):
                        continue
                    fol_loc = fol_.replace('${SCRATCH_FOLDER}', '')
                    fol_sto = fol_loc.replace(dirname(self.config.dir), STORAGE)
                    if not isdir(fol_sto):
                        os.makedirs(fol_sto)
                    if not isdir(dirname(fol_loc)):
                        os.makedirs(dirname(fol_loc))
                    if not isdir(fol_loc):
                        os.symlink(fol_sto, fol_loc)

    def get_cmds(self, soft):
        self.cmds = {}
        for sam_or_pool, cmds in soft.cmds.items():
            if sum(map(bool, cmds)):
                self.scratch(soft, sam_or_pool, cmds)
                self.avail.setdefault(soft.name, []).append(sam_or_pool)

    def print_status(self, m, n, name, h, soft):
        gap = (m - len(name) - len(str(n))) + 1
        print('\t%s [%s: %s] %s%s' % (n, name, h, ('.'*gap), ('.'*8)), end=' ')
        print_status_table(soft, self.config.show_status)
        return 1

    def write_logs(self):
        if not isdir(self.log_dir):
            os.makedirs(self.log_dir)
        fp = '%s/%s.txt' % (self.log_dir, compute_hash(self.hash))
        logs = self.get_logs(get_dates(fp))
        with open(fp, 'w') as o:
            for log in logs:
                o.write('%s\n' % log)
            if self.link_script:
                o.write(self.link_script)
        print('\nLog file: %s' % fp)

    def get_logs(self, dates):
        logs = self.scripts
        logs.extend(self.get_date_logs(dates))
        logs.extend(self.get_input_logs())
        return logs

    def get_date_logs(self, dates) -> list:
        date_logs = ['\n------------------']
        date_logs.extend(['Date: %s' % date for date in dates])
        date_logs.append('Date: %s' % self.time)
        date_logs.append('------------------')
        return date_logs

    def get_input_logs(self):
        input_logs = ['\n* Input files:']
        for key, arg in [
            ('illumina_dirs', '--fastq-dir-illumina'),
            ('pacbio_dirs', '--fastq-dir-pacbio'),
            ('nanopore_dirs', '--fastq-dir-nanopore'),
            ('output_dir', '--output-dir'),
            ('meta_fp', '--metadata'),
            ('pipeline_tsv', '--pipeline'),
            ('databases_yml', '--databases'),
            ('user_params_yml', '--user-params'),
            ('coassembly_yml', '--co-assembly'),
            ('strains_yml', '--strains')
        ]:
            vals = self.config.__dict__[key]
            if not vals:
                continue
            if not isinstance(vals, tuple):
                vals = (vals,)
            for val in vals:
                if key == 'meta_fp':
                    input_logs.append('* Input configuration files:')
                input_logs.append('%s%s%s\t%s' % (arg, ' ' * (25-len(arg)),
                                                  get_md5(val), val))
        return input_logs

    def get_links_chunks(self):
        cmds = {}
        if self.config.links_chunks:
            n_chunks = self.config.links_chunks
            links = set([y for x in self.soft_links.values() for y in x])
            if n_chunks and 1 < n_chunks <= len(links):
                cmds_ = dict(enumerate(list(links)))
                chunk = [list(x) for x in np.array_split(list(cmds_), n_chunks)]
                for cdx, c in enumerate(chunk):
                    if len(c):
                        cmds['-%s' % str(cdx + 1)] = [cmds_[x] for x in c]
            elif links:
                cmds[''] = list(links)
        else:
            for name, links in self.soft_links.items():
                if links:
                    cmds['-%s' % name] = links
        return cmds

    def get_links_dir(self):
        links_dir = '%s/%s' % (self.log_dir, compute_hash(self.hash))
        if not isdir(links_dir):
            os.makedirs(links_dir)
        if not isdir('%s/scripts' % links_dir):
            os.makedirs('%s/scripts' % links_dir)
        return links_dir

    def write_screen_jobs(self, links_dir, scripts):
        self.move_sh = '%s/move.sh' % links_dir
        with open(self.move_sh, 'w') as o:
            for name, (screen_sh, out) in sorted(scripts.items()):
                o.write('sh %s\n' % screen_sh)
            o.write('screen -ls\n')
            o.write('echo "To list running screen session(s): screen -ls"\n')
            o.write('echo "To get into a screen session: screen -r <ID>"\n')
            o.write('echo "To detach when within screen session: <ctrl-d>"\n')
            o.write('echo "To kill a screen session from within: <ctrl-k>"\n')

    def write_links(self):
        links_dir = self.get_links_dir()
        scripts = self.get_bring_links_scripts(links_dir)
        if scripts:
            self.write_screen_jobs(links_dir, scripts)
            self.print_links(scripts)

    def print_links(self, scripts):
        s = '\n%s\n' % ('*' * 46)
        m = '%sSome data needed for these jobs is stored away%s' % (s, s)
        tab = [x for x in self.links_stats.items() if x[1]]
        tab = [tuple(list(x) + self.move_shs.get(x[0], [''])) for x in tab]
        t = tabulate(tab, headers='', tablefmt='presto', showindex=False)
        m += '\t- %s' % t.replace(
            '\n', '\n\t- '
        ).replace(
            '|', ' :'
        ).replace(
            ' : sh', ' files  -> sh'
        ).replace(
            ' :\n', ' files\n'
        )
        m += '\n\t-> all softwares '
        n = len(scripts)
        m += '(%s screen session' % n
        if n > 1:
            m += 's'
        m += ')  -> sh %s\n' % self.move_sh
        self.link_script = m
        print(m)

    def get_bring_links_scripts(self, links_dir):
        scripts = {}
        chunks = self.get_links_chunks()
        for cdx, (chunk, links) in enumerate(chunks.items()):

            base = '%s/scripts/move%s' % (links_dir, chunk)
            out = '%s_out.txt' % base
            screen_sh = '%s_screen.sh' % base
            sh = '%s.sh' % base

            part, name = '', 'move'
            if chunk[1:].isdigit():
                part += ' [ %s / %s ]' % (chunk[1:], len(chunks))
                name += '_%s_of_%s' % (chunk[1:], len(chunks))
            else:
                part += ' [ %s ]' % chunk[1:]
                name += '_%s' % chunk[1:]
                self.move_shs[chunk[1:]] = ['sh %s' % screen_sh]

            scripts[name] = (screen_sh, out)
            with open(screen_sh, 'w') as screen_o:
                echo = 'Started screen in detached mode (ID: %s)' % name
                echo += '\nCheck whether some moves went wrong: %s' % out
                screen_o.write('screen -dmS %s /bin/bash "%s"\n' % (name, sh))
                screen_o.write('echo "%s"\n' % echo)

            with open(sh, 'w') as o:
                message = 'Bringing data from storage%s' % part
                o.write('touch %s\n' % out)
                o.write('echo "%s"\n' % message)
                o.write('sleep %s\n' % (cdx + 2))
                for link_ in links:
                    dest = link_.replace('${SCRATCH_FOLDER}', '')
                    src = self.links[dest]
                    o.write('if [ -L %s ]; then\n' % dest)
                    # o.write("    m0=`md5sum %s | cut -d' ' -f 1`\n" % src)
                    o.write('    rm -rf %s\n' % dest)
                    o.write('    rsync %s %s\n' % (src, dest))
                    # o.write("    m1=`md5sum %s | cut -d' ' -f 1`\n" % dest)
                    # o.write('    if [ "$m0" != "$m1" ]; then echo ')
                    # o.write('"$m0 $m1 %s %s" >> %s; fi\n' % (src, dest, out))
                    o.write('fi\n\n')
                o.write('echo "done"\n')

        return scripts

    def software_cmds(self):
        m = max(len(x) for x in self.commands.softs) + 1
        n = 0
        for name, hashes_softs in self.commands.softs.items():
            for hashed, soft in hashes_softs.items():
                n += self.print_status(m, n, name, hashed, soft)
                self.write_bash(name, hashed, soft)
                if not len(soft.cmds):
                    continue
                self.hash += str(soft.hash)
                self.get_links(soft)
                self.get_modules(name)
                self.get_cmds(soft)
                self.get_main_sh(name, hashed, soft)
                if self.cmds:
                    self.get_chunks(soft.params['chunks'])
                    self.write_jobs(name, soft)
                    self.write_main()
                self.write_provenance(name, soft)

    def get_links(self, soft):
        self.links_stats[soft.name] = len(soft.links)
        self.soft_links[soft.name] = soft.links
        self.links.update(soft.links)

    def get_sh(self, name: str, chunk: str, soft=None) -> None:
        """

        Parameters
        ----------
        name : str
            Name of the current software of the pipeline workflow
        chunk : str
            Name of the current chunk of commands
        soft
        """
        if soft:
            self.sh = '%s/jobs/run_%s.sh' % (
                soft.dir, chunk.replace('/', '_'))
        else:
            self.sh = '%s/databases/jobs/build_%s_%s.sh' % (
                self.config.dir, name, chunk)
        mkdr(dirname(self.sh))

    def prep_script(self, params: dict) -> None:
        job_script = '%s.slm' % splitext(self.sh)[0]
        if self.config.torque:
            job_script = '%s.pbs' % splitext(self.sh)[0]
        self.job_fps.append(job_script)
        with open(job_script, 'w') as o:
            o.write('\n'.join(set_directives(self, params)))
            o.write('\n# ------ directives END ------\n\n')
            o.write('\n'.join(set_preamble(self, params, job_script)))
            o.write('\n# ------ preamble END ------\n\n')
            o.write('\n'.join(set_scratching(self, params)))
            o.write('\n# ------ scratching END ------\n\n')
            if 'TMPDIR' in os.environ:
                o.write('\n'.join(set_tmpdir(self)))
                o.write('\n# ------ tmpdir END ------\n\n')
            with open(self.sh) as f:
                for line in f:
                    o.write(line)
            o.write('\n# ------ commands END ------\n\n')
            if 'TMPDIR' in os.environ:
                o.write('\n\nrm -rf ${TMPDIR}\n')
                o.write('# ------ clear END ------\n\n')
            o.write('\necho "Done!"\n')
        os.remove(self.sh)


    # def prep_script(self, params: dict) -> None:
    #     # mandatory options
    #     self.cmd = [
    #         'Xhpc',
    #         '-i', self.sh,
    #         '-j', self.job_name,
    #         '-t', str(params['time']),
    #         '-c', str(params['cpus']),
    #         '-M', str(params['mem']), params['mem_dim'],
    #         '--no-stat', '--no-abspath']
    #     # whether the cpus requests is per node
    #     if params['nodes']:
    #         self.cmd.extend(['-n', str(params['nodes'])])
    #     # always provide an account
    #     if self.config.account:
    #         self.cmd.extend(['-a', self.config.account])
    #     # get the job script file path and use it
    #     job_script = '%s.slm' % splitext(self.sh)[0]
    #     if self.config.torque:
    #         self.cmd.append('--torque')
    #         job_script = '%s.pbs' % splitext(self.sh)[0]
    #     self.job_fps.append(job_script)
    #     self.cmd.extend(['-o', job_script])
    #     # machine-specific setup: env activating and slurm partition
    #     if params['machine']:
    #         self.cmd.append('--%s' % params['machine'])
    #     if params['partition']:
    #         self.cmd.append('--p-partition %s' % params['partition'])
    #     # whether an environment must be used for the current software
    #     if params['env']:
    #         self.cmd.extend(['-e', params['env']])
    #     # setup the scratch location to be used for the current software
    #     if isinstance(params['scratch'], int):
    #         self.cmd.extend(['--localscratch', str(params['scratch'])])
    #     elif params['scratch'] == 'scratch':
    #         self.cmd.append('--scratch')
    #     elif params['scratch'] == 'userscratch':
    #         self.cmd.append('--userscratch')
    #     if not self.config.verbose:
    #         self.cmd.append('--quiet')
    #
    # def call_cmd(self):
    #     cmd = ' '.join(self.cmd)
    #     subprocess.call(cmd.split())
    #     os.remove(self.sh)

    def write_chunks(self, chunk_keys: list, soft):
        with open(self.sh, 'w') as sh:
            if self.modules:
                sh.write('module purge\n')
            for module in self.modules:
                sh.write('module load %s\n' % module)
            if self.config.cleanup:
                cleanup = 'cleanup rm -rf ${TMPDIR}'
                if soft.params['scratch'] and self.config.jobs:
                    cleanup += ' ${SCRATCH_FOLDER}/*'
                sh.write('%s\n' % cleanup)
            for chunk_key in chunk_keys:
                for cmd in self.cmds[chunk_key]:
                    if soft.params['scratch'] and self.config.jobs:
                        sh.write('%s\n' % cmd)
                    else:
                        sh.write('%s\n' % cmd.replace('${SCRATCH_FOLDER}', ''))

    def get_job_name(self, name: str, chunk: str, soft=None):
        self.job_name = name + '.' + self.pjct
        if soft:
            self.job_name += '.' + ''.join([
                x for x in soft.dir.split('after_')[-1] if x not in 'aeiuoy'])
        self.job_name += '.' + chunk.replace('/', '_')

    def write_dummy_oe(self):
        outro = open('%s/outro.o' % TESTS).readlines()
        outro_err = open('%s/outro_err.o' % TESTS).readlines()
        job_dir = '%s/output' % dirname(self.sh)
        if not isdir(job_dir):
            os.makedirs(job_dir)
        for job_fp in glob.glob('%s/#%s*' % (job_dir, self.job_name)):
            os.remove(job_fp)
        jobs = range(1, random.randrange(2, 4))
        for job in jobs:
            job_fp = 'slurm-%s_000000%s' % (self.job_name, job)
            lines = ['random %s' % r for r in range(5, random.randrange(6, 20))]
            if job == jobs[-1] and random.choice([0, 1]):
                lines.append('Done!')
            for oe in ['.o', '.e']:
                job_fp_oe = '%s/%s%s' % (job_dir, job_fp, oe)
                with open(job_fp_oe, 'w') as o_dev:
                    o_dev.write('\n'.join(lines))
                    if random.choice([0, 1]):
                        o_dev.write('\n%s' % ''.join(outro))
                    else:
                        o_dev.write('\n%s' % ''.join(outro_err))

    def write_script(self, soft=None):
        if self.config.jobs:
            if soft:
                self.prep_script(soft.params)
            else:
                self.prep_script(self.config.params)
            # self.call_cmd()
        else:
            self.job_fps.append('%s.sh' % splitext(self.sh)[0])

    def write_jobs(self, name: str, soft=None):
        for chunk, chunk_keys in self.chunks.items():
            self.get_sh(name, chunk, soft)
            self.write_chunks(chunk_keys, soft)
            self.get_job_name(name, chunk, soft)
            if self.config.dev:
                self.write_dummy_oe()
            self.write_script(soft)

    def display(self):
        for db_soft, name_main in self.run.items():
            if len(self.run[db_soft]):
                print()
                soft_print = '========== %s(s) scripts ========== ' % db_soft
                print(soft_print)
                self.scripts.append(soft_print)
            for name, main in name_main.items():
                if main.startswith('Input'):
                    main_print = '>%s\n%s' % (name, main)
                else:
                    main_print = '>%s\nsh %s' % (name, main)
                print(main_print)
                self.scripts.append(main_print)
