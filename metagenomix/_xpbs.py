# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
from os.path import isfile
from typing import TextIO


def run_xpbs(out_sh: str, out_pbs: str, job_name: str,
             qiime_env: str, time: str, n_nodes: str,
             n_procs: str, mem_num: str, mem_dim: str, chmod: str,
             written: int, single: str, o: TextIO = None,
             noloc: bool = True, slurm: bool = False, jobs: bool = True,
             tmp: str = None) -> None:
    """
    Run the Xpbs script assorted with print or writing in higher-level command.
    """
    if written:
        if os.getcwd().startswith('/panfs'):
            out_sh_lines = open(out_sh).readlines()
            with open(out_sh, 'w') as sh:
                if not jobs:
                    sh.write('conda activate %s\n' % qiime_env)
                for out_sh_line in out_sh_lines:
                    sh.write(out_sh_line.replace(os.getcwd(), ''))
        if jobs:
            xpbs_call(out_sh, out_pbs, job_name, qiime_env,
                      time, n_nodes, n_procs, mem_num,
                      mem_dim, chmod, slurm, noloc, tmp)
        if single:
            if slurm:
                launcher = 'sbatch'
            else:
                launcher = 'qsub'
            if os.getcwd().startswith('/panfs'):
                out_pbs = out_pbs.replace(os.getcwd(), '')
            if not o:
                print(single)
                if jobs:
                    print('[TO RUN]', launcher, out_pbs)
                else:
                    print('[TO RUN] sh', out_sh)
            else:
                if jobs:
                    o.write('%s %s\n' % (launcher, out_pbs))
                else:
                    o.write('sh %s\n' % out_sh)
    else:
        os.remove(out_sh)
        if isfile(out_pbs):
            os.remove(out_pbs)


def xpbs_call(out_sh: str, out_pbs: str, prjct_nm: str,
              qiime_env: str, time: str, n_nodes: str,
              n_procs: str, mem_num: str, mem_dim: str,
              chmod: str, slurm: bool, noloc: bool,
              tmp: str = None) -> None:
    """
    Call the subprocess to run Xpbs on the current bash script.
    """
    cmd = [
        'Xpbs',
        '-i', out_sh,
        '-j', prjct_nm,
        '-o', out_pbs,
        '-e', qiime_env,
        '-t', time,
        '-n', n_nodes,
        '-p', n_procs,
        '-M', mem_num, mem_dim,
        '-c', chmod,
        '--noq'
    ]
    if tmp:
        cmd.extend(['-T', tmp])
    if not noloc:
        cmd.append('--no-loc')
    if slurm:
        cmd.append('--slurm')
    subprocess.call(cmd)

    if os.getcwd().startswith('/panfs'):
        out_pbs_lines = open(out_pbs).readlines()
        with open(out_pbs, 'w') as pbs:
            for out_pbs_line in out_pbs_lines:
                pbs.write(out_pbs_line.replace(os.getcwd(), ''))


def print_message(message: str, sh_pbs: str, to_run: str, jobs: bool) -> None:
    if message:
        print(message)
    if os.getcwd().startswith('/panfs'):
        to_run = to_run.replace(os.getcwd(), '')
    if jobs:
        print('[TO RUN]', sh_pbs, to_run)
    else:
        print('[TO RUN]', 'sh', to_run.replace('.pbs', '.sh'))
