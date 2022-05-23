# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import subprocess
from os.path import dirname, splitext

from metagenomix._io_utils import mkdr, get_chunks


class CreateScripts(object):

    def __init__(self, config):
        self.config = config
        self.jobs_folders = {}
        self.cmd = ''
        self.soft = ''
        self.shs = {}
        self.non_empty_shs = []
        self.all_sh_chunks = []
        self.to_chunk = {}
        self.job_name = ''
        self.run_main = ''

    def get_job_name(
            self,
            sam: str,
            adx: int
    ) -> None:
        """

        Parameters
        ----------
        sam : str
            Current sample (could be None).
        adx : int
            Chunk number
        """
        nickname = [x for x in self.config.project if x.lower() not in 'aeiouy']
        if sam:
            self.job_name = '.'.join([self.soft, ''.join(nickname), sam, adx])
        else:
            self.job_name = '.'.join([self.soft, ''.join(nickname), adx])

    def get_run_sh_fp(
            self,
            soft_prev: str,
            sam: str,
            adx: int
    ) -> tuple:
        """

        Parameters
        ----------
        soft_prev : str
        sam : str
        adx : int

        Returns
        -------
        run_sh: str
        run_top_sh: str
        run_pbs: str
        run_slm: str
        previous_sh: str
        """
        if sam:
            run_sh = '%s/%s/jobs/run_%s_after_%s_%s_%s.sh' % (
                self.config.dir, self.soft, self.soft, soft_prev, sam, adx)
        else:
            run_sh = '%s/%s/jobs/run_%s_after_%s_%s.sh' % (
                self.config.dir, self.soft, self.soft, soft_prev, adx)
        run_top_sh = '%s_top.sh' % splitext(run_sh)[0]
        previous_sh = glob.glob('%s/run_%s_after_%s_*.sh' % (
            dirname(run_sh), self.soft, soft_prev))
        run_dir = dirname(run_sh)
        mkdr(run_dir)
        run_pbs = run_sh.replace('.sh', '.pbs')
        run_slm = run_sh.replace('.sh', '.slm')
        return run_sh, run_top_sh, run_pbs, run_slm, previous_sh

    def get_non_empty_shs(
            self,
            all_sh: str
    ) -> None:
        """Reduce a list of bash scripts to those that are not empty."""
        self.non_empty_shs = [sh for sh in all_sh if len([
            x for x in open(sh).readlines() if len(x.strip())])]

    def get_sh_chunks(
            self,
            all_sh: list,
            n_chunks: int
    ) -> None:
        """Get the chunk of scripts to run, one chunk will be one job.

        Parameters
        ----------
        all_sh : list
            List of bash scripts.
        n_chunks : int
            Number of chunk of scripts to prepare, i.e., to run in queue.
        """
        if n_chunks:
            if len(all_sh) < n_chunks:
                self.all_sh_chunks = [[x] for x in all_sh]
            else:
                self.all_sh_chunks = get_chunks(all_sh, 0, n_chunks)
        else:
            if len(all_sh) <= 4:
                self.all_sh_chunks = [[x] for x in all_sh]
            if len(all_sh) > 8:
                self.all_sh_chunks = get_chunks(all_sh, 0, 8)
            else:
                self.all_sh_chunks = get_chunks(all_sh, 2)

    def prep_script(
            self,
            run_sh: str,
            run_pbs: str,
            cnd: str,
            wall_time: str,
            procs: int,
            nodes: int,
            mem_n: str,
            mem_u: str,
    ) -> None:
        """Get the Torque or Slurm script.

        Parameters
        ----------
        run_sh : str
        run_pbs : str
        cnd : str
        wall_time : str
        nodes : int
        procs : int
        mem_n : int
        mem_u : str
        """
        self.cmd = [
            'Xhpc',
            '-i', run_sh,
            '-o', run_pbs,
            '-j', self.job_name,
            '-e', cnd,
            '-t', wall_time,
            '-c', str(procs),
            '-n', str(nodes),
            '-M', str(mem_n), mem_u,
            '--no-stat'
        ]

    def prep_slm_script(
            self,
    ) -> None:
        """
        Get the Slurm script.
        """
        self.cmd.extend(['--slurm',
                         '-T', '/panfs/flejzerowicz',
                         '-l', '/panfs/flejzerowicz'])

    def call_cmd(self):
        if self.config.show:
            if self.soft in self.config.which_to_show:
                print(' '.join(self.cmd))
        else:
            print(' '.join(self.cmd))
        subprocess.call(self.cmd)

    def get_main_sh(
            self,
            soft_prev: str,
            all_pbs: str,
            all_slm: str
    ) -> None:
        # write the file that will launch all the software jobs
        self.run_main = '%s/%s/run_main_%s_after_%s.sh' % (
            self.config.dir, self.soft, self.soft, soft_prev)
        with open(self.run_main, 'w') as o:
            if self.config.slurm:
                for slm in all_slm:
                    o.write('sbatch %s\n' % slm)
            else:
                for pbs in all_pbs:
                    o.write('qsub %s\n' % pbs)

    def get_pbs_slm(
            self,
            all_sh,
            soft,
            soft_prev,
            mem_n_u,
            procs,
            cnd,
            sam=None,
            n_chunks=None,
    ):
        """

        Parameters
        ----------
        all_sh
        soft
        soft_prev
        mem_n_u
        procs
        cnd
        sam
        n_chunks

        Returns
        -------

        """
        # get pbs parameters
        nodes, procs = self.get_n_procs(procs)

        all_sh = self.get_non_empty_shs(all_sh)
        all_sh_chunks = self.get_sh_chunks(all_sh, n_chunks)
        mem_n, mem_u = mem_n_u

        all_pbs = []
        all_slm = []
        for adx, all_sh_chunk in enumerate(all_sh_chunks):
            run_sh, run_top_sh, run_pbs, run_slm, previous_sh = self.get_run_sh_fp(
                soft, soft_prev, sam, adx)

            # actual .sh command writing
            with open(run_sh, 'w') as sh:
                for chunk in all_sh_chunk:
                    sh.write('\n')
                    with open(chunk) as f:
                        for line in f:
                            sh.write(line)
            self.get_job_name(sam, adx)

            if self.config.slurm and soft in ['instrain', 'graspx']:
                with open(run_top_sh, 'w') as sh_top:
                    sh_top.write('sh %s\n' % run_sh)
                self.prep_slm_script(
                    run_top_sh, run_slm, cnd, nodes, wall_time, procs,
                    mem_n, mem_u)
                all_slm.append(run_slm)
            else:
                if self.config.slurm:
                    self.prep_slm_script(
                        run_sh, run_slm, cnd, nodes, wall_time, procs,
                        mem_n, mem_u)
                    all_slm.append(run_slm)
                else:
                    self.prep_pbs_script(
                        run_sh, run_pbs, cnd, nodes, wall_time,
                        procs, mem_n, mem_u)
                    all_pbs.append(run_pbs)
            self.call_cmd()

        return all_pbs, all_slm

    # def get_n_procs(self) -> None:
    #     """
    #     Get the number of nodes and CPUs.
    #     """
    #     self.nodes = 1
    #     if self.procs > 12:
    #         self.nodes = (((self.procs - 1) // 12) + 1)
    #         p = self.procs // self.nodes
    #         if self.nodes * p < self.procs:
    #             p += 1
    #         self.procs = p

    # def get_sh(self, sam: str = None) -> str:
    #     """Get the path to the job file for one sample
    #     (or one all-samples analysis).
    #
    #     Parameters
    #     ----------
    #     sam : str
    #         Current sample name (could be None).
    #
    #     Returns
    #     -------
    #     run_sh : str
    #         Path to the job file for one sample (or one all-samples analysis).
    #     """
    #     run_dir = '%s/%s/jobs/staged' % (self.config.dir, self.name)
    #     mkdr(run_dir)
    #     run_sh = '%s/staged_%s_after_%s' % (run_dir, self.name, self.prev)
    #     if sam:
    #         run_sh += '_%s' % sam.replace(' ', '_')
    #     run_sh += '.sh'
    #     return run_sh
    #
    # def get_outdir(self) -> str:
    #     """Get the path to the output folder for the current analysis.
    #
    #     Returns
    #     -------
    #     out_dir : str
    #         Path to the output folder for the current analysis.
    #     """
    #     out_dir = '%s/%s/after_%s' % (self.config.dir, self.name, self.prev)
    #     mkdr(out_dir)
    #     return out_dir
