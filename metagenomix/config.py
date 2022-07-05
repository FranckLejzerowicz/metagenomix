# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import subprocess
import pkg_resources
from os.path import exists, isdir, isfile, splitext
from metagenomix._io_utils import read_yaml, fill_fastq
from metagenomix._metadata import read_metadata
import pandas as pd
import numpy as np

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources")


class AnalysesConfig(object):
    """Collect the data associated with each dataset passed but the user."""
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.meta = pd.DataFrame()
        self.pooling_groups = []
        self.conda_envs = {}
        self.modules = {}
        self.soft_paths = []
        self.fastq_scratch = {}
        self.fastq = {}
        self.params = {}
        self.dir = ''
        self.r = {}
        # self.cazy_focus_dbs = {}  <---- see databases "self.cazys"
        self.instrain = {'refs': {}, 'bams': {}}

    def run(self):
        self.check_xhpc_install()
        self.get_conda_envs()
        self.set_metadata()
        self.set_fastq()
        self.get_r()
        self.set_output()
        self.parse_yamls()
        self.set_coassembly()
        self.update_metadata()
        self.get_default_params()

    def check_xhpc_install(self):
        """Try to get the install path of third party tool
        [Xhpc](https://github.com/FranckLejzerowicz/Xhpc).
            If it exists, nothing happens and the code proceeds.
            Otherwise, the code ends and tells what to do.
        """
        if self.jobs:
            ret, _ = subprocess.getstatusoutput('which Xhpc')
            if ret:
                raise IOError('Xhpc not installed')
            else:
                xpbs_config = subprocess.getoutput('Xhpc --show-config')
                if not xpbs_config.startswith('* email address config file'):
                    raise IOError('Need to run Xhpc to setup config file')

    def get_conda_envs(self):
        """Get the names of the conda environments."""
        for env in subprocess.getoutput('conda env list').split('\n'):
            if env.startswith('#') or '*' in env or not env.strip():
                continue
            name, path = env.split()
            self.conda_envs[name] = path

    def set_metadata(self):
        """Read metadata with first column as index."""
        if not isfile(self.meta_fp):
            raise IOError('No file "%s"' % self.meta_fp)
        self.meta = read_metadata(self.meta_fp)

    def set_coassembly(self):
        """Create a metadata variable for the groups on which to co-assemble."""
        if self.coassembly:
            for name, coassembly in self.coassembly.items():
                self.get_pooling_groups(name, coassembly)
        else:
            self.pooling_groups = ['assembly_per_sample']
            self.meta['assembly_per_sample'] = self.meta.sample_name

    def get_pooling_groups(
            self,
            name: str,
            coassembly: dict
    ) -> None:
        """Get the column(s) corresponding to the factors groups
        identifying the samples to co-assemble.

        Parameters
        ----------
        name : str
            Name of the co-assembly.
        coassembly : dict
            Factors groups identifying the samples to co-assemble
        """
        cols = {}
        for col, factors in coassembly.items():
            if col not in self.meta.columns:
                raise IOError(
                    'Co-assembly variable "%s" not in %s' % (col, self.meta_fp))
            factors_flat = [str(f) for factor in factors for f in factor]
            if not set(factors_flat).issubset(set(self.meta[col])):
                raise IOError('Co-assembly factors for variable "%s" not'
                              ' in %s' % (col, self.meta_fp))
            d = {str(y): '_'.join(map(str, x)) for x in factors for y in x}
            cols[col] = [d[x] if x in d else np.nan for x in self.meta[col]]
        ser = pd.DataFrame(cols).fillna('')
        col = ser.apply(lambda x: '-'.join([str(i) for i in x if i]), axis=1)
        self.meta[name] = col.replace('', np.nan)
        self.pooling_groups.append(name)

    def get_fastq_paths(self) -> list:
        fastqs = []
        for fastq_dir in self.fastq_dirs:
            fastqs.extend(glob.glob(fastq_dir + '/*.fastq*'))
        return fastqs

    def get_fastq_samples(self):
        fastq_paths = self.get_fastq_paths()
        fastq = fill_fastq(fastq_paths, set(self.meta.sample_name))
        # keep only the `.fastq.gz` files (if `.fastq` files are also present)
        for sam, fastqs_ in fastq.items():
            if len([x for x in fastqs_ if '.gz' in x]):
                fqs = [x for x in fastqs_ if '.gz' in x]
            else:
                fqs = fastqs_
            self.fastq[sam] = fqs
            self.fastq_scratch[sam] = ['${SCRATCH_FOLDER}%s' % x for x in fqs]

    def set_fastq(self):
        """
        Check that fastq folder exists and that it contains fastq files.
        """
        if sum([not exists(x) for x in self.fastq_dirs]):
            raise IOError('Fastq folder(s) do not exist')
        self.get_fastq_samples()

    def get_r(self):
        self.r = {}
        for sam, r1_r2 in self.fastq.items():
            self.r[sam] = []
            if len(r1_r2) == 2:
                r1, r2 = r1_r2
                diffs = [i for i in range(len(r1)) if r1[i] != r2[i]]
                if len(diffs) > 1:
                    raise IOError('Too different fastq names for sample "%s":\n'
                                  '- %s\n- %s' % (sam, r1, r2))
                self.r[sam] = [r1[diffs[0] - 1:].split('.fastq')[0],
                               r2[diffs[0] - 1:].split('.fastq')[0]]

    def update_metadata(self):
        self.meta.set_index('sample_name', inplace=True)
        self.meta = self.meta.loc[list(self.fastq.keys())]
        self.meta.fillna('Unspecified', inplace=True)
        self.meta_fp = '%s_pipeline.tsv' % splitext(self.meta_fp)[0]
        self.meta.to_csv(self.meta_fp, sep='\t')

    def set_output(self):
        """Check if main output folder exists or create it."""
        if not isdir(self.output_dir):
            os.makedirs(self.output_dir)
        self.dir = os.path.abspath(self.output_dir)

    def parse_yamls(self):
        for arg in list(self.__dict__):
            if arg.endswith('_yml'):
                yaml = read_yaml(self.__dict__[arg])
                setattr(self, arg[:-4], yaml)
            if arg == 'pipeline_tsv':
                with open(self.__dict__[arg]) as f:
                    read = [['edit']] + [
                        x.strip().split() for x in f.readlines()
                        if x[0] != '#' and len(x.strip().split())]
                    setattr(self, arg[:-4], read)

    def get_default_params(self):
        """
        Get the run parameters based on the default
        values, that are possibly updated by the user.
        """
        if isfile('%s/run_params.yml' % RESOURCES):
            self.params = read_yaml('%s/run_params.yml' % RESOURCES)['default']
