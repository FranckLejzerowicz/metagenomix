# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import sys
import subprocess
import pkg_resources
from os.path import abspath, basename, isdir, isfile, splitext
from metagenomix._io_utils import read_yaml, get_fastq_files, get_fastq_paths
from metagenomix._metadata import read_metadata
import pandas as pd
import numpy as np

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources")


class AnalysesConfig(object):
    """Collect the data associated with each dataset passed but the user."""
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.meta = pd.DataFrame()
        self.tools = {}
        self.pooling_groups = []
        self.conda_envs = {}
        self.modules = {}
        self.soft_paths = []
        self.techs = ['illumina', 'pacbio', 'nanopore']
        self.techs_fastqs = {}
        self.fastq_mv = {}
        self.fastq = {}
        self.params = {}
        self.params_dbs = {}
        self.read = []
        self.dir = ''
        self.r = {}
        # self.cazy_focus_dbs = {}  <---- see databases "self.cazys"
        self.instrain = {'refs': {}, 'bams': {}}
        self.merge = {'midas2'}

    def run(self):
        self.check_xhpc_install()
        self.get_conda_envs()
        self.set_metadata()
        self.get_tools()
        self.get_techs()
        self.set_fastqs()
        self.show_fastqs()
        self.get_r()
        self.set_output()
        self.parse_yamls()
        self.get_params_dbs()
        self.set_coassembly()
        self.update_metadata()
        self.get_default_params()

    def check_xhpc_install(self):
        """Try to get the install path of third party tool
        [Xhpc](https://github.com/FranckLejzerowicz/Xhpc).
            If it exists, nothing happens and the code proceeds.
            Otherwise, the code ends and tells what to do.
        """
        if self.jobs and self.command == 'create':
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
            if env.startswith('#') or '*' in env or len(env.split()) != 2:
                continue
            name, path = env.split()
            self.conda_envs[name] = path

    def set_metadata(self):
        """Read metadata with first column as index."""
        if not isfile(self.meta_fp):
            raise IOError('No file "%s"' % self.meta_fp)
        self.meta = read_metadata(self.meta_fp)

    def get_tools(self):
        self.tools['fastq'] = 'raw data'
        with open('%s/softwares.txt' % RESOURCES) as f:
            for line in f:
                tool, category = line.strip().split('\t')
                self.tools[tool] = category
                self.tools.setdefault(category.split(' (')[0], []).append(tool)

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
            t = (col, self.meta_fp)
            if col not in self.meta.columns:
                sys.exit('Co-assembly variable "%s" not in %s' % t)
            factors_flat = [str(f) for factor in factors for f in factor]
            if not set(factors_flat).issubset(set(self.meta[col])):
                sys.exit('Co-assembly factors for variable "%s" not in %s' % t)
            d = {str(y): '_'.join(map(str, x)) for x in factors for y in x}
            cols[col] = [d[y] if y in d else self.meta.sample_name[x]
                         for x, y in enumerate(self.meta[col])]
        ser = pd.DataFrame(cols).fillna('')
        col = ser.apply(lambda x: '-'.join([str(i) for i in x if i]), axis=1)
        self.meta[name] = col.replace('', np.nan)
        self.pooling_groups.append(name)

    def init_fastq(self, sam):
        def fill_fq_dict(techs):
            return dict(((tech, sam), []) for tech in techs)
        if sam not in self.fastq:
            self.fastq[sam] = fill_fq_dict(self.techs_fastqs)
            self.fastq_mv[sam] = fill_fq_dict(self.techs_fastqs)

    def fill_fastq(self, fastqs: list):
        """Populate a `fastq` dict with for each sample (keys) the list of
        fastq file paths (values), which would be either of length 1 if there
        is only one fastq file for the sample, or of length 2 if there are
        two, paired fastq file paths.

        Notes
        -----
        In fact, the values could have a length 4 since all *fastq* files are
        considered and there could be both `.fastq` and `.fastq.gz` file in the
        input folder. Only the `.fastq.gz` files would be used in a later stage.

        Parameters
        ----------
        fastqs : list
            All fastq files in the input folder

        Returns
        -------
        fastq : dict
            Fastq file path(s) per sample
        """
        fastq = {}
        sams = set(self.meta.sample_name)
        for fq in sorted(fastqs):
            for sam in sams:
                if re.match('%s_R?[1-2].fastq(.gz)?' % sam, basename(fq)):
                    if sam in fastq:
                        fastq[sam].append(abspath(fq))
                    else:
                        fastq[sam] = [abspath(fq)]
                    break
        return fastq

    def get_fastq_samples(self, tech):
        # keep only the `.fastq.gz` files (if `.fastq` files are also present)
        for sam, fastqs in self.fill_fastq(self.techs_fastqs[tech]).items():
            self.init_fastq(sam)
            fqs = get_fastq_files(fastqs)
            key = (tech, sam)
            self.fastq[sam][key] = fqs
            self.fastq_mv[sam][key] = ['${SCRATCH_FOLDER}%s' % x for x in fqs]

    def get_techs(self):
        for tech_dir in [x for x in self.__dict__.keys() if x.endswith('dirs')]:
            tech = tech_dir.split('_')[0]
            fastq_paths = get_fastq_paths(self.__dict__[tech_dir])
            if fastq_paths:
                self.techs_fastqs[tech] = fastq_paths

    def set_fastqs(self):
        """
        Check that fastq folder exists and that it contains fastq files.
        """
        for tech in self.techs_fastqs:
            self.get_fastq_samples(tech)
        if not sum([len(y) for _, x in self.fastq.items() for y in x.values()]):
            sys.exit('Input fastq folder(s) do not exist')

    def show_fastqs(self):
        if self.verbose:
            max_sam_len = max([len(x) for x in self.fastq])
            print('\n========\n Inputs\n========\n')
            print('sample%s %s' % (' '*(max_sam_len - 6), ' '.join(self.techs)))
            for sam in self.fastq:
                print('%s%s' % (sam, ' '*(max_sam_len - len(sam))), end='')
                for tech in self.techs:
                    n = 0
                    if self.fastq[sam].get((tech, sam)):
                        n = len(self.fastq[sam][(tech, sam)])
                    print(' %s%s' % (n, ' '*(len(tech) - len(str(n)))), end='')
                print()

    def get_r(self):
        self.r = {}
        for sam in self.fastq:
            if 'illumina' not in self.fastq[sam]:
                continue
            r1_r2 = self.fastq[sam][('illumina', sam)]
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
        self.meta = self.meta.loc[list(self.fastq)]
        self.meta.fillna('Unspecified', inplace=True)
        self.meta_fp = '%s_pipeline.tsv' % splitext(self.meta_fp)[0]
        self.meta.to_csv(self.meta_fp, sep='\t')

    def set_output(self):
        """Check if main output folder exists or create it."""
        if not isdir(self.output_dir):
            os.makedirs(self.output_dir)
        self.dir = os.path.abspath(self.output_dir)

    def add_mergers(self):
        mergers = {}
        for softs in self.read:
            soft = softs[-1].split('_')[0]
            if soft in self.merge:
                mergers[soft] = [softs[-1], '%s_merge' % soft]
        self.read.extend(list(mergers.values()))

    def parse_yamls(self):
        for arg in list(self.__dict__):
            if arg.endswith('_yml'):
                yaml = read_yaml(self.__dict__[arg])
                setattr(self, arg[:-4], yaml)
            if arg == 'pipeline_tsv':
                with open(self.__dict__[arg]) as f:
                    self.read = [
                        x.strip().split() for x in f.readlines()
                        if x[0] != '#' and len(x.strip().split())]
                    self.add_mergers()
                    setattr(self, arg[:-4], self.read)

    def get_params_dbs(self):
        for soft, params in self.user_params.items():
            if "databases" in params:
                self.params_dbs[soft] = params.get("databases")

    def get_default_params(self):
        """
        Get the run parameters based on the default
        values, that are possibly updated by the user.
        """
        if isfile('%s/run_params.yml' % RESOURCES):
            self.params = read_yaml('%s/run_params.yml' % RESOURCES)['default']
        if self.chunks:
            self.params['chunks'] = self.chunks
