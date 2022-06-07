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
        self.fastq = {}
        self.fastq_scratch = {}
        self.r = {}
        self.dir = ''
        self.params = {}
        # self.cazy_focus_dbs = {}  <---- see databases "self.cazys"
        self.instrain = {'refs': {}, 'bams': {}}
        self.plass_type = ''  # 'nucl'
        self.midas_foci = {'all': ('', '')}
        # self.midas_foci = {'CLA': ('', 'CLAmicrobes_20species.txt')}
        self.midas_strain_tracking = []
        self.simka_params = {}
        self.humann_profile = {}
        # self.humann_profile = {'cla': 'cla_profile_rep82_7_filt.txt'}
        self.stand_alone_groupings = {'metamarker': ('metawrap_ref', 'is_ccz')}

    def init(self):
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
                fastqs = [x for x in fastqs_ if '.gz' in x]
            else:
                fastqs = fastqs_
            self.fastq[sam] = fastqs
            self.fastq_scratch[sam] = ['${SCRATCH_DIR}%s' % x for x in fastqs]

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
        self.dir = self.output_dir

    def parse_yamls(self):
        for arg in list(self.__dict__):
            if arg.endswith('_yml'):
                yaml = read_yaml(self.__dict__[arg])
                setattr(self, arg[:-4], yaml)
            if arg == 'pipeline_tsv':
                with open(self.__dict__[arg]) as f:
                    read = [['edit_fastqs']] + [
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

    # def get_soft_paths(self) -> None:
    #     self.soft_paths = [
    #         'prepare_ioncom_inputs.py',
    #         '/opt/infernal/1.0.2/bin/cmsearch',
    #         '/opt/trimmomatic/0.36/trimmomatic-0.36.jar',
    #         '/home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208',
    #         '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/bowtie2',
    #         '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/diamond',
    #         '/home/flejzerowicz/softs',
    #         '/home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin',
    #         '/home/flejzerowicz/usr/local/genometools/bin',
    #         '/home/flejzerowicz/usr/bin/hmmsearch',
    #         '/home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin/spades.py',
    #         '/home/flejzerowicz/softs/SPAdes-3.14.1-Linux/bin/spades.py',
    #         '/home/flejzerowicz/softs/SPAdes-3.15.0-corona-2020-07-15/spades.py',
    #         '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/metaphlan_databases/mpa_v20_m200.pkl',
    #         '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/metaphlan_databases',
    #         '/home/flejzerowicz/databases/dbCAN/db',
    #         '/home/flejzerowicz/databases/midas_db_v1.2',
    #         '/home/flejzerowicz/databases/checkM',
    #         '/home/flejzerowicz/databases/metaphlan3',
    #         '/home/flejzerowicz/databases/PFAM/Pfam-A.hmm',
    #         '/home/flejzerowicz/databases/PFAM/Pfam-A.hmm.dat',
    #         '/home/flejzerowicz/databases/PFAM/Pfam-A.fasta',
    #         '/databases/humann2_data/full_chocophlan.v0.1.1/chocophlan',
    #         '/databases/humann2_data/uniref90/uniref',
    #         '/databases/bowtie/Human/Human',
    #         '/databases/uniref/uniref50/uniref50_2018_02.dmnd',
    #         '/databases/db-shogun-20181114/functions-uniprot',
    #         '/databases/db-shogun-20181114/functions-refseq',
    #         '/databases/db-shogun-20181114/functions-kegg',
    #         '/databases/genome/rep82/shogun',
    #         '/databases/genome/rep200/shogun',
    #         '/databases/db-shogun-20181114',
    #         '/projects/wol/release/proteins/all.faa'
    #         '/projects/wol/profiling/dbs/wol/shogun',
    #         '/projects/wol/20170307/diamond/prok.dmnd',
    #         '/home/flejzerowicz/databases/wol/lineage.txt',
    #         '/home/flejzerowicz/databases/wol/genome_sizes.txt',
    #         '/home/flejzerowicz/softs/python_scripts/kegg_query.py',
    #         '/home/flejzerowicz/databases/checkM/genome_tree/genome_tree.taxonomy.tsv',
    #         '/databases/gtdb/release202/taxonomy/gtdb_taxonomy.tsv',
    #         #         '/home/flejzerowicz/databases/wol/g2tid.txt',
    #         #         '/home/flejzerowicz/databases/wol/g2gg.tax',
    #         #         '/home/flejzerowicz/databases/wol/gg2tid.txt',
    #         #         '/home/flejzerowicz/databases/wol/nucl2g.txt',
    #         #         '/home/flejzerowicz/databases/wol/nucl2gg.tax',
    #         #         '/home/flejzerowicz/databases/wol/tid2gg.txt',
    #         #         '/home/flejzerowicz/databases/wol/nodes.dmp',
    #         #         '/home/flejzerowicz/databases/wol/names.dmp',
    #         #         '/projects/wol/release/annotation/uniref.map.xz',
    #         #         '/projects/wol/release/annotation/metacyc/enzrxn2reaction.map',
    #         #         '/projects/wol/release/annotation/metacyc/pathway2class.map',
    #         #         '/projects/wol/release/annotation/metacyc/protein.map',
    #         #         '/projects/wol/release/annotation/metacyc/protein2enzrxn.map',
    #         #         '/projects/wol/release/annotation/metacyc/protein2gene.map',
    #         #         '/projects/wol/release/annotation/metacyc/reaction2ec.map',
    #         #         '/projects/wol/release/annotation/metacyc/reaction2pathway.map',
    #         #         '/projects/wol/release/annotation/coords.txt.xz',
    #         #         '/projects/wol/release/annotation/go/process.tsv.xz',
    #         #         '/projects/wol/release/annotation/go/function.tsv.xz',
    #         #         '/projects/wol/release/annotation/go/component.tsv.xz',
    #         #         '/projects/wol/release/annotation/refseq/refseq.map.xz',
    #         #         '/projects/wol/release/annotation/eggnog/eggnog.map.xz',
    #     ]
