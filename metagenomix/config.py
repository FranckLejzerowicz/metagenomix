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
from os.path import exists, isdir, isfile
from metagenomix._io_utils import read_yaml
from metagenomix._metadata import read_metadata
import pandas as pd
import numpy as np

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources")


class AnalysesConfig(object):
    """Collect the data associated with each dataset passed but the user
    """
    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)
        self.meta = pd.DataFrame()
        self.co_assembly_var = None
        self.pooling_groups = []
        self.conda_envs = {}
        self.conda_path = ''
        self.fastq = {}
        self.r = {}
        self.dir = ''
        self.db_type = 'genomes'
        # self.db_type = 'uniref'
        self.params = {'time': 48, 'mem_num': 10, 'mem_dim': "gb",
                       'env': "mg", 'nodes': 1, 'cpus': 1, 'chunks': 6}
        # self.cazy_focus_dbs = {}  <---- see databases "self.cazys"
        self.instrain = {'refs': {}, 'bams': {}}
        self.plass_assembly = ''  # 'nucl'
        self.midas_focus = {'all': ('', '')}
        # self.midas_focus = {'CLAmicrobes': ('', 'CLAmicrobes_20species.txt'), 'all': ('', '')}
        self.midas_strain_tracking = []
        self.integron_focus = { # COULD BE THE SAME AS THE HMMs COLLECTED IN database()
            # 'PFAM_Fatty_acid_cistrans_isomerase_PF06934': 'PFAM_Fatty_acid_cistrans_isomerase_PF06934_seed_profile.hmm',
            # 'PFAM_MCRA_PF06100': 'PFAM_MCRA_PF06100_seed_profile.hmm',
            # 'PFAM_Short_chain_fatty_acid_transporter_PF026677': 'PFAM_Short_chain_fatty_acid_transporter_PF026677_seed_profile.hmm',
            # 'UniProt_map_IDs2NCBI_ali_kept_final': 'UniProt_map_IDs2NCBI_ali_kept_final_profile.hmm',
            'dbCAN': '/Users/franck/databases/dbCAN/db/dbCAN.txt'}
        self.humann2_taxo_profile_glob = {}
        # self.humann2_taxo_profile_glob = {'cla_papers': 'cla_papers_taxonomic_profile_rep82_7_filt.txt'}
        self.flash_params = {'min_overlap': 30, 'max_mismatch_density': 0}
        self.stand_alone_groupings = {'metamarker': ('metawrap_ref', 'is_ccz')}
        self.highmem_nodes = ['35', '36', '72', '73']
        self.highmem_node_idx = 0
        self.soft_paths = []
        # self.workflow = {}

    def init(self):
        self.check_xpbs_install()
        self.get_conda_envs()
        self.set_metadata()
        self.set_fastq()
        self.get_r()
        self.set_output()
        self.parse_yamls()
        self.set_coassembly()
        self.update_metadata()
        self.get_default_params()
        self.get_soft_paths()

    def check_xpbs_install(self):
        """Try to get the install path of third party tool
        [Xpbs](https://github.com/FranckLejzerowicz/Xpbs).
            If it exists, nothing happens and the code proceeds.
            Otherwise, the code ends and tells what to do.
        """
        if self.jobs:
            ret, _ = subprocess.getstatusoutput('which Xpbs')
            if ret:
                raise IOError('Xpbs not installed')
            else:
                xpbs_config = subprocess.getoutput('Xpbs --show-config')
                if xpbs_config.startswith('$HOME'):
                    raise IOError('Xpbs installed but config.txt need editing')

    def get_conda_envs(self):
        """Get the names of the conda environments."""
        for env in subprocess.getoutput('conda env list').split('\n'):
            if '/envs/' in env:
                self.conda_envs[env.split('/')[-1]] = env.split()[-1]
                self.conda_path = env.split()[-1].split('/envs/')[0]

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
            self.pooling_groups = ['per_sample_no_coassembly']
            self.meta['per_sample_no_coassembly'] = self.meta.sample_name

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
        cols = []
        for col, factors in coassembly.items():
            if col not in self.meta.columns:
                raise IOError('Co-assembly variable not in %s' % self.meta_fp)
            factors_flat = [f for factor in factors for f in factor]
            if not set(factors_flat).issubset(set(self.meta[col])):
                raise IOError('Co-assembly factors not in %s' % self.meta[col])
            d = {y: '_'.join(map(str, x)) for x in factors for y in x}
            cols.append([d[x] if x in d else np.nan for x in self.meta[col]])
        col = ['_'.join(map(str, x)) for x in zip(*cols)]
        self.meta[name] = col
        self.pooling_groups.append(name)

    def get_fastq_paths(self) -> list:
        fastqs = []
        for fastq_dir in self.fastq_dirs:
            fastqs.extend(glob.glob(fastq_dir + '/*.fastq*'))
        return fastqs

    def get_fastq_samples(self):
        fastqs = self.get_fastq_paths()
        sams = set(self.meta.sample_name)
        self.fastq = {}
        for file in sorted(fastqs):
            for sam in sams:
                if sam in file:
                    if sam in self.fastq:
                        self.fastq[sam].append(file)
                    else:
                        self.fastq[sam] = [file]
                    break

    def set_fastq(self):
        """
        Check that fastq folder exists and that it contains fastq files.
        """
        if sum([not exists(x) for x in self.fastq_dirs]):
            raise IOError('Fastq folder(s) do not exist')
        self.get_fastq_samples()

    def get_r(self):
        self.r = {}
        for sam, (r1, r2) in self.fastq.items():
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
        self.meta_fp = '%s_pipeline.tsv' % os.path.splitext(self.meta_fp)[0]
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
            if arg.endswith('_tsv'):
                with open(self.__dict__[arg]) as f:
                    read = [x.strip().split() for x
                            in f.readlines() if x[0] != '#']
                    setattr(self, arg[:-4], read)

    def get_default_params(self):
        """
        Get the run parameters based on the default
        values, that are possibly updated by the user.
        """
        if isfile('%s/run_params.yml' % RESOURCES):
            self.params = read_yaml('%s/run_params.yml' % RESOURCES)['default']

    def get_soft_paths(self) -> None:
        self.soft_paths = [
            'prepare_ioncom_inputs.py',
            '/opt/infernal/1.0.2/bin/cmsearch',
            '/opt/trimmomatic/0.36/trimmomatic-0.36.jar',
            '/home/flejzerowicz/softs/cd-hit-v4.6.8-2017-1208',
            '%s/bin' % self.conda_path,
            '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/bowtie2',
            '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/diamond',
            '/home/flejzerowicz/softs',
            '/home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin',
            '/home/flejzerowicz/softs/gOTU/gOTU_from_maps.py',
            '/home/flejzerowicz/usr/local/genometools/bin',
            '/home/flejzerowicz/usr/bin/hmmsearch',
            '/home/flejzerowicz/softs/SPAdes-3.13.0-Linux/bin/spades.py',
            '/home/flejzerowicz/softs/SPAdes-3.14.1-Linux/bin/spades.py',
            '/home/flejzerowicz/softs/SPAdes-3.15.0-corona-2020-07-15/spades.py',
            '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/metaphlan_databases/mpa_v20_m200.pkl',
            '/home/flejzerowicz/usr/miniconda3/envs/humann2/bin/metaphlan_databases',
            '/home/flejzerowicz/databases/dbCAN/db',
            '/home/flejzerowicz/databases/midas_db_v1.2',
            '/home/flejzerowicz/databases/checkM',
            '/home/flejzerowicz/databases/metaphlan3',
            '/home/flejzerowicz/databases/PFAM/Pfam-A.hmm',
            '/home/flejzerowicz/databases/PFAM/Pfam-A.hmm.dat',
            '/home/flejzerowicz/databases/PFAM/Pfam-A.fasta',
            '/databases/humann2_data/full_chocophlan.v0.1.1/chocophlan',
            '/databases/humann2_data/uniref90/uniref',
            '/databases/bowtie/Human/Human',
            '/databases/uniref/uniref50/uniref50_2018_02.dmnd',
            '/databases/db-shogun-20181114/functions-uniprot',
            '/databases/db-shogun-20181114/functions-refseq',
            '/databases/db-shogun-20181114/functions-kegg',
            '/databases/genome/rep82/shogun',
            '/databases/genome/rep200/shogun',
            '/databases/db-shogun-20181114',
            '/projects/wol/release/proteins/all.faa'
            '/projects/wol/profiling/dbs/wol/shogun',
            '/projects/wol/20170307/diamond/prok.dmnd',
            '/home/flejzerowicz/databases/wol/lineage.txt',
            '/home/flejzerowicz/databases/wol/genome_sizes.txt',
            '/home/flejzerowicz/softs/python_scripts/kegg_query.py',
            '/home/flejzerowicz/databases/checkM/genome_tree/genome_tree.taxonomy.tsv',
            '/databases/gtdb/release202/taxonomy/gtdb_taxonomy.tsv',
            #         '/home/flejzerowicz/databases/wol/g2tid.txt',
            #         '/home/flejzerowicz/databases/wol/g2gg.tax',
            #         '/home/flejzerowicz/databases/wol/gg2tid.txt',
            #         '/home/flejzerowicz/databases/wol/nucl2g.txt',
            #         '/home/flejzerowicz/databases/wol/nucl2gg.tax',
            #         '/home/flejzerowicz/databases/wol/tid2gg.txt',
            #         '/home/flejzerowicz/databases/wol/nodes.dmp',
            #         '/home/flejzerowicz/databases/wol/names.dmp',
            #         '/projects/wol/release/annotation/uniref.map.xz',
            #         '/projects/wol/release/annotation/metacyc/enzrxn2reaction.map',
            #         '/projects/wol/release/annotation/metacyc/pathway2class.map',
            #         '/projects/wol/release/annotation/metacyc/protein.map',
            #         '/projects/wol/release/annotation/metacyc/protein2enzrxn.map',
            #         '/projects/wol/release/annotation/metacyc/protein2gene.map',
            #         '/projects/wol/release/annotation/metacyc/reaction2ec.map',
            #         '/projects/wol/release/annotation/metacyc/reaction2pathway.map',
            #         '/projects/wol/release/annotation/coords.txt.xz',
            #         '/projects/wol/release/annotation/go/process.tsv.xz',
            #         '/projects/wol/release/annotation/go/function.tsv.xz',
            #         '/projects/wol/release/annotation/go/component.tsv.xz',
            #         '/projects/wol/release/annotation/refseq/refseq.map.xz',
            #         '/projects/wol/release/annotation/eggnog/eggnog.map.xz',
        ]
