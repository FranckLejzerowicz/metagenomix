# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import unittest
import subprocess
import pkg_resources
from os.path import isdir

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from metagenomix.config import AnalysesConfig

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')
FOLDER = pkg_resources.resource_filename('metagenomix', 'tests/unittests')


class TestIO(unittest.TestCase):

    def setUp(self) -> None:

        self.config_input = {
            'meta_fp': '',
            'fastq_dirs': '%s/fastqs' % FOLDER,
            'output_dir': '%s/output' % FOLDER,
            'project': 'unittest',
            'coassembly': '',
            'user_params_yml': None,
            'force': False,
            'jobs': False,
            'slurm': False
        }
        self.conda_envs = [x.split('/')[-1] for x in glob.glob(
            os.environ['CONDA_EXE'].replace('bin/conda', 'envs/*'))]

        self.nofile_fp = 'nofile.txt'

        self.empty_meta_fp = '%s/metadata/empty_meta.txt' % FOLDER
        with open(self.empty_meta_fp, 'w'):
            pass

        self.only_header_fp = '%s/metadata/only_header.txt' % FOLDER
        with open(self.only_header_fp, 'w') as o:
            o.write('sample_name\tdepth\tlight')

        self.one_sample_meta_fp = '%s/metadata/meta_one_sample.txt' % FOLDER
        with open(self.one_sample_meta_fp, 'w') as o:
            o.write('some_name\n')
            o.write('A\nX')

        self.meta_fp = '%s/metadata/meta.txt' % FOLDER
        self.meta_pipeline_fp = '%s/metadata/meta_pipeline.tsv' % FOLDER
        with open(self.meta_fp, 'w') as o:
            o.write('sample_name\tdepth\tlight\n')
            o.write('A\t10\tyes\nB\t2000\tno\nX\t0\tyes')

        self.alt_header_fp = '%s/metadata/alt_header.txt' % FOLDER
        with open(self.alt_header_fp, 'w') as o:
            o.write('some_name\tdepth\tlight\n')
            o.write('A\t10\tyes\nB\t2000\tno\nX\t0\tyes')

        self.no_coassembly_fp = '%s/metadata/no_coassembly.txt' % FOLDER
        with open(self.no_coassembly_fp, 'w') as o:
            o.write('sample_name\tdepth\tassembly_per_sample\n')
            o.write('A\t10\tgp1\nB\t2000\tgp2\nX\t0\tgp2\n')

        self.meta = pd.DataFrame(
            {'sample_name': ['A', 'B', 'X'], 'depth': ['10', '2000', '0'],
             'light': ['yes', 'no', 'yes']})

        self.meta_no_coassembly = pd.DataFrame({
            'sample_name': ['A', 'B', 'X'], 'depth': ['10', '2000', '0'],
            'assembly_per_sample': ['gp1', 'gp2', 'gp2']})

        self.meta_fastqed = pd.DataFrame({
            'depth': ['10', '2000'], 'light': ['yes', 'no']}, index=['A', 'B'])
        self.meta_fastqed.index.name = 'sample_name'

        self.fastqs = {}
        self.fastqs1 = []
        for fastq in ['A_1.fastq', 'A_2.fastq', 'B_1.fastq', 'B_2.fastq']:
            fastq_path = '%s/fastqs/%s' % (FOLDER, fastq)
            self.fastqs.setdefault(fastq[0], []).append(fastq_path)
            self.fastqs1.append(fastq_path)
            with open(fastq_path, 'w'): pass

        self.fastqs2 = []
        for fastq in ['C_1.fastq', 'C_2.fastq', 'D_1.fastq', 'D_2.fastq']:
            fastq_path = '%s/output/%s' % (FOLDER, fastq)
            self.fastqs2.append(fastq_path)
            with open(fastq_path, 'w'): pass

        self.fastqs3 = []
        self.fastqs3_only_gz = []
        for fastq in ['A_1.fastq', 'A_2.fastq', 'A_1.fastq.gz', 'A_2.fastq.gz']:
            fastq_path = '%s/fastqs_gz/%s' % (FOLDER, fastq)
            self.fastqs3.append(fastq_path)
            if fastq.endswith('gz'):
                self.fastqs3_only_gz.append(fastq_path)
            with open(fastq_path, 'w'): pass

        self.yaml1 = {'a': 1, 'b': 2}
        self.yaml2 = {'x': 1, 'y': 2}
        self.yaml1_yml = '%s/yamls/yaml1.yml' % FOLDER
        self.yaml2_yml = '%s/yamls/yaml2.yml' % FOLDER
        with open(self.yaml1_yml, 'w') as o:
            o.write('a: 1\nb: 2')
        with open(self.yaml2_yml, 'w') as o:
            o.write('x: 1\ny: 2')

        self.params = {'time': 48, 'mem_num': 10, 'mem_dim': "gb", 'scratch': 0,
                       'env': "mg", 'nodes': 1, 'cpus': 1, 'chunks': 6}
        self.user_params_yml = '%s/yamls/shogun.yml' % FOLDER
        with open(self.user_params_yml, 'w') as o:
            o.write('shogun:\n  time: 1\n')

    def test_check_xpbs_install(self):

        conf = AnalysesConfig(**self.config_input)

        if subprocess.getstatusoutput('which Xpbs')[0]:
            with self.assertRaises(IOError) as e:
                conf.check_xpbs_install()
            self.assertEqual('Xpbs not installed', str(e.exception))
        else:
            conf.check_xpbs_install()
            conf.jobs = True
            conf.check_xpbs_install()

    def test_get_conda_envs(self):

        conf = AnalysesConfig(**self.config_input)
        conf.get_conda_envs()
        self.assertEqual(sorted([x for x in conf.conda_envs if x != 'base']),
                         sorted(self.conda_envs))

    def test_set_metadata(self):

        conf = AnalysesConfig(**self.config_input)

        conf.meta_fp = self.meta_fp
        conf.set_metadata()
        assert_frame_equal(conf.meta, self.meta)

        conf.meta_fp = self.alt_header_fp
        conf.set_metadata()
        assert_frame_equal(conf.meta, self.meta)

        conf.meta_fp = self.nofile_fp
        with self.assertRaises(IOError) as e:
            conf.set_metadata()
        self.assertEqual('No file "%s"' % conf.meta_fp, str(e.exception))

        conf.meta_fp = self.only_header_fp
        with self.assertRaises(IOError) as e:
            conf.set_metadata()
        self.assertEqual('No sample in "%s"' % conf.meta_fp, str(e.exception))

        conf.meta_fp = self.empty_meta_fp
        with self.assertRaises(IOError) as e:
            conf.set_metadata()
        self.assertEqual('File "%s" is empty' % conf.meta_fp, str(e.exception))

    def test_set_coassembly(self):

        conf = AnalysesConfig(**self.config_input)

        conf.meta_fp = self.meta_fp

        conf.set_metadata()
        conf.set_coassembly()
        meta_coassembly = self.meta.copy()
        meta_coassembly['assembly_per_sample'] = meta_coassembly['sample_name']
        assert_frame_equal(conf.meta, meta_coassembly)

        conf.coassembly = {'wrong': {'no_a_var': []}}
        with self.assertRaises(IOError) as e:
            conf.set_coassembly()
        self.assertEqual(
            'Co-assembly variable "no_a_var" not in %s' % conf.meta_fp,
            str(e.exception))

        conf.coassembly = {'depth_coa': {'depth': [['a']]}}
        with self.assertRaises(IOError) as e:
            conf.set_coassembly()
        self.assertEqual(
            'Co-assembly factors for variable "depth" not in %s' % conf.meta_fp,
            str(e.exception))

        # conf.coassembly = {'depth_coa': {'depth': [['10'], ['2000']]}}
        conf.coassembly = {'depth_coa': {'depth': [['10'], ['2000']]}}
        conf.set_metadata()
        conf.set_coassembly()
        meta_coassembly = self.meta.copy()
        meta_coassembly['depth_coa'] = ['10', '2000', np.nan]
        assert_frame_equal(conf.meta, meta_coassembly)

        conf.coassembly = {'depth_coa': {'depth': [['10'], ['2000']],
                                         'light': [['yes'], ['no']]}}
        conf.set_metadata()
        conf.set_coassembly()
        meta_coassembly = self.meta.copy()
        meta_coassembly['depth_coa'] = ['10-yes', '2000-no', 'yes']
        assert_frame_equal(conf.meta, meta_coassembly)

        conf.coassembly = {'depth_coa': {'depth': [['2000']],
                                         'light': [['no']]}}
        conf.set_metadata()
        conf.set_coassembly()
        meta_coassembly = self.meta.copy()
        meta_coassembly['depth_coa'] = [np.nan, '2000-no', np.nan]
        assert_frame_equal(conf.meta, meta_coassembly)

        conf.meta_fp = self.no_coassembly_fp

        conf.coassembly = ''
        conf.set_metadata()
        conf.set_coassembly()
        meta_coassembly = self.meta_no_coassembly.copy()
        meta_coassembly['assembly_per_sample'] = meta_coassembly['sample_name']
        assert_frame_equal(conf.meta, meta_coassembly)

    def test_get_fastq_paths(self):

        conf = AnalysesConfig(**self.config_input)

        conf.meta_fp = self.meta_fp
        fastq_fps = conf.get_fastq_paths()
        self.assertEqual(fastq_fps, [])

        conf.fastq_dirs = ('%s/fastqs' % FOLDER,)
        fastq_fps = conf.get_fastq_paths()
        self.assertEqual(sorted(fastq_fps), sorted(self.fastqs1))

        conf.fastq_dirs = ('%s/fastqs' % FOLDER, '%s/output' % FOLDER)
        fastq_fps = sorted(conf.get_fastq_paths())
        self.assertEqual(fastq_fps, sorted(self.fastqs1 + self.fastqs2))

        conf.fastq_dirs = ('%s/fastqs_gz' % FOLDER,)
        fastq_fps = sorted(conf.get_fastq_paths())
        self.assertEqual(fastq_fps, sorted(self.fastqs3))

    def test_get_fastq_samples(self):

        conf = AnalysesConfig(**self.config_input)
        conf.meta_fp = self.meta_fp

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, {})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, self.fastqs)

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER, '%s/output' % FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, self.fastqs)

        conf.meta_fp = self.one_sample_meta_fp
        conf.set_metadata()
        conf.fastq_dirs = (FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, {})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, {'A': self.fastqs['A']})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER, '%s/output' % FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, {'A': self.fastqs['A']})

        conf.set_metadata()
        conf.fastq_dirs = ('%s/fastqs_gz' % FOLDER,)
        conf.get_fastq_samples()
        self.assertEqual(conf.fastq, {'A': self.fastqs3_only_gz})

    def test_set_fastq(self):

        conf = AnalysesConfig(**self.config_input)

        conf.meta_fp = self.meta_fp
        conf.set_metadata()
        conf.fastq_dirs = (FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, {})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, self.fastqs)

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER, '%s/output' % FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, self.fastqs)

        conf.meta_fp = self.one_sample_meta_fp
        conf.set_metadata()
        conf.fastq_dirs = (FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, {})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, {'A': self.fastqs['A']})

        conf.set_metadata()
        conf.fastq_dirs = (FOLDER, '%s/fastqs' % FOLDER, '%s/output' % FOLDER,)
        conf.set_fastq()
        self.assertEqual(conf.fastq, {'A': self.fastqs['A']})

        conf.set_metadata()
        conf.fastq_dirs = ('%s/inexistant' % FOLDER, '%s/nofolder' % FOLDER,)
        with self.assertRaises(IOError) as e:
            conf.set_fastq()
        self.assertEqual('Fastq folder(s) do not exist', str(e.exception))

    def test_get_r(self):

        conf = AnalysesConfig(**self.config_input)
        conf.fastq = {'A': ['A_1.fastq', 'A_2.fastq'],
                      'B': ['B_1.fastq.gz', 'B_2.fastq.gz']}
        conf.get_r()
        self.assertEqual(conf.r, {'A': ['_1', '_2'], 'B': ['_1', '_2']})

        conf = AnalysesConfig(**self.config_input)
        conf.fastq = {'A': ['A_1.fastq', 'A_R2.fastq'],
                      'B': ['B_1.fastq.gz', 'B_2.fastq.gz']}
        with self.assertRaises(IOError) as e:
            conf.get_r()
        self.assertEqual('Too different fastq names for sample "A":\n'
                         '- A_1.fastq\n- A_R2.fastq', str(e.exception))

    def test_update_metadata(self):

        conf = AnalysesConfig(**self.config_input)

        conf.fastq = {'A': ['A_1.fastq', 'A_2.fastq'],
                      'B': ['B_1.fastq.gz', 'B_2.fastq.gz']}
        conf.meta_fp = self.meta_fp
        conf.set_metadata()
        conf.update_metadata()
        assert_frame_equal(conf.meta, self.meta_fastqed)
        self.assertEqual(conf.meta_fp, self.meta_pipeline_fp)
        os.remove(self.meta_pipeline_fp)

        conf.fastq = {'A': ['A_1.fastq', 'A_2.fastq']}
        conf.meta_fp = self.meta_fp
        conf.set_metadata()
        conf.update_metadata()
        assert_frame_equal(conf.meta, self.meta_fastqed.loc[['A'], :])
        self.assertEqual(conf.meta_fp, self.meta_pipeline_fp)
        os.remove(self.meta_pipeline_fp)

    def test_set_output(self):

        conf = AnalysesConfig(**self.config_input)

        conf.output_dir = '%s/folder' % FOLDER
        conf.set_output()
        self.assertEqual(conf.output, '%s/folder' % FOLDER)
        self.assertTrue(isdir(conf.output))

        conf.output_dir = FOLDER
        conf.set_output()
        self.assertEqual(conf.output, FOLDER)
        self.assertTrue(isdir(conf.output))

    def test_parse_yamls(self):

        config_input = self.config_input
        config_input.update({'yaml1_yml': self.yaml1_yml})
        config_input.update({'yaml2_yml': self.yaml2_yml})
        conf = AnalysesConfig(**config_input)
        conf.parse_yamls()
        self.assertEqual(self.yaml1, conf.yaml1)
        self.assertEqual(self.yaml2, conf.yaml2)

        config_input = self.config_input
        config_input['pipeline_tsv'] = self.only_header_fp
        conf = AnalysesConfig(**config_input)
        conf.parse_yamls()
        self.assertEqual([['sample_name', 'depth', 'light']], conf.pipeline)

        config_input = self.config_input
        config_input['whatever_tsv'] = self.no_coassembly_fp
        conf = AnalysesConfig(**config_input)
        conf.parse_yamls()
        self.assertEqual([['sample_name', 'depth', 'no_coassembly'],
                          ['A', '10', 'gp1'],
                          ['B', '2000', 'gp2'],
                          ['X', '0', 'gp2']], conf.whatever)

    def test_get_params(self):

        conf = AnalysesConfig(**self.config_input)
        conf.parse_yamls()
        conf.get_default_params()
        self.assertEqual(self.params, conf.params)
        self.assertEqual({}, conf.user_params)

        conf = AnalysesConfig(**self.config_input)
        conf.user_params_yml = self.user_params_yml
        conf.parse_yamls()
        conf.get_default_params()
        self.assertEqual(self.params, conf.params)
        self.assertEqual({'shogun': {'time': 1}}, conf.user_params)

        conf = AnalysesConfig(**self.config_input)
        conf.user_params_yml = 'not_a_file.yml'
        with self.assertRaises(IOError) as e:
            conf.parse_yamls()
        self.assertEqual('No yaml file "not_a_file.yml"', str(e.exception))

    # def test_get_pooling_groups(self):
    #
    #     conf = AnalysesConfig(**self.config_input)
    #     conf.get_pooling_groups

    def tearDown(self) -> None:
        os.remove(self.meta_fp)
        os.remove(self.only_header_fp)
        os.remove(self.alt_header_fp)
        os.remove(self.empty_meta_fp)
        os.remove(self.no_coassembly_fp)
        os.remove(self.one_sample_meta_fp)
        os.remove(self.yaml1_yml)
        os.remove(self.yaml2_yml)
        for fastq in self.fastqs1:
            os.remove(fastq)
        for fastq in self.fastqs2:
            os.remove(fastq)
        for fastq in self.fastqs3:
            os.remove(fastq)


if __name__ == '__main__':
    unittest.main()
