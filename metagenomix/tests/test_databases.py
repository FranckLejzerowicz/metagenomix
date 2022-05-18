# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import yaml
import unittest
import pkg_resources
from os.path import isfile

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import io
from unittest.mock import patch

from metagenomix.databases import ReferenceDatabases
from metagenomix.config import AnalysesConfig

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')
FOLDER = pkg_resources.resource_filename('metagenomix', 'tests/unittests')


class TestDatabases(unittest.TestCase):

    def setUp(self) -> None:

        self.mar = '%s/databases/MAR' % FOLDER
        self.pfam = '%s/databases/PFAM' % FOLDER
        self.dbcan = '%s/databases/dbCAN' % FOLDER
        self.species = '%s/databases/dbCAN/species.txt' % FOLDER
        self.squeezemeta = '%s/databases/squeezemeta' % FOLDER
        self.wol = '%s/databases/wol' % FOLDER
        self.midas = '%s/databases/midas' % FOLDER
        self.metaphlan = '%s/databases/metaphlan' % FOLDER

        self.yaml = {'mar': {'path': self.mar},
                     'pfam': {'path': self.pfam, 'terms': ['oxin', 'anti']},
                     'dbcan': {'path': self.dbcan, 'taxa': {'A': self.species}},
                     'squeezemeta': {'path': self.squeezemeta},
                     'wol': {'path': self.wol},
                     'midas': {'path': self.midas},
                     'metaphlan': {'path': self.metaphlan}}
        self.databases_fp = '%s/yamls/databases.yml' % FOLDER
        with open(self.databases_fp, "w") as o:
            yaml.dump(self.yaml, o)
        self.config = {'databases_yml': self.databases_fp, 'midas_foci': {}}

        self.yaml_simple = {'A': {'path': 'path/to/A'},
                            'B': {'path': 'path/to/B'}}
        self.databases_simple_fp = '%s/yamls/databases_simple.yml' % FOLDER
        with open(self.databases_simple_fp, "w") as o:
            yaml.dump(self.yaml_simple, o)
        self.config_simple = {'databases_yml': self.databases_simple_fp}

        self.no_path = {'example1': {}, 'example2': {'something': 'here'}}
        self.no_path_fp = '%s/yamls/yaml_no_path.yml' % FOLDER
        with open(self.no_path_fp, "w") as o:
            yaml.dump(self.no_path, o)
        self.config_no_path = {'databases_yml': self.no_path_fp}

        self.no_config = {'databases_yml': None}

        self.paths = []
        self.configs = {}
        for database, database_d in self.yaml.items():
            self.cur_fp = '%s/yamls/databases_%s.yml' % (FOLDER, database)
            with open(self.cur_fp, "w") as o:
                yaml.dump({database: self.yaml[database]}, o)
            self.paths.append(self.cur_fp)
            self.configs[database] = {'databases_yml': self.cur_fp}

        self.meta_dir = '%s/dbCAN-seq' % self.dbcan
        self.cazyme_seq_list = '%s/CAZyme_seq_list' % self.meta_dir
        self.meta_fp = '%s/metadata.txt' % self.meta_dir
        self.meta = pd.DataFrame(
            [['GCF_000007365.1', 'Buchnera'], ['GCF_000021405.1', 'Borrelia'],
             ['GCF_001047655.1', 'Leptospira']], columns=['gen', 'genome_name'])
        self.meta.set_index('gen', inplace=True)

    @patch('sys.stdout', new_callable=io.StringIO)
    def assert_stdout(self, databases, expected, mock_stdout):
        databases.show_loaded_databases()
        self.assertEqual(mock_stdout.getvalue(), expected)

    def test_show_loaded_databases_prints(self):
        config = AnalysesConfig(**self.no_config)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        expected = 'No database passed using option `-d`\n'
        self.assert_stdout(databases, expected)

        config = AnalysesConfig(**self.config_no_path)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        expected = 'Databases in %s\n' % self.no_path_fp
        expected += '[Database] "example1" folder not found\n'
        expected += '[Database] "example2" folder not found\n'
        self.assert_stdout(databases, expected)

        config = AnalysesConfig(**self.config_simple)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        expected = 'Databases in %s\n' % self.databases_simple_fp
        expected += '- A: path/to/A\n'
        expected += '- B: path/to/B\n'
        self.assert_stdout(databases, expected)

        config = AnalysesConfig(**self.config)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        expected = 'Databases in %s\n' % self.databases_fp
        for k, v in sorted(self.yaml.items()):
            expected += '- %s: %s\n' % (k, v['path'])
        self.assert_stdout(databases, expected)

    def test_show_loaded_databases(self):

        config = AnalysesConfig(**self.no_config)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        res = databases.show_loaded_databases()
        self.assertFalse(res)

        config = AnalysesConfig(**self.config)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        res = databases.show_loaded_databases()
        self.assertTrue(res)

    def test_set_squeezemeta(self):
        config = AnalysesConfig(**self.config)
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.set_squeezemeta()
        self.assertEqual(databases.database, 'squeezemeta')

    def test_set_mar(self):
        pass

    def test_set_pfam(self):
        anti = 'PF12574.11__120_KDa_Rickettsia_surface_antigen'
        oxin = 'PF10417.12__C_terminal_domain_of_1_Cys_peroxiredoxin'
        hmms = {'oxin': {oxin: ['%s/oxin/%s.hmm' % (self.pfam, oxin),
                                '%s/oxin/%s.dmnd' % (self.pfam, oxin)]},
                'anti': {anti: ['%s/anti/%s.hmm' % (self.pfam, anti),
                                '%s/anti/%s.dmnd' % (self.pfam, anti)]}}
        hmms_pd = pd.DataFrame([
            ['1-cysPrx_C', 'PF10417.12', 'C-terminal domain of 1-Cys peroxiredoxin',
             '21.1; 21.1;', 'Domain', 40, np.nan, np.nan],
            ['120_Rick_ant', 'PF12574.11', '120 KDa Rickettsia surface antigen',
             '25; 25;', 'Family', 238, np.nan, np.nan]],
            columns=['ID', 'AC', 'DE', 'GA', 'TP', 'ML', 'CL', 'NE'])
        config = AnalysesConfig(**self.configs['pfam'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.set_pfam()
        self.assertEqual(databases.database, 'pfam')
        self.assertEqual(databases.hmms, hmms)
        assert_frame_equal(databases.hmms_pd, hmms_pd)

    def test_get_dbcan_hmms(self):
        expected = {}
        for hmm in ['a', 'b', 'c']:
            hmm_fp = '%s/%s.hmm' % (self.dbcan, hmm)
            with open(hmm_fp, 'w') as o: pass
            expected[hmm] = hmm_fp

        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.get_dbcan_hmms()
        self.assertEqual(databases.hmms, expected)

        for _, hmm_fp in expected.items():
            os.remove(hmm_fp)

        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.get_dbcan_hmms()
        self.assertEqual(databases.hmms, {})

    def test_write_dbcan_subset(self):
        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.dbcan_meta = self.meta

        folder = '%s/subset' % config.databases['dbcan']['path']
        os.makedirs(folder)
        os.makedirs(self.cazyme_seq_list)
        gcf_fas = '%s/GCF_000007365.1.fasta' % self.cazyme_seq_list
        with open(gcf_fas, 'w') as o:
            o.write('>header\nACGT\n')
        fas, dia = '%s/A.fa' % folder, '%s/A.dmnd' % folder
        expected = "diamond makedb --in %s -d %s\n" % (fas, dia)
        observed = databases.write_dbcan_subset(['Buchnera'], folder, 'A')
        self.assertEqual(observed, expected)
        self.assertTrue(isfile(fas))

        os.remove(fas)
        os.remove(gcf_fas)
        os.rmdir(self.cazyme_seq_list)
        os.rmdir(self.meta_dir)

        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        with open(dia, 'w') as o: pass
        observed = databases.write_dbcan_subset(['Buchnera'], folder, 'A')
        self.assertEqual(observed, '')
        os.remove(dia)
        os.rmdir(folder)

    def test_set_dbcan_taxa(self):
        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.set_dbcan_taxa()
        self.assertEqual(databases.cazys, {})

        # folder = '%s/subset' % config.databases['dbcan']['path']
        # with open(self.species, 'w') as o:
        #     o.write('a_species\n')
        # databases.set_dbcan_taxa()
        # self.assertEqual(databases.cazys, {'A': folder})
        # os.remove(self.species)
        # os.rmdir(folder)

    def test_set_dbcan(self):
        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.set_dbcan()
        self.assertEqual(databases.database, 'dbcan')
        self.assertEqual(databases.hmms, {})
        assert_frame_equal(databases.dbcan_meta, pd.DataFrame())

        os.makedirs(self.meta_dir)
        self.meta.to_csv(self.meta_fp, sep='\t')
        config = AnalysesConfig(**self.configs['dbcan'])
        config.parse_yamls()
        databases = ReferenceDatabases(config)
        databases.set_dbcan()
        assert_frame_equal(databases.dbcan_meta, self.meta)
        os.remove(self.meta_fp)
        os.rmdir(self.meta_dir)

        expected = {}
        for hmm in ['a', 'b', 'c']:
            hmm_fp = '%s/%s.hmm' % (self.dbcan, hmm)
            with open(hmm_fp, 'w') as o: pass
            expected[hmm] = hmm_fp
        databases.set_dbcan()
        self.assertEqual(databases.hmms, expected)
        for _, hmm_fp in expected.items():
            os.remove(hmm_fp)

    def test_set_wol(self):
        pass

    def test_set_midas(self):
        pass

    def test_set_humann(self):
        pass

    def test_set_metaphlan(self):
        pass

    def tearDown(self) -> None:
        os.remove(self.databases_fp)
        os.remove(self.databases_simple_fp)
        os.remove(self.no_path_fp)


if __name__ == '__main__':
    unittest.main()
