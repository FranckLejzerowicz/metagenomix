# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest

import pandas as pd
import pkg_resources
from pandas.testing import assert_frame_equal, assert_series_equal

from metagenomix._metadata import read_metadata, get_cur_type, get_first_column

FOLDER = pkg_resources.resource_filename('metagenomix', 'tests/unittests')


class ReadMetadata(unittest.TestCase):

    def setUp(self) -> None:

        self.empty_meta_fp = '%s/metadata/empty_meta.txt' % FOLDER
        with open(self.empty_meta_fp, 'w'):
            pass

        self.only_header_fp = '%s/metadata/only_header.txt' % FOLDER
        with open(self.only_header_fp, 'w') as o:
            o.write('sample_name\tdepth\tlight')
        self.only_header = pd.DataFrame(
            {'sample_name': [], 'depth': [], 'light': []})

        self.meta_fp = '%s/metadata/meta.txt' % FOLDER
        with open(self.meta_fp, 'w') as o:
            o.write('sample_name\tdepth\tlight\n')
            o.write('A\t10\tyes\nB\t2000\tno\nX\t0\tyes')
        self.meta = pd.DataFrame({'sample_name': ['A', 'B', 'X'],
                                  'depth': ['10', '2000', '0'],
                                  'light': ['yes', 'no', 'yes']})

        self.alt_header_fp = '%s/metadata/alt_header.txt' % FOLDER
        with open(self.alt_header_fp, 'w') as o:
            o.write('some_name\tdepth\tlight\n')
            o.write('A\t10\tyes\nB\t2000\tno\nX\t0\tyes')
        self.alt_header = self.meta.copy().rename(
            columns={'sample_name': 'some_name'})

    def test_read_metadata(self):
        with self.assertRaises(IOError) as e:
            read_metadata(self.empty_meta_fp)
        self.assertEqual(
            'File "%s" is empty' % self.empty_meta_fp, str(e.exception))

        meta = read_metadata(self.only_header_fp)
        self.assertEqual(str(meta), str(self.only_header))
        self.assertEqual(meta.values.tolist(), self.only_header.values.tolist())

        meta = read_metadata(self.meta_fp)
        assert_frame_equal(meta, self.meta)

        meta = read_metadata(self.alt_header_fp)
        assert_frame_equal(meta, self.meta)

    def test_get_first_column(self):
        first_column = get_first_column(self.meta_fp)
        self.assertEqual(first_column, 'sample_name')
        first_column = get_first_column(self.alt_header_fp)
        self.assertEqual(first_column, 'some_name')
        first_column = get_first_column(self.only_header_fp)
        self.assertEqual(first_column, 'sample_name')
        with self.assertRaises(IOError) as e:
            get_first_column(self.empty_meta_fp)
        self.assertEqual(
            'File "%s" is empty' % self.empty_meta_fp, str(e.exception))

    def tearDown(self) -> None:
        os.remove(self.meta_fp)
        os.remove(self.only_header_fp)
        os.remove(self.alt_header_fp)
        os.remove(self.empty_meta_fp)


class CheckType(unittest.TestCase):

    def setUp(self) -> None:

        self.short_fp = '%s/short.fa' % FOLDER
        self.nucl_fp = '%s/nucl.fa' % FOLDER
        self.prot_fp = '%s/prot.fa' % FOLDER
        self.neither_fp = '%s/neither.fa' % FOLDER
        with open(self.short_fp, 'w') as o:
            o.write('>header\n')
            o.write('ACTGTGATTG\n')

        with open(self.nucl_fp, 'w') as o:
            o.write('>header\n')
            o.write('ACTGTGATTGCANTAN\n')

        with open(self.prot_fp, 'w') as o:
            o.write('>header\n')
            o.write('POSYYRTWBNHCCRCG\n')

        with open(self.neither_fp, 'w') as o:
            o.write('Nothing\n')
            o.write('Special\n')

    def test_get_cur_type(self):
        ret = get_cur_type(self.short_fp)
        self.assertEqual(ret, 'nucl')
        ret = get_cur_type(self.nucl_fp)
        self.assertEqual(ret, 'nucl')
        ret = get_cur_type(self.prot_fp)
        self.assertEqual(ret, 'prot')
        ret = get_cur_type(self.neither_fp)
        self.assertEqual(ret, 'prot')

    def tearDown(self) -> None:
        os.remove(self.short_fp)
        os.remove(self.nucl_fp)
        os.remove(self.prot_fp)
        os.remove(self.neither_fp)


if __name__ == '__main__':
    unittest.main()
