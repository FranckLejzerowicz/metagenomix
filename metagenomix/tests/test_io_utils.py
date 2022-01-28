# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from unittest.mock import patch

import os
import subprocess

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

import pkg_resources
from os.path import isdir, isfile, splitext
from metagenomix._io_utils import (
    read_yaml, get_fastq_header, get_cat_zcat, edit_fastq_cmd, mkdr, get_chunks,
    get_pfam_file, get_hmm_dat, get_pfams_cmd, reads_lines)

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')
FOLDER = pkg_resources.resource_filename('metagenomix', 'tests/unittests')


class TestIOUtils(unittest.TestCase):

    def setUp(self) -> None:

        self.fastq_fp = '%s/fastqs/file.fastq' % FOLDER
        self.fastq_gz_fp = '%s.fastq.gz' % splitext(self.fastq_fp)[0]
        o = open(self.fastq_fp, 'w')
        o.write("@A.1.1 1 length=10\n")
        o.write("CTCCTGTGCC\n")
        o.write("+A.1.1 1 length=10\n")
        o.write("GGGGGGIIII\n")
        o.write("@A.2.1 2 length=10\n")
        o.write("GTCAGAGATA\n")
        o.write("+A.2.1 2 length=10\n")
        o.write("AGAAAGGGGG\n")
        o.close()
        subprocess.call(('gzip --keep %s' % self.fastq_fp).split())

        self.to_format = '%s/to_format.dat' % FOLDER
        self.formatted = '%s/formatted.tsv' % FOLDER
        o = open(self.to_format, 'w')
        o.write('# WHATEVER\n')
        o.write('#=GF AA  X\n')
        o.write('#=GF BB  X\n')
        o.write('//\n')
        o.write('# SOEVER\n')
        o.write('#=GF BB  Y\n')
        o.write('#=GF CC  Y\n')
        o.write('//\n')
        o.close()
        subprocess.call(('gzip --keep %s' % self.to_format).split())

        self.empty_fp = '%s/empty.tsv' % FOLDER
        o = open(self.empty_fp, 'w')
        o.close()

    def test_read_yaml(self):
        self.assertEqual({}, read_yaml(None))
        self.assertEqual({}, read_yaml(''))

        res = {'default': {'time': "48", 'nodes': "1", 'cpus': "1",
                           'mem_num': "5", 'mem_dim': "gb", 'env': "mg"}}
        file_path = '%s/run_params.yml' % RESOURCES
        self.assertEqual(res, read_yaml(file_path))

        file_path = 'not_a_file.yml'
        with self.assertRaises(IOError) as exc:
            read_yaml(file_path)
        self.assertEqual('No yaml file "%s"' % file_path, str(exc.exception))

    def test_get_fastq_header(self):
        first_line = get_fastq_header(self.fastq_fp)
        self.assertEqual(first_line, "@A.1.1 1 length=10")
        first_line = get_fastq_header(self.fastq_gz_fp)
        self.assertEqual(first_line, "@A.1.1 1 length=10")

    def test_get_cat_zcat(self):
        cat = get_cat_zcat(self.fastq_fp)
        self.assertEqual(cat, 'cat')
        cat = get_cat_zcat(self.fastq_gz_fp)
        self.assertEqual(cat, 'gzcat')

    def test_edit_fastq_cmd(self):
        ref_cmd = 'cat %s | ' % self.fastq_fp
        ref_cmd += "awk '{ if (NR%4==1) { print $1\"/1\" } "
        ref_cmd += "else if (NR%2 == 1) { print \"+\" } "
        ref_cmd += "else { print } }' | gzip > %s_renamed\n" % self.fastq_fp
        ref_cmd += "mv %s_renamed %s.gz\n" % (self.fastq_fp, self.fastq_fp)
        cmd = edit_fastq_cmd(self.fastq_fp, 1, 'illumina')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = ref_cmd.replace("/1", "/2")
        cmd = edit_fastq_cmd(self.fastq_fp, 2, 'illumina')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = 'gzcat %s | ' % self.fastq_gz_fp
        ref_cmd += "awk '{ if (NR%4==1) { print $1\"/1\" } "
        ref_cmd += "else if (NR%2 == 1) { print \"+\" } "
        ref_cmd += "else { print } }' | gzip > %s_renamed\n" % self.fastq_gz_fp
        ref_cmd += "mv %s_renamed %s\n" % (self.fastq_gz_fp, self.fastq_gz_fp)
        cmd = edit_fastq_cmd(self.fastq_gz_fp, 1, 'illumina')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = ref_cmd.replace("/1", "/2")
        cmd = edit_fastq_cmd(self.fastq_gz_fp, 2, 'illumina')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = 'cat %s | ' % self.fastq_fp
        ref_cmd += "awk '{ if (NR%4==1) { gsub(\".1$\",\"/1\",$1) "
        ref_cmd += "else if (NR%2 == 1) { print \"+\" } "
        ref_cmd += "else { print } }' | gzip > %s_renamed\n" % self.fastq_fp
        ref_cmd += "mv %s_renamed %s.gz\n" % (self.fastq_fp, self.fastq_fp)
        cmd = edit_fastq_cmd(self.fastq_fp, 1, 'ebi')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = ref_cmd.replace(".1$", ".2$").replace("/1", "/2")
        cmd = edit_fastq_cmd(self.fastq_fp, 2, 'ebi')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = 'gzcat %s | ' % self.fastq_gz_fp
        ref_cmd += "awk '{ if (NR%4==1) { gsub(\".1$\",\"/1\",$1) "
        ref_cmd += "else if (NR%2 == 1) { print \"+\" } "
        ref_cmd += "else { print } }' | gzip > %s_renamed\n" % self.fastq_gz_fp
        ref_cmd += "mv %s_renamed %s\n" % (self.fastq_gz_fp, self.fastq_gz_fp)
        cmd = edit_fastq_cmd(self.fastq_gz_fp, 1, 'ebi')
        self.assertEqual(ref_cmd, cmd)

        ref_cmd = ref_cmd.replace(".1$", ".2$").replace("/1", "/2")
        cmd = edit_fastq_cmd(self.fastq_gz_fp, 2, 'ebi')
        self.assertEqual(ref_cmd, cmd)

    def test_mkdr(self):
        folder = '%s/test_folder' % FOLDER
        self.assertFalse(isdir(folder))
        mkdr(folder)
        self.assertTrue(isdir(folder))
        os.rmdir(folder)
        self.assertFalse(isdir(folder))

    def test_get_chunks(self):
        ref_chunks = [[1, 2], [3, 4], [5]]
        chunks = get_chunks([1, 2, 3, 4, 5], 2)
        self.assertEqual(chunks, ref_chunks)

        ref_chunks = [[1, 2, 3], [4, 5]]
        chunks = get_chunks([1, 2, 3, 4, 5], 0, 3)
        self.assertEqual(chunks, ref_chunks)

        ref_chunks = [[1, 2, 3, 4], [5]]
        chunks = get_chunks([1, 2, 3, 4, 5], 0, 4)
        self.assertEqual(chunks, ref_chunks)

        ref_chunks = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        chunks = get_chunks([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 4)
        self.assertEqual(chunks, ref_chunks)

        ref_chunks = [[1, 2, 3, 4], [5, 6, 7, 8], [9]]
        chunks = get_chunks([1, 2, 3, 4, 5, 6, 7, 8, 9], 0, 4)
        self.assertEqual(chunks, ref_chunks)

        ref_chunks = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        chunks = get_chunks([1, 2, 3, 4, 5, 6, 7, 8, 9], 0, 3)
        self.assertEqual(chunks, ref_chunks)

    @patch('builtins.input', return_value='no')
    def test_get_pfam_file(self, input):
        ret = get_pfam_file(self.fastq_fp)
        self.assertTrue(ret)

    def test_get_hmm_dat(self):
        ret = pd.DataFrame({'AA': ['X', np.nan],
                            'BB': ['X', 'Y'],
                            'CC': [np.nan, 'Y']})
        ret_pd = get_hmm_dat('%s.gz' % self.to_format, self.formatted)
        assert_frame_equal(ret, ret_pd)
        self.assertTrue(isfile(self.formatted))

        ret_pd = get_hmm_dat('%s.gz' % self.to_format, self.formatted)
        assert_frame_equal(ret, ret_pd)
        os.remove(self.formatted)

    def test_get_pfams_cmd(self):

        self.assertFalse(isdir('./n'))
        pfam_key_pd = pd.DataFrame()
        ret_d, ret_cmd = get_pfams_cmd('HMM', pfam_key_pd, '.', 'n')
        self.assertTrue(isdir('./n'))
        self.assertEqual({}, ret_d)
        self.assertEqual('', ret_cmd)

        pfam_key_pd = pd.DataFrame({'AC': ['A', 'B'],
                                    'DE': ['1', '2']})
        d = {'A__1': ['./n/A__1.hmm', './n/A__1.dmnd'],
             'B__2': ['./n/B__2.hmm', './n/B__2.dmnd']}
        hmm1 = 'hmmfetch HMM A >> ./n/A__1.hmm\n'
        dmnd1 = 'cp ./fastas/A.fa ./fastas/A__1.fa\n'
        dmnd1 += 'diamond makedb --in ./fastas/A__1.fa -d ./n/A__1.dmnd\n'
        hmm2 = 'hmmfetch HMM B >> ./n/B__2.hmm\n'
        dmnd2 = 'cp ./fastas/B.fa ./fastas/B__2.fa\n'
        dmnd2 += 'diamond makedb --in ./fastas/B__2.fa -d ./n/B__2.dmnd\n'
        ret_d, ret_cmd = get_pfams_cmd('HMM', pfam_key_pd, '.', 'n')
        self.assertEqual(d, ret_d)
        cmd = hmm1 + dmnd1 + hmm2 + dmnd2
        self.assertEqual(cmd, ret_cmd)

        arefiles = [('./n/A__1.hmm', hmm1), ('./n/B__2.hmm', hmm2),
                    ('./n/A__1.dmnd', dmnd1), ('./n/B__2.dmnd', dmnd2)]
        for (fp, c) in arefiles:
            o = open(fp, 'w')
            o.close()
            cmd = cmd.replace(c, '')
            ret_d, ret_cmd = get_pfams_cmd('HMM', pfam_key_pd, '.', 'n')
            self.assertEqual(cmd, ret_cmd)
        for (fp, _) in arefiles:
            os.remove(fp)
        os.rmdir('./n')

    def test_reads_lines(self):

        ret = reads_lines(self.empty_fp)
        self.assertEqual(ret, [])
        ret = reads_lines(self.to_format)
        self.assertEqual(ret, ['# WHATEVER', '#=GF AA  X', '#=GF BB  X', '//',
                               '# SOEVER', '#=GF BB  Y', '#=GF CC  Y', '//'])

    def tearDown(self) -> None:
        os.remove(self.fastq_fp)
        os.remove(self.fastq_gz_fp)
        os.remove(self.to_format)
        os.remove('%s.gz' % self.to_format)


if __name__ == '__main__':
    unittest.main()
