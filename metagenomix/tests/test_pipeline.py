# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from metagenomix.pipeline import Graph, Soft, Workflow


class TestGraph(unittest.TestCase):

    def test_add_edge(self):
        graph = Graph(3)
        pass

    def test_print_paths_util(self):
        graph = Graph(3)
        pass

    def test_print_paths(self):
        pass


class TestSoft(unittest.TestCase):

    def test_get_softs(self):

        soft = Soft()
        self.assertEqual(soft.prev, 'fastq')
        self.assertEqual(soft.name, '')

        soft.get_softs(['A', 'B'])
        self.assertEqual(soft.prev, 'A')
        self.assertEqual(soft.name, 'B')

        soft = Soft()
        soft.get_softs(['B'])
        self.assertEqual(soft.prev, 'fastq')
        self.assertEqual(soft.name, 'B')

        soft = Soft()
        soft.get_softs([[], []])
        self.assertEqual(soft.prev, [])
        self.assertEqual(soft.name, [])


class TestPipeline(unittest.TestCase):

    def setUp(self) -> None:
        # pipeline = Pipeline()
        pass

    def test_collect_soft_name(self):
        pass

    def test_get_names_idx(self):
        pass

    def test_make_graph(self):
        pass

    def test_get_paths(self):
        pass

    def test_init(self):
        pass

    def test_get_params(self):
        pass

    def test_generic_command(self):
        pass

    def test_get_input_path(self):
        pass

    def test_get_scratch(self):
        pass

    def test_run(self):
        pass


if __name__ == '__main__':
    unittest.main()
