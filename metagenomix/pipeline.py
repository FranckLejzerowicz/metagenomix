# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from collections import defaultdict
from metagenomix import parameters
from metagenomix.parameters import *


class Soft(object):

    def __init__(self, config):
        self.name = ''
        self.prev = None
        self.nodes = 1
        self.procs = 1
        self.scratch = 0
        self.io = {'I': {'f': set(), 'd': set(), 'D': set()},
                   'O': {'f': set(), 'd': set(), 'D': set()}}
        self.params = dict(config.params)
        self.inputs = {}
        self.outputs = {}
        self.cmds = {}
        self.dirs = set()

    def get_softs(self, softs):
        if len(softs) == 1:
            self.name = softs[0]
        else:
            self.prev, self.name = softs


class Graph(object):
    """https://www.geeksforgeeks.org/find-paths-given-source-destination/"""
    def __init__(self, vertices):
        self.V = vertices
        self.graph = defaultdict(list)
        self.paths = defaultdict(list)

    def add_edge(self, u, v):
        self.graph[u].append(v)

    def print_paths_util(self, u, d, visited, path):
        visited[u] = True
        path.append(u)
        if u == d:
            self.paths[path[-1]].append(list(path))
        else:
            for i in self.graph[u]:
                if not visited[i]:
                    self.print_paths_util(i, d, visited, path)
        path.pop()
        visited[u] = False

    def print_paths(self, s, d):
        visited = [False] * self.V
        path = []
        self.print_paths_util(s, d, visited, path)


class Workflow(object):
    """Collect the data associated with each dataset passed but the user
    """
    def __init__(self, config, databases) -> None:
        self.config = config
        self.databases = databases
        self.graph = None
        self.names = ['fastq']
        self.softs = {}
        self.names_idx = {}
        self.names_idx_rev = {}

    def init(self) -> None:
        """
        Parse the list of softwares, collect their sequence in the
        attribute `self.softs`, and collect all the paths of their
        sequences in the object `self.graph.paths`.
        """
        for softs in self.config.pipeline:
            self.collect_soft_name(softs)
            soft = Soft(self.config)
            soft.get_softs(softs)
            self.softs[softs[-1]] = soft
        self.get_names_idx()
        self.make_graph()
        self.get_paths()
        self.set_params()

    def collect_soft_name(self, softs: list) -> None:
        """Collect the sequential list of softwares.

        Parameters
        ----------
        softs : list
            One or two names of softwares that are run on
            the fastq files or after one another, respectively.
        """
        if softs[-1] not in self.names:
            if len(softs) > 1 and softs[0] not in self.names:
                raise IOError('"%s" not planned before "%s"' % tuple(softs))
            self.names.append(softs[-1])

    def get_names_idx(self) -> None:
        """
        Get a sequential numeric index (key) per software (value) in invoked
        order and the reversed, i.e., software (key) per numeric index (value).
        """
        self.names_idx = {y: x for x, y in enumerate(self.names)}
        self.names_idx_rev = {x: y for x, y in enumerate(self.names)}

    def make_graph(self) -> None:
        """
        Make the graph of all possible sequences of the
        list of softwares.
        """
        self.graph = Graph(len(self.names))
        for idx, soft in self.softs.items():
            if soft.prev:
                self.graph.add_edge(self.names_idx[soft.prev],
                                    self.names_idx[soft.name])
            else:
                self.graph.add_edge(0, self.names_idx[soft.name])

    def get_paths(self) -> None:
        """
        Collect all the paths of their sequences in the
        object `self.graph.paths`.
        """
        for soft in self.names[1:]:
            self.graph.print_paths(0, self.names_idx[soft])
        self.graph.paths = {self.names_idx_rev[idx]: [
            [self.names_idx_rev[p] for p in path] for path in paths]
            for idx, paths in self.graph.paths.items()}

    def make_dirs(self):
        for name, soft in self.softs.items():
            for directory in sorted(soft.dirs):
                if not isdir(directory):
                    os.makedirs(directory)

    def check_basic_params(self, soft_params, soft):
        ints = ['time', 'procs', 'mem_num', 'chunks']
        for param, value in soft_params.items():
            soft.params[param] = value
            if param in ints:
                check_ints(param, value, soft)
            elif param == 'mem_dim':
                check_mems(param, value, soft)
            elif param == 'env':
                check_env(self.config, value, soft)
            elif param == 'path':
                check_path(value)

    def set_scratch(self, soft):
        if self.config.localscratch:
            soft.params['scratch'] = self.config.localscratch
        elif self.config.scratch:
            soft.params['scratch'] = 'scratch'
        elif self.config.userscratch:
            soft.params['scratch'] = 'userscratch'

    def set_user_params(self, soft):
        user_params = self.config.user_params.get(soft.name, {})
        print()
        print("software:\t:\t", soft.name)
        # print("user_params\t:\t", user_params)
        # print("soft.params\t:\t", soft.params)
        func = 'check_%s' % soft.name
        if hasattr(parameters, func) and callable(getattr(parameters, func)):
            check_ = getattr(parameters, func)
            check_(user_params, soft, self.databases, self.config)
        self.check_basic_params(user_params, soft)

    def set_params(self) -> None:
        """
        Update the default params assigned to each software with the
        params passed by the user for each of the software.
        """
        for _, soft in self.softs.items():
            self.set_user_params(soft)
            self.set_scratch(soft)
            print("soft.params")
            print(soft.params)


    # def collect_paths(self):
    #
    #     if soft.soft_prev == 'qiita_wol':
    #         soft.output_paths[(None, 'qiita_wol')] = qiita_wol
    #         if soft == 'instrain':
    #             instrain_refs, instrain_bams = get_instrain_refs_bams(soft,
    #                                                                   soft_prev)
