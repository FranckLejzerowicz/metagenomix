# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict


class Soft(object):

    def __init__(self, config):
        self.name = ''
        self.prev = None
        self.nodes = 1
        self.procs = 1
        self.scratch = 0
        self.io = {'I': {'f': [], 'd': [], 'D': []},
                   'O': {'f': [], 'd': [], 'D': []}}
        self.params = dict(config.params)
        self.inputs = {}
        self.outputs = {}

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
    def __init__(self, config) -> None:
        self.config = config
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
        self.get_params()

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
                raise IOError('Software "%s" not planned before "%s"' % softs)
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

    def get_params(self) -> None:
        """
        Update the default params assigned to each software with the
        params passed by the user for each of the software.
        """
        mem_dim = ['kb', 'mb', 'gb']
        integer_params = ['time', 'procs', 'mem_num', 'chunks']
        for _, soft in self.softs.items():
            soft_params = self.config.user_params.get(soft.name, {})
            for param, value in soft_params.items():
                soft.params[param] = value
                if param in integer_params:
                    if not str(value).isdigit():
                        raise IOError('Param "%s" for "%s" must be integer' % (
                            param, soft.name))
                elif param == 'mem_dim':
                    if value not in mem_dim:
                        raise IOError('Param "%s" for "%s" must be of %s' % (
                            param, soft.name, str(mem_dim)))
                elif param == 'env':
                    if value in self.config.conda_envs:
                        raise EnvironmentError('"%s" not a conda env' % value)
                else:
                    raise IOError('"%s" is an unknown runparameter' % param)
                soft.params[param] = value

    # def collect_paths(self):
    #
    #     if soft.soft_prev == 'qiita_wol':
    #         soft.output_paths[(None, 'qiita_wol')] = qiita_wol
    #         if soft == 'instrain':
    #             instrain_refs, instrain_bams = get_instrain_refs_bams(soft,
    #                                                                   soft_prev)
