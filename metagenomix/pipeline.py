# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import yaml
from collections import defaultdict
from metagenomix import parameters
from metagenomix.parameters import *


class Soft(object):

    def __init__(self, config):
        self.name = ''
        self.prev = None
        self.scratch = None   # no use of the scratch file system by default
        self.params = dict(config.params)  # init with default params
        # key attributes to be filled by each tool-specific code
        self.io = {}
        self.inputs = {}
        self.outputs = {}
        self.cmds = {}
        self.dirs = set()
        self.defaults = {}

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
        self.skip = {}

    def run(self) -> None:
        """
        Parse the list of softwares, collect their sequence in the
        attribute `self.softs`, and collect all the paths of their
        sequences in the object `self.graph.paths`.
        """
        for softs in self.config.pipeline:
            self.collect_soft_name(softs)
            soft = Soft(self.config)
            soft.get_softs(softs)
            self.check_workflow(softs)
            self.softs[softs[-1]] = soft
        self.get_names_idx()
        self.make_graph()
        self.get_paths()
        self.set_params()

    def check_workflow(self, softs):
        if softs[-1] in self.softs:
            already_used = self.softs[softs[-1]]
            prev, name = already_used.prev, already_used.name
            message = 'Error in workflow "%s":\n' % self.config.pipeline_tsv
            message += "\tCan't run \"%s\" after \"%s\": " % tuple(softs[::-1])
            message += '"%s" already planned to use after "%s"\n' % (name, prev)
            message += "Each tool must yield a single output re-used as input\n"
            message += '-> please re-run for each different workflow "graph"\n'
            message += "Exiting\n"
            sys.exit(message)

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

    def check_basic_params(self, user_params, soft):
        ints = ['time', 'procs', 'mem_num', 'chunks']
        for param, value in user_params.items():
            if param in ints:
                check_int(param, value, soft.name)
            elif param == 'mem_dim':
                check_mems(param, value, soft.name)
            elif param == 'env':
                check_env(self.config, value, soft.name)
            elif param == 'path':
                check_path(self.config, value, soft.name)
            elif param == 'scratch':
                check_scratch(value, soft.name)
            soft.params[param] = value

    def ignored_params(self, user_params, soft):
        valid_params = set(soft.defaults) | set(self.config.params)
        for param in sorted(user_params):
            if param not in valid_params:
                name = soft.name
                print('[%s] Param "%s" unknown and ignored' % (name, param))
                del user_params[param]

    def set_scratch(self, soft):
        """scratch set on command line overrides per-software scratches"""
        if self.config.localscratch:
            soft.params['scratch'] = self.config.localscratch
        elif self.config.scratch:
            soft.params['scratch'] = 'scratch'
        elif self.config.userscratch:
            soft.params['scratch'] = 'userscratch'

    def set_user_params(self, soft):
        name = soft.name.split('_')[0]
        user_params = dict(self.config.user_params.get(name, {}))
        func = 'check_%s' % name
        if hasattr(parameters, func) and callable(getattr(parameters, func)):
            check_ = getattr(parameters, func)
            soft.defaults = check_(self, user_params, soft)
            self.ignored_params(user_params, soft)

        self.check_basic_params(user_params, soft)
        if isinstance(soft.params['scratch'], int):
            soft.params['mem_num'] = soft.params['scratch']
            soft.params['mem_dim'] = 'gb'

    def skip_params(self, soft) -> bool:
        if '_' in soft.name:
            if self.skip.get(soft.name.split('_')[0], False):
                return True
            self.skip[soft.name.split('_')[0]] = True

    def show_params(self, soft):
        params_show = dict(
            x for x in soft.params.items() if x[0] not in self.config.params)
        if params_show:
            print('[Parameters] %s' % soft.name)
            print('=' * 30)
            if soft.name.startswith('search'):
                databases = params_show['databases']
                del params_show['databases']
                params = {'search': {soft.name.split('_')[-1]: params_show,
                                     'databases': databases}}
                print(yaml.dump(params))
            else:
                print(yaml.dump(params_show))
            print('%s defaults %s' % ('-' * 10, '-' * 10))
            print(yaml.dump(soft.defaults))
            print('=' * 30)

    def set_params(self) -> None:
        """
        Update the default params assigned to each software with the
        params passed by the user for each of the software.
        """
        for _, soft in self.softs.items():
            # print()
            # print()
            # print('-----------------')
            # print(soft.name)
            self.set_scratch(soft)
            self.set_user_params(soft)
            # print('soft.params:', soft.params)
            # print('-----------------')
        if self.config.show_params:
            for _, soft in self.softs.items():
                self.show_params(soft)

    # def collect_paths(self):
    #
    #     if soft.soft_prev == 'qiita_wol':
    #         soft.output_paths[(None, 'qiita_wol')] = qiita_wol
    #         if soft == 'instrain':
    #             instrain_refs, instrain_bams = get_instrain_refs_bams(soft,
    #                                                                   soft_prev)
