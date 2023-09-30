# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import yaml
from collections import Counter
from metagenomix.core import parameters
from metagenomix.core.parameters import *
from metagenomix.core.graph import Graph
from metagenomix.core.software import Soft


class Workflow(object):
    """Collect the data associated with each dataset passed but the user
    """
    def __init__(self, config, databases) -> None:
        self.config = config
        self.databases = databases
        self.graph = None
        self.names = ['None']
        self.workflow = []
        self.softs = {}
        self.hashes = {}
        self.params = {}
        self.defaults = {}
        self.names_idx = {}
        self.names_idx_rev = {}
        self.skip = {}
        self.tools = {}
        self.name = None

    def check_workflow(self):
        counts = Counter([tuple(softs) for softs in self.config.pipeline])
        workflow_issue = False
        for softs, count in counts.items():
            if count > 1:
                workflow_issue = True
                print('[pipeline] step "%s" duplicated' % ' '.join(softs))
        if workflow_issue:
            sys.exit('[pipeline] Please fix duplicates... Exiting')

    def fill_workflow(self, step):
        if len(step) > 2:
            pipeline_tsv = self.config.pipeline_tsv
            sys.exit('[config] Max 2 names per step in "%s"' % pipeline_tsv)
        if len(step) == 1:
            self.workflow.append(('None', step[0]))
        else:
            self.workflow.append(tuple(step))

    def visit(self) -> None:
        """
        Parse the list of softwares, collect their sequence in the
        attribute `self.softs`, and collect all the paths of their
        sequences in the object `self.graph.paths`.
        """
        self.check_workflow()
        for step in self.config.pipeline:
            self.collect_step_names(step)
            self.fill_workflow(step)

    def setup(self) -> None:
        self.get_names_idx()
        self.make_graph()
        self.get_paths()
        self.show_graph()

    def collect_step_names(self, step: list) -> None:
        """Collect the sequential list of softwares.

        Parameters
        ----------
        step : list
            One or two names of softwares that are run on
            the fastq files or after one another, respectively.
        """
        if step[-1] not in self.names:
            if len(step) > 1 and step[0] not in self.names:
                raise IOError('"%s" not planned before "%s"' % tuple(step))
            self.names.append(step[-1])

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
        for prev, name in self.workflow:
            if prev == 'None':
                self.graph.add_edge(0, self.names_idx[name])
            else:
                self.graph.add_edge(self.names_idx[prev],
                                    self.names_idx[name])

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

    def show_graph(self):
        # if self.config.verbose:
        for i, js in self.graph.paths.items():
            print('\t%s:' % i)
            for j in js:
                print('\t -', '->'.join(j))

    def check_basic_params(self, user_params):
        ints = ['time', 'procs', 'mem', 'chunks']
        for param, value in user_params.items():
            if param in ints:
                check_int(param, value, self.name)
            elif param == 'mem_dim':
                check_mems(param, value, self.name)
            elif param == 'env':
                check_env(self.config, value, self.name)
            elif param == 'path':
                check_path(self.config, value, self.name)
            elif param == 'scratch':
                check_scratch(value, self.name)
            self.params[self.name][param] = value

    def ignored_params(self, user_params):
        # take union of defaults params known from software check and run_params
        valid_params = set(self.defaults[self.name]) | set(self.config.params)
        # then, for each if the user params
        for param in sorted(user_params):
            # if param never returned by checks (default) or in run_params
            if param not in valid_params:
                # print that it is unknown and stop
                sys.exit('[%s] Param "%s" unknown' % (self.name, param))

    def set_scratch(self):
        """scratch set on command line overrides per-software scratches"""
        if self.config.localscratch:
            self.params[self.name]['scratch'] = self.config.localscratch
        elif self.config.scratch:
            self.params[self.name]['scratch'] = 'scratch'
        elif self.config.userscratch:
            self.params[self.name]['scratch'] = 'userscratch'

    def get_user_params(self, name_, name):
        """The params set by user for a tool's subcommand take precedence over
        the params of the general tool. For example, if params are set both for:
         - `metawrap`
         - `metawrap_refine`
        then the values given under `metawrap_refine` take precedence.
        """
        # params of the general tool
        user_params = dict(self.config.user_params.get(name, {}))
        # params of the specific tool subcommand/module
        user_params.update(dict(self.config.user_params.get(name_, {})))
        return user_params

    def set_user_params(self):
        # get the name of the software ("metawrap_binning" would be "metawrap")
        nam = self.name.split('_')[0]
        # get the user parameters (looking up to class attribute value)
        user_params = self.get_user_params(self.name, nam)
        # get the name of the parameter-checking function for the software
        func = 'check_%s' % nam
        # if this parameter-checking function exists in the "parameters" module
        if hasattr(parameters, func) and callable(getattr(parameters, func)):
            # get the function as an object
            check_ = getattr(parameters, func)
            # run this function to check that all parameters are valid
            # and get the returned value as default (for --show-params etc)
            self.defaults[self.name] = check_(self, user_params)
            # if here, params were valids, then check that the params are known
            self.ignored_params(user_params)
        self.check_basic_params(user_params)

    @staticmethod
    def get_params_dict(name, params_show):
        if name.startswith('search'):
            databases = params_show['databases']
            del params_show['databases']
            params = {'search': {name.split('_')[-1]: params_show,
                                 'databases': databases}}
        else:
            params = params_show
        return params

    def print_params(self, soft):
        params_show = dict(x for x in self.params[self.name].items())
        if params_show:
            x = '=' * (13 + len(self.name))
            print('\n%s\n[%s] Parameters\n%s' % (x, self.name, x))
            params = self.get_params_dict(self.name, params_show)
            print(yaml.dump(params))
            print('%s defaults %s' % ('-' * 10, '-' * 10))
            print(yaml.dump(self.defaults[self.name]))
            print('=' * 30)

    def write_params(self):
        params_show = dict(x for x in self.params[self.name].items())
        params_dict = self.get_params_dict(self.name, params_show)
        params = {self.name: params_dict}
        p = '/Users/franck/programs/metagenomix/metagenomix/resources/params'
        with open('%s/%s.yml' % (p, self.name), 'w') as o:
            yaml.dump(params, o)

    def parametrize(self) -> None:
        """
        Update the default params assigned to each software with the
        params passed by the user for each of the software.
        """
        for name in self.graph.paths.keys():
            self.name = name
            self.params[name] = {}
            self.defaults[name] = {}
            self.set_scratch()
            self.set_user_params()
            # self.write_params()

    def show_params(self) -> None:
        print('* Showing parameters (user-defined and defaults):')
        for _, soft in self.softs.items():
            if 'all' in self.config.show_params:
                self.print_params(soft)
            else:
                for s in self.config.show_params:
                    if re.search(s, soft.name):
                        self.print_params(soft)
                        break

    def prepare(self) -> None:
        """
        Parse the list of softwares, collect their sequence in the
        attribute `self.softs`, and collect all the paths of their
        sequences in the object `self.graph.paths`.
        """
        softs = {}
        for name, paths in self.graph.paths.items():
            if name not in self.softs:
                self.softs[name] = {}
            for path in paths:
                soft = Soft(self.config)
                soft.set_soft(self.params, path)
                soft.get_hash(self.params[name], softs)
                self.softs[name][soft.hashed] = soft
                self.hashes[tuple(path)] = soft.hashed
