# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from os.path import abspath

from metagenomix._io_utils import compute_hash
from metagenomix._inputs import show_inputs
from metagenomix.softwares.alignment import *
from metagenomix.softwares.annotation import *
from metagenomix.softwares.args import *
from metagenomix.softwares.assembly import *
from metagenomix.softwares.binning import *
from metagenomix.softwares.metamarker import *
from metagenomix.softwares.genomics import *
from metagenomix.softwares.phlans import *
from metagenomix.softwares.plasmids import *
from metagenomix.softwares.pooling import pooling
from metagenomix.softwares.preprocess import *
from metagenomix.softwares.profiling import *
from metagenomix.softwares.simka import *
from metagenomix.softwares.strains import *
from metagenomix.softwares.viruses import *
# from metagenomix.softwares.mapping import mapping


class Commands(object):

    def __init__(self, config, databases, workflow):
        self.config = config
        self.databases = databases
        self.graph = workflow.graph
        self.softs = workflow.softs
        self.outputs = {}
        self.cmds = {}
        self.args = {}
        self.pools = {}
        self.longs = None
        self.inputs = None
        self.method = None
        self.soft = None
        self.sam_pool = None
        self.dir = ''
        self.out = []
        self.struc = list
        self.links = {}
        self.links_stats = {}
        self.holistics = [
            'antismash',
            'simka',
            'quast',
            'qiita_wol',
            'mag_data',
            'metamarker',
            'woltka',
            'drep',
            'strainphlan'
        ]

    def collect(self):
        for sdx, softs in enumerate(self.config.pipeline):
            print(' - %s' % softs[-1])
            self.soft = self.softs[softs[-1]]
            self.soft.add_soft_path(self.softs)
            self.get_inputs()
            self.get_hash()
            self.get_dir()
            self.generic_command()
            self.show_messages()

    def get_inputs(self):
        """Update the `inputs` attribute of the software object."""
        if not self.soft.prev or self.soft.name == 'map__drep':
            if self.soft.params['scratch'] and self.config.jobs:
                self.inputs = self.config.fastq_mv
            else:
                self.inputs = self.config.fastq
        else:
            self.inputs = self.softs[self.soft.prev].outputs
        show_inputs(self)

    def get_hash(self):
        avoid = {'time', 'nodes', 'mem', 'mem_dim', 'env', 'chunks',
                 'scratch', 'machine', 'partition'}
        params = dict(x for x in self.soft.params.items() if x[0] not in avoid)
        path = self.soft.path
        hashes = [self.softs[x].hash for x in self.soft.path[1:-1]]
        self.soft.hash = (params, path, hashes)
        self.soft.hashed = compute_hash(self.soft.hash)

    def get_dir(self):
        self.dir = abspath('%s/%s/after_%s_%s' % (
            self.config.dir, self.soft.name, self.soft.prev, self.soft.hashed))
        self.soft.dir = self.dir
        if self.soft.params['scratch'] and self.config.jobs:
            self.dir = '${SCRATCH_FOLDER}%s' % self.dir

    def is_pool(self):
        self.struc = list
        if set(self.inputs) == set(self.pools) or self.soft.name == 'pooling':
            self.struc = dict

    def generic_command(self):
        self.sam_pool = ''
        self.is_pool()
        self.soft.io = {}
        if self.soft.name in self.holistics:
            self.prep_job()
        elif self.soft.name == 'pooling':
            self.pooling()
        else:
            for sam_or_pool in sorted(self.inputs):
                self.sam_pool = sam_or_pool
                self.prep_job()
        self.register_command()

    def show_messages(self):
        for message in self.soft.messages:
            print('[%s] %s' % (self.soft.name, message))

    def update_dirs(self):
        self.soft.dirs.update(set([x.replace('${SCRATCH_FOLDER}/', '/')
                                   for x in self.outputs['dirs']]))

    def init_outputs(self):
        self.outputs = {'cmds': {}, 'outs': {}, 'dirs': [], 'bash': [],
                        'io': {('I', 'd'): {}, ('I', 'f'): {},
                               ('O', 'd'): {}, ('O', 'f'): {}}}

    def init_io(self, key):
        if key not in self.soft.io:
            self.soft.io[key] = {}

    def fill_soft_io(self):
        for i, j in itertools.product(*[['I', 'O'], ['d', 'f']]):
            for key, io in self.outputs['io'].get((i, j), {}).items():
                self.init_io((self.sam_pool, key))
                self.soft.io[(self.sam_pool, key)][(i, j)] = io

    def unpack_cmds(self):
        for tech, cmds in self.outputs['cmds'].items():
            self.cmds[(self.sam_pool, tech)] = cmds

    def extract_data(self):
        if self.outputs.get('cmds'):
            self.unpack_cmds()
            self.fill_soft_io()
        if self.soft.name in self.holistics:
            self.soft.outputs = self.outputs['outs']
        else:
            self.soft.outputs[self.sam_pool] = self.outputs['outs']
        self.soft.bash = self.outputs['bash']

    def prep_job(self):
        self.init_outputs()     # initialize data structure collection
        self.call_method()      # collect commands, outputs, io, dirs
        self.extract_data()     # fill the useful self.soft attributes
        self.update_dirs()

    def pooling(self):
        for pool in self.config.pooling_groups:
            self.pools[pool] = {}
            self.soft.outputs[pool] = {}
            pooling(self, pool)

    def call_method(self):
        """Call the command-preparing method from this class (for the
        softwares that are easy to deal with), or from auxillary modules
        located in the softwares submodules path."""
        names = self.soft.name.split('_')
        name = names[0]
        if name in globals():
            globals()[name](self)
        else:
            raise ValueError('No method for software %s' % self.soft.name)

    def register_command(self):
        self.softs[self.soft.name].cmds = dict(self.cmds)
        self.cmds = {}
