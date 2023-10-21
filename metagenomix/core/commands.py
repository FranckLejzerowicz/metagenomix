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
from metagenomix.softwares.midas2 import *
from metagenomix.softwares.simka import *
from metagenomix.softwares.strains import *
from metagenomix.softwares.viruses import *
from metagenomix.softwares.mapping import *
from metagenomix.softwares.anvio import *
from metagenomix.softwares.squeezemeta import *


class Commands(object):

    def __init__(self, config, databases, workflow):
        self.config = config
        self.databases = databases
        self.graph = workflow.graph
        self.pipeline = workflow.workflow
        self.hashes = workflow.hashes
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
        self.path = []
        self.status = {}
        self.links = {}
        self.links_stats = {}
        self.holistics = {
            'simka',
            'quast',
            'qiita_wol',
            'mag_data',
            'metamarker',
            'woltka',
            'drep',
            'strainphlan',
            'midas2_merge'
        }

    def collect(self, log):
        for sdx, (name, paths) in enumerate(self.graph.paths.items()):
            log.info(' [%s] %s' % (sdx, name))
            for path in paths:
                self.path = path
                hashed = self.hashes[tuple(path)]
                self.soft = self.softs[name][hashed]
                if self.config.verbose:
                    log.info('   %s > %s' % (' ' * len(str(sdx)), hashed))
                self.get_inputs(log)  # Update self.inputs to previous output
                self.get_dir()
                self.make_holistic()
                self.generic_command()
                self.show_messages(log)

    def get_inputs(self, log):
        """Update the `inputs` attribute of the software object."""
        if self.soft.prev == 'None':
            # if running on raw data: use the fastq files (possibly on scratch)
            self.inputs = self.config.fastq
            if self.soft.params['scratch'] and self.config.jobs:
                self.inputs = self.config.fastq_mv
        else:
            # if the previous software is not None (i.e., not on the raw data)
            prev_hash = self.hashes[tuple(self.path[:-1])]
            self.inputs = self.softs[self.soft.prev][prev_hash].outputs
        show_inputs(self, log)

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

    def make_holistic(self):
        if self.soft.name in ['abritamr', 'diting']:
            if self.soft.params['samples'] == 'all':
                self.holistics.add(self.soft.name)

    def generic_command(self):
        self.sam_pool = ''
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

    def show_messages(self, log):
        for message in self.soft.messages:
            log.info('[%s] %s' % (self.soft.name, message))

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

    def unpack_outputs(self):
        if self.soft.name in self.holistics:
            self.soft.outputs = self.outputs['outs']
        else:
            outputs = dict(x for x in self.outputs['outs'].items() if x[1])
            self.soft.outputs[self.sam_pool] = outputs

    def extract_data(self):
        if self.outputs.get('cmds'):
            self.unpack_cmds()
            self.fill_soft_io()
        self.unpack_outputs()
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
            raise ValueError('No method for software "%s"' % self.soft.name)

    def register_command(self):
        self.softs[self.soft.name][self.soft.hashed].cmds = dict(self.cmds)
        self.cmds = {}
