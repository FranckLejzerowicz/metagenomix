# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import yaml
import itertools
from os.path import abspath

from metagenomix._inputs import show_inputs
from metagenomix.tools.pooling import pooling
from metagenomix.tools.preprocess import *
from metagenomix.tools.alignment import *
from metagenomix.tools.simka import *
from metagenomix.tools.profiling import *
from metagenomix.tools.phlans import *
from metagenomix.tools.assembly import *
from metagenomix.tools.annotation import *
from metagenomix.tools.binning import *
from metagenomix.tools.plasmids import *
from metagenomix.tools.viruses import *
from metagenomix.tools.args import *
from metagenomix.tools.genomics import *
from metagenomix.tools.metamarker import *


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
        # self.pool = None
        self.longs = None
        self.inputs = None
        self.method = None
        self.soft = None
        # self.sam = None
        self.sam_pool = None
        self.dir = ''
        self.out = []
        self.struc = list
        self.holistics = [
            'antismash'
            'multiqc',
            'simka',
            'quast'
            'qiita_wol',
            'mag_data',
            'metamarker',
            'woltka',
            'drep',
            'strainphlan'
        ]

    def run(self):
        for sdx, softs in enumerate(self.config.pipeline):
            # print()
            # print('*' * 30)
            # print('>>>', softs)
            # print('*' * 30)
            # print()
            self.soft = self.softs[softs[-1]]
            self.get_inputs()
            # print()
            # print('-' * 100)
            # print(self.inputs)
            # print('-' * 100)
            self.get_dir()
            self.generic_command()
            # print(' * ' * 75)
            # print(yaml.dump(self.softs[self.soft.name].cmds))
            # print(' * ' * 75)
            # print('- ' * 50)
            # print(self.soft.outputs)
            # print('- ' * 50)

    def make_dirs(self):
        for name, soft in self.softs.items():
            for directory in sorted(soft.dirs):
                if not isdir(directory):
                    os.makedirs(directory)

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

    def get_dir(self):
        self.dir = abspath('%s/%s/after_%s' % (
            self.config.dir, self.soft.name, self.soft.prev))
        self.soft.dirs.add(self.dir)
        if self.soft.params['scratch'] and self.config.jobs:
            self.dir = '${SCRATCH_FOLDER}%s' % self.dir

    def is_pool(self):
        self.struc = list
        if set(self.inputs) == set(self.pools) or self.soft.name == 'pooling':
            self.struc = dict

    def generic_command(self):
        self.sam_pool = ''
        # self.sam = ''
        self.is_pool()
        self.soft.io = {}
        if self.soft.name in self.holistics:
            self.prep_job()
        elif self.soft.name == 'pooling':
            self.pooling()
        else:
            for sam_or_pool in sorted(self.inputs):
                self.sam_pool = sam_or_pool
                # self.sam = sam_or_pool
                # self.pool = sam_or_pool
                self.prep_job()
        self.register_command()

    def update_dirs(self):
        self.soft.dirs.update(set([x.replace('${SCRATCH_FOLDER}/', '/')
                                   for x in self.outputs['dirs']]))

    def init_outputs(self):
        self.outputs = {'cmds': {}, 'outs': {}, 'dirs': [],
                        'io': {('I', 'd'): {}, ('I', 'f'): {},
                               ('O', 'd'): {}, ('O', 'f'): {}}}

    def call_method(self):
        """Call the command-preparing method from this class (for the tools that
        are easy to deal with), or from auxillary modules located in the `tools`
        submodules path."""
        func = self.soft.name.split('_')[0]
        if func in globals():
            globals()[func](self)
        # ---------------------- TEMPORARY -----------------------
        # this is for the functions that might be defined here
        elif func in dir(self) and callable(getattr(self, func)):
            getattr(self, func)()
        else:
            raise ValueError('No method for software %s' % self.soft.name)

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

    def prep_job(self):
        self.init_outputs()     # initialize data structure collection
        self.call_method()      # collect commands, outputs, io, dirs
        self.extract_data()     # fill the useful self.soft attributes
        self.update_dirs()

    def pooling(self):
        for pool in self.config.pooling_groups:
            self.soft.io = {}
            self.pools[pool] = {}
            self.soft.outputs[pool] = {}
            pooling(self, pool)

    # IN PROGRESS
    def mapping(self):
        refs = {}
        fastqs, group_fps = self.config.fastq, self.inputs[self.pool]
        if self.soft.prev == 'drep':
            for algo, fp in self.inputs[self.pool].items():
                ref_fasta = '%s.fa' % fp
                if not isfile(ref_fasta):
                    cmd = 'if [ ! -f %s ]; then cat %s/*.fa > %s; fi' % (
                        ref_fasta, fp, ref_fasta)
                    self.outputs['cmds'].setdefault(
                        (self.pool, algo), []).append(cmd)
                refs[algo] = ref_fasta
        elif self.soft.prev == 'metawrap_reassemble':
            # needs revision...
            refs = {'': glob.glob(self.inputs[self.pool][self.sam][1])}
        elif self.soft.prev == 'spades':
            refs = {group: fps[1] for group, fps in group_fps.items()}

        self.outputs['outs'] = {}
        for group, fasta in refs.items():
            self.outputs['outs'][group] = {}
            out_dir = '%s/%s' % (self.dir, self.pool)
            if group:
                out_dir += '/%s' % group
            self.outputs['dirs'].append(out_dir)
            for sam in self.pools[self.pool][group]:
                bam_out = '%s/%s.bam' % (out_dir, sam)
                self.outputs['outs'][group][sam] = bam_out
                cmd = 'minimap2'
                cmd += ' -a -x sr'
                cmd += ' %s' % fasta
                cmd += ' %s' % fastqs[sam][0]
                cmd += ' %s' % fastqs[sam][1]
                cmd += ' | samtools view -F 0x104 -b -'
                cmd += ' | samtools sort -o %s - ' % bam_out
                if not isfile(bam_out) or self.config.force:
                    self.outputs['cmds'].setdefault(
                        (self.pool, group), []).append(cmd)
                bam_out_bai = '%s.bai' % bam_out
                cmd = 'samtools index %s' % bam_out
                if not isfile(bam_out_bai) or self.config.force:
                    self.outputs['cmds'].setdefault(
                        (self.pool, group), []).append(cmd)
                io_update(self, i_f=([fasta] + fastqs[sam]),
                          o_f=[bam_out, bam_out_bai])

    def prep_map__spades_prodigal(self):
        if 'prodigal' not in self.softs or 'mapping' not in self.softs:
            return None
        if self.softs['prodigal'].prev != 'spades':
            return None
        if self.softs['mapping'].prev != 'spades':
            return None
        prodigals_fps = self.softs['prodigal'].outputs
        sams_fps = self.softs['mapping'].outputs
        group_fps = self.inputs[self.pool]
        self.outputs['outs'] = {}
        for group, fps in group_fps.items():
            self.outputs['outs'][group] = {}
            for sam in self.pools[self.pool][group]:
                bam = sams_fps[self.pool][group][sam]
                prot = prodigals_fps[self.pool][group][1]
                out_dir = '%s/%s/%s' % (self.dir, self.pool, sam)
                out = '%s/reads.txt' % out_dir
                if not isfile(out):
                    cmd = 'pysam_reads_to_prodigal.py \\\n'
                    cmd += '-prodigal %s \\\n' % prot
                    cmd += '-bam %s \\\n' % bam
                    cmd += '-out %s\n' % out
                    self.outputs['cmds'].setdefault(
                        (self.pool, group), []).append(cmd)
                self.outputs['outs'][group][sam] = out
                io_update(self, i_f=[prot, bam, out])

    def prep_mocat(self):
        self.cmds = {}

    def prep_mapDamage(self):
        self.cmds = {}

    def prep_metaclade(self):
        self.cmds = {}

    def prep_ezTree(self):
        self.cmds = {}

    def prep_instrain(self):
        self.cmds = {}

    def prep_map__cazy_spades(self):
        self.cmds = {}

    def prep_map__cazy_macsyfinder(self):
        self.cmds = {}

    def prep_map__spades_bins(self):
        self.cmds = {}

    def prep_map__drep(self):
        self.cmds = {}

    def prep_mag_data(self):
        self.cmds = {}

    def prep_yamb(self):
        self.cmds = {}

    def prep_summarize_yamb(self):
        self.cmds = {}

    def prep_mycc(self):
        self.cmds = {}

    def prep_graspx(self):
        self.cmds = {}

    def prep_deepvirfinder(self):
        self.cmds = {}

    def prep_deepvirfinder_extractVirs(self):
        self.cmds = {}

    def prep_closed_ref_qiime1(self):
        self.cmds = {}

    def prep_phyloflash(self):
        self.cmds = {}

    def prep_wish(self):
        self.cmds = {}

    def prep_thmm(self):
        self.cmds = {}

    def prep_cazy(self):
        self.cmds = {}

    def prep_itasser(self):
        self.cmds = {}

    def register_command(self):
        self.softs[self.soft.name].cmds = dict(self.cmds)
        self.cmds = {}
