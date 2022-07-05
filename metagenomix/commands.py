# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import inspect
import itertools
from os.path import abspath, basename, isdir, isfile, splitext

from metagenomix._io_utils import get_out_dir
from metagenomix.tools.preprocess import (
    edit, count, fastqc, fastp, cutadapt, atropos, kneaddata, filtering)
from metagenomix.tools.simka import simka
from metagenomix.tools.alignment import bowtie2, flash
from metagenomix.tools.profiling import shogun, woltka, kraken2, midas, metaxa2
from metagenomix.tools.phlans import metaphlan, humann, strainphlan
from metagenomix.tools.assembly import (
    pooling, spades, quast, plass, viralverify)
from metagenomix.tools.annotation import (
    prodigal, integron_finder, macsyfinder, ioncom, search, antismash, prokka)
from metagenomix.tools.binning import metawrap
from metagenomix.tools.drep import get_drep_bins, get_drep_inputs
from metagenomix.tools.metamarker import metamarker


class Commands(object):

    def __init__(self, config, databases, workflow):
        self.config = config
        self.databases = databases
        self.softs = workflow.softs
        self.outputs = {}
        self.cmds = {}
        self.args = {}
        self.pools = {}
        self.pool = None
        self.pool_bool = False
        self.inputs = None
        self.method = None
        self.soft = None
        self.sam = None
        self.dir = ''
        self.out = []
        self.io = set
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
            self.soft = self.softs[softs[-1]]
            print('[Collecting commands] #%s: %s' % (sdx, self.soft.name))
            self.get_inputs()
            self.get_dir()
            self.generic_command()
            self.update_dirs()

    def get_inputs(self):
        """Update the `inputs` attribute of the software object."""
        if not self.soft.prev or self.soft.name == 'map__drep':
            if self.soft.params['scratch']:
                self.inputs = self.config.fastq_scratch
            else:
                self.inputs = self.config.fastq
        else:
            self.inputs = self.softs[self.soft.prev].outputs

    def get_dir(self):
        self.dir = abspath('%s/%s/after_%s' % (
            self.config.dir, self.soft.name, self.soft.prev))
        self.soft.dirs.add(self.dir)
        if self.soft.params['scratch']:
            self.dir = '${SCRATCH_FOLDER}%s' % self.dir

    def is_pool(self):
        self.pool_bool = False
        self.struc = list
        self.io = set
        if set(self.inputs) == set(self.pools):
            self.pool_bool = True
            self.struc = dict
            self.io = dict

    def generic_command(self):
        self.sam = ''
        self.is_pool()
        self.soft.io = {}
        if self.soft.name in self.holistics:
            self.prep_job()
        elif self.soft.name == 'pooling':
            self.pooling()
        else:
            for sam_or_pool in sorted(self.inputs):
                self.sam = sam_or_pool
                self.pool = sam_or_pool
                self.prep_job()
        self.register_command()

    def fill_soft_io(self):
        for i, j in itertools.product(*[['I', 'O'], ['d', 'f']]):
            if self.io == set:
                if self.sam not in self.soft.io:
                    self.soft.io[self.sam] = {}
                self.soft.io[self.sam][(i, j)] = self.outputs[
                    'io'].get((i, j), self.io())
            else:
                for k, v in self.outputs['io'].get((i, j), self.io()).items():
                    if (self.sam, k) not in self.soft.io:
                        self.soft.io[(self.sam, k)] = {}
                    self.soft.io[(self.sam, k)][(i, j)] = v

    def update_dirs(self):
        self.soft.dirs = set([
            # x.replace('${SCRATCH_FOLDER}%s' % self.dir, self.dir)
            x.replace('${SCRATCH_FOLDER}/', '/')
            for x in self.soft.dirs])

    def init_outputs(self):
        self.outputs = {'cmds': self.struc(), 'outs': self.struc(), 'dirs': [],
                        'io': {('I', 'd'): self.io(), ('I', 'f'): self.io(),
                               ('O', 'd'): self.io(), ('O', 'f'): self.io()}}

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

    def extract_data(self):
        if self.outputs.get('cmds', []):
            self.cmds[self.sam] = self.outputs['cmds']
            self.fill_soft_io()
        if self.soft.name in self.holistics:
            self.soft.outputs[self.soft.name] = self.outputs['outs']
        else:
            self.soft.outputs[self.pool] = self.outputs['outs']

    def prep_job(self):
        self.init_outputs()     # initialize data structure collection
        self.call_method()      # collect commands, outputs, io, dirs
        self.extract_data()     # fill the useful self.soft attributes

    def pooling(self):
        for pool in self.config.pooling_groups:
            self.pools[pool] = {}
            self.soft.outputs[pool] = {}
            self.soft.io[pool] = {('I', 'f'): dict(), ('I', 'd'): dict(),
                                  ('O', 'f'): dict(), ('O', 'd'): dict()}
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
        print(sdfkvjjb)

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

    def prep_drep(self):
        self.outputs['outs'] = {}
        genomes = get_drep_bins(self.soft.prev, self.pools, self.inputs)
        for pool, group_paths in genomes.items():
            self.outputs['outs'][pool] = {}
            for stringency, sam_paths in group_paths.items():
                drep_dir = '%s/%s' % (self.dir, pool)
                if stringency:
                    drep_dir += '/%s' % stringency
                drep_in, paths = get_drep_inputs(drep_dir, sam_paths)
                io_update(self, i_f=([drep_in] + list(paths)))
                for algorithm in ['fastANI', 'ANIn']:
                    drep_out = '%s/%s' % (drep_dir, algorithm)
                    self.outputs['dirs'].append(drep_out)
                    log = '%s/log' % drep_out
                    figure = '%s/figure' % drep_out
                    data_table = '%s/data_tables' % drep_out
                    dereplicated_genomes = '%s/dereplicated_genomes' % drep_out
                    outputs = [data_table, log, figure, dereplicated_genomes]
                    self.outputs['outs'][pool][algorithm] = outputs
                    if len(glob.glob('%s/*.fa' % dereplicated_genomes)):
                        continue
                    if isdir(drep_out):
                        io_update(self, i_d=drep_out)
                    cmd = 'dRep dereplicate'
                    cmd += ' %s' % drep_out
                    cmd += ' --S_algorithm %s' % algorithm
                    cmd += ' --ignoreGenomeQuality'
                    if len(list(paths)) > 5000 and algorithm == 'fastANI':
                        cmd += ' --multiround_primary_clustering'
                        cmd += ' --greedy_secondary_clustering'
                        cmd += ' --run_tertiary_clustering'
                    cmd += ' -pa 0.9 -sa 0.98 -nc 0.3 -cm larger'
                    cmd += ' -p %s' % self.soft.params['cpus']
                    cmd += ' -g %s' % drep_in
                    self.outputs['cmds'].setdefault(pool, []).append(cmd)
                    io_update(self, o_d=outputs)

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

    def prep_gtdbtk(self):
        self.cmds = {}

    def prep_graspx(self):
        self.cmds = {}

    def prep_checkm(self):
        self.cmds = {}

    def prep_deepvirfinder(self):
        self.cmds = {}

    def prep_deepvirfinder_extractVirs(self):
        self.cmds = {}

    def prep_metaxa2(self):
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

    # def add_transfer_localscratch(self):
    #     pass
    #     # if self.scratch:
    #     #     self.add_transfer_localscratch()
