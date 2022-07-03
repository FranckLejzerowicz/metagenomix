# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob
from os.path import isdir, isfile

from metagenomix._cmds import caller
from metagenomix._io_utils import io_update


def get_sams_fastqs(mode: str, fastqs: dict):
    """Get the paths of the fastq files per sample or altogether.

    Parameters
    ----------
    mode : str
        Mode of processing for blobology:
        - "coassembly": all co-assembled samples.
        - "sample": each sample independently.
    fastqs : dict
        Path to the raw fastq files per sample.

    Returns
    -------
    sams_fastqs : dict
        Either same as `fastqs` of all fastqs under an empty key.
    """
    sams_fastqs = fastqs
    if mode == 'coassembly':
        sams_fastqs = {'': [y for x in fastqs.values() for y in x]}
    return sams_fastqs


def quantify(self):
    for group in self.pools[self.pool]:

        bins = self.inputs[self.pool][group][-1]
        contigs = self.softs['spades'].outputs[self.pool][group][1]
        fastqs = {sam: self.config.fastq[sam] for sam
                  in self.pools[self.pool][group]}

        self.outputs['outs'][group] = {}
        io_update(self, i_f=contigs, i_d=bins, key=group)

        for mode in self.soft.params['blobology']:
            self.outputs['outs'][group][mode] = {}

            out = '%s/%s/%s/%s' % (self.dir, self.pool, group, mode)
            sams_fastqs = get_sams_fastqs(mode, fastqs)
            for sam, fastq in sams_fastqs.items():
                if sam:
                    out += '/%s' % sam
                cmd = 'metawrap quant_bins'
                cmd += ' -b %s' % bins
                cmd += ' -o %s' % out
                cmd += ' -a %s' % contigs
                cmd += ' -t %s' % self.soft.params['cpus']
                cmd += ' %s' % ' '.join(fastq)

                self.outputs['dirs'].append(out)
                self.outputs['outs'][group][mode][sam] = out
                self.outputs['cmds'].setdefault(group, []).append(cmd)
                io_update(self, i_f=fastq, o_d=out, key=group)


def classify_or_annotate(self, command):
    for group in self.pools[self.pool]:

        if self.soft.prev == 'metawrap_refine':
            bin_dirs = {'': [self.inputs[self.pool][group][1]]}
        elif self.soft.prev == 'metawrap_reassemble':
            bin_dirs = self.inputs[self.pool][group]
        else:
            sys.exit('No metawrap classify_bins after %s' % self.soft.prev)

        self.outputs['outs'][group] = {}

        for sam, bin_dir in bin_dirs.items():
            for bin in bin_dir:
                out = '%s/%s/%s' % (self.dir, self.pool, group)
                if 'permissive' in bin:
                    out += '/permissive'
                elif 'strict' in bin:
                    out += '/strict'
                if sam:
                    out += '/%s' % sam
                cmd = 'metawrap %s' % command
                cmd += ' -b %s' % bin
                cmd += ' -o %s' % out
                cmd += ' -t %s' % self.soft.params['cpus']

                self.outputs['dirs'].append(out)
                self.outputs['outs'][group].setdefault(sam, []).append(out)
                io_update(self, i_d=bin, o_d=out, key=group)

                if command == 'annotate_bins':
                    tab = glob.glob('%s/bin_funct_annotations/*.tab' % out)
                    if self.config.force or not tab:
                        self.outputs['cmds'].setdefault(group, []).append(cmd)
                elif command == 'classify_bins':
                    tab = '%s/bin_taxonomy.tab' % out
                    if self.config.force or not isfile(tab):
                        self.outputs['cmds'].setdefault(group, []).append(cmd)


def classify(self):
    classify_or_annotate(self, 'classify_bins')


def annotate(self):
    classify_or_annotate(self, 'annotate_bins')


def blobology(self):
    for group in self.pools[self.pool]:
        self.outputs['outs'][group] = {}
        bins = self.inputs[self.pool][group][-1]
        contigs = self.softs['spades'].outputs[self.pool][group][1]
        fastqs = {
            sam: self.config.fastq[sam] for sam in self.pools[self.pool][group]}
        io_update(self, i_f=contigs, i_d=bins, key=group)
        for mode in self.soft.params['blobology']:
            self.outputs['outs'][group][mode] = {}
            out = '%s/%s/%s/%s' % (self.dir, self.pool, group, mode)
            sams_fastqs = get_sams_fastqs(mode, fastqs)
            for sam, fastq in sams_fastqs.items():
                if sam:
                    out += '/%s' % sam
                self.outputs['dirs'].append(out)
                io_update(self, i_f=fastq, o_d=out, key=group)
                plot = '%s/contigs.binned.blobplot' % out
                if self.config.force or not isfile(plot):
                    cmd = 'metawrap blobology'
                    cmd += ' -o %s' % out
                    cmd += ' -a %s' % contigs
                    cmd += ' -t %s' % self.soft.params['cpus']
                    cmd += ' --bins %s' % bins
                    cmd += ' %s' % ' '.join(fastq)
                    self.outputs['cmds'].setdefault(group, []).append(cmd)
                    self.outputs['outs'][group][mode][sam] = out


def reassembly_bins_cmd(self, sam, out, bins):
    cmd = 'metawrap reassemble_bins'
    cmd += ' -o %s' % out
    cmd += ' -1 %s' % self.config.fastq[sam][0]
    cmd += ' -2 %s' % self.config.fastq[sam][1]
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -m %s' % self.soft.params['mem_num']
    cmd += ' -c %s' % self.soft.params['min_completion_reassembly']
    cmd += ' -x %s' % self.soft.params['min_contamination_reassembly']
    cmd += ' -b %s' % bins
    cmd += ' --skip-checkm'
    if self.soft.params['cpus'] > 1:
        cmd += ' --parallel\n'
    return cmd


def reassembly_cmd(self, bins, group, sam, out):
    cmd = ''
    for mode in self.soft.params['reassembly']:
        mode_dir = '%s/reassembled_bins_%s' % (out, mode)
        io_update(self, i_d=mode_dir, key=group)
        if not self.config.force and isdir(mode_dir):
            continue
        cmd += 'mkdir %s\n' % mode_dir
        cmd += 'mv %s/reassembled_bins/*%s.fa %s/.\n' % (out, mode, mode_dir)
    if cmd:
        cmd = reassembly_bins_cmd(self, sam, out, bins) + cmd
    return cmd


def reassemble(self):
    for group in self.pools[self.pool]:
        bins = self.inputs[self.pool][group][-1]
        io_update(self, i_d=bins, key=group)
        for sam in self.pools[self.pool][group]:
            out = '%s/%s/%s/%s' % (self.dir, self.pool, group, sam)
            self.outputs['outs'].setdefault(group, []).append(out)
            self.outputs['dirs'].append(out)
            io_update(self, i_f=self.config.fastq[sam], o_d=out, key=group)
            cmd = reassembly_cmd(self, bins, group, sam, out)
            if cmd:
                self.outputs['cmds'].setdefault(group, []).append(cmd)


def refine_cmd(self, out_dir, bin_folders) -> tuple:
    n_bins = 0
    cmd = 'metawrap bin_refinement'
    cmd += ' -o %s' % out_dir
    for fdx, folder in enumerate(bin_folders):
        if len(glob.glob('%s/bin.*.fa' % folder)) or self.config.force:
            n_bins += 1
        cmd += ' -%s %s' % (['A', 'B', 'C'][fdx], folder)
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -c %s' % self.soft.params['min_completion']
    cmd += ' -x %s' % self.soft.params['min_contamination']
    return cmd, n_bins


def refine(self):
    for group in self.pools[self.pool]:
        out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
        bin_folders = self.inputs[self.pool][group][:-1]

        cmd, n_bins = refine_cmd(self, out_dir, bin_folders)
        if n_bins == len(bin_folders):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
        out = '%s/metawrap_%s_%s' % (out_dir,
                                     self.soft.params['min_completion'],
                                     self.soft.params['min_contamination'])
        stats, bins = '%s.stats' % out, '%s_bins' % out
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = [stats, bins]
        io_update(self, i_d=bin_folders, o_d=[stats, bins], key=group)


def get_binners(self, out, binned):
    binners = []
    work_files = '%s/work_files' % out
    for binner, bins in binned.items():
        if self.config.force or not isdir(bins) or not isdir(work_files):
            binners.append(binner)
    return binners


def binning_cmd(self, fastqs, out, contigs, binned):
    binners = get_binners(self, out, binned)
    if binners:
        cmd = 'metawrap binning'
        cmd += ' -o %s' % out
        cmd += ' -a %s' % contigs
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -l %s' % self.soft.params['min_contig_length']
        cmd += ' --universal'
        for binner in binners:
            cmd += ' --%s' % binner
        cmd += ' %s' % ' '.join(fastqs)
        return cmd


def binning(self):
    for group in self.pools[self.pool]:
        tmp = '$TMPDIR/mtwrp_%s_%s' % (self.pool, group)
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        self.outputs['dirs'].append(out)
        binned = {binner: '%s/%s_bins' % (out, binner)
                  for binner in self.soft.params['binners']}
        # outputs
        bin_dirs = sorted(binned.values()) + ['work_files']
        self.outputs['outs'][group] = bin_dirs
        # inputs
        contigs = self.inputs[self.pool][group][1]
        fastqs = [fastq for sam in self.pools[self.pool][group]
                  for fastq in self.config.fastq[sam]]
        cmd = binning_cmd(self, fastqs, out, contigs, binned)
        if cmd:
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=([contigs] + fastqs), i_d=tmp,
                      o_d=bin_dirs, key=group)


def metawrap(self) -> None:
    """Create command lines for metaWRAP

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    # This function splits the name of the software and calls as function
    # the last underscore-separated field (which is in this module)
    caller(self, __name__)
