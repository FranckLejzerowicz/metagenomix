# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob

from metagenomix._io_utils import caller, io_update, to_do


def get_sams_fastqs(
        mode: str,
        fastqs: dict
) -> dict:
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


def get_fqs(
        fastq: list
) -> tuple:
    """Get the fastq files for the current sample to reassemble

    Parameters
    ----------
    fastq : list
        Paths to fastq files for the current sample to reassemble

    Returns
    -------
    fqs : list
        Paths to gunzipped fastq files for the current sample to reassemble
    cmd : str
        gunzip commands
    """
    fqs = []
    cmd = ''
    for f in fastq:
        if f.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (f, f.replace('.gz', ''))
            fqs.append(f.replace('.gz', ''))
        else:
            fqs.append(f)
    return fqs, cmd


def quantify(self):
    """Quantify the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
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
                fqs, fqs_cmd = get_fqs(fastq)
                cmd = 'metawrap quant_bins'
                cmd += ' -b %s' % bins
                cmd += ' -o %s' % out
                cmd += ' -a %s' % contigs
                cmd += ' -t %s' % self.soft.params['cpus']
                cmd += ' %s\n' % ' '.join(fqs)
                if fqs_cmd:
                    cmd = fqs_cmd + cmd + 'rm %s\n' % ' '.join(fqs)
                self.outputs['dirs'].append(out)
                self.outputs['outs'][group][mode][sam] = out
                self.outputs['cmds'].setdefault(group, []).append(cmd)
                io_update(self, i_f=fastq, o_d=out, key=group)


def classify_or_annotate(self, command):
    """Classify or annotate the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
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
                    tab = '%s/bin_funct_annotations/*.tab' % out.replace(
                        '${SCRATCH_FOLDER}', '')
                    if self.config.force or not glob.glob(tab):
                        self.outputs['cmds'].setdefault(group, []).append(cmd)
                elif command == 'classify_bins':
                    tab = '%s/bin_taxonomy.tab' % out
                    if self.config.force or to_do(tab):
                        self.outputs['cmds'].setdefault(group, []).append(cmd)


def classify(self):
    """Classify the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
    classify_or_annotate(self, 'classify_bins')


def annotate(self):
    """Annotate the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
    classify_or_annotate(self, 'annotate_bins')


def get_blobology_cmd(
        self,
        fastq_fps: list,
        out: str,
        contigs: str,
        bins: str
) -> str:
    """Classify or annotate the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    fastq_fps : list
        Path the input fastq files
    out : str
        Path to the metaWRAP bin reassembly output folder
    contigs : str
        Path to the input contigs file
    bins : str
        Path to the input bins folder
    """
    fqs, fqs_cmd = get_fqs(fastq)
    cmd = 'metawrap blobology'
    cmd += ' -o %s' % out
    cmd += ' -a %s' % contigs
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' --bins %s' % bins
    cmd += ' %s\n' % ' '.join(fqs)
    if fqs_cmd:
        cmd = fqs_cmd + cmd + 'rm %s\n' % ' '.join(fqs)
    return cmd


def blobology(self):
    """Make blobology plots for the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
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
            for sam, fastq_fps in sams_fastqs.items():
                if sam:
                    out += '/%s' % sam
                self.outputs['dirs'].append(out)
                io_update(self, i_f=fastq_fps, o_d=out, key=group)
                plot = '%s/contigs.binned.blobplot' % out
                if self.config.force or to_do(plot):
                    cmd = get_blobology_cmd(self, fastq_fps, out, contigs, bins)
                    self.outputs['cmds'].setdefault(group, []).append(cmd)
                    self.outputs['outs'][group][mode][sam] = out


def reassembly_bins_cmd(
        self,
        sam: str,
        out: str,
        bins: str
) -> str:
    """Collect metaWRAP bin reassembly command

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    sam : str
        Sample name
    out : str
        Path to the metaWRAP bin reassembly output folder
    bins : str
        Path to the input bins folder

    Returns
    -------
    cmd : str
        metaWRAP bin reassembly command
    """
    fqs, fqs_cmd = get_fqs(self.config.fastq[sam])
    cmd = 'metawrap reassemble_bins'
    cmd += ' -o %s' % out
    for idx, fq in enumerate(fqs):
        cmd += ' -%s %s' % ((idx + 1), fq)
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -m %s' % self.soft.params['mem']
    cmd += ' -c %s' % self.soft.params['min_completion_reassembly']
    cmd += ' -x %s' % self.soft.params['min_contamination_reassembly']
    cmd += ' -b %s' % bins
    cmd += ' --skip-checkm'
    if self.soft.params['cpus'] > 1:
        cmd += ' --parallel'
    cmd += '\n'
    if fqs_cmd:
        cmd = fqs_cmd + cmd + 'rm %s\n' % ' '.join(fqs)
    return cmd


def reassembly_cmd(
        self,
        bins: str,
        group: str,
        sam: str,
        out: str
) -> str:
    """Collect metaWRAP bin reassembly command

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    bins : str
        Path to the input bins folder
    group : str
        Name of the sample group within the pool
    sam : str
        Sample name
    out : str
        Path to the metaWRAP bin reassembly output folder

    Returns
    -------
    cmd : str
        metaWRAP bin reassembly command
    """
    cmd = ''
    for mode in self.soft.params['reassembly']:
        mode_dir = '%s/reassembled_bins_%s' % (out, mode)
        self.outputs['dirs'].append(mode_dir)
        self.outputs['outs'].setdefault(group, []).append(mode_dir)
        io_update(self, i_d=mode_dir, key=group)
        if not self.config.force and not to_do(folder=mode_dir):
            self.soft.status.add('Done')
            continue
        cmd += 'mkdir %s\n' % mode_dir
        cmd += 'mv %s/reassembled_bins/*%s.fa %s/.\n' % (out, mode, mode_dir)
    if cmd:
        cmd = reassembly_bins_cmd(self, sam, out, bins) + cmd
    return cmd


def reassemble(self):
    """Reassemble the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
        .status
            Tool status
    """
    for group in self.pools[self.pool]:
        bins = self.inputs[self.pool][group][-1]
        if to_do(folder=bins):
            self.soft.status.add('Must run %s (%s)' % (self.soft.prev, group))
        io_update(self, i_d=bins, key=group)
        for sam in self.pools[self.pool][group]:
            out = '%s/%s/%s/%s' % (self.dir, self.pool, group, sam)
            self.outputs['dirs'].append(out)
            io_update(self, i_f=self.config.fastq[sam], o_d=out, key=group)
            cmd = reassembly_cmd(self, bins, group, sam, out)
            if cmd:
                self.outputs['cmds'].setdefault(group, []).append(cmd)


def refine_cmd(
        self,
        out_dir: str,
        bin_folders: list
) -> tuple:
    """Collect metaWRAP bin refinement command

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
        .config
            Configurations
    out_dir : str
        Path to the output folder
    bin_folders : list
        Path to the input bin folder

    Returns
    -------
    cmd : str
        metaWRAP bin refinement command
    n_bins : int
        Number of binners used
    """
    n_bins = 0
    cmd = 'metawrap bin_refinement'
    cmd += ' -o %s' % out_dir
    for fdx, folder in enumerate(bin_folders):
        folders = '%s/bin.*.fa' % folder.replace('${SCRATCH_FOLDER}', '')
        if len(glob.glob(folders)) or self.config.force:
            n_bins += 1
        cmd += ' -%s %s' % (['A', 'B', 'C'][fdx], folder)
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -c %s' % self.soft.params['min_completion']
    cmd += ' -x %s' % self.soft.params['min_contamination']
    return cmd, n_bins


def refine(self):
    """Refine the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    for group in self.pools[self.pool]:
        out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
        bin_folders = self.inputs[self.pool][group][:-1]
        if sum([to_do(folder=x) for x in bin_folders]):
            self.soft.status.add('Must run %s (%s)' % (self.soft.prev, group))

        cmd, n_bins = refine_cmd(self, out_dir, bin_folders)
        if n_bins == len(bin_folders):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
        out = '%s/metawrap_%s_%s' % (out_dir,
                                     self.soft.params['min_completion'],
                                     self.soft.params['min_contamination'])
        stats, bins = '%s.stats' % out, '%s_bins' % out
        self.outputs['outs'][group] = [stats, bins]
        if not self.config.force and not to_do(folder=bins):
            self.soft.status.add('Done')
            continue
        self.outputs['dirs'].append(out_dir)
        io_update(self, i_d=bin_folders, o_d=[stats, bins], key=group)


def get_binners(
        self,
        out: str,
        binned: dict
) -> list:
    """Get the binning tools that were not run yet.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
        .config
            Configurations
    out : str
        Path to the output folder
    binned : dict
        Path to the output folder per binner tool

    Returns
    -------
    binners : list
        Binner tools used for the metaWRAP binning
    """
    binners = []
    work_files = '%s/work_files' % out
    for binner, bins in binned.items():
        if self.config.force or to_do(folder=bins) or to_do(folder=work_files):
            binners.append(binner)
    return binners


def binning_cmd(
        self,
        fastq: list,
        out: str,
        contigs: str,
        binned: dict
) -> str:
    """Collect metaWRAP binning command

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastq : list
        Paths to the input fastq files
    out : str
        Path to the output folder
    contigs : str
        Path to the input contigs file
    binned : dict
        Path to the output folder per binner tool

    Returns
    -------
    cmd : str
        metaWRAP binning comand
    """
    binners = get_binners(self, out, binned)
    if binners:
        fqs, fqs_cmd = get_fqs(fastq)
        cmd = 'metawrap binning'
        cmd += ' -o %s' % out
        cmd += ' -a %s' % contigs
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -l %s' % self.soft.params['min_contig_length']
        cmd += ' --universal'
        for binner in binners:
            cmd += ' --%s' % binner
        cmd += ' %s\n' % ' '.join(fqs)
        if fqs_cmd:
            cmd = fqs_cmd + cmd + 'rm %s\n' % ' '.join(fqs)
        return cmd


def binning(self):
    """Bin the assembled contigs using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    for group in self.pools[self.pool]:
        tmp = '$TMPDIR/mtwrp_%s_%s' % (self.pool, group)
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        self.outputs['dirs'].append(out)
        binned = {binner: '%s/%s_bins' % (out, binner)
                  for binner in self.soft.params['binners']}
        bin_dirs = sorted(binned.values()) + ['work_files']
        self.outputs['outs'][group] = bin_dirs
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
    module_call = caller(self, __name__)
    module_call(self)
