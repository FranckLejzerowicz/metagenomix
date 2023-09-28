# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
import glob
import pkg_resources
from os.path import basename, dirname, splitext
from metagenomix._io_utils import (
    caller, status_update, io_update, to_do, get_assembly, get_assembly_graph)
from metagenomix.core.parameters import tech_params

SCRIPTS = pkg_resources.resource_filename("metagenomix", "resources/scripts")


def process_outputs(
        self,
        key: tuple,
        group: str,
        out: list,
) -> bool:
    """

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    key : tuple
        (tech, group)
    group : str
        Current pool group
    out : list
        Path(s) to metaWRAP output folder(s)

    Returns
    -------
    process : bool
        Whether to process the sample command or not
    """
    if group in self.soft.params['skip_samples']:
        self.outputs['outs'][key] = []
        process = True
    else:
        self.outputs['outs'][key] = out
        process = False
    return process


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
        fastq_fps: list
) -> tuple:
    """Get the fastq files for the current sample to reassemble

    Parameters
    ----------
    fastq_fps : list
        Paths to fastq files for the current sample to reassemble

    Returns
    -------
    gz : bool
        Whether the fastq were gzipped
    fqs : list
        Paths to gunzipped fastq files for the current sample to reassemble
    cmd : str
        gunzip commands
    """
    gz = False
    fqs = []
    cmd = ''
    for f in fastq_fps:
        if f.endswith('.gz'):
            f_q = f.replace('.gz', '')
            cmd += 'gunzip -c %s > %s\n' % (f, f_q)
            gz = True
        else:
            f_q = f
        if f_q.endswith('_R1.fastq'):
            f_q_ = f_q.replace('_R1.fastq', '_1.fastq')
            cmd += 'mv %s %s\n' % (f_q, f_q_)
        elif f_q.endswith('_R2.fastq'):
            f_q_ = f_q.replace('_R2.fastq', '_2.fastq')
            cmd += 'mv %s %s\n' % (f_q, f_q_)
        else:
            f_q_ = f_q
        fqs.append(f_q_)
    return gz, fqs, cmd


def quantify_cmd(
        self,
        bins: str,
        fastq: list,
        out: str,
        contigs: str
) -> str:
    """Collect metaWRAP quantify command.

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
    bins : str
        Path to the metaWRAP bins folder
    fastq : list
        Paths to the input fastq files
    out : str
        Paths to the output folder
    contigs : str
        Paths to the input contigs fasta file

    Returns
    -------
    cmd : str
        metaWRAP quantify command
    """
    gz, fqs, fqs_cmd = get_fqs(fastq)
    cmd, cmd_rm = '', ''
    if contigs.endswith('.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        contigs = contigs.rstrip('.gz')
        cmd_rm += 'rm %s\n' % contigs

    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

    cmd += 'metawrap quant_bins'
    cmd += ' -b %s' % bins
    cmd += ' -o %s' % out
    cmd += ' -a %s' % contigs
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' %s\n' % ' '.join(fqs)
    if fqs_cmd:
        cmd = fqs_cmd + cmd
    if gz:
        cmd += 'rm %s\n' % ' '.join(fqs)
    cmd += cmd_rm
    return cmd


def get_fastqs(
        self,
        group: str,
) -> dict:
    """Collect the illumina reads for the samples of the current pool group.

    Parameters
    ----------
    self : Commands class instance
        .pool : str
            Pool name
        .pools : str
            All Pools and their samples
        .status
            Tool status
    group : str
        Current pool group

    Returns
    -------
    fastqs : dict
        Paths to fastqs files per sample
    """

    for software in self.soft.path[::-1]:
        if self.config.tools[software] == 'preprocessing':
            break
    else:
        sys.exit('[metawrap_bining] No "preprocessing" fastqs found')

    fastqs = {}
    for sam in self.pools[self.sam_pool][group]:
        if self.softs[software].outputs[sam].get(('illumina', sam), []):
            fastqs[sam] = self.softs[software].outputs[sam][('illumina', sam)]
    return fastqs


def quantify(self):
    """Quantify the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)
        self.outputs['outs'][key] = {}

        bins = inputs[0]
        assembler = self.softs['metawrap_binning'].prev
        contigs = self.softs[assembler].outputs[self.sam_pool][key][0]
        fastqs = get_fastqs(self, group)
        if not fastqs:
            self.soft.add_status(
                tech, self.sam_pool, group=group, message='no illumina data')
            print('[metawrap_blobology] No illumina reads for group %s' % group)
            continue

        io_update(self, i_f=contigs, i_d=bins, key=key)

        for mode in self.soft.params['blobology']:
            self.outputs['outs'][key][mode] = {}

            out_ = '%s/%s/%s/%s/%s' % (
                self.dir, tech, self.sam_pool, group, mode)
            sams_fastqs = get_sams_fastqs(mode, fastqs)
            for sam, fastq in sams_fastqs.items():
                to_dos = status_update(self, tech, fastq, group=group)
                to_dos.extend(
                    status_update(self, tech, [bins], group=group, folder=True))
                if sam:
                    o = out_ + '/%s' % sam
                else:
                    o = out_
                self.outputs['dirs'].append(o)
                cmd = quantify_cmd(self, bins, fastq, o, contigs)
                if self.config.force or to_do('%s/bin_abundance_table.tab' % o):
                    self.outputs['outs'][key][mode][sam] = o
                    if to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, i_f=fastq, o_d=o, key=key)
                    self.soft.add_status(
                        tech, self.sam_pool, 1, group=group, genome=sam)
                else:
                    self.soft.add_status(
                        tech, self.sam_pool, 0, group=group, genome=sam)


def classify_or_annotate_dirs(
        self,
        inputs: list
) -> dict:
    """

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
    inputs : list
        Paths to the input files

    Returns
    -------
    bin_dirs : dict
        Path to the bins input folder per assembly solution
    """
    if self.soft.prev == 'metawrap_refine':
        bin_dirs = {'': [inputs[0]]}
    elif self.soft.prev == 'metawrap_reassemble':
        bin_dirs = inputs[self.sam_pool]
    else:
        sys.exit('No metawrap classify_bins after %s' % self.soft.prev)
    return bin_dirs


def classify_or_annotate_out(
        self,
        tech: str,
        bin_dir: str,
        group: str,
        sam: str
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP
        .pool : str
            Pool name.
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    bin_dir : str
        Path to the bins input folder
    group : str
        Group for the current co-assembly and
    sam : str
        Name of a sample

    Returns
    -------
    out : str
        Path to the output folder
    """
    out = '/'.join([self.dir, tech, self.sam_pool, group])
    if 'permissive' in bin_dir:
        out += '/permissive'
    elif 'strict' in bin_dir:
        out += '/strict'
    if sam:
        out += '/%s' % sam
    return out


def classify_or_annotate_cmd(
        self,
        command: str,
        bin_dir: str,
        out: str
) -> str:
    """Collect metaWRAP classify_bins or annotate_bins command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    command : str
        "classify_bins" or "annotate_bins"
    bin_dir : str
        Path to the bins input folder
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        metaWRAP classify_bins or annotate command_bins
    """
    cmd = ''
    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

    cmd += 'metawrap %s' % command
    cmd += ' -b %s' % bin_dir
    cmd += ' -o %s' % out
    cmd += ' -t %s' % self.soft.params['cpus']
    return cmd


def classify_or_annotate(
        self,
        command: str
) -> None:
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
    command : str
        "classify_bins" or "annotate_bins"
    """
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)
        self.outputs['outs'][key] = {}

        for sam, bins_dir in classify_or_annotate_dirs(self, inputs).items():
            for bin_dir in bins_dir:

                out = classify_or_annotate_out(self, tech, bin_dir, group, sam)
                self.outputs['dirs'].append(out)

                to_dos = status_update(
                    self, tech, [bin_dir], group=group, genome=sam, folder=True)
                cmd = classify_or_annotate_cmd(self, command, bin_dir, out)
                self.outputs['outs'][key][sam] = [out]

                condition = classify_or_annotate_condition(command, out)
                if self.config.force or condition:
                    if to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, i_d=bin_dir, o_d=out, key=key)
                    self.soft.add_status(
                        tech, self.sam_pool, 1, group=group, genome=sam)
                else:
                    self.soft.add_status(
                        tech, self.sam_pool, 0, group=group, genome=sam)


def classify_or_annotate_condition(
        command: str,
        out: str
) -> bool:
    """Get the boolean telling whether the softwares was already run.

    Parameters
    ----------
    command : str
        "classify_bins" or "annotate_bins"
    out : str
        Path to the output folder

    Returns
    -------
    condition : bool
        True if classify_bins or annotate_bins was not already run
    """
    out = out.replace('${SCRATCH_FOLDER}', '')
    if command == 'annotate_bins':
        tab = '%s/bin_funct_annotations/*.tab' % out
        condition = not glob.glob(tab)
    elif command == 'classify_bins':
        tab = '%s/bin_taxonomy.tab' % out
        condition = to_do(tab)
    return condition


def classify(self):
    """Classify the bins obtained using metaWRAP.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for metaWRAP classify_bins
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
            Path to pipeline output folder for metaWRAP annotate_bins
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


def blobology_cmd(
        self,
        fastq_fps: list,
        out: str,
        contigs: str,
        bins: str
) -> str:
    """Collect metaWRAP blobology command.

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

    Returns
    -------
    cmd : str
        metaWRAP blobology command
    """
    gz, fqs, fqs_cmd = get_fqs(fastq_fps)
    cmd, cmd_rm = '', ''

    if contigs.endswith('.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        contigs = contigs.rstrip('.gz')
        cmd_rm += 'rm %s\n' % contigs

    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

    cmd += 'metawrap blobology'
    cmd += ' -o %s' % out
    cmd += ' -a %s' % contigs
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' --bins %s' % bins
    cmd += ' %s\n' % ' '.join(fqs)
    if fqs_cmd:
        cmd = fqs_cmd + cmd
    if gz:
        cmd += 'rm %s\n' % ' '.join(fqs)
    cmd += cmd_rm
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)

        self.outputs['outs'][key] = {}
        bins = inputs[0]

        assembler = self.softs['metawrap_binning'].prev
        contigs = self.softs[assembler].outputs[self.sam_pool][key][0]
        fastqs = get_fastqs(self, group)
        if not fastqs:
            # self.soft.add_status(tech, self.sam_pool, [fastqs], group=group)
            print('[metawrap_blobology] No illumina reads for group %s' % group)
            continue

        io_update(self, i_f=contigs, i_d=bins, key=key)
        for mode in self.soft.params['blobology']:
            self.outputs['outs'][key][mode] = {}
            out_ = '/'.join([self.dir, tech, self.sam_pool, group, mode])
            sams_fastqs = get_sams_fastqs(mode, fastqs)
            for sam, fastq_fps in sams_fastqs.items():
                if mode == 'coassembly' and len(fastq_fps) <= 2:
                    continue
                to_dos = status_update(
                    self, tech, fastq_fps, group=group, genome=sam)
                to_dos.extend(status_update(
                    self, tech, [bins], group=group, genome=sam, folder=True))
                if sam:
                    out = out_ + '/%s' % sam
                else:
                    out = out_

                self.outputs['dirs'].append(out)
                base = basename(splitext(contigs.rstrip('.gz'))[0])
                blobplot = '%s/%s.binned.blobplot' % (out, base)
                self.outputs['outs'][key][mode][sam] = out
                if self.config.force or to_do(blobplot):
                    cmd = blobology_cmd(self, fastq_fps, out, contigs, bins)
                    if to_dos:
                        self.outputs['cmds'].setdefault(key, []).append(False)
                    else:
                        self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, i_f=fastq_fps, o_d=out, key=key)
                    self.soft.add_status(
                        tech, self.sam_pool, 1, group=group, genome=sam)
                else:
                    self.soft.add_status(
                        tech, self.sam_pool, 0, group=group, genome=sam)


def reassembly_bins_cmd(
        self,
        fastq: list,
        out: str,
        bins: str
) -> str:
    """Collect metaWRAP bin reassembly command

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    fastq : list
        Paths to the input fastq files
    out : str
        Path to the metaWRAP bin reassembly output folder
    bins : str
        Path to the input bins folder

    Returns
    -------
    cmd : str
        metaWRAP bin reassembly command
    """
    gz, fqs, fqs_cmd = get_fqs(fastq)
    cmd = ''
    if self.soft.params['path']:
        cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

    cmd += 'metawrap reassemble_bins'
    cmd += ' -o %s' % out
    for idx, fq in enumerate(fqs):
        cmd += ' -%s %s' % ((idx + 1), fq)
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -m %s' % self.soft.params['mem']
    cmd += ' -c %s' % self.soft.params['min_completion_reassembly']
    cmd += ' -x %s' % self.soft.params['max_contamination_reassembly']
    cmd += ' -b %s' % bins
    cmd += ' --skip-checkm'
    if self.soft.params['cpus'] > 1:
        cmd += ' --parallel'
    cmd += '\n'
    if fqs_cmd:
        cmd = fqs_cmd + cmd
    if gz:
        cmd += 'rm %s\n' % ' '.join(fqs)
    return cmd


def reassembly_cmd(
        self,
        fastq: list,
        bins: str,
        tech: str,
        group: str,
        out: str,
        key: tuple
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
    fastq : list
        Paths to the input fastq files
    bins : str
        Path to the input bins folder
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    group : str
        Name of the current co-assembly group
    out : str
        Path to the metaWRAP bin reassembly output folder
    key : tuple
        (tech, group)

    Returns
    -------
    cmd : str
        metaWRAP bin reassembly command
    """
    cmd = ''
    for mode in self.soft.params['reassembly']:
        mode_dir = '%s/reassembled_bins_%s' % (out, mode)
        self.outputs['dirs'].append(mode_dir)
        self.outputs['outs'].setdefault(key, []).append(mode_dir)
        io_update(self, i_d=mode_dir, key=key)
        if not self.config.force and not to_do(mode_dir):
            continue
        cmd += 'mkdir %s\n' % mode_dir
        cmd += 'mv %s/reassembled_bins/*%s.fa %s/.\n' % (out, mode, mode_dir)
    if cmd:
        cmd = reassembly_bins_cmd(self, fastq, out, bins) + cmd
        self.soft.add_status(tech, self.sam_pool, 1, group=group)
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group=group)
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)

        bins = inputs[0]
        # if to_do(folder=bins):
        #     self.soft.add_status(tech, self.sam_pool, [bins], group=group)

        fastqs = get_fastqs(self, group)
        for sam, fastq in fastqs.items():
            if not fastq:
                print('[metawrap_reassemble] No illumina reads for %s' % sam)
                continue
            to_dos = status_update(self, tech, fastq, group=group, genome=sam)
            to_dos.extend(status_update(
                self, tech, [bins], group=group, genome=sam, folder=True))

            out = '/'.join([self.dir, tech, self.sam_pool, group, sam])
            self.outputs['dirs'].append(out)

            cmd = reassembly_cmd(self, fastq, bins, tech, group, out, key)
            if cmd:
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=self.config.fastq[sam], i_d=bins,
                          o_d=out, key=key)
                self.soft.add_status(
                    tech, self.sam_pool, 1, group=group, genome=sam)
            else:
                self.soft.add_status(
                    tech, self.sam_pool, 0, group=group, genome=sam)


def refine_cmd(
        self,
        out_dir: str,
        bin_folders: list,
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
    bins : str
        Path to the output folder
    names : str
        Path to the output names file
    """
    out = '%s/metawrap_%s_%s' % (out_dir, self.soft.params['min_completion'],
                                 self.soft.params['max_contamination'])
    bins = '%s_bins' % out
    names = '%s.names' % bins
    stats = '%s.stats' % bins
    cmd = ''
    if to_do(bins) or to_do(stats):

        if self.soft.params['path']:
            cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

        cmd += 'rm -rf %s\n' % out_dir
        cmd += 'metawrap bin_refinement'
        cmd += ' -o %s' % out_dir
        for fdx, folder in enumerate(bin_folders):
            cmd += ' -%s %s' % (['A', 'B', 'C'][fdx], folder)
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -c %s' % self.soft.params['min_completion']
        cmd += ' -x %s\n' % self.soft.params['max_contamination']
        cmd += 'rm -rf %s/work_files\n' % out_dir
        cmd += 'for i in %s/{ma,metab,c}*_bins; do rm -rf $i; done\n' % out_dir
    return cmd, bins, names


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
    for (tech, group), bin_dirs in self.inputs[self.sam_pool].items():

        key = (tech, group)

        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        to_dos = status_update(self, tech, bin_dirs, group=group, folder=True)

        cmd, bins, names = refine_cmd(self, out_dir, bin_dirs)
        if process_outputs(self, key, group, [bins]):
            continue
        # self.outputs['dirs'].append(dirname(out_dir))
        self.outputs['dirs'].append(out_dir)

        stats = '%s.stats' % bins
        if self.config.force or to_do(bins) or to_do(stats):
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_d=bin_dirs, o_d=out_dir, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def get_binners(
        self,
        binned: dict
) -> list:
    """Get the binning softwares that were not run yet.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
        .config
            Configurations
    binned : dict
        Path to the output folder per binner tool

    Returns
    -------
    binners : list
        Binner softwares used for the metaWRAP binning
    """
    binners = []
    for binner, bins in binned.items():
        fastas = glob.glob('%s/*.fa' % bins.replace('${SCRATCH_FOLDER}', ''))
        if self.config.force or not fastas:
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
    binners = get_binners(self, binned)
    cmd, cmd_rm = '', ''
    if binners:
        if contigs.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
            contigs = contigs.rstrip('.gz')
            cmd_rm += 'rm %s\n' % contigs

        if self.soft.params['path']:
            cmd += 'export PATH=$PATH:%s/bin\n' % self.soft.params['path']

        gz, fqs, fqs_cmd = get_fqs(fastq)
        cmd += 'metawrap binning'
        cmd += ' -o %s' % out
        cmd += ' -a %s' % contigs
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -l %s' % self.soft.params['min_contig_length']
        cmd += ' --universal'
        for binner in binners:
            cmd += ' --%s' % binner
        cmd += ' %s\n' % ' '.join(fqs)
        if fqs_cmd:
            cmd = fqs_cmd + cmd
        if gz:
            cmd += 'rm %s\n' % ' '.join(fqs)
        cmd += cmd_rm

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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():

        for software in self.soft.path[::-1]:
            if self.config.tools[software] == 'preprocessing':
                break
        else:
            sys.exit('[metawrap_bining] No "preprocessing" fastqs found')
        fastqs = [fq for sam in self.pools[self.sam_pool][group] for fq in
                  self.softs[software].outputs[sam].get(('illumina', sam), [])]

        to_dos = []
        if not fastqs:
            self.soft.add_status(tech, self.sam_pool, [fastqs], group=group,
                                 message='no illumina data')
            print('[metawrap_binning] No illumina reads for group %s' % group)
            continue
        else:
            to_dos.extend(status_update(self, tech, fastqs, group=group))

        contigs = inputs[0]
        to_dos.extend(status_update(self, tech, [contigs], group=group))

        tmp = '$TMPDIR/mtwrp_%s_%s_%s' % (self.sam_pool, tech, group)
        out = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out)

        binned = {binner: '%s/%s_bins' % (out, binner)
                  for binner in self.soft.params['binners']}
        bin_dirs = sorted(binned.values())  # + ['%s/work_files' % out]

        key = (tech, group)
        if process_outputs(self, key, group, bin_dirs):
            continue

        cmd = binning_cmd(self, fastqs, out, contigs, binned)
        if cmd:
            if to_dos:
                self.outputs['cmds'].setdefault(key, []).append(False)
            else:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=([contigs] + fastqs), i_d=tmp,
                      o_d=bin_dirs, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


def metawrap_(self):
    # for step in [binning, refine, ..]:
    #     step(...)
    pass


def metawrap(self) -> None:
    """MetaWRAP - a flexible pipeline for genome-resolved metagenomic data
    analysis.
    MetaWRAP aims to be an easy-to-use metagenomic wrapper suite that
    accomplishes the core tasks of metagenomic analysis from start to finish:
    read quality control, assembly, visualization, taxonomic profiling,
    extracting draft genomes (binning), and functional annotation. Additionally,
    metaWRAP takes bin extraction and analysis to the next level (see module
    overview below). While there is no single best approach for processing
    metagenomic data, metaWRAP is meant to be a fast and simple approach
    before you delve deeper into parameterization of your analysis. MetaWRAP
    can be applied to a variety of environments, including gut, water, and soil
    microbiomes (see metaWRAP paper for benchmarks). Each individual module of
    metaWRAP is a standalone program, which means you can use only the modules
    you are interested in for your data.

    References
    ----------
    Uritskiy, Gherman V., Jocelyne DiRuggiero, and James Taylor. "MetaWRAP—a
    flexible pipeline for genome-resolved metagenomic data analysis."
    Microbiome 6.1 (2018): 1-13.

    Notes
    -----
    GitHub  : https://github.com/bxlab/metaWRAP
    Paper   : https://doi.org/10.1186/s40168-018-0541-1

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    # This function splits the name of the software and calls as function
    # the last underscore-separated field (which is in this module)
    module_call = caller(self, __name__)
    module_call(self)


def yamb_cmd(
        self,
        tech: str,
        out_dir: str,
        fastxs: list,
        contigs: str
) -> str:
    """Collect the YAMB command line.

    Parameters
    ----------
    self
    tech : str
        Technology
    out_dir : str
        Path to the output folder
    fastxs : list
        Path(s) to the input files
    contigs : str
        Path to the input contigs file

    Returns
    -------
    cmd : str
        Command line for YAMB
    """
    params = tech_params(self, tech)
    cmd, cmd_rm = '', ''
    if contigs.endswith('.gz'):
        cmd += 'gunzip -c %s > %s\n' % (contigs, contigs.rstrip('.gz'))
        contigs = contigs.rstrip('.gz')
        cmd_rm += 'rm %s\n' % contigs

    if params['path']:
        cmd += '%s/yamb.sh' % params['path']
    else:
        cmd += 'yamb.sh'
    cmd += ' -c %s' % contigs
    if len(fastxs) == 3:
        cmd += ' -s %s' % fastxs[0]
        cmd += ' -f %s' % fastxs[1]
        cmd += ' -r %s' % fastxs[2]
    if len(fastxs) == 2:
        cmd += ' -f %s' % fastxs[0]
        cmd += ' -r %s' % fastxs[1]
    if len(fastxs) == 1:
        cmd += ' -f %s' % fastxs[0]
    cmd += ' -o %s' % out_dir
    cmd += ' -t %s' % params['cpus']
    cmd += ' -m %s' % params['m']
    cmd += ' -l %s\n' % params['l']
    cmd += cmd_rm
    cmd += 'for i in %s/*; do gzip -q $i; done\n' % out_dir
    return cmd


def get_yamb(self, tech, group, contigs, fastxs):
    to_dos = []
    if not fastxs:
        self.soft.add_status(tech, self.sam_pool, [fastxs], group=group)
        sys.exit('[YAMB] No input reads for group %s' % group)
    else:
        to_dos.extend(status_update(self, tech, fastxs, group=group))
    to_dos.extend(status_update(self, tech, [contigs], group=group))

    key = (tech, group)
    out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
    self.outputs['dirs'].append(out_dir)
    self.outputs['outs'][key] = [out_dir]

    out = '%s/yamb-pp-*-hdbscan.csv.gz' % out_dir
    if self.config.force or not glob.glob(out):
        if to_dos:
            self.outputs['cmds'].setdefault(key, []).append(False)
        else:
            cmd = yamb_cmd(self, tech, out_dir, fastxs, contigs)
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=([contigs] + fastxs), o_d=out_dir, key=key)
        self.soft.add_status(tech, self.sam_pool, 1, group=group)
    else:
        self.soft.add_status(tech, self.sam_pool, 0, group=group)


def yamb(self):
    """
    YAMB (Yet Another Metagenome Binner) - semi-automatic pipeline for
    metagenomic contigs binning. It's based on t-Distributed Stochastic
    Neighbor Embedding (t-SNE) technique for dimensionality reduction,
    that uses tetranucleotide frequency distribution in metagenomic contigs
    and mean contig coverage, and HDBSCAN clustering method.

    References
    ----------
    Korzhenkov, A., 2019. YAMB: metagenome binning using nonlinear
    dimensionality reduction and density-based clustering. bioRxiv, p.521286.

    Notes
    -----
    GitHub  : https://github.com/laxeye/YAMB
    Paper   : https://doi.org/10.1101/521286

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    reads = self.config.fastq_mv
    if '_' in self.soft.name:
        reads_tool = self.soft.name.split('_')[-1]
        reads = self.softs[reads_tool].outputs

    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        fastxs = []
        if tech in self.config.techs:
            fastxs += [fastq for sam in self.pools[self.sam_pool][group]
                       for fastq in reads[sam].get((tech, sam), [])]
        else:
            for t in self.config.techs:
                if t in tech:
                    fastxs += [fastq for sam in self.pools[self.sam_pool][group]
                               for fastq in reads[sam].get((t, sam), [])]
        get_yamb(self, tech, group, inputs[0], fastxs)


def binspreader_cmd(
        self,
        tech: str,
        graph: str,
        inputs: dict,
        group: str,
        out_dir: str
):
    """Collect the command line for BinSPreader.

    Parameters
    ----------
    self
    tech : str
    graph : str
    inputs : dict
    group : str
    out_dir : str

    Returns
    -------
    cmd : str
        BinSPreader command line
    """
    params = tech_params(self, tech)
    contigs_per_bin = '%s/contigs_per_bin.tsv' % out_dir
    cmd = '%s/contigs_per_bin.py -i %s -m %s -o %s\n' % (
        SCRIPTS, inputs, self.soft.prev, contigs_per_bin)

    if params['path']:
        cmd += '%s/bin-refine' % params['path']
    else:
        cmd += 'bin-refine'
    cmd += ' %s' % graph
    cmd += ' %s' % contigs_per_bin
    cmd += ' %s' % out_dir
    cmd += ' -t %s' % params['cpus']
    # cmd += ' --paths %s' % paths
    if params['dataset']:
        cmd += ' --dataset %s' % params['dataset']
    for param in ['l', 'e', 'n', 'la']:
        cmd += ' -%s %s' % (param, params[param])
    for param in ['metaalpha', 'bin_weight']:
        cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
    for boolean in ['m', 'cami', 'zero_bin', 'tall_multi', 'bin_dist',
                    'sparse_propagation', 'no_unbinned_bin', 'reads',
                    'length_threshold', 'distance_bound']:
        if params[boolean]:
            if len(boolean) > 1:
                cmd += ' --%s %s' % (param.replace('_', '-'), params[param])
            else:
                cmd += ' -%s %s' % (param, params[param])
    if params['Smax']:
        cmd += ' -Smax'
    else:
        cmd += ' -Smle'
    if params['Rcorr']:
        cmd += ' -Rcorr'
    else:
        cmd += ' -Rprop'

    tmp_dir = '$TMPDIR/bnsprdr_%s_%s_%s' % (self.sam_pool, tech, group)
    cmd += ' --tmp-dir %s' % tmp_dir
    return cmd


def binspreader(self):
    """
    BinSPreader is a novel tool that attempts to refine metagenome-assembled
    genomes (MAGs) obtained from existing tools. BinSPreader exploits the
    assembly graph topology and other connectivity information, such as
    paired-end and Hi-C reads, to refine the existing binning, correct
    binning errors, propagate binning from longer contigs to shorter contigs
    and infer contigs belonging to multiple bins.

    References
    ----------
    Tolstoganov, I., Kamenev, Y., Kruglikov, R., Ochkalova, S. and
    Korobeynikov, A., 2022. Binspreader: refine binning results for fuller
    mag reconstruction. Iscience, 25(8), p.104770.

    Notes
    -----
    Code    : https://github.com/ablab/spades/releases/tag/binspreader-recombseq
    Docs    : https://cab.spbu.ru/software/binspreader/
    Paper   : https://doi.org/10.1016/j.isci.2022.104770

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    if self.config.tools[self.soft.prev] == 'binning':
        assembly = get_assembly(self)
        for (tech, group), inputs in self.inputs[self.sam_pool].items():

            graph = get_assembly_graph(self, tech, group, assembly)

            inp = inputs[0]
            to_dos = status_update(self, tech, [inp], group=group, folder=True)
            to_dos.extend(status_update(self, tech, [graph], group=group))

            key = (tech, group)
            out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
            self.outputs['dirs'].append(out_dir)
            self.outputs['outs'][key] = [out_dir]

            if self.config.force or to_do('%s/binning.tsv' % out_dir):
                if to_dos:
                    self.outputs['cmds'].setdefault(key, []).append(False)
                else:
                    cmd = binspreader_cmd(self, tech, graph, inputs,
                                          group, out_dir)
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, i_f=graph, i_d=inp, o_d=out_dir, key=key)
                self.soft.add_status(tech, self.sam_pool, 1, group=group)
            else:
                self.soft.add_status(tech, self.sam_pool, 0, group=group)


def mycc(self):
    """
    Automated binning tool that combines genomic signatures, marker genes and
    optional contig coverages within one or multiple samples, in order to
    visualize the metagenomes and to identify the reconstructed genomic
    fragments.

    References
    ----------
    Lin, H.H. and Liao, Y.C., 2016. Accurate binning of metagenomic
    contigs via automated clustering sequences using information of genomic
    signatures and marker genes. Scientific reports, 6(1), pp.1-8.

    Notes
    -----
    Code    : https://sourceforge.net/projects/sb2nhri/files/MyCC/
    Docs    : http://sourceforge.net/projects/sb2nhri/files/MyCC/.
    Paper   : https://doi.org/10.1038/srep24175

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    pass


def metabinner(self):
    """MetaBinner consists of two modules: 1) “Component module” includes
    steps 1-4, developed for generating high-quality, diverse component
    binning results; and 2) “Ensemble module” includes step 5, developed for
    recovering individual genomes from the component binning results.
    MetaBinner is an ensemble binning method, but it does not need the
    outputs of other individual binners. Instead, MetaBinner generates
    multiple high-quality component binning results based on the proposed
    “partial seed” method, for further integration. Please see our manuscript
    for the details.

    References
    ----------
    Wang, Z., Huang, P., You, R., Sun, F. and Zhu, S., 2021. MetaBinner: a
    high-performance and stand-alone ensemble binning method to recover
    individual genomes from complex microbial communities. bioRxiv.

    Notes
    -----
    GitHub  : https://github.com/ziyewang/MetaBinner
    Paper   : https://doi.org/10.1101/2021.07.25.453671
    Docs    : http://www.cecill.info/

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def semibin(self):
    """Command tool for metagenomic binning with semi-supervised deep
    learning using information from reference genomes in Linux and MacOS.

    References
    ----------
    Pan, S., Zhu, C., Zhao, X.M. and Coelho, L.P., 2022. A deep siamese
    neural network improves metagenome-assembled genomes in microbiome
    datasets across different environments. Nature communications, 13(1),
    pp.1-12.

    Notes
    -----
    GitHub  : https://github.com/BigDataBiology/SemiBin/
    Paper   : https://doi.org/10.1038/s41467-022-29843-y

    Parameters
    ----------
    self : Commands class instance
    """
    pass


def vamb(self):
    """Vamb is a metagenomic binner which feeds sequence composition
    information from a contig catalogue and co-abundance information from BAM
    files into a variational autoencoder and clusters the latent representation.
    It performs excellently with multiple samples, and pretty good on
    single-sample data. Vamb is implemented purely in Python (with a little
    bit of Cython) and can be used both from command line and from within a
    Python interpreter.

    References
    ----------
    Nissen, J.N., Johansen, J., Allesøe, R.L., Sønderby, C.K., Armenteros,
    J.J.A., Grønbech, C.H., Jensen, L.J., Nielsen, H.B., Petersen, T.N.,
    Winther, O. and Rasmussen, S., 2021. Improved metagenome binning and
    assembly using deep variational autoencoders. Nature biotechnology,
    39(5), pp.555-560.

    Notes
    -----
    GitHub  : https://github.com/RasmussenLab/vamb
    Paper   : https://doi.org/10.1038/s41587-020-00777-4

    Parameters
    ----------
    self : Commands class instance
    """
    pass
