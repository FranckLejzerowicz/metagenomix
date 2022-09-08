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
        fastq_fps: list
) -> tuple:
    """Get the fastq files for the current sample to reassemble

    Parameters
    ----------
    fastq_fps : list
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
    for f in fastq_fps:
        if f.endswith('.gz'):
            cmd += 'gunzip -c %s > %s\n' % (f, f.replace('.gz', ''))
            fqs.append(f.replace('.gz', ''))
        else:
            fqs.append(f)
    return fqs, cmd


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
        Paths to the input contifs fasta file

    Returns
    -------
    cmd : str
        metaWRAP quantify command
    """
    fqs, fqs_cmd = get_fqs(fastq)
    cmd = 'metawrap quant_bins'
    cmd += ' -b %s' % bins
    cmd += ' -o %s' % out
    cmd += ' -a %s' % contigs
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' %s\n' % ' '.join(fqs)
    if fqs_cmd:
        cmd = fqs_cmd + cmd + 'rm %s\n' % ' '.join(fqs)
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
    fastqs = {}
    for sam in self.pools[self.sam_pool][group]:
        if self.config.fastq_mv[sam].get(('illumina', sam)):
            fastqs[sam] = self.config.fastq_mv[sam][('illumina', sam)]
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
                if sam:
                    out = out_ + '/%s' % sam
                else:
                    out = out_
                self.outputs['dirs'].append(out)

                cmd = quantify_cmd(self, bins, fastq, out, contigs)
                if self.config.force or glob.glob('%s/*' % out):
                    self.outputs['outs'][key][mode][sam] = out
                    self.outputs['cmds'].setdefault(key, []).append(cmd)
                    io_update(self, i_f=fastq, o_d=out, key=key)
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
    cmd = 'metawrap %s' % command
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

                if to_do(folder=bin_dir):
                    self.soft.add_status(
                        tech, self.sam_pool, [bin_dir], group=group, genome=sam)

                cmd = classify_or_annotate_cmd(self, command, bin_dir, out)
                self.outputs['outs'][key][sam] = [out]

                condition = classify_or_annotate_condition(command, out)
                if self.config.force or condition:
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
    if command == 'annotate_bins':
        tab = '%s/bin_funct_annotations/*.tab' % out.replace(
            '${SCRATCH_FOLDER}', '')
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
    fqs, fqs_cmd = get_fqs(fastq_fps)
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)

        self.outputs['outs'][key] = {}
        bins = inputs[0]

        assembler = self.softs['metawrap_binning'].prev
        contigs = self.softs[assembler].outputs[self.sam_pool][key][0]
        fastqs = get_fastqs(self, group)
        if not fastqs:
            self.soft.add_status(tech, self.sam_pool, [fastqs], group=group)
            print('[metawrap_blobology] No illumina reads for group %s' % group)
            continue

        io_update(self, i_f=contigs, i_d=bins, key=key)
        for mode in self.soft.params['blobology']:
            self.outputs['outs'][key][mode] = {}
            out = '/'.join([self.dir, tech, self.sam_pool, group, mode])
            sams_fastqs = get_sams_fastqs(mode, fastqs)
            for sam, fastq_fps in sams_fastqs.items():
                if sam:
                    out += '/%s' % sam
                self.outputs['dirs'].append(out)
                plot = '%s/contigs.binned.blobplot' % out
                self.outputs['outs'][key][mode][sam] = out
                if self.config.force or to_do(plot):
                    cmd = blobology_cmd(self, fastq_fps, out, contigs, bins)
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
    fqs, fqs_cmd = get_fqs(fastq)
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
        if not self.config.force and not to_do(folder=mode_dir):
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
        if to_do(folder=bins):
            self.soft.add_status(tech, self.sam_pool, [bins], group=group)

        fastqs = get_fastqs(self, group)
        for sam, fastq in fastqs.items():
            if not fastq:
                self.soft.add_status(tech, self.sam_pool, [fastq], group=group,
                                     message='no illumina data')
                print('[metawrap_reassemble] No illumina reads for %s' % sam)
                continue

            out = '/'.join([self.dir, tech, self.sam_pool, group, sam])
            self.outputs['dirs'].append(out)

            cmd = reassembly_cmd(self, fastq, bins, tech, group, out, key)
            if cmd:
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)

        out_dir = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out_dir)

        bin_folders = inputs[:-1]
        to_dos = [x for x in bin_folders if to_do(folder=x)]
        if to_dos:
            self.soft.add_status(tech, self.sam_pool, to_dos, group=group)

        cmd, n_bins = refine_cmd(self, out_dir, bin_folders)
        out = '%s/metawrap_%s_%s' % (out_dir,
                                     self.soft.params['min_completion'],
                                     self.soft.params['min_contamination'])
        stats, bins = '%s.stats' % out, '%s_bins' % out
        self.outputs['outs'][key] = [bins, stats]
        if not self.config.force and not to_do(folder=bins):
            self.soft.add_status(tech, self.sam_pool, 0, group=group)
            continue
        if n_bins == len(bin_folders):
            self.outputs['cmds'].setdefault(key, []).append(cmd)
        io_update(self, i_d=bin_folders, o_f=stats, o_d=bins, key=key)
        self.soft.add_status(tech, self.sam_pool, 1, group=group)


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
        bin_0 = '%s/bin.0.fa' % bins
        bin_1 = '%s/bin.1.fa' % bins
        if self.config.force or (to_do(bin_0) and to_do(bin_1)):
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
    cmd = ''
    if binners:
        fqs, fqs_cmd = get_fqs(fastq)
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
    for (tech, group), inputs in self.inputs[self.sam_pool].items():
        key = (tech, group)

        fastqs = [fastq for sam in self.pools[self.sam_pool][group] for fastq
                  in self.config.fastq_mv[sam].get(('illumina', sam), [])]
        if not fastqs:
            self.soft.add_status(tech, self.sam_pool, [fastqs], group=group,
                                 message='no illumina data')
            print('[metawrap_binning] No illumina reads for group %s' % group)
            continue

        tmp = '$TMPDIR/mtwrp_%s_%s_%s' % (self.sam_pool, tech, group)
        out = '/'.join([self.dir, tech, self.sam_pool, group])
        self.outputs['dirs'].append(out)

        binned = {binner: '%s/%s_bins' % (out, binner)
                  for binner in self.soft.params['binners']}
        bin_dirs = sorted(binned.values()) + ['work_files']
        self.outputs['outs'][key] = bin_dirs

        cmd = binning_cmd(self, fastqs, out, inputs[0], binned)
        if cmd:
            self.outputs['cmds'].setdefault(key, []).append(cmd)
            io_update(self, i_f=([inputs[0]] + fastqs), i_d=tmp,
                      o_d=bin_dirs, key=key)
            self.soft.add_status(tech, self.sam_pool, 1, group=group)
        else:
            self.soft.add_status(tech, self.sam_pool, 0, group=group)


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
