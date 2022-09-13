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
from os.path import basename, splitext

RESOURCES = pkg_resources.resource_filename("metagenomix", "resources")


def files_to_show_list(
        self,
        s: str,
        tech: str,
        t: str,
        seps: list,
        fs: set,
        ins
) -> None:
    fastq_mv = self.config.fastq_mv
    seps.append('│%s' % (' ' * len(t)))
    for i in ins:
        if '/%s/' % tech in i:
            fs.add(i)
        elif s in fastq_mv and i in fastq_mv[s][(tech, s)]:
            fs.add(i)
        elif isinstance(i, tuple) and '/%s/' % tech in i[0]:
            fs.add('%s (%s)' % i)


def files_to_show_dict(
        t: str,
        ins: dict,
        seps: list,
        fs: set
) -> None:
    seps.append('│%s' % (' ' * len(t)))
    fs.update(['%s (%s)' % (y, x) for x, y in ins.items()])


def fill_seps_and_fs(
        self,
        inputs: dict,
        sam: str,
        tech: str,
        t: str,
        seps: list,
        fs: set
) -> None:
    if (t, sam) in inputs[sam]:
        ins = inputs[sam][(t, sam)]
    elif t in inputs[sam]:
        ins = inputs[sam][t]
    if isinstance(ins, list) and len(ins):
        files_to_show_list(self, sam, tech, t, seps, fs, ins)
    elif isinstance(ins, dict) and [x for x in ins.values()
                                    if '/%s/' % tech in x]:
        files_to_show_dict(t, ins, seps, fs)
    else:
        no_file_to_show(t, seps)


def no_file_to_show(
        t: str,
        seps: list
) -> None:
    seps.append('X%s' % (' ' * len(t)))


def files_to_show(
        self,
        inputs: dict,
        techs: list,
        sam: str,
        tech: str
) -> tuple:
    seps, fs = [], set()
    techs = techs[:1 + techs.index(tech)]
    for t in techs:
        if (t, sam) in inputs[sam]:
            fill_seps_and_fs(self, inputs, sam, tech, t, seps, fs)
        elif t in inputs[sam]:
            fill_seps_and_fs(self, inputs, sam, tech, t, seps, fs)
        else:
            no_file_to_show(t, seps)
    sep = ''.join(seps)
    return sep, fs


def get_inputs_to_show(self) -> tuple:
    sams = []
    inputs = {}
    techs = set()
    if set(self.inputs).issubset(set(self.pools)):
        col = 'pool (group)'
        for pool, group_inputs in self.inputs.items():
            for group, files in group_inputs.items():
                if set(files).issubset(self.config.techs):
                    pool_group = '%s (%s)' % (pool, group)
                    sams.append(pool_group)
                    inputs[pool_group] = files
                    techs.update(self.config.techs)
                elif isinstance(group, tuple):
                    pool_group = '%s (%s)' % (pool, group[1])
                    sams.append(pool_group)
                    if pool_group not in inputs:
                        inputs[pool_group] = {}
                    inputs[pool_group][group[0]] = files
                    techs.add(group[0])
        if not inputs:
            col = 'co-assembly'
            techs.update([y for x in self.inputs.values() for y in x])
            sams = sorted(self.inputs)
            inputs = dict(self.inputs)
    else:
        col = 'sample'
        sams = sorted(self.inputs.keys())
        inputs = dict(self.inputs)
        techs.update(self.config.techs)
    techs = sorted(techs)
    return inputs, techs, sams, col


def show_inputs(self):
    if self.config.verbose:
        inputs, techs, sams, c = get_inputs_to_show(self)
        mlen = max([len(x) for x in sams])
        print('\n%s\n[%s] inputs (after %s)\n%s\n' % (
            ('-' * 30), self.soft.name, self.soft.prev, ('-' * 30)))
        print('%s%s %s' % (c, ' ' * (mlen-len(c)), ' '.join(techs)))
        for sam in sams:
            if sam not in inputs:
                continue
            show_sam = True
            s = '%s%s' % (sam, ' ' * (mlen-len(sam)))
            for tdx, tech in enumerate(techs[::-1]):
                sep, fs = files_to_show(self, inputs, techs, sam, tech)
                if sep.rstrip()[-1] != 'X':
                    sep = '%s └─' % sep.rstrip()[:-2]
                for f in fs:
                    if show_sam:
                        print('%s %s %s' % (s, sep.strip(), f))
                        show_sam = False
                    else:
                        print('%s %s %s' % ((' ' * len(s)), sep.strip(), f))
            print()


def get_post_prot_dir(
        self,
        input_fp: str,
) -> str:
    """Edit the "after_TOOL" to "after_TOOL1_TOOL2" with "TOOL1" and "TOOL2"
    being the two last tool run one after the other.

    Notes
    -----
    This function is (currently) only called for the softwares that run after
    Prodigal (and in the future, after other protein-prediction softwares),
    since the protein prediction throughout the pipeline can be called after
    different softwares (e.g., after metawrap_refine, after drep, or, for long
    read data, after error correction). Hence, editing the output directory
    to reflect the provenance and thus obtain different output folders if
    the pipeline is re-run for different paths at the Prodigal step.

    In the future, the provenance as labelled by the output folder's "after_"
    nomenclature, will be replaced by hashing key for the concatenation of
    the names of the software (and their parameters) that led to the output.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
    input_fp : str
        Path to the input file or folder

    Returns
    -------
    out_dir : str
        Path to the main output folder
    """
    if 'after_' in input_fp:
        after = input_fp.split('/after_')[-1].split('/')[0]
        out_dir = self.dir.replace('after', 'after_%s' % after)
    else:
        out_dir = self.dir
    return out_dir


def genome_key(
        tech: str,
        sam_group: str,
        genome: str = ''
) -> tuple:
    """Get a concatenation of the variables fields to split the commands as
    separate jobs later when writing out.

    Parameters
    ----------
    tech : str
        Technology: 'illumina', 'pacbio', 'nanopore', or hybrif naming
    sam_group : str
        Sample name or group for the current co-assembly
    genome : str
        Name of the genome/MAG (or empty for single-file outputs)

    Returns
    -------
    key : tuple
        Variables names for the current analytic level
    """
    key = (tech, sam_group)
    if genome:
        key += (genome,)
    return key


def genome_out_dir(
        self,
        tech: str,
        fasta: str,
        sam_group: str,
        genome: str = ''
) -> str:
    """Get a concatenation of the variables fields to split the commands as
    separate jobs later when writing out.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for MacSyFinder
        .sam_pool : str
            Sample or co-assembly name
    tech : str
        Technology: 'illumina', 'pacbio', 'nanopore', or hybrif naming
    fasta : str
        Path to the input fasta file
    sam_group : str
        Sample name or group for the current co-assembly
    genome : str
        Name of the genome/MAG (or empty for single-file outputs)

    Returns
    -------
    key : str
        Concatenation of variables names for the current analytic level
    """
    # post_prot_dir = get_post_prot_dir(self, fasta)
    # self.soft.dir = post_prot_dir.replace('${SCRATCH_FOLDER}', '')

    # out_dir = '/'.join([post_prot_dir, tech])
    out_dir = '/'.join([self.dir, tech])
    if self.sam_pool:
        out_dir += '/%s' % self.sam_pool
    if sam_group != self.sam_pool:
        out_dir += '/%s' % sam_group
    if genome:
        out_dir += '/' + genome
    return out_dir


def check_input(
        self,
        tech: str,
        raw: bool
) -> bool:
    """Checks whether to skip annotation for the current per-sample illumina
    data.

    Parameters
    ----------
    self : Commands class instance
        .soft.name : str
            Current software in the pipeline
        .soft.prev : str
            Previous software in the pipeline
        .softs : dict
            All softwares' Soft() class instances
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    raw : bool
        Whether illumina data cna be processed as raw/per-sample paired reads

    Returns
    -------
    bool
        Whether to skip annotation for the current per-sample illumina data
    """
    if tech == 'illumina':
        prev, name = self.soft.prev, self.soft.name
        tools = ['plass']
        if prev not in tools and not raw:
            print('[%s] illumina annotation only after %s' % (name, tools))
            return True
        typ = 'nuclassemble'
        if name in ['search_diamond', 'search_hmmer', 'macsyfinder']:
            typ = 'assemble'
        is_type = (self.softs['plass'].params['type'] != typ)
        if prev == 'plass' and is_type and not raw:
            print('[%s] illumina annotation if plass type: %s' % (name, typ))
            return True
    return False


def sample_inputs(
        self,
        techs=None,
        raw: bool = False,
        index: int = 1
) -> dict:
    """Get per-sample input files per technology, depending on the current
    softwares.

    Notes
    -----
    - For Prodigal's gene prediction, the input must be nucleotides and could
      be either the nucleotide assembly of Plass "nuclassemble" (if run after
      Plass), or the per-sample long read data (i.e., pacbio or nanopore).
    - For any other tool, the inputs would be one of the outputs (the first
      of the list by default) form the previous tool for all three technologies.
    The returned files are in a dict since for inputs after co-assembly,
    the mechanism is to iterate over a dict (genome/MAG (key) to fasta (value)).
    Hence, this mechanism can be used: here it iterates over one file.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to software's output folder
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    techs : list
        Technologies possibly used per sample for the current softwares
    raw : bool
        Whether illumina data cna be processed as raw/per-sample paired reads
    index : int
        Index in the list of input files (i.e. outputs from the previous step)

    Returns
    -------
    inputs : dict
        Paths to inputs per-sample input files per technology
    """
    if techs is None:
        techs = ['illumina', 'pacbio', 'nanopore']

    inputs = {}
    for tech in techs:
        if self.inputs[self.sam_pool].get((tech, self.sam_pool)):
            if check_input(self, tech, raw):
                continue
            fastxs = self.inputs[self.sam_pool][(tech, self.sam_pool)]
            if fastxs:
                inputs[tech] = {'': fastxs}
    return inputs


def group_inputs(
        self,
        inputs: list,
        folder: bool = False,
        index: int = 0
) -> dict:
    """Get the input fasta files for a step that takes as input
    either the contigs of an assembly, the genomes of a binning or after
    dereplication, or either nucleotide (protein-coding) or protein sequences.

    Parameters
    ----------
    self : Commands class instance
        .prev : str
            Previous software in the pipeline
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
    inputs: list
        Path(s) to the input files(s) or folder(s)
    folder : bool
        Whether the returned paths should be the folders (not the files)
    index : int
        The index of the fasta file for the current input

    Returns
    -------
    fastas : dict
        Path(s) to the fasta file(s) for an assembly or a series of MAGs
    """
    fastas = {}
    if self.soft.prev in self.config.tools['assembling']:
        fastas[''] = [inputs[index]]
    elif self.soft.prev in ['metawrap_refine', 'drep']:
        for fa in genomes_fastas(self, inputs[index], folder):
            fastas[splitext(basename(fa))[0]] = [fa]
    elif self.soft.prev == 'metawrap_reassemble':
        for inp in inputs:
            sam, stringency = inp.rsplit('/')[-2:]
            for fa in genomes_fastas(self, inp, folder):
                key = '%s/%s/%s' % (sam, stringency, splitext(basename(fa))[0])
                if folder:
                    key = '%s/%s' % (sam, stringency)
                fastas[key] = [fa]
    elif self.soft.prev == 'prodigal':
        if len(inputs) == 1:
            fastas[''] = [inputs[index]]
        else:
            for inp in inputs:
                if 'after_metawrap_reassemble' in inp:
                    fastas['/'.join(inp.rsplit('/', 3)[-3:])] = [inp]
                else:
                    fastas[basename(inp)] = [inp]
    elif self.soft.name.startswith('checkm_'):
        # fastas[''] = inputs
        fastas = inputs
    else:
        print(self.soft.prev)
        print(inputs)
        sys.exit('[check group_inputs()]')
    return fastas


def genomes_fastas(
        self,
        path: str,
        folder: bool,
        ext: str = '.fa'
) -> list:
    """Collect the paths to the input fasta files contained in the folder
    of a binning or dereplication result.

    Parameters
    ----------
    self : Commands class instance
        .config
            Configurations
    path : str
        Path to the input folder for the fasta files of multiple genomes
    folder : bool
        Whether the returned paths should be the folders (not the files)
    ext : str
        Extension to the genomes fasta files

    Returns
    -------
    genome_paths : list
        Paths to the input fasta files
    """
    if folder:
        return [path]
    if self.config.dev:
        genome_paths = ['%s/a%s' % (path, ext), '%s/b%s' % (path, ext)]
    else:
        genome_paths = glob.glob('%s/*%s' % (path, ext))
    return genome_paths


def get_extension(self) -> str:
    """Get the `.fa` or `.fna` fasta filename extension.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev : str
            Previous software in the pipeline

    Returns
    -------
    ext : str
        Fasta filename extension
    """
    ext = 'fa'
    if self.soft.prev == 'yamb':
        ext = 'fna'
    return ext


def get_assembler(self) -> str:
    """Get the name of the assembler used for "assembling" by parsing the
    pipeline steps that led to this software.

    Parameters
    ----------
    self : Commands class instance
        .graph : dict
            Pipeline paths
        .soft.name
            Current software in the pipeline
        .config
            Configurations

    Returns
    -------
    assembler : str
        Name of the assembler that produced the genomes/MAGs
    """
    steps = self.graph.paths[self.soft.name][0]
    for step in steps:
        if 'assembling' in self.config.tools.get(step, []):
            return step
    tools = '%s/softwares.txt' % RESOURCES
    sys.exit('[%s] None of "%s" is an "assembling" tool (see %s)' % (
        self.soft.name, '", "'.join(steps), tools))


def get_reads(self) -> dict:
    """Get the output dict structure containing the per-sample reads for the
    last step performed before assembly.

    Parameters
    ----------
    self : Commands class instance
        .graph : dict
            Path of the pipeline steps to the current software
        .soft.name : str
            Current software in the pipeline
        .softs : dict
            Software class instances
        .config
            Configurations

    Returns
    -------
    reads : dict
        per-sample reads for the last step performed before assembly
    """
    steps = self.graph.paths[self.soft.name][0]
    steps_before_pooling = steps[:steps.index('pooling')]
    step = 'fastq'
    for step_before_pooling in steps_before_pooling:
        if self.config.tools[step_before_pooling] == 'paired read merging':
            break
        step = step_before_pooling
    reads = self.softs[step].outputs
    return reads


def add_folder(
        self,
        name: str,
        out_dir: str,
        step: str = ''
) -> str:
    """Appended the lorikeet output folder name with the check module name
    if the module is not called explicitly in the pipeline. That is, if
    the pipeline specifies "drep lorikeet" and not "drep lorikeet_call".

    Parameters
    ----------
    self
    name : str
        Name of the software
    out_dir : str
        Path to the output folder
    step : str
        Name of the current check module

    Returns
    -------
    out : str
        Path to the output folder possibly appended with module name
    """
    out = out_dir
    if self.soft.name == name:
        out += '/%s' % step
    return out
