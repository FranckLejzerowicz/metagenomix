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
        elif s in fastq_mv and i in fastq_mv[s].get((tech, s)):
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
        if (t, sam) in inputs[sam] or t in inputs[sam]:
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


def show_inputs(self, log):
    if self.config.verbose:
        inputs, techs, sams, c = get_inputs_to_show(self)
        mlen = max([len(x) for x in sams])
        log.info('\n%s\n[%s] inputs (after %s)\n%s\n' % (
            ('-' * 30), self.soft.name, self.soft.prev, ('-' * 30)))
        log.info('%s%s %s' % (c, ' ' * (mlen-len(c)), ' '.join(techs)))
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
                        log.info('%s %s %s' % (s, sep.strip(), f))
                        show_sam = False
                    else:
                        log.info('%s %s %s' % ((' ' * len(s)), sep.strip(), f))
            log.info()


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
        genome: str = '',
        array: int = 0
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
    array : int
        Number of array jobs to spawn

    Returns
    -------
    key : tuple
        ((Variables names for the current analytic level), number of arrays)
    """
    key = (tech, sam_group)
    if genome:
        key += (genome,)
    key = (key, array)
    return key


def genome_out_dir(
        self,
        tech: str,
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
    sam_group : str
        Sample name or group for the current co-assembly
    genome : str
        Name of the genome/MAG (or empty for single-file outputs)

    Returns
    -------
    key : str
        Concatenation of variables names for the current analytic level
    """
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
        Whether illumina data can be processed as raw/per-sample paired reads

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
        if not raw:
            prev_h = self.hashes[tuple(self.soft.path[:-1])]
            is_type = (self.softs['plass'][prev_h].params['type'] != typ)
            if prev == 'plass' and is_type:
                message = 'Illumina annotation if plass type: %s' % typ
                self.soft.messages.add(message)
                return True
    return False


def sample_inputs(
        self,
        techs=None,
        raw: bool = False,
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

    Returns
    -------
    inputs : dict
        Paths to inputs per-sample input files and per technology
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
        index: int = 0,
        target: str = None
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
    target : str
        Previous tool that is the actual data target

    Returns
    -------
    fastas : dict
        Path(s) to the fasta file(s) for an assembly or a series of MAGs
    """
    fastas = {}
    if target is None:
        prev = self.soft.prev
    else:
        prev = target

    if prev in self.config.tools['assembling']:
        fastas[''] = [inputs[index]]
    elif prev in ['metawrap_refine', 'drep']:
        for fa in genomes_fastas(self, inputs[index], folder):
            fastas[splitext(basename(fa))[0]] = [fa]
    elif prev == 'metawrap_reassemble':
        for inp in inputs:
            sam, stringency = inp.rsplit('/')[-2:]
            for fa in genomes_fastas(self, inp, folder):
                key = '%s/%s/%s' % (sam, stringency, splitext(basename(fa))[0])
                if folder:
                    key = '%s/%s' % (sam, stringency)
                fastas[key] = [fa]
    elif prev == 'prodigal':
        if len(inputs) == 1:
            fastas[''] = [inputs[index]]
        else:
            for inp in inputs:
                if 'after_metawrap_reassemble' in inp:
                    fastas['/'.join(inp.rsplit('/', 3)[-3:])] = [inp]
                else:
                    fastas[basename(inp)] = [inp]
    # elif self.soft.name.startswith('checkm_'):
    #     fastas = inputs
    else:
        fastas = inputs
        # print(prev)
        # print(inputs)
        # sys.exit('[check group_inputs()]')
    return fastas


def glob_scratched(
        path: str,
        ext: str
) -> list:
    """

    Parameters
    ----------
    path : str
    ext : str

    Returns
    -------
    globbed : list
    """
    globbed = []
    path_ = path.replace('${SCRATCH_FOLDER}', '')
    for fp in glob.glob('%s/*%s' % (path_, ext)):
        fpo = '%s/%s' % (path, basename(fp))
        globbed.append(fpo)
    return globbed


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
        genome_paths = glob_scratched(path, ext)
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


def get_assembler(self) -> tuple:
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
    sdx : int
        Index in the path that corresponds to the assembler
    assembler : str
        Name of the assembler that produced the genomes/MAGs
    """
    for sdx, step in enumerate(self.soft.path):
        if 'assembling' in self.config.tools.get(step, []):
            return sdx, step
    tools = '%s/softwares.txt' % RESOURCES
    sys.exit('[%s] None of "%s" is an "assembling" tool (see %s)' % (
        self.soft.name, '", "'.join(self.soft.path), tools))


def get_contigs_from_path(self, tech, group, ret_list=False):
    """Return the path to the contigs file that was generated by the workflow.

    Parameters
    ----------
    self : Commands class instance
        .graph : dict
            Pipeline paths
        .soft.name
            Current software in the pipeline
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    group : str
        Name of a co-assembly pool's group
    ret_list : bool
        Whether to return a list of single or all (tech, group) contigs

    Returns
    -------
    contigs : str
        Path to the contigs file
    """
    adx, assembler = get_assembler(self)
    assembler_hash = self.hashes[tuple(self.path[:(adx + 1)])]
    assembly = self.softs[assembler][assembler_hash].outputs[self.sam_pool]
    if ret_list:
        if (tech, group) in assembly:
            contigs = assembly[(tech, group)][:1]
        else:
            groups = self.pools[self.sam_pool]
            contigs = [assembly[(tech, x)][0] for x in groups]
    else:
        contigs = assembly[(tech, group)][0]
    return contigs


def get_assembly_graph(
        self,
        tech: str,
        group: str,
) -> str:
    """Get path to an assembly graph.

    Parameters
    ----------
    self
    tech : str
        Technology/ies iterated over
    group : str
        Group for the current co-assembly and

    Returns
    -------
    graph : str
        Paths to the assembly graph file
    """
    _, assembler = get_assembler(self)
    contigs = get_contigs_from_path(self, tech, group)
    if assembler == 'spades':
        graphs = '%s/assembly_graph_with_scaffolds.gfa.gz' % dirname(contigs)
    elif assembler == 'megahit':
        graphs = '%s/*.contigs.fg.gz' % dirname(contigs)
    else:
        sys.exit('[%s] Not avail for "%s" assembler' % (self.soft.name,
                                                        assembler))
    if not self.config.dev and not glob.glob(graphs):
        sys.exit('[%s] No assembly graph for "%s"' % (self.soft.name,
                                                      assembler))
    return graphs[-1]


def find_software_path(self, soft):
    if soft not in self.graph.paths:
        sys.exit('[%s] "%s" not in the pipeline' % soft)
    ref = '>'.join(self.path)
    for path in self.graph.paths[soft]:
        j = '>'.join(path[:-1])
        if ref.startswith(j) and self.hashes[tuple(path)] in self.softs[soft]:
            return path
    sys.exit('[%s] "%s" not in a workflow related path' % (soft, self.name))


def get_reads(self, step: str = 'preprocessing', soft: str = None) -> dict:
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
    step : str
        Category of the step (default: preprocessed reads)
    soft: str
        Software to get the reads from (default: None)

    Returns
    -------
    reads : dict
        per-sample reads for the last step performed before assembly
    """
    reads = self.config.fastq_mv
    hdx, has_soft, has_step = None, '', ''
    for sdx, software in enumerate(self.soft.path[::-1]):
        if software == soft and not has_soft:
            hdx, has_soft = sdx, software
        if self.config.tools.get(software) == step and not has_step:
            hdx, has_step = sdx, software
    if has_soft:
        softs = tuple(self.soft.path[:(self.soft.path.index(has_soft) + 1)])
        reads = self.softs[has_soft][self.hashes[softs]].outputs
    elif has_step:
        softs = tuple(self.soft.path[:(self.soft.path.index(has_step) + 1)])
        reads = self.softs[has_step][self.hashes[softs]].outputs
    return reads


def get_contigs(
        self,
        tech: str,
        pool: str
) -> dict:
    """Collect the paths to the contigs for all the co-assembly groups of the
    current co-assembly pool name, and those that are not yet generated.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
        .soft.params
            Parameters
        .config
            Configurations
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pool : str
        Name of the co-assembly pool

    Returns
    -------
    contigs : dict
        Paths to the contigs assembly files per co-assembly group
    """
    contigs = {}
    for group, assembly_outputs in sorted(self.inputs[pool].items()):
        if group[0] != tech:
            continue
        contigs[group[1]] = assembly_outputs[0]
    return contigs


def get_sams(self, group: str, pool: str = None) -> list:
    """Get the sample name(s).

    Parameters
    ----------
    self
    group : str
        Name of a co-assembly pool's group
    pool : str
        Name of the co-assembly

    Returns
    -------
    sams : list
        Sample name(s)
    """
    if pool:
        pools = self.pools[pool]
    else:
        pools = self.pools[self.sam_pool]
    if group in pools:
        sams = pools[group]
    else:
        sams = sorted(set(sam for pool in pools.values() for sam in pool))
    return sams


def get_group_reads(
        self,
        tech: str,
        group: str,
        reads: dict,
        pool: str = None
) -> dict:
    """Subset to the path(s) of input fastq file(s) per sample and
    tech/sample for the samples in the current co-assembly group(s).

    Parameters
    ----------
    self
    tech : str
        Technology/ies
    group : str
        Name of a co-assembly pool's group
    reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    pool : str
        Name of the co-assembly

    Returns
    -------
    group_reads : dict
        Path(s) to the input fastq file(s) per sample and tech/sample
    """
    sams = get_sams(self, group, pool)
    per_tech = self.soft.params.get('per_tech', False)
    group_reads = {}
    for sam in sams:
        for tech_sam, files in reads[sam].items():
            if not files:
                continue
            if per_tech:
                if tech == tech_sam[0]:
                    group_reads.setdefault(sam, {}).update({tech_sam: files})
            else:
                group_reads.setdefault(sam, {}).update({tech_sam: files})
    return group_reads


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


def get_args(fastas, contigs) -> tuple:
    selection = fastas + '/selection.tsv.gz'
    fasta = fastas + '/selection.fasta'
    cmds = 'python3 %s/scripts/subset_fasta.py -i %s -s %s -o %s\n' % (
        RESOURCES, contigs, selection, fasta)
    cmd_rm = 'rm %s\n' % fasta
    return selection, fasta, cmds, cmd_rm


def get_plasmids(self, fastas, contigs) -> tuple:
    if self.soft.prev == 'plasmidfinder':
        plasmids = fastas + '/results_tab.tsv.gz'
        cmd = 'tail -n +2 | cut -f 5'
    elif self.soft.prev == 'plasforest':
        plasmids = fastas + '/plasmids.csv.gz'
        cmd = 'grep Plasmid | cut -d"," -f 1'
    elif self.soft.prev == 'viralverify':
        plasmids = fastas + '/contigs_result_table.csv.gz'
        cmd = 'cut -d "," -f 1,2 | grep Plasmid | cut -d"," -f 1'
    else:
        sys.exit('[%s] Cannot run after %s' % (self.soft.name, self.soft.prev))

    subset = plasmids.replace('.gz', '.ref')
    fasta = plasmids.replace('.gz', '.fa')
    cmds = 'zcat %s | %s | cut -d" " -f 1 > %s\n' % (plasmids, cmd, subset)
    cmds += 'python3 %s/scripts/subset_fasta.py -i %s -s %s -o %s\n' % (
        RESOURCES, contigs, subset, fasta)
    cmd_rm = 'rm %s %s\n' % (subset, fasta)
    return plasmids, fasta, cmds, cmd_rm


def get_arg_inputs(self, inputs):
    reports = {}
    cmd, cmd_rm = '', ''
    module = self.soft.prev
    if module.startswith('deeparg'):
        module = 'deeparg'
        inp = "%s/nucleotide.sequences.mapping.ARG.gz" % inputs
        rep = splitext(inp)[0]
        reports[rep] = ('%s.out' % module, 'nucleotide.sequences.mapping.ARG')
        cmd += "gunzip -c %s > %s\n" % (inp, rep)
        cmd_rm += "rm %s\n" % rep
    elif module == 'abricate':
        for db, inp in inputs.items():
            rep = splitext(inp)[0]
            cmd += 'gunzip -c %s > %s\n' % (inp, rep)
            cmd_rm += 'rm %s\n' % rep
            reports[rep] = ('%s-%s.out' % (module, db), '')
    elif module == 'abritamr':
        module = 'amrfinderplus'
        inp = "%s/amrfinder.out.gz" % inputs
        rep = splitext(inp)[0]
        reports[rep] = ('%s.out' % module, inputs.split('/')[-1])
        cmd += "gunzip -c %s > %s\n" % (inp, rep)
        cmd_rm += "rm %s\n" % rep
    elif module == 'staramr':
        inp = "%s/resfinder.tsv.gz" % inputs
        rep = splitext(inp)[0]
        reports[rep] = ('%s.out' % module, '')
        cmd += "gunzip -c %s > %s\n" % (inp, rep)
        cmd_rm += "rm %s\n" % rep
    elif module == 'amrplusplus':
        # out = "gene.tsv"
        pass
    elif module == 'ariba':
        pass
    elif module == 'csstar':
        # out = "output.tsv"
        pass
    elif module == 'fargene':
        pass
    elif module == 'groot':
        # out = "output.tsv"
        pass
    elif module == 'kmerresistance':
        pass
    elif module == 'resfams':
        # out = "resfams.tblout"
        pass
    elif module == 'resfinder':
        pass
    elif module == 'mykrobe':
        # out = "out.json"
        pass
    elif module == 'kmerresistance':
        pass
    elif module == 'pointfinder':
        pass
    elif module == 'rgi':
        # out = "out.txt"
        pass
    elif module == 'srax':
        pass
    elif module == 'srst2':
        pass
        # fp = "resfinder.tsv"
    elif module == 'tbprofiler':
        pass
    return cmd, cmd_rm, module, reports