# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from metagenomix._io_utils import to_do, status_update


def pool_cmd(
        self,
        tech: str,
        pool: str,
        paths: list,
        fasta: str,
        group: str
) -> None:
    """Write the pooling command and collect the output and io for FLASh.

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .config
            Configurations
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pool : str
        Name of the pool
    paths : list
        Path the input fasta files
    fasta : str
        Path to an output fasta file
    group : str
        Name of the sample group within the pool
    """
    if self.config.force or to_do(fasta):
        cmd = ''
        for pdx, path in enumerate(paths):
            if pdx:
                cmd += 'cat %s >> %s\n' % (path, fasta)
            else:
                cmd += 'cat %s > %s\n' % (path, fasta)
        if cmd:
            self.cmds.setdefault((tech, (pool, group)), []).append(cmd)


def extension_paths(
        paths_to_merge: dict,
        fastqs: list
) -> None:
    """Fill the `paths_to_merge` dict with the fastq files to merge based on
    the matching of their extension which are the extensions to be expected
    based on the number of files available for the sample/technology and the
    previous steps of the pipeline.

    Parameters
    ----------
    paths_to_merge : dict
        Paths to the input files to merge per type of output
    fastqs : list
        Path the input files for the current sample and technology
    """
    extensions = {1: ['extendedFrags.fastq.gz', 'fastq.gz'],
                  2: ['notCombined_1.fastq.gz', 'notCombined_2.fastq.gz',
                      'R1.fastq.gz', 'R2.fastq.gz'],
                  3: ['extendedFrags.fastq.gz',
                      'notCombined_1.fastq.gz',
                      'notCombined_2.fastq.gz']}
    for fastq in fastqs:
        exts = extensions[len(fastqs)]
        for ext in exts:
            if fastq.endswith(ext):
                paths_to_merge.setdefault(ext, []).append(fastq)
                break
        else:
            sys.exit('[pooling] File "%s" has no extension %s' % (fastq, exts))


def collect_paths_to_merge(
        self,
        tech: str,
        sams: list,
) -> dict:
    """Get the list of input file paths to merge in order to constitute one
    file representing the current pool, and this for the current orientations,
    that can be 1 file (usually for long read data), 2 files (for paired-end
    illumina reads that have not been merged), or 3 files (for paired-end
    illimina reads that have been merged).

    Notes
    -----
    * "orientation" means all the per-sample files that are read 1 or
      read2, or in the case of merged reads, the "extended" or "notCombined"
      outputs.
    * the returned dict could be empty if there is no input for the
      current technology, in which case the pooling simply won't happen.

    Parameters
    ----------
    self : Commands class instance
        .inputs : dict
            Input files
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sams : list
        List of samples for the current pooling group

    Returns
    -------
    paths_to_merge : dict
        Paths to the input files to merge per type of output
    """
    paths_to_merge = {}
    for sam in sams:
        if (tech, sam) in self.inputs[sam] and self.inputs[sam][(tech, sam)]:
            fastqs = self.inputs[sam][(tech, sam)]
            extension_paths(paths_to_merge, fastqs)
    return paths_to_merge


def combine_single(
        paths_to_merge: dict
) -> None:
    """Combine the 'fastq.gz' file with the 'extendedFrags.fastq.gz' file
    if it happens that the samples' fastqs to combine are from samples that
    were both single and paired-end.

    Parameters
    ----------
    paths_to_merge : dict
        Paths to the input files to merge per type of output
    """
    if 'fastq.gz' in paths_to_merge and len(paths_to_merge) > 1:
        ext = 'extendedFrags.fastq.gz'
        paths_to_merge.setdefault(ext, []).extend(paths_to_merge['fastq.gz'])
        del paths_to_merge['fastq.gz']


def get_fasta_pools(
        self,
        tech: str,
        out: str,
        sams: list,
        pool: str,
        group: str
) -> list:
    """Collect the paths to the fasta file resulting from the pooling
    for the current group, including the preparation of the command for
    making this pool.

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    out : str
        Path the output folder
    sams : list
        Samples to pool
    pool : str
        Name of the pool
    group : str
        Name of the sample group within the pool

    Returns
    -------
    fasta_fps : list
        Paths to the fasta file resulting from the pooling for the current group
    """
    paths_to_merge = collect_paths_to_merge(self, tech, sams)
    if self.soft.params['pool_single_and_merged']:
        combine_single(paths_to_merge)
    fasta_fps = []
    for extension, paths in sorted(paths_to_merge.items()):
        status_update(self, tech, paths, group=group)
        fasta = pool_fasta(self, tech, out, extension, paths, pool, group)
        fasta_fps.append(fasta)
    return fasta_fps


def fasta_fp(
        out: str,
        group: str,
        extension: str
) -> str:
    """Get the path to the fasta file that results from merging.

    Parameters
    ----------
    out : str
        Path the output folder
    group : str
        Name of the sample group within the pool
    extension : str
        Extension to the files to merge

    Returns
    -------
    fasta : str
        Fasta file name for the merging result
    """
    fasta = '%s/%s.%s' % (out, group, extension)
    fasta = fasta.replace(' ', '_').replace('..', '.')
    return fasta


def make_pool(
        self,
        tech: str,
        pool: str,
        group: str,
        paths: list,
        fasta: str
) -> None:
    """

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    pool : str
        Name of the pool
    group : str
        Name of the sample group within the pool
    paths : list
        Path to the input files to merge
    fasta : str
        Fasta file name for the merging result
    """
    # collect the inputs to pool as to move to scratch
    add_to_pool_io(self, ('I', 'f'), tech, pool, group, paths)
    # collect the pooling command line
    pool_cmd(self, tech, pool, paths, fasta, group)


def pool_fasta(
        self,
        tech: str,
        out: str,
        extension: str,
        paths: list,
        pool: str,
        group: str
) -> str:
    """

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    tech: str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    out : str
        Path the output folder
    extension : str
        Extension to the files to merge
    paths : list
        Paths to the input files to merge
    pool : str
        Name of the pool
    group : str
        Name of the sample group within the pool

    Returns
    -------
    fasta : str
        Path to the file resulting from the pooling
    """
    # only pool if there is min 2 samples being merged
    if len(paths) > 1:
        fasta = fasta_fp(out, group, extension)
        make_pool(self, tech, pool, group, paths, fasta)
    else:
        fasta = paths[0]
    return fasta


def get_pools(
        self,
        pool: str,
        group: str,
        sams: list
) -> None:
    """Get the output fasta files returned by FLASh and write the commands.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for pools
        .inputs : dict
            Input files
        .cmds
            Command lines
        .soft.dirs : set
            All directories to create
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    pool : str
        Name of the pool
    group : str
        Name of the sample group within the pool
    sams : list
        Samples to pool
    """
    # not possible to iterate over samples since pooling considers ALL samples
    for tech in self.config.techs:
        # initialize the data structure directly for the Commands instance
        out = '%s/%s/%s' % (self.dir, tech, pool)
        self.soft.dirs.add(out.replace('${SCRATCH_FOLDER}', ''))
        # get the fasta files of the pools (per orientation)
        fasta_fps = get_fasta_pools(self, tech, out, sams, pool, group)
        # fill the data structures
        add_to_pool_io(self, ('O', 'd'), tech, pool, group, [out])
        add_to_pool_io(self, ('O', 'f'), tech, pool, group, fasta_fps)
        self.soft.outputs[pool][group][tech] = fasta_fps


def add_to_pool_io(
        self,
        io: tuple,
        tech: str,
        pool: str,
        group: str,
        values: list
) -> None:
    """Add to the self.soft.io data structure the path to move to/from
    scratch for the pooling step specifically (which is pivotal between
    per sample processing to per pool/pool group).

    Parameters
    ----------
    self : Commands class instance
        .soft.io : dict
            All files to I/O is scratch is used
    io : tuple
        ('I', 'd'), ('I', 'f'), ('O', 'd'), or ('O', 'f')
    tech : str
        Name pf the technology
    pool : str
        Name of the pool
    group : str
        Name of the pool group
    values: list
        Paths to move to/from scratch for the pooling
    """
    key = (tech, (pool, group))
    if key not in self.soft.io:
        self.soft.io[key] = {}
    if io not in self.soft.io[key]:
        self.soft.io[key].setdefault(io, set()).update(values)


def pooling(
        self,
        pool: str
) -> None:
    """Create command lines for pooling samples' reads using FLASh.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    pool : str
        Name of the pool
    """
    for group, group_pd in self.config.meta.groupby(pool):
        # get the full list of samples and the samples per pooling group
        sams = group_pd.index.tolist()
        self.pools[pool][group] = sams
        self.soft.outputs[pool][group] = {}
        # get the outputs for the current group and collect pooling commands
        get_pools(self, pool, group, sams)
