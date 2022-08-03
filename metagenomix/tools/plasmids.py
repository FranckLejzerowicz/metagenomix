# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from os.path import basename, dirname
from metagenomix._io_utils import io_update, to_do


def plasforest_cmd(
        self,
        out_dir: str,
        out_fp: str,
        in_fp: str
) -> str:
    """Collect plasforest command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    out_dir : str
        Path to the output folder
    out_fp : str
        Path to the output file
    in_fp : str
        Path to the input file

    Returns
    -------
    cmd : str
        plasforest command
    """
    binary = self.soft.params['binary']
    cmd = 'cd %s\n' % out_dir
    cmd += 'cp %s/plasforest.sav %s/.\n' % (dirname(binary), out_dir)
    cmd += 'cp %s/*.fasta* %s/.\n' % (dirname(binary), out_dir)
    cmd += 'python3 %s' % binary
    cmd += ' -i %s' % in_fp
    cmd += ' -o %s' % out_fp
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    cmd += '\nrm %s/plasforest.sav\n' % out_dir
    cmd += 'rm %s/*.fasta*\n' % out_dir
    return cmd


def plasforest(self) -> None:
    """Classify assembly contigs as plasmids or not.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for plasforest
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
    """
    if self.soft.prev not in self.config.tools['assembling']:
        sys.exit('[plasforest] can only be run after assembly (on contigs)')

    for (tech, group), spades_outs in self.inputs[self.pool].items():
        out_dir = '%s/%s/%s/%s' % (self.dir, tech, self.pool, group)
        self.outputs['dirs'].append(out_dir)

        self.outputs['outs'][(tech, group)] = out_dir
        out_fp = '%s/plasmids.csv' % out_dir

        if self.config.force or to_do(out_fp):
            tech_group = '_'.join([tech, group])
            cmd = plasforest_cmd(self, out_dir, out_fp, spades_outs[1])
            self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], i_d=out_dir,
                      o_f=out_fp, key=tech_group)


def plasmidfinder_cmd(
        self,
        key: str,
        fp: str,
        base: str,
        base_out: str
) -> str:
    """Collect PlasmidFinder command.

    Parameters
    ----------
    self : Commands class instance
        .pool : str
            Pool name.
        .soft.params
            PlasmidFinder parameters
    key : str
        Concatenation of the technology and/or co-assembly pool group name
    base : str
        Basename for the current sample/MAG
    base_out : str
        Path to the output folder for the current sample/MAG
    fp : str
        Path to the input file

    Returns
    -------
    cmd : str
        PlasmidFinder command
    """
    tmp_dir = '$TMPDIR/plasmidfinder_%s_%s_%s' % (self.pool, key, base)
    cmd = 'mkdir -p %s\n' % tmp_dir
    cmd += 'plasmidfinder.py'
    if len(fp) == 2:
        cmd += ' --infile %s' % ' '.join(fp)
    else:
        cmd += ' --infile %s' % fp[0]
    cmd += ' --outputPath %s' % base_out
    cmd += ' --tmp_dir %s' % tmp_dir
    cmd += ' --methodPath %s' % self.soft.params['methodPath']
    cmd += ' --databasePath %s' % self.soft.params['databasePath']
    if 'databases' in self.soft.params:
        cmd += ' --databases %s' % self.soft.params['databases']
    cmd += ' --mincov %s' % self.soft.params['mincov']
    cmd += ' --threshold %s' % self.soft.params['threshold']
    if self.soft.params['extented_output']:
        cmd += ' --extented_output'
    cmd += '\nrm -rf %s\n' % tmp_dir
    # cmd += ' --speciesinfo_json %s' % speciesinfo_json
    return cmd


def plasmidfinder_io(
        self,
        tech: str,
        sam_group: str,
        path: str,
        ext: str
) -> tuple:
    """

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for PlasmidFinder
        .prev : str
            Previous software in the pipeline
        .config
            Configurations
        .sam
            Sample name
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam_group : str
        Sample name or name of a co-assembly pool's group
    path : str or list
        Path(s) to the input file/folder(s)
    ext : str
        Extension of the input fasta file

    Returns
    -------
    fps : dict
        Path(s) to the input file(s) per basename
    out : str
        Path to the output folder
    """
    outs = [self.dir, tech]
    if sam_group == self.sam:
        fps = {sam_group: path}
    elif self.soft.prev in self.config.tools['assembling']:
        fps = {sam_group: path}
    else:
        outs.append(sam_group)
        if self.config.dev:
            fps = ['%s/a%s' % (path, ext), '%s/b%s' % (path, ext)]
        else:
            fps = glob.glob('%s/*%s' % (path, ext))
        if 'reassembled_bins' in path:
            outs.extend(path.split('/')[-2:])
        fps = {basename(x).split(ext)[0].split('.fastq')[0]: [x] for x in fps}
    out = '/'.join(outs)
    return fps, out


def genome_dirs(
        self,
        inputs: dict,
        sam_group: str
) -> tuple:
    """

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .sam
            Sample name
    inputs : dict
        Path(s) to the input file(s)/folder(s) per tech or (tech/hybrid, group)
    sam_group : str
        Sample name or pool group

    Returns
    -------
    dirs : list
    ext : str
    """
    ext = '.fa'
    if sam_group == self.sam:
        dirs = [inputs]
    elif self.soft.prev in self.config.tools['assembling']:
        dirs = [[inputs[1]]]
    elif self.soft.prev == 'drep':
        dirs = [inputs]
    elif self.soft.prev == 'metawrap_refine':
        dirs = [inputs[-1]]
    elif self.soft.prev == 'metawrap_reassemble':
        dirs = inputs
    elif self.soft.prev == 'yamb':
        ext = '.fna'
        dirs = glob.glob('%s/bins-yamb-pp-*' % inputs)
    else:
        sys.exit('[%s] Not after "%s"' % (self.soft.name, self.soft.prev))
    return dirs, ext


def plasmidfinder(self) -> None:
    """Detect plasmids in fasta files of contigs or binned/dereplicated MAGs
    using PlasmidFinder.

    Parameters
    ----------
    self : Commands class instance
        .name : str
            Name of the current software in the pipeline
        .dir : str
            Path to pipeline output folder for PlasmidFinder
        .prev : str
            Previous software in the pipeline
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
        .sam
            Sample name
    """
    for (tech, sam_group), inputs in self.inputs[self.pool].items():
        key = '_'.join([tech, sam_group])
        paths, ext = genome_dirs(self, inputs, sam_group)
        for path in paths:
            fps, out = plasmidfinder_io(self, tech, sam_group, path, ext)
            self.outputs['outs'].setdefault((tech, sam_group), []).append(out)
            cmd = ''
            for base, fp in fps.items():
                base_out = out + '/%s' % base
                json = '%s/data.json' % base_out
                self.outputs['dirs'].append(base_out)
                if self.config.force or not to_do(json):
                    cmd += plasmidfinder_cmd(self, key, fp, base, base_out)
                    io_update(self, i_f=fp, o_d=base_out, key=key)
            if cmd:
                self.outputs['cmds'].setdefault(key, []).append(cmd)
