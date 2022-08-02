# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from metagenomix._io_utils import io_update, to_do


def plasforest_cmd(
        self,
        out_fp: str,
        in_fp: str
) -> str:
    """Collect plasforest command.

    Parameters
    ----------
    self : Commands class instance
        .soft.params
            Parameters
    out_fp : str
        Path to the output file
    in_fp : str
        Path to the input file

    Returns
    -------
    cmd : str
        plasforest command
    """
    cmd = 'python3 %s' % self.soft.params['binary']
    cmd += ' -i %s' % in_fp
    cmd += ' -o %s' % out_fp
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    return cmd


def plasforest(self) -> None:
    """Classify assembly contigs as plasmids or not.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for plasforest
        .prev : str
            Previous software in
        .pool : str
            Pool name.
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .config
            Configurations
    """
    assemblers = [
        'plass', 'spades', 'spades_metaviral', 'spades_plasmid', 'spades_bio',
        'flye', 'canu', 'necat', 'megahit', 'unicycler']
    if self.soft.prev not in assemblers:
        sys.exit('[plasforest] can only be run after assembly (on contigs)')

    for (tech, group), spades_outs in self.inputs[self.pool].items():
        out_dir = '%s/%s/%s/%s' % (self.dir, tech, self.pool, group)
        self.outputs['dirs'].append(out_dir)

        self.outputs['outs'][(tech, group)] = out_dir
        out_fp = '%s/plasmids.csv' % out_dir

        if self.config.force or to_do(out_fp):
            tech_group = '_'.join([tech, group])
            cmd = plasforest_cmd(self, out_fp, spades_outs[1])
            self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out_dir, key=tech_group)



def plasmidfinder_cmd(
        self,
        key: str,
        out: str,
        inputs
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
    out : str
        Path to the output folder
    inputs : str or list
        Path(s) to the input file(s)

    Returns
    -------
    cmd : str
        PlasmidFinder command
    """
    tmp_dir = '$TMPDIR/plasmidfinder_%s_%s' % (self.pool, key)
    cmd = 'plasmidfinder.py'
    if isinstance(inputs, list):
        cmd += ' --infile %s' % ' '.join(inputs)
    else:
        cmd += ' --infile %s' % inputs
    cmd += ' --outputPath %s' % out
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
        tech_group: str,
        path: str,
        ext: str,
        assemblers: list
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
    tech_group : str or tuple
        Technology and/or co-assembly pool group name
    path : str or list
        Path(s) to the input file/folder(s)
    ext : str
        Extension of the input fasta file
    assemblers : list
        All the assembler software in usr in the pipeline

    Returns
    -------
    fps : list
        Path(s) to the input file(s)
    out : str
        Path to the output folder
    """
    if tech_group in self.config.techs or self.soft.prev in assemblers:
        fps = path
        out = '%s/%s/%s' % (self.dir, tech_group, self.sam)
    else:
        if self.config.dev:
            fps = ['%s/a.%s' % (path, ext), '%s/b.%s' % (path, ext)]
        else:
            fps = glob.glob('%s/*.%s' % (path, ext))
        if isinstance(tech_group, str):
            out = '%s/%s' % (self.dir, tech_group)
            if 'reassembled_bins' in path:
                out += '/' + path.split('/')[-1]
        else:
            out = '%s/%s' % (self.dir, '/'.join(tech_group))
    return fps, out


def genome_dirs(
        self,
        inputs: dict,
        sam_group: str,
        assemblers: list
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
    assemblers : list
        All the assembler software in use in the pipeline

    Returns
    -------
    dirs : list
    ext : str
    """
    ext = 'fa'
    if sam_group == self.sam:
        dirs = [inputs]
        print(1)
    elif self.soft.prev in assemblers:
        dirs = [[inputs[sam_group][1]]]
        print(2)
    elif self.soft.prev == 'drep':
        dirs = [inputs[sam_group]]
        print(3)
    elif self.soft.prev == 'metawrap_refine':
        dirs = [inputs[sam_group][-1]]
        print(4)
    elif self.soft.prev == 'metawrap_reassemble':
        dirs = inputs[sam_group]
        print(5)
    elif self.soft.prev == 'yamb':
        ext = 'fna'
        dirs = glob.glob('%s/bins-yamb-pp-*' % inputs[sam_group])
        print(6)
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
    print()
    print('>cutadapt', list(self.softs['cutadapt'].outputs)[-1])
    print(self.softs['cutadapt'].outputs[list(self.softs['cutadapt'].outputs)[-1]])
    print()
    print('>spades', list(self.softs['spades'].outputs)[-1])
    print(self.softs['spades'].outputs[list(self.softs['spades'].outputs)[-1]])
    print()
    print('>metawrap_refine', list(self.softs['metawrap_refine'].outputs)[-1])
    print(self.softs['metawrap_refine'].outputs[list(self.softs['metawrap_refine'].outputs)[-1]])
    print()
    print('>metawrap_reassemble', list(self.softs['metawrap_reassemble'].outputs)[-1])
    print(self.softs['metawrap_reassemble'].outputs[list(self.softs['metawrap_reassemble'].outputs)[-1]])
    print()
    print('>drep', list(self.softs['drep'].outputs)[-1])
    print(self.softs['drep'].outputs[list(self.softs['drep'].outputs)[-1]])

    assemblers = [
        'plass', 'spades', 'spades_metaviral', 'spades_plasmid', 'spades_bio',
        'flye', 'canu', 'necat', 'megahit', 'unicycler']
    print()
    print()
    print()
    for (tech, group), inputs in self.inputs[self.pool].items():
        if tech == 'illumina':
            continue
        tech_group = '_'.join([tech, group])
        print()
        print('tech, group:', (tech, group), tech_group)
        paths, ext = genome_dirs(self, inputs, group, assemblers)
        for path in paths:
            print('*** path\t:', path)
            fps, out = plasmidfinder_io(self, tech_group, path, ext, assemblers)
            self.outputs['outs'].setdefault(tech_group, []).append(out)
            self.outputs['dirs'].append(out)
            print("*** fps\t\t:", fps)
            print("*** out\t\t:", out)
            cmd = ''
            if group == self.sam or self.soft.prev in assemblers:
                if self.config.force or not glob.glob('%s/*' % out):
                    cmd += plasmidfinder_cmd(self, tech_group, out, fps)
                    self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
            else:
                for fp in fps:
                    if self.config.force or not glob.glob('%s/*' % out):
                        cmd += plasmidfinder_cmd(self, tech_group, out, fp)
            if cmd:
                self.outputs['cmds'].setdefault(tech_group, []).append(cmd)
                io_update(self, i_f=fps, o_d=out, key=tech_group)
