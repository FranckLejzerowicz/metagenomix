# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys
from os.path import splitext
from metagenomix._io_utils import (io_update, to_do)
from metagenomix.parameters import tech_params


def viralverify_cmd(
        self,
        contigs: str,
        out: str
) -> str:
    """Collect ViralVerify command.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Path to the reference databases
        .outputs : dict
            All outputs
    contigs : str
        Path to the input contigs fasta file
    out : str
        Path to the output folder

    Returns
    -------
    cmd : str
        ViralVerify command
    """
    cmd = 'export PATH=$PATH:%s' % self.soft.params['path']
    cmd += '\nviralverify'
    cmd += ' -f %s' % contigs
    cmd += ' -o %s' % out
    if self.soft.params['db']:
        cmd += ' --db'
    if self.soft.params['p']:
        cmd += ' --p'
    cmd += ' -t %s' % self.soft.params['cpus']
    cmd += ' -thr %s' % self.soft.params['thr']
    cmd += ' --hmm %s' % '%s/Pfam-A.hmm' % self.databases.paths.get('pfam')
    return cmd


def viralverify(self) -> None:
    """Classify the contigs as viral or non-viral.

    Parameters
    ----------
    self : Commands class instance
        .databases : dict
            Path to the reference databases
        .dir : str
            Path to pipeline output folder for filtering
        .pool : str
            Co-assembly pool name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params : dict
            Parameters
        .config
            Configurations
    """
    assemblers = [
        'plass', 'spades', 'spades_metaviral', 'spades_plasmid', 'spades_bio',
        'flye', 'canu', 'necat', 'megahit', 'unicycler']
    if self.soft.prev not in assemblers:
        sys.exit('[viralVerify] can only be run on assembly output (contigs)')

    pfam = '%s/Pfam-A.hmm' % self.databases.paths.get('pfam')
    if not self.config.dev and to_do(pfam):
        sys.exit('[viralVerify] Needs the Pfam .hmm database in database yaml')

    for (group, techs), inputs in self.inputs[self.pool].items():

        contigs = inputs[1]

        out = '%s/%s/%s/%s' % (self.dir, techs, self.pool, group)
        self.outputs['dirs'].append(out)

        out_fp = '%s_result_table.csv' % splitext(contigs)[0]
        self.outputs['outs'][group] = out_fp

        cmd = viralverify_cmd(self, contigs, out)
        if self.config.force or to_do(out_fp):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=contigs, o_d=out, key=group)
