# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import basename, isfile, splitext
from metagenomix._io_utils import io_update


def get_quast_contigs_labels(self, pool, group_samples):
    contigs = []
    labels = []
    for group, fps in sorted(self.inputs[pool].items()):
        contigs.append(fps[1])
        if self.soft.params['label']:
            samples = group_samples[group]
            vals = self.config.meta.loc[samples, self.soft.params['label']]
            labels.append('__'.join(sorted(set(vals))))
    return contigs, labels


def quast(self):
    for pool, group_samples in self.pools.items():
        contigs, labels = get_quast_contigs_labels(self, pool, group_samples)

        out_dir = '%s/%s' % (self.dir, pool)
        if self.soft.params['label']:
            out_dir += '_%s' % self.soft.params['label']

        cmd = 'metaquast.py'
        cmd += ' %s' % ' '.join(contigs)
        if labels:
            cmd += ' --labels %s' % ','.join(labels)
        cmd += ' --output-dir %s' % out_dir
        cmd += ' --min-contig %s' % self.soft.params['min_contig']
        cmd += ' --ambiguity-usage %s' % self.soft.params['ambiguity_usage']
        for boolean in ['circos', 'glimmer', 'rna_finding',
                        'conserved_genes_finding', 'space_efficient']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean.replace('_', '-')
        if self.soft.params['fast']:
            cmd += ' --fast'
        else:
            for boolean in ['no_check', 'no_plots', 'no_html', 'no_icarus',
                            'no_snps', 'no_gc', 'no_sv', 'no_read_stats']:
                if self.soft.params[boolean]:
                    cmd += ' --%s' % boolean.replace('_', '-')

        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][pool] = out_dir
        self.outputs['cmds'].setdefault(pool, []).append(cmd)
        io_update(self, i_f=contigs, o_d=out_dir, key=pool)


def spades_cmd(fastas: list, params: dict, tmp: str, out: str):
    """Create command lines for spades.

    Parameters
    ----------
    fastas : list
        Input fasta files
    params : dict
        Parameters for Spades
    tmp : str
        Temporary folder for Spades
    out : str
        Output folder for Spades
    """
    cmd = 'spades.py'
    cmd += ' -m %s' % params['mem_num']
    cmd += ' -k %s' % ','.join(map(str, params['k']))
    cmd += ' -t %s' % params['cpus']
    if params['bio']:
        cmd += ' --bio'
    if params['meta']:
        cmd += ' --meta'
    if params['plasmid']:
        cmd += ' --plasmid'
    if params['only_assembler']:
        cmd += ' --only-assembler'
    cmd += ' --tmp-dir %s -o %s' % (tmp, out)
    for fasta in fastas:
        if 'extendedFrags' in fasta:
            cmd += ' --merge %s' % fasta
        elif 'notCombined_1' in fasta:
            cmd += ' -1 %s' % fasta
        elif 'notCombined_2' in fasta:
            cmd += ' -2 %s' % fasta
    return cmd


def spades(self) -> None:
    """Create command lines for spades.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for spades
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
    for group, fastas in self.inputs[self.pool].items():
        tmp = '$TMPDIR/spades_%s_%s' % (self.pool, group)
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = 'mkdir -p %s\n' % tmp
        cmd += '%s\n' % spades_cmd(fastas, self.soft.params, tmp, out)
        cmd += 'rm -rf %s\n' % tmp
        before_rr = '%s/before_rr.fasta' % out
        contigs = '%s/contigs.fasta' % out
        first_pe = '%s/first_pe_contigs.fasta' % out
        scaffolds = '%s/scaffolds.fasta' % out
        log = '%s/spades.log' % out
        outs = [before_rr, contigs, first_pe, scaffolds, log]
        if self.config.force or not isfile(contigs):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=fastas, i_d=tmp, o_f=outs, o_d=out, key=group)
        self.outputs['outs'][group] = outs


def viralverify(self):
    if self.soft.prev != 'spades':
        sys.exit('[viralVerify] can only be run on assembly output')
    pfam = '%s/Pfam-A.hmm' % self.databases.pfams.get('dir')
    if not self.config.dev and not isfile(pfam):
        sys.exit('[viralVerify] Needs the Pfam .hmm database in database yaml')
    for group, spades_outs in self.inputs[self.pool].items():
        out = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = 'export PATH=$PATH:%s' % self.soft.params['path']
        cmd += '\nviralverify'
        cmd += ' -f %s' % spades_outs[1]
        cmd += ' -o %s' % out
        if self.soft.params['db']:
            cmd += ' --db'
        if self.soft.params['p']:
            cmd += ' --p'
        cmd += ' -t %s' % self.soft.params['cpus']
        cmd += ' -thr %s' % self.soft.params['thr']
        cmd += ' --hmm %s' % pfam
        out_fp = '%s_result_table.csv' % splitext(spades_outs[1])[0]
        if self.config.force or not isfile(out_fp):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out, key=group)
        self.outputs['outs'][group] = out_fp


def plass(self):
    """Create command lines for Plass.

    self : Commands class instance
        .dir : str
            Path to pipeline output folder for plass
        .sam : str
            Sample name
        .inputs : dict
            Input files
        .outputs : dict
            All outputs
        .soft.params
            Parameters
        .config
            Configurations
    """
    out = '%s/%s' % (self.dir, self.sam)
    self.outputs['dirs'].append(out)
    inputs = self.inputs[self.sam]
    if self.soft.prev == 'flash':
        inputs = [x.replace('.fasta', '.fastq') for x in inputs if
                  'notCombined_' in basename(x)]
    out_fp = '%s/plass_%sassembly.fasta' % (out, self.soft.params.type)
    self.outputs['outs'].append(out_fp)
    tmp_dir = '$TMPDIR/plass_%s' % self.sam
    cmd = 'plass %sassemble' % self.soft.params.type
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --num-iterations 6'
    cmd += ' --min-seq-id 0.8'
    cmd += ' -c 0.8'
    cmd += ' %s' % ' '.join(sorted(inputs))
    cmd += ' %s' % out_fp
    cmd += ' %s' % tmp_dir
    if self.config.force or not isfile(out_fp):
        cmd = 'mkdir -p %s\n%s\nrm -rf %s' % (tmp_dir, cmd, tmp_dir)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=inputs, i_d=tmp_dir, o_f=out_fp)


def get_merges(self) -> list:
    """Get the extension of the file to merge to make the pools.

    Parameters
    ----------
    self : Commands class instance
        .soft.prev
            Software run just before pooling

    Returns
    -------
    merges : list
        Extension of the outputs to merge
    """
    merges = ['']
    if self.soft.prev == 'flash':
        merges = ['extendedFrags', 'notCombined_1', 'notCombined_2']
    return merges


def pool_cmd(self, pool: str, out_dir: str, paths: list,
             fasta: str, group: str) -> None:
    """Write the pooling command and collect the output and io for FLASh.

    Parameters
    ----------
    self : Commands class instance
        .cmds
            Command lines
        .soft.io : dict
            All files to I/O is scratch is used
        .config
            Configurations
    pool : str
        Name of the pool
    out_dir : str
        Path to the output directory
    paths : list
        Path the input fasta files
    fasta : str
        Path to an output fasta file
    group : str
        Name of the sample group within the pool
    """
    if self.config.force or not isfile(fasta):
        for pdx, path in enumerate(paths):
            if pdx:
                cmd = 'cat %s >> %s' % (path, fasta)
            else:
                cmd = 'cat %s > %s' % (path, fasta)
            self.cmds.setdefault(pool, []).append(cmd)
        io_update(self, i_f=paths, o_f=fasta, o_d=out_dir, key=group)


def get_pools(self, pool: str, group: str, sams: list) -> list:
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

    Returns
    -------
    fasta_fps : list
        Paths to the output fasta files
    """
    out_dir = self.dir + '/' + pool
    self.soft.dirs.add(out_dir.replace('${SCRATCH_FOLDER}', ''))
    fasta_fps = []
    for merge in get_merges(self):
        paths = [fp for sam in sams for fp in self.inputs[sam] if merge in fp]
        # only pool if there is min 2 samples being merged
        if len(sams) > 1:
            fasta = '%s/%s.%s.fasta' % (out_dir, group, merge)
            fasta = fasta.replace(' ', '_').replace('..', '.')
            pool_cmd(self, pool, out_dir, paths, fasta, group)
        else:
            if len(paths) > 1:
                raise ValueError('Error in pooling group...')
            fasta = paths[0]
        fasta_fps.append(fasta)
    return fasta_fps


def pooling(self, pool):
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
        # get the outputs for the current group and collect pooling commands
        self.soft.outputs[pool][group] = get_pools(self, pool, group, sams)
