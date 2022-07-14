# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from os.path import dirname, isdir, isfile
from metagenomix._io_utils import mkdr, io_update


def get_drep_bins(self) -> dict:
    """Get the genomes in an iterable format accommodating
    both the output of `metawrap_refine()` and of `metawrap_reassemble()`.

    Parameters
    ----------
    self : Commands class instance
        Contains all attributes needed for dereplicating on the current sample:
            prev : str
                Name of the software preceding drep.
            pools : dict
                Samples per group as dispatched according to the pooling design.
            inputs : dict
                Outputs from the software preceding drep.

    Returns
    -------
    bins : dict
        Bins to dereplicate.
    """
    genomes = {}
    for pool in self.pools:
        genomes[pool] = {}
        for group, group_paths in self.inputs[pool].items():
            if self.soft.prev == 'metawrap_refine':
                genomes[pool].setdefault('', []).append((group, group_paths[1]))
            elif self.soft.prev == 'metawrap_reassemble':
                for path in group_paths:
                    suf = path.split('_')[-1]
                    genomes[pool].setdefault(suf, []).append((group, path))
            else:
                sys.exit('[drep] Not possible after "%s"' % self.soft.prev)
    return genomes


def get_drep_inputs(drep_dir: str, sam_paths: list):
    """Write the file containing the inputs bins to drep
    and the list of paths to these bins.

    Parameters
    ----------
    drep_dir :
        Path to the drep output folder.
    sam_paths : list
        List containing for each (group, bins_folder).

    Returns
    -------
    drep_in : str
        File containing the paths corresponding to each bin.
    paths : list
        List of paths corresponding to each bin.
    """
    paths = []
    mkdr(drep_dir)
    drep_in = '%s/input_genomes.txt' % drep_dir
    with open(drep_in, 'w') as o:
        for group, bins_folder in sam_paths:
            for path in glob.glob('%s/*fa' % bins_folder):
                paths.append(path)
                o.write('%s\n' % path)
    return drep_in, paths


def drep(self):
    genomes = get_drep_bins(self)
    for group, group_paths in genomes.items():
        for stringency, sam_paths in group_paths.items():
            drep_dir = '%s/%s' % (self.dir, group)
            if stringency:
                drep_dir += '/%s' % stringency
            drep_in, paths = get_drep_inputs(drep_dir, sam_paths)
            for algorithm in self.soft.params['S_algorithm']:
                key = (group, stringency, algorithm)
                io_update(self, i_f=([drep_in] + paths), key=key)
                drep_out = '%s/%s' % (drep_dir, algorithm)
                dereps = '%s/dereplicated_genomes' % drep_out
                if not self.config.force and len(glob.glob('%s/*.fa' % dereps)):
                    continue
                if isdir(drep_out):
                    io_update(self, i_d=drep_out, key=key)
                cmd = 'dRep dereplicate'
                cmd += ' %s' % drep_out
                cmd += ' --S_algorithm %s' % algorithm
                cmd += ' --ignoreGenomeQuality'
                chunk_size = self.soft.params['primary_chunksize']
                if len(paths) > chunk_size:
                    cmd += ' --multiround_primary_clustering'
                    cmd += ' --run_tertiary_clustering'
                    if algorithm == 'fastANI':
                        cmd += ' --greedy_secondary_clustering'
                else:
                    if self.soft.params['multiround_primary_clustering']:
                        cmd += ' --multiround_primary_clustering'
                    if self.soft.params['greedy_secondary_clustering']:
                        cmd += ' --greedy_secondary_clustering'
                    if self.soft.params['run_tertiary_clustering']:
                        cmd += ' --run_tertiary_clustering'
                if self.soft.params['SkipMash']:
                    cmd += ' --SkipMash'
                if self.soft.params['SkipSecondary']:
                    cmd += ' --SkipSecondary'
                coverage_method = self.soft.params['coverage_method']
                cmd += ' --MASH_sketch %s' % self.soft.params['MASH_sketch']
                cmd += ' --n_PRESET %s' % self.soft.params['n_PRESET']
                cmd += ' --coverage_method %s' % coverage_method
                cmd += ' --cov_thresh %s' % self.soft.params['cov_thresh']
                cmd += ' --S_ani %s' % self.soft.params['S_ani']
                cmd += ' --P_ani %s' % self.soft.params['P_ani']
                cmd += ' --warn_dist %s' % self.soft.params['warn_dist']
                cmd += ' --warn_sim %s' % self.soft.params['warn_sim']
                cmd += ' --warn_aln %s' % self.soft.params['warn_aln']
                cmd += ' --processors %s' % self.soft.params['cpus']
                cmd += ' --clusterAlg %s' % self.soft.params['clusterAlg']
                cmd += ' --genomes %s' % drep_in

                log = '%s/log' % drep_out
                figure = '%s/figure' % drep_out
                tables = '%s/data_tables' % drep_out
                self.outputs['dirs'].append(drep_out)
                self.outputs['outs'][key] = [tables, log, figure, dereps]
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                io_update(self, o_d=[tables, log, figure, dereps], key=key)


def get_genome_dirs(self, group):
    ext = 'fa'
    if self.soft.prev == 'drep':
        dirs = [self.inputs[group][-1]]
    elif self.soft.prev == 'metawrap_refine':
        dirs = [self.inputs[self.pool][group][-1]]
    elif self.soft.prev == 'metawrap_reassemble':
        dirs = [self.inputs[self.pool][group][0]]
    elif self.soft.prev == 'yamb':
        ext = 'fna'
        dirs = glob.glob('%s/bins-yamb-pp-*' % self.inputs[self.pool][group])
    else:
        sys.exit('[checkm] Not possible after "%s"' % self.soft.prev)
    return dirs, ext


def get_groups(self):
    if self.soft.prev == 'drep':
        groups = self.inputs.keys()
    else:
        groups = self.pools[self.pool].keys()
    return groups


def checkm(self):
    if self.soft.params['coverage']:
        if 'mapping_spades' not in self.softs:
            sys.exit('Run "spades mapping_spades" before checkm coverage')
        bams = self.softs['mapping_spades'].outputs
    # print()
    # print("spades")
    # print(spades)
    # print()
    # print("self.inputs")
    # print(self.inputs)
    groups = get_groups(self)
    self.outputs['outs'] = {}
    for group in groups:
        genome_dirs, ext = get_genome_dirs(self, group)
        if len(genome_dirs) > 1:
            print(group)
            print(genome_dirs)
            print(genome_dirsfdsa)
        # print()
        # print("self.pool:", self.pool)
        # print("group:", group)
        # print("genome_dirs:", genome_dirs)
        # print("ext:", ext)
        cmd = checkm_tetra(self, self.dir, group)
        for genome_dir in genome_dirs:
            out = self.dir
            if self.soft.prev == 'drep':
                out += '/%s' % '/'.join(group)
            else:
                if self.soft.prev == 'yamb':
                    out += '/' + genome_dir.split('/')[-1]
                out += '/%s/%s' % (self.pool, group)
            self.outputs['outs'].setdefault(group, []).append(out)
            if not self.config.force and not glob.glob('%s/*' % genome_dir):
                continue
            self.outputs['dirs'].append(out)
            io_update(self, i_d=genome_dir, key=group)

            cov = checkm_coverage(self, out, group, ext, genome_dir, cmd)
            # tree - lineage_wf (workflow)
            tree, cmd_tree = checkm_tree(self, out, group, ext, genome_dir, cmd)
            tree_qa = checkm_tree_qa(self, out, group, tree, cmd_tree, cmd)
            lineage = checkm_lineage_set(self, out, group, tree, tree_qa, cmd)
            analyze, cmd_analyze = checkm_analyze(self, out, group, ext,
                                                  lineage, genome_dir, cmd)
            checkm_qa(self, out, group, cov, lineage, analyze, cmd_analyze, cmd)
            # other controls
            checkm_unbinned(self, out, group, ext, genome_dir, cmd)
            if cmd:
                cmd_set = 'checkm data setRoot %s\n' % self.soft.params['data']
                cmd = cmd_set + cmd
            self.outputs['cmds'].setdefault(group, []).append(cmd)


def checkm_coverage(self, out_dir, key, ext, genome_dir, cmd):
    cov = ''
    if self.soft.params['coverage'] and self.soft.prev != 'drep':
        cov = '%s/coverage.txt' % out_dir
        if self.config.force or not isfile(cov):
            io_update(self, i_f=bams[sam][pool][sam][''],
                      o_f=coverage, key=key)
            cmd += '\ncheckm coverage'
            cmd += ' --extension %s' % ext
            if self.soft.params['all_reads']:
                cmd += ' --all_reads    '
            cmd += ' --threads %s' % self.soft.params['cpus']
            cmd += ' %s' % genome_dir
            cmd += ' %s' % cov
            cmd += ' %s\n' % bams[sam][pool][sam]['']
        else:
            io_update(self, i_f=cov, key=key)
    return cov


def checkm_tree(self, out_dir, key, ext, genome_dir, cmd):
    tree = '%s/lineage/tree' % out_dir
    self.outputs['dirs'].append(tree)
    cmd_tree = '\ncheckm tree'
    cmd_tree += ' --extension %s' % ext
    cmd_tree += ' --pplacer_threads %s' % self.soft.params['cpus']
    cmd_tree += ' %s' % genome_dir
    cmd_tree += ' %s\n' % tree
    if self.config.force or not len(glob.glob('%s/*' % tree)):
        io_update(self, o_d=tree, key=key)
        cmd += cmd_tree
        return tree, ''
    else:
        io_update(self, i_d=tree, key=key)
        return tree, cmd_tree


def checkm_tree_qa(self, out, key, tree, cmd_tree, cmd):
    tree_qa = '%s/lineage/tree_qa' % out
    io_update(self, o_d=tree_qa, key=key)
    self.outputs['dirs'].append(tree_qa)
    for (ext, out_format) in [('txt', '2'), ('nwk', '4')]:
        tree_fpo = '%s/tree_qa.%s' % (tree_qa, ext)
        if self.config.force or not isfile(tree_fpo):
            if cmd_tree:
                cmd += cmd_tree
                cmd_tree = ''
            cmd += '\ncheckm tree_qa'
            cmd += ' --out_format %s' % out_format
            cmd += ' --file %s' % tree_fpo
            if self.soft.params['ali']:
                cmd += ' --ali'
            if self.soft.params['nt']:
                cmd += ' --nt'
            if self.soft.params['tab_table']:
                cmd += ' --tab_table'
            cmd += ' %s\n' % tree
    return tree_qa


def checkm_lineage_set(self, out, key, tree, tree_qa, cmd):
    lineage = '%s/lineage/lineage.ms' % out
    self.outputs['dirs'].append(dirname(tree_qa))
    if self.config.force or not isfile(lineage):
        cmd += '\ncheckm lineage_set'
        cmd += ' %s' % tree
        cmd += ' %s\n' % lineage
    else:
        io_update(self, i_f=lineage, key=key)
    return lineage


def checkm_analyze(self, out, key, ext, lineage, genome_dir, cmd):
    analyze = '%s/lineage/analyze' % out
    cmd_analyze = '\ncheckm analyze'
    cmd_analyze += ' --extension %s' % ext
    cmd_analyze += ' --threads %s' % self.soft.params['cpus']
    if self.soft.params['ali']:
        cmd_analyze += ' --ali'
    if self.soft.params['nt']:
        cmd_analyze += ' --nt'
    cmd_analyze += ' %s' % lineage
    cmd_analyze += ' %s' % genome_dir
    cmd_analyze += ' %s\n' % analyze
    if not len(glob.glob('%s/*' % analyze)):
        io_update(self, o_d=analyze, key=key)
        self.outputs['dirs'].append(analyze)
        cmd += cmd_analyze
        return analyze, 0
    else:
        io_update(self, i_d=analyze, key=key)
        return analyze, cmd_analyze


def checkm_qa(self, out, key, cov, lineage, analyze, cmd_analyze, cmd):
    dir_qa = '%s/lineage/qa' % out
    io_update(self, o_d=dir_qa, key=key)
    self.outputs['dirs'].append(dir_qa)
    for (qa, out_format) in [
        ('table.txt', 2),
        ('marker_genes.fas', 9),
        ('marker_genes.pos', 8)
    ]:
        qa_fp = '%s/qa_%s' % (dir_qa, qa)
        if self.config.force or not isfile(qa_fp):
            if cmd_analyze:
                cmd += cmd_analyze
                cmd_analyze = ''
            cmd += '\ncheckm qa'
            cmd += ' --out_format %s' % out_format
            cmd += ' --file %s' % qa_fp
            cmd += ' --threads %s' % self.soft.params['cpus']
            if self.soft.params['tab_table']:
                cmd += ' --tab_table'
            if self.soft.params['individual_markers']:
                cmd += ' --individual_markers'
            if self.soft.params['skip_adj_correction']:
                cmd += ' --skip_adj_correction'
            if self.soft.params['skip_pseudogene_correction']:
                cmd += ' --skip_pseudogene_correction'
            if self.soft.params['ignore_thresholds']:
                cmd += ' --ignore_thresholds'
            if cov and self.soft.prev != 'drep':
                cmd += ' --coverage_file %s' % cov
            cmd += ' --aai_strain %s' % self.soft.params['aai_strain']
            cmd += ' --alignment_file %s' % self.soft.params['alignment_file']
            cmd += ' --e_value %s' % self.soft.params['e_value']
            cmd += ' --length %s' % self.soft.params['length']
            cmd += ' %s' % lineage
            cmd += ' %s\n' % analyze


def checkm_unbinned(self, out_dir, key, ext, genome_dir, cmd):
    if self.soft.prev not in ['drep']:
        contigs = self.softs['spades'].outputs[self.pool][key][0]
        unbinned = '%s/unbinned' % out_dir
        unbinned_fa = '%s/unbinned.fa' % unbinned
        unbinned_stats = '%s/unbinned_stats.tsv' % unbinned
        self.outputs['dirs'].append(unbinned)
        io_update(self, i_f=contigs, o_d=unbinned, key=key)
        if self.config.force or not isfile(unbinned_fa):
            cmd += '\ncheckm unbinned'
            cmd += ' --extension %s' % ext
            cmd += ' --min_seq_len %s' % self.soft.params['min_seq_len']
            cmd += ' %s' % genome_dir
            cmd += ' %s' % contigs
            cmd += ' %s' % unbinned_fa
            cmd += ' %s\n' % unbinned_stats


def checkm_tetra(self, out_dir, key):
    cmd = ''
    if self.soft.prev not in ['drep']:
        contigs = self.softs['spades'].outputs[self.pool][key]
        print(contigs)
        print(contigsdsa)
        tetra = '%s/tetra' % out_dir
        tetra_fpo = '%s/tetra.txt' % tetra
        self.outputs['dirs'].append(tetra)
        io_update(self, o_d=tetra, key=key)
        if self.config.force or not isfile(tetra_fpo):
            cmd += '\nif [ ! -f "%s" ]; then\n' % tetra_fpo
            cmd += 'checkm tetra'
            cmd += ' --threads %s' % self.soft.params['cpus']
            cmd += ' %s' % contigs
            cmd += ' %s\n' % tetra_fpo
            cmd += 'fi\n ' % tetra_fpo
    return cmd
