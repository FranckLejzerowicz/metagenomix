# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys
from os.path import dirname
from metagenomix._io_utils import io_update, to_do


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
    genomes : dict
        Genomes bins to dereplicate.
    """
    genomes = {}
    for pool in self.pools:
        genomes[pool] = {}
        for group, group_paths in self.inputs[self.pool].items():
            if self.soft.prev == 'metawrap_refine':
                genomes[pool].setdefault('', []).append((group, group_paths[1]))
            elif self.soft.prev == 'metawrap_reassemble':
                for path in group_paths:
                    suf = path.split('_')[-1]
                    genomes[pool].setdefault(suf, []).append((group, path))
            else:
                sys.exit('[drep] Not possible after "%s"' % self.soft.prev)
    return genomes


def get_drep_inputs(self, drep_dir: str, sam_paths: list):
    """Write the file containing the inputs bins to drep
    and the list of paths to these bins.

    Parameters
    ----------
    drep_dir :
        Path to the drep output folder
    sam_paths : list
        List containing for each (group, bins_folder)

    Returns
    -------
    cmd : str
        Command to create the input
    drep_in : str
        File containing the paths corresponding to each bin
    paths : list
        List of paths corresponding to each bin
    """
    cmd = ''
    paths = []
    drep_in = '%s/input_genomes.txt' % drep_dir
    for group, bins in sam_paths:
        bin_paths = glob.glob('%s/*fa' % bins.replace('${SCRATCH_FOLDER}', ''))
        if self.config.dev:
            bin_paths = ['%s/a.fa' % bins.replace('${SCRATCH_FOLDER}', ''),
                         '%s/b.fa' % bins.replace('${SCRATCH_FOLDER}', '')]
        for path in bin_paths:
            paths.append(path)
            if cmd:
                cmd += 'echo "%s" >> %s\n' % (path, drep_in)
            else:
                cmd += 'echo "%s" > %s\n' % (path, drep_in)
    if cmd:
        cmd += 'envsubst < %s > %s.tmp\n' % (drep_in, drep_in)
        cmd += 'mv %s.tmp %s\n' % (drep_in, drep_in)
    return cmd, drep_in, paths


def drep_cmd(self, algorithm: str, drep_in: str, drep_out: str,
             paths: list, cmd: str) -> str:
    cmd = '\ndRep dereplicate'
    cmd += ' %s' % drep_out
    cmd += ' --S_algorithm %s' % algorithm
    cmd += ' --ignoreGenomeQuality'
    cmd += ' --processors %s' % self.soft.params['cpus']
    chunk_size = self.soft.params['primary_chunksize']
    if len(paths) > chunk_size:
        cmd += ' --multiround_primary_clustering'
        cmd += ' --run_tertiary_clustering'
        if algorithm == 'fastANI':
            cmd += ' --greedy_secondary_clustering'
    else:
        for boolean in ['multiround_primary_clustering',
                        'greedy_secondary_clustering',
                        'run_tertiary_clustering']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean
    for boolean in ['SkipMash', 'SkipSecondary']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean
    cmd += ' --coverage_method %s' % self.soft.params['coverage_method']
    for param in [
        'S_ani', 'P_ani', 'warn_sim', 'warn_aln', 'n_PRESET',
        'warn_dist', 'cov_thresh', 'clusterAlg', 'MASH_sketch'
    ]:
        cmd += ' --%s %s' % (param, self.soft.params[param])
    cmd += ' --genomes %s' % drep_in
    return cmd

def get_key(stringency, algorithm):
    if stringency:
        key = (stringency, algorithm)
    else:
        key = (algorithm,)
    return key


def get_drep_out(self, pool, stringency, algorithm):
    drep_out = '%s/%s' % (self.dir, pool)
    if stringency:
        drep_out += '/%s' % stringency
    drep_out += '/%s' % algorithm
    self.outputs['dirs'].append(drep_out)
    return drep_out


def drep(self):
    prev = self.soft.prev
    genomes = get_drep_bins(self)
    for pool, pool_paths in genomes.items():
        for stringency, sam_paths in pool_paths.items():
            for algorithm in self.soft.params['S_algorithm']:
                key = get_key(stringency, algorithm)
                drep_out = get_drep_out(self, pool, stringency, algorithm)
                cmd, drep_in, paths = get_drep_inputs(self, drep_out, sam_paths)
                dereps = '%s/dereplicated_genomes' % drep_out
                if not paths:
                    self.soft.status.add('Must run %s (%s)' % (prev, self.pool))
                    if self.config.dev:
                        paths = ['a.fa', 'b.fa', 'c.fa']
                out_dereps = '%s/*.fa' % dereps.replace('${SCRATCH_FOLDER}', '')
                if not self.config.force and glob.glob(out_dereps):
                    self.soft.status.add('Done')
                    continue
                io_update(self, i_f=([drep_in] + paths), o_d=drep_out, key=key)
                cmd += drep_cmd(self, algorithm, drep_in, drep_out, paths, cmd)
                self.outputs['cmds'].setdefault(key, []).append(cmd)
                if pool not in self.outputs['outs']:
                    self.outputs['outs'][pool] = {}
                self.outputs['outs'][pool][key] = dereps


def get_genomes_out(self, group, genome_dir):
    out = self.dir
    if self.soft.prev == 'drep':
        out += '/%s' % '/'.join(group)
    else:
        if self.soft.prev == 'yamb':
            out += '/' + genome_dir.split('/')[-1]
        out += '/%s/%s' % (self.pool, group)
    return out


def checkm(self):
    if self.soft.params['coverage']:
        if 'mapping_spades' not in self.softs:
            sys.exit('Run "spades mapping_spades" before checkm coverage')
        bams = self.softs['mapping_spades'].outputs
    groups = get_groups(self)
    self.outputs['outs'] = {}
    for group in groups:
        genome_dirs, ext = get_genome_dirs(self, group)
        print(group)
        print(genome_dirs)
        if len(genome_dirs) > 1:
            print(genome_dirsfdsa)
        # print()
        # print("self.pool:", self.pool)
        # print("group:", group)
        # print("genome_dirs:", genome_dirs)
        # print("ext:", ext)
        cmd = checkm_tetra(self, self.dir, group)
        for genome_dir in genome_dirs:
            out = get_genomes_out(self, group, genome_dir)
            self.outputs['outs'].setdefault(group, []).append(out)
            out_dirs = '%s/*' % genome_dir.replace('${SCRATCH_FOLDER}', '')
            if not self.config.force and not glob.glob(out_dirs):
                continue
            self.outputs['dirs'].append(out)
            io_update(self, i_d=genome_dir, key=group)

            cov = checkm_coverage(self, out, group, ext, genome_dir, cmd)
            tree, cmd_tree = checkm_tree(self, out, group, ext, genome_dir, cmd)
            tree_qa = checkm_tree_qa(self, out, group, tree, cmd_tree, cmd)
            lineage = checkm_lineage_set(self, out, group, tree, tree_qa, cmd)
            analyze, cmd_analyze = checkm_analyze(self, out, group, ext,
                                                  lineage, genome_dir, cmd)
            checkm_qa(self, out, group, cov, lineage, analyze, cmd_analyze, cmd)
            checkm_unbinned(self, out, group, ext, genome_dir, cmd)
            if cmd:
                cmd_set = 'checkm data setRoot %s\n' % self.soft.params['data']
                cmd = cmd_set + cmd
            self.outputs['cmds'].setdefault(group, []).append(cmd)


def checkm_coverage(self, out_dir, key, ext, genome_dir, cmd):
    cov = ''
    if self.soft.params['coverage'] and self.soft.prev != 'drep':
        cov = '%s/coverage.txt' % out_dir
        if self.config.force or to_do(cov):
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
    trees = glob.glob('%s/*' % tree.replace('${SCRATCH_FOLDER}', ''))
    if self.config.force or not len(trees):
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
        if self.config.force or to_do(tree_fpo):
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
    if self.config.force or to_do(lineage):
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
    if not len(glob.glob('%s/*' % analyze.replace('${SCRATCH_FOLDER}', ''))):
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
        if self.config.force or to_do(qa_fp):
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
        if self.config.force or to_do(unbinned_fa):
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
        tetra = '%s/tetra' % out_dir
        tetra_fpo = '%s/tetra.txt' % tetra
        self.outputs['dirs'].append(tetra)
        io_update(self, o_d=tetra, key=key)
        if self.config.force or to_do(tetra_fpo):
            cmd += '\nif [ ! -f "%s" ]; then\n' % tetra_fpo
            cmd += 'checkm tetra'
            cmd += ' --threads %s' % self.soft.params['cpus']
            cmd += ' %s' % contigs
            cmd += ' %s\n' % tetra_fpo
            cmd += 'fi\n ' % tetra_fpo
    return cmd


def coconet(self):
    if self.soft.prev != 'spades':
        sys.exit('[coconet] can only be run on assembly output')
    bams = {}
    if 'mapping_spades' in self.softs:
        bams = self.softs['mapping_spades'].outputs
    for group, spades_outs in self.inputs[self.pool].items():
        out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
        cmd = '\ncoconet run'
        cmd += ' --fasta %s' % spades_outs[1]
        cmd += ' --output %s' % out_dir
        cmd += ' --threads %s' % self.soft.params['cpus']
        for boolean in ['quiet', 'no_rc', 'silent', 'continue',
                        'patience', 'recruit_small_contigs']:
            if self.soft.params[boolean]:
                cmd += ' --%s' % boolean
        for p in [
            'debug', 'flag', 'min_ctg_len', 'min_prevalence', 'min_dtr_size',
            'min_mapping_quality', 'min_aln_coverage', 'fragment_step',
            'n_train', 'n_test', 'batch_size', 'test_batch', 'load_batch',
            'cover_filters', 'cover_kernel', 'cover_stride', 'merge_neurons',
            'kmer', 'wsize', 'wstep', 'n_frags', 'max_neighbors', 'n_clusters',
            'vote_threshold', 'gamma1', 'gamma2', 'fragment_length', 'theta',
            'test_ratio', 'learning_rate'
        ]:
            cmd += ' --%s %s' % (p.replace('_', '-'), self.soft.params[p])
        for p in ['tlen_range', 'compo_neurons', 'cover_neurons']:
            cmd += ' --%s %s' % (p.replace('_', '-'),
                                 ' '.join(map(str, self.soft.params[p])))
        if bams:
            for bam in bams:
                print(bam)
                print(kjsrbf)
            cmd += ' --bam %s\n' % ' '.join(bams)
        log_fp = '%s/coconet.log' % out_dir
        if self.config.force or to_do(log_fp):
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out_dir, key=group)
        self.outputs['outs'][group] = out_dir


def tiara(self):
    if self.soft.prev != 'spades':
        sys.exit('[tiara] can only be run on assembly output')
    for group, spades_outs in self.inputs[self.pool].items():
        out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir
        out_fp = '%s/classifications.txt' % out_dir
        if self.config.force or to_do(out_fp):
            cmd = 'cd %s\n' % out_dir
            cmd += 'tiara'
            cmd += ' --input %s' % spades_outs[1]
            cmd += ' --output %s' % out_fp
            cmd += ' --threads %s' % self.soft.params['cpus']
            for boolean in ['probabilities', 'verbose', 'gzip']:
                if self.soft.params[boolean]:
                    cmd += ' --%s' % boolean
            for param in ['min_len', 'first_stage_kmer', 'second_stage_kmer']:
                cmd += ' --%s %s' % (param, self.soft.params[param])
            for param in ['prob_cutoff', 'to_fasta']:
                cmd += ' --%s %s' % (param, ' '.join(self.soft.params[param]))
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out_dir, key=group)



def plasforest_cmd(self, out_fp: str, in_fp: str):
    cmd = 'python3 %s' % self.soft.params['path']
    cmd += ' -i %s' % in_fp
    cmd += ' -o %s' % out_fp
    if self.soft.params['size_of_batch']:
        cmd += ' --size_of_batch %s' % self.soft.params['size_of_batch']
    cmd += ' --threads %s' % self.soft.params['cpus']
    for boolean in ['b', 'f', 'r']:
        if self.soft.params[boolean]:
            cmd += ' -%s' % boolean
    return cmd


def plasforest(self):
    if self.soft.prev != 'spades':
        sys.exit('[plasforest] can only be run on assembly output')
    for group, spades_outs in self.inputs[self.pool].items():
        out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
        self.outputs['dirs'].append(out_dir)
        self.outputs['outs'][group] = out_dir
        out_fp = '%s/plasmids.csv' % out_dir
        if self.config.force or to_do(out_fp):
            cmd = plasforest_cmd(self, out_fp, spades_outs[1])
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=spades_outs[1], o_d=out_dir, key=group)


def get_genome_dirs(self, group):
    ext = 'fa'
    if not group:
        dirs = [self.inputs[self.sam]]
    elif self.soft.prev == 'spades':
        dirs = [[self.inputs[self.pool][group][1]]]
    elif self.soft.prev == 'drep':
        dirs = [self.inputs[self.pool][group]]
    elif self.soft.prev == 'metawrap_refine':
        dirs = [self.inputs[self.pool][group][-1]]
    elif self.soft.prev == 'metawrap_reassemble':
        dirs = self.inputs[self.pool][group]
    elif self.soft.prev == 'yamb':
        ext = 'fna'
        dirs = glob.glob('%s/bins-yamb-pp-*' % self.inputs[self.pool][group])
    else:
        sys.exit('[%s] Not possible after "%s"' % (self.soft.name,
                                                   self.soft.prev))
    return dirs, ext


def get_groups(self):
    if self.soft.prev == 'drep':
        groups = self.inputs[self.pool].keys()
    elif self.pool in self.pools:
        groups = self.pools[self.pool].keys()
    else:
        groups = ['']
    return groups


def plasmidfinder_cmd(self, group: str, out_dir: str, inp):
    tmp_dir = '$TMPDIR/plasmidfinder_%s_%s' % (self.pool, group)
    cmd = 'plasmidfinder.py'
    if isinstance(inp, list):
        cmd += ' --infile %s' % ' '.join(inp)
    else:
        cmd += ' --infile %s' % inp
    cmd += ' --outputPath %s' % out_dir
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


def get_plasmidfinder_io(self, group, path, ext):
    if not group or self.soft.prev == 'spades':
        fps = path
        out_dir = '%s/%s' % (self.dir, self.sam)
    else:
        if self.config.dev:
            fps = ['%s/a.%s' % (path, ext), '%s/b.%s' % (path, ext)]
        else:
            fps = glob.glob('%s/*.%s' % (path, ext))
        if isinstance(group, str):
            out_dir = '%s/%s' % (self.dir, group)
            if 'reassembled_bins' in path:
                out_dir += '/' + path.split('/')[-1]
        else:
            out_dir = '%s/%s' % (self.dir, '/'.join(group))
    return fps, out_dir


def plasmidfinder(self):
    groups = get_groups(self)
    for group in groups:
        paths, ext = get_genome_dirs(self, group)
        for path in paths:
            fps, out_dir = get_plasmidfinder_io(self, group, path, ext)
            if group:
                io_update(self, i_f=fps, o_d=out_dir, key=group)
                self.outputs['outs'].setdefault(group, []).append(out_dir)
            else:
                io_update(self, i_f=fps, o_d=out_dir)
                self.outputs['outs'].append(out_dir)
            self.outputs['dirs'].append(out_dir)
            if not group or self.soft.prev == 'spades':
                if self.config.force or not glob.glob('%s/*' % out_dir):
                    cmd = plasmidfinder_cmd(self, group, out_dir, fps)
                    if group:
                        self.outputs['cmds'].setdefault(group, []).append(cmd)
                    else:
                        self.outputs['cmds'].append(cmd)
            else:
                for fp in fps:
                    if self.config.force or not glob.glob('%s/*' % out_dir):
                        cmd = plasmidfinder_cmd(self, group, out_dir, fp)
                        self.outputs['cmds'].setdefault(group, []).append(cmd)
