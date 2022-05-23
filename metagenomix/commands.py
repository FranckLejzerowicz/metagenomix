# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import numpy as np
from os.path import basename, dirname, isdir, isfile, splitext
from metagenomix._io_utils import (
    get_edit_fastq_cmd, count_reads_cmd, get_out_dir, write_hmms)
from metagenomix.tools.simka import (
    check_simka_params, get_simka_input, simka_cmd, simka_pcoa_cmd)
from metagenomix.tools.midas import midas
from metagenomix.tools.phlans import metaphlan, humann, strainphlan
from metagenomix.tools.shogun import shogun
from metagenomix.tools.woltka import woltka
from metagenomix.tools.metawrap import get_sams_fastqs
from metagenomix.tools.drep import get_drep_bins, get_drep_inputs
# from metagenomix.tools.metamarker import metamarker


class Commands(object):

    def __init__(self, config, databases, workflow):
        self.config = config
        self.databases = databases
        self.softs = workflow.softs
        self.analysis = ''
        self.dir = ''
        self.cmds = {}
        self.args = {}
        self.pools = {}
        self.sam = None
        self.pool = None
        self.inputs = None
        self.method = None
        self.soft = None
        self.out = []
        self.holistics = [
            'multiqc',
            'simka',
            'simka_pcoa',
            'qiita_wol',
            'mag_data',
            'metamarker',
            'woltka',
            'drep',
            'strainphlan'
        ]

    def collect(self):
        for sdx, softs in enumerate(self.config.pipeline):
            self.analysis = softs[-1]
            self.soft = self.softs[self.analysis]
            self.get_inputs()
            self.get_method()
            self.get_dir()
            self.generic_command()

    def get_method(self):
        """Get the command-preparing method from this class (for the tools that
        are easy to deal with), or from auxillary modules located in the `tools`
        submodules path."""
        func = 'prep_%s' % self.analysis
        if func in dir(self) and callable(getattr(self, func)):
            self.method = getattr(self, func)
        else:
            raise ValueError('No method for software %s' % self.analysis)

    def get_inputs(self):
        """Update the `inputs` attribute of the software object."""
        if not self.soft.prev or self.soft.name == 'map__drep':
            self.inputs = self.config.fastq
        else:
            self.inputs = self.softs[self.soft.prev].outputs

    def get_dir(self):
        self.dir = '%s/%s/after_%s' % (
            self.config.dir, self.soft.name, self.soft.prev)
        self.softs[self.soft.name].dirs.add(self.dir)

    def generic_command(self):
        self.sam = None
        # print()
        # print('name:', self.soft.name)
        if self.soft.name in self.holistics:
            self.prep_job()
        elif self.soft.name == 'pooling':
            self.pooling()
        else:
            for sam in sorted(self.inputs):
                self.sam = sam
                self.pool = sam
                self.prep_job()
        # print('---- cmds:')
        # print(self.cmds)
        # print('---- outputs:')
        # print(self.soft.outputs)
        # print('---- io:')
        # print(self.soft.io)
        self.register_command()

    def prep_job(self):
        self.out = []
        self.method()
        if 'simka' in self.soft.name or self.soft.name in ['drep', 'mag_data']:
            self.soft.outputs[self.soft.name] = self.out
        else:
            self.soft.outputs[self.pool] = self.out

    def prep_edit_fastqs(self):
        cmd1 = get_edit_fastq_cmd(self.config.fastq[self.sam], 1)
        cmd2 = get_edit_fastq_cmd(self.config.fastq[self.sam], 2)
        self.cmds[self.sam] = [cmd1, cmd2]

    def pooling(self):
        for pool in self.config.pooling_groups:
            self.pools[pool] = {}
            self.softs[self.soft.name].outputs[pool] = {}
            for group, group_pd in self.config.meta.groupby(pool):
                sams = group_pd.index.tolist()
                self.pools[pool][group] = sams
                self.prep_pooling(group, sams, pool)

    def prep_pooling(self, group, sams, pool):
        out_dir = self.dir + '/' + pool
        self.soft.dirs.add(out_dir)
        merges = ['']
        if self.soft.prev == 'flash':
            merges = ['extendedFrags', 'notCombined_1', 'notCombined_2']
        for merge in merges:
            paths = [x for s in sams for x in self.inputs[s] if merge in x]
            if len(sams) > 1:
                fas = '%s/%s.%s.fasta' % (out_dir, group, merge)
                fas = fas.replace(' ', '_').replace('..', '.')
                for pdx, path in enumerate(paths):
                    if pdx:
                        cmd = 'cat %s >> %s' % (path, fas)
                    else:
                        cmd = 'cat %s > %s' % (path, fas)
                    self.cmds.setdefault(pool, []).append(cmd)
                self.soft.io['I']['f'].update(paths)
                self.soft.io['O']['f'].update([fas, out_dir])
            else:
                if len(paths) > 1:
                    raise ValueError('Error in pooling group...')
                fas = paths[0]
            if group in self.soft.outputs[pool]:
                self.soft.outputs[pool][group].append(fas)
            else:
                self.soft.outputs[pool][group] = [fas]

    def prep_count_reads_grep(self):
        out = '%s/%s_read_count.tsv' % (self.dir, self.sam)
        self.out.append(out)
        if not isfile(out):
            self.cmds[self.sam] = []
            inputs = self.inputs[self.sam]
            for idx, input_path in enumerate(inputs):
                cmd = count_reads_cmd(idx, input_path, out, self.sam)
                self.cmds[self.sam].append(cmd)
            self.soft.io['I']['f'].update(inputs)
            self.soft.io['O']['f'].add(out)

    def prep_fastqc(self):
        out_dir = '%s/%s' % (self.dir, self.sam)
        inputs = self.inputs[self.sam]
        cmd = 'fastqc %s %s -o %s' % (inputs[0], inputs[1], out_dir)
        self.soft.io['I']['f'].update(inputs)
        self.soft.io['O']['d'].add(out_dir)
        self.soft.dirs.add(out_dir)
        self.cmds[self.sam] = list([cmd])
        self.out = out_dir

    def prep_simka(self):
        inp = get_simka_input(self.dir, self.inputs)
        self.soft.io['I']['f'].add(inp)
        smin = True
        tmp_dir = '%s/simka' % self.config.scratch
        self.soft.dirs.add(tmp_dir)
        k_space, n_space = check_simka_params(self.soft.params)
        for k in map(int, k_space):
            self.cmds[k] = []
            for n in map(int, n_space):
                out_dir = '%s/k%s/n%s' % (self.dir, k, n)
                cmd = simka_cmd(self.soft, smin, inp, out_dir, k, n, tmp_dir)
                self.out.append(out_dir)
                self.soft.dirs.add(out_dir)
                self.cmds[k].append(cmd)
                self.soft.io['O']['d'].add(out_dir)

    def prep_simka_pcoa(self):
        for idx, input_path in enumerate(self.inputs['simka']):
            self.cmds[idx] = []
            for mdx, mat in enumerate(glob.glob('%s/mat_*.csv*' % input_path)):
                cmd = simka_pcoa_cmd(mat, self.config.meta_fp)
                if cmd:
                    self.cmds[idx].append(cmd)

    def prep_cutadapt(self):
        r1_o = '%s/%s.R1.fastq.gz' % (self.dir, self.sam)
        r2_o = '%s/%s.R2.fastq.gz' % (self.dir, self.sam)
        cmd = 'cutadapt'
        cmd += ' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
        cmd += ' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
        cmd += ' -m 10 -o %s -p %s ' % (r1_o, r2_o)
        cmd += ' '.join(self.inputs[self.sam])
        self.out = [r1_o, r2_o]
        self.cmds[self.sam] = list([cmd])
        self.soft.io['I']['f'].update(self.inputs)
        self.soft.io['O']['f'].update([r1_o, r2_o])

    def prep_midas(self):
        self.soft.io['I']['f'].update(self.inputs)
        for focus, db_species in self.config.midas_foci.items():
            io, cmds, outputs = midas(
                self.dir, self.sam, self.inputs, self.databases.midas['path'],
                self.soft.params['cpus'], focus, db_species)
            if outputs:
                self.out = outputs
            self.softs[self.soft.name].dirs.update(outputs)
            self.cmds[self.sam] = cmds
            self.soft.io['I']['d'].update(io['I'])
            self.soft.io['O']['d'].update(io['O'])

    def prep_metaphlan(self):
        io, cmd, outputs = metaphlan(
            self.dir, self.sam, self.inputs, self.databases.metaphlan['path'],
            self.soft.params, self.softs['count_reads_grep'].outputs,
            self.config.strains)
        for output in outputs:
            self.soft.dirs.add(dirname(output))
        self.out = outputs
        self.cmds[self.sam] = cmd
        self.soft.io['I']['f'].update(io['I'])
        self.soft.io['O']['f'].update(io['O'])

    def prep_humann(self):
        io, cmd, outputs = humann(
            self.dir, self.sam, self.inputs, self.databases.humann['path'],
            self.soft.params, self.config.humann_profile)
        self.out = outputs
        self.cmds[self.sam] = cmd
        self.soft.dirs.update([x[0] for x in list(outputs)])
        self.soft.io['I']['f'].update(io['I'])
        self.soft.io['O']['f'].update(io['O']['f'])
        self.soft.io['O']['d'].update(io['O']['d'])

    def prep_phylophlan(self):
        self.cmds = {}

    def prep_strainphlan(self):
        io, cmd, outputs, dirs = strainphlan(
            self.dir, self.inputs, self.soft.params, self.databases.wol,
            self.config.strains)
        self.cmds[''] = cmd
        self.out = outputs
        self.soft.dirs.update(dirs)
        self.soft.io['I']['f'].update(io['I']['f'])
        self.soft.io['I']['d'].update(io['I']['d'])
        self.soft.io['O']['d'].update(io['O'])

    def prep_flash(self):
        min_overlap = self.soft.params['min_overlap']
        max_mismatch_density = self.soft.params['max_mismatch_density']
        out = self.dir + '/%s_%s' % (min_overlap, max_mismatch_density)
        self.softs[self.soft.name].dirs.add(out)
        rad = out + '/' + self.sam
        ext = '%s.extendedFrags.fastq' % rad
        ext_fa = ext.replace('.fastq', '.fasta')
        nc1 = '%s.notCombined_1.fastq' % rad
        nc1_fa = nc1.replace('.fastq', '.fasta')
        nc2 = '%s.notCombined_2.fastq' % rad
        nc2_fa = nc2.replace('.fastq', '.fasta')
        if not isfile(ext_fa) or not isfile(nc1_fa) or not isfile(nc2_fa):
            cur_cmd = 'flash %s %s -m %s -x %s -d %s -o %s -t %s' % (
                self.inputs[self.sam][0], self.inputs[self.sam][1], min_overlap,
                max_mismatch_density, out, self.sam, self.soft.params['cpus'])
            self.cmds[self.sam] = [
                cur_cmd, 'seqtk seq -A %s > %s' % (ext, ext_fa),
                         'seqtk seq -A %s > %s' % (nc1, nc1_fa),
                         'seqtk seq -A %s > %s' % (nc2, nc2_fa)]
        outputs = [ext_fa, nc1_fa, nc2_fa]
        self.out = outputs
        self.soft.io['I']['f'].update([ext, nc1, nc2])
        self.soft.io['O']['f'].update(outputs)

    def prep_shogun(self):
        io, cmds, outputs, dirs = shogun(
            self.dir, self.sam, self.inputs, self.soft.params,
            self.databases.shogun, self.config)
        self.out = outputs
        self.cmds[self.sam] = cmds
        self.soft.dirs.update(dirs)
        self.soft.io['I']['f'].update(io['I'])
        self.soft.io['O']['f'].update(io['O']['f'])
        self.soft.io['O']['d'].update(io['O']['d'])

    def prep_woltka(self):
        io, cmd, outputs, dirs = woltka(
            self.dir, self.inputs, self.databases.wol['path'])
        self.cmds[''] = cmd
        self.out = outputs
        self.soft.dirs.update(dirs)
        self.soft.io['I']['f'].update(io['I']['f'])
        self.soft.io['I']['d'].update(io['I']['d'])
        self.soft.io['O']['f'].update(io['O']['f'])
        self.soft.io['O']['d'].update(io['O']['d'])

    def prep_plass(self):
        cmds = []
        out = '%s/%s' % (self.dir, self.sam)
        self.soft.dirs.add(out)
        inputs = self.inputs[self.sam]
        if self.soft.prev == 'flash':
            inputs = [x.replace('.fasta', '.fastq') for x in inputs if
                      'notCombined_' in basename(x)]
        out_fp = '%s/plass_%sassembly.fasta' % (out, self.config.plass_type)
        self.out = out_fp
        tmp_dir = '$TMPDIR/plass_%s' % self.sam
        cmd = 'plass %sassemble' % self.config.plass_type
        cmd += ' --threads %s' % self.soft.params['cpus']
        cmd += ' --num-iterations 6'
        cmd += ' --min-seq-id 0.8'
        cmd += ' -c 0.8'
        cmd += ' %s' % ' '.join(sorted(inputs))
        cmd += ' %s' % out_fp
        cmd += ' %s' % tmp_dir
        if not isfile(out_fp):
            cmds.append('mkdir -p %s' % tmp_dir)
            cmds.append(cmd)
            cmds.append('rm -rf %s' % tmp_dir)
            self.soft.io['I']['f'].update(inputs)
            self.soft.io['I']['d'].add(tmp_dir)
            self.soft.io['O']['f'].add(out_fp)
        self.cmds[self.sam] = cmds

    def prodigal_cmd(self, contigs_fp: str, out_dir: str) -> list:
        gbk_out_fp = '%s/gene.coords.gbk' % out_dir
        prot_translated_fp = '%s/protein.translations.fasta' % out_dir
        potential_genes_fp = '%s/potential.starts.fasta' % out_dir
        cmd = 'prodigal'
        cmd += ' -i %s' % contigs_fp
        cmd += ' -o %s' % gbk_out_fp
        cmd += ' -a %s' % prot_translated_fp
        cmd += ' -s %s' % potential_genes_fp
        cmd += ' -p meta'
        outputs = [gbk_out_fp, prot_translated_fp, potential_genes_fp]
        if not isfile(gbk_out_fp):
            self.cmds[self.sam] = cmd
            self.soft.io['I']['f'].add(contigs_fp)
            self.soft.io['O']['f'].update(outputs)
        return outputs

    def prep_prodigal(self):
        inputs = self.inputs[self.sam]
        if self.pool in self.pools:
            self.out = {}
            for group in self.pools[self.pool]:
                contigs_fp = inputs[group][1]
                out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
                outputs = self.prodigal_cmd(contigs_fp, out_dir)
                self.out[group] = outputs
        else:
            contigs_fp = self.inputs[self.sam]
            out_dir = '%s/%s' % (self.dir, self.sam)
            outputs = self.prodigal_cmd(contigs_fp, out_dir)
            self.out = outputs

    def prep_spades(self):
        self.out = {}
        self.cmds[self.pool] = []
        for group, fastas in self.inputs[self.pool].items():
            tmp_dir = '%s/spades_%s' % (self.config.scratch, self.pool)
            out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
            cmd = 'spades.py'
            cmd += ' -m 240 -k 33,55,77,99,127'
            cmd += ' -t %s' % self.soft.params['cpus']
            cmd += ' --meta --only-assembler'
            cmd += ' --tmp-dir %s -o %s' % (tmp_dir, out_dir)
            for fasta in fastas:
                if 'extendedFrags' in fasta:
                    cmd += ' --merge %s' % fasta
                elif 'notCombined_1' in fasta:
                    cmd += ' -1 %s' % fasta
                elif 'notCombined_2' in fasta:
                    cmd += ' -2 %s' % fasta
            before_rr = '%s/before_rr.fasta' % out_dir
            contigs = '%s/contigs.fasta' % out_dir
            first_pe = '%s/first_pe_contigs.fasta' % out_dir
            scaffolds = '%s/scaffolds.fasta' % out_dir
            log = '%s/spades.log' % out_dir
            outputs = [before_rr, contigs, first_pe, scaffolds, log]
            if not isfile(contigs):
                self.cmds[self.pool].extend(['mkdir -p %s' % tmp_dir, cmd])
                self.soft.io['I']['d'].add(tmp_dir)
                self.soft.io['I']['f'].update(fastas)
                self.soft.io['O']['d'].add(out_dir)
                self.soft.io['O']['f'].update(outputs)
            self.out[group] = outputs

    def prep_read_mapping(self):
        refs = {}
        fastqs, group_fps = self.config.fastq, self.inputs[self.pool]
        if self.soft.prev == 'drep':
            for algo, fp in input_paths['drep'][pool].items():
                ref_fasta = '%s.fa' % fp
                if not isfile(ref_fasta):
                    cmd = 'if [ ! -f %s ]; then cat %s/*.fa > %s; fi' % (
                        ref_fasta, fp, ref_fasta)
                    read_mapping_cmds.append(cmd)
                refs[algo] = ref_fasta
        elif self.soft.prev == 'metawrap_analysis_reassemble_bins':
            refs = {'': glob.glob(input_paths[sam][pool][sam][1])}
        elif self.soft.prev == 'spades':
            refs = {group: fps[1] for group, fps in group_fps.items()}

        self.out = {}
        for group, fasta in refs.items():
            self.out[group] = {}
            out_dir = '%s/%s' % (self.dir, self.pool)
            if group:
                out_dir += '/%s' % group
            self.soft.dirs.add(out_dir)
            for sam in self.pools[self.pool][group]:
                bam_out = '%s/%s.bam' % (out_dir, sam)
                self.out[group][sam] = bam_out
                cmd = 'minimap2'
                cmd += ' -a -x sr'
                cmd += ' %s' % fasta
                cmd += ' %s' % fastqs[sam][0]
                cmd += ' %s' % fastqs[sam][1]
                cmd += ' | samtools view -F 0x104 -b -'
                cmd += ' | samtools sort -o %s - ' % bam_out
                if not isfile(bam_out) or self.config.force:
                    self.cmds.setdefault(self.pool, []).append(cmd)
                bam_out_bai = '%s.bai' % bam_out
                cmd = 'samtools index %s' % bam_out
                if not isfile(bam_out_bai) or self.config.force:
                    self.cmds.setdefault(self.pool, []).append(cmd)
                self.soft.io['I']['f'].add(fasta)
                self.soft.io['I']['f'].update(fastqs[sam])
                self.soft.io['O']['f'].update([bam_out, bam_out_bai])

    def prep_map__spades_prodigal(self):
        if 'prodigal' not in self.softs or 'read_mapping' not in self.softs:
            return None
        if self.softs['prodigal'].prev != 'spades':
            return None
        if self.softs['read_mapping'].prev != 'spades':
            return None
        prodigals_fps = self.softs['prodigal'].outputs
        sams_fps = self.softs['read_mapping'].outputs
        group_fps = self.inputs[self.pool]
        self.out = {}
        for group, fps in group_fps.items():
            self.out[group] = {}
            for sam in self.pools[self.pool][group]:
                bam_fp = sams_fps[self.pool][group][sam]
                prot_fp = prodigals_fps[self.pool][group][1]
                out_dir = '%s/%s/%s' % (self.dir, self.pool, sam)
                out_fp = '%s/reads.txt' % out_dir
                if not isfile(out_fp):
                    cmd = 'pysam_reads_to_prodigal.py \\\n'
                    cmd += '-prodigal %s \\\n' % prot_fp
                    cmd += '-bam %s \\\n' % bam_fp
                    cmd += '-out %s\n' % out_fp
                    self.cmds.setdefault(self.pool, []).append(cmd)
                self.out[group][sam] = out_fp
                self.soft.io['I']['f'].update([prot_fp, bam_fp, out_fp])

    def custom_cmd(self, input_file, out_dir, dia_hmm):
        for target, gene_hmm_dia in self.databases.hmms_dias.items():
            if target == 'cazy':
                continue
            for gene, (hmm, dia) in gene_hmm_dia.items():
                sam_dir = '%s/%s' % (out_dir, gene)
                self.soft.dirs.add(sam_dir)
                out = '%s/%s.tsv' % (sam_dir, self.sam)
                if not isfile(out):
                    if dia_hmm == 'hmm':
                        stdout = '%s_hmmer.out' % splitext(out)[0]
                        self.soft.io['O']['f'].update([out, stdout])
                        self.soft.io['I']['f'].add(hmm)
                        cmd = 'hmmsearch'
                        cmd += ' --cut_tc'
                        cmd += ' --tblout %s' % out
                        cmd += ' -o %s' % stdout
                        cmd += ' --cpu %s' % self.soft.params['cpus']
                        cmd += ' %s %s' % (hmm, input_file)
                    else:
                        self.soft.io['O']['f'].add(out)
                        self.soft.io['I']['f'].add(dia)
                        tmp = '%s/%s_tmp' % (sam_dir, self.sam)
                        cmd = 'mkdir -p %s\n' % tmp
                        cmd += 'diamond blastp'
                        cmd += ' -d %s' % dia
                        cmd += ' -q %s' % input_file
                        cmd += ' -o %s' % out
                        cmd += ' -k 1 --id 80'
                        cmd += ' -p %s' % self.soft.params['cpus']
                        cmd += ' -t %s' % tmp
                    self.cmds.setdefault(target, []).append(cmd)
                    self.soft.io['O']['f'].add(out)
                else:
                    self.soft.io['I']['f'].add(out)
                out2 = out.replace('.tsv', '_contigs.tsv')
                if not isfile(out2):
                    self.soft.io['O']['f'].add(out2)
                    cmd = 'extract_custom_searched_contigs.py'
                    cmd += ' -i %s' % out
                    cmd += ' -o %s' % out2
                    self.cmds.setdefault(target, []).append(cmd)
                self.out.setdefault(target, []).extend([out, out2])

    def prep_custom(self, dia_hmm):
        self.out = {}
        if self.pool in self.pools:
            for group in self.pools[self.pool]:
                self.sam = group
                o_dir, fp = get_out_dir(self.dir, self.inputs, self.pool, group)
                self.soft.io['I']['f'].add(fp)
                self.custom_cmd(fp, o_dir, dia_hmm)
        else:
            o_dir, fp = get_out_dir(self.dir, self.inputs, self.sam)
            self.soft.io['I']['f'].add(fp)
            self.custom_cmd(fp, o_dir, dia_hmm)

    def prep_diamond_custom(self):
        self.prep_custom('dia')

    def prep_hmmer_custom(self):
        self.prep_custom('hmm')

    def macsyfinder_cmd(self, fp, o_dir, folder):
        self.soft.io['I']['f'].add(fp)
        self.soft.dirs.add(o_dir)
        models = ['CAS', 'TFF-SF', 'TXSS']
        for model in models:
            out_dir = '%s/%s' % (o_dir, model)
            self.soft.io['O']['d'].add(out_dir)
            self.soft.dirs.add(out_dir)
            cmd = 'macsyfinder'
            cmd += ' --db-type unordered'
            cmd += ' --replicon-topology linear'
            cmd += ' --e-value-search 0.1'
            cmd += ' --coverage-profile 0.5'
            cmd += ' -o %s' % out_dir
            cmd += ' --res-search-suffix _hmm.tsv'
            cmd += ' --res-extract-suffix _out.tsv'
            cmd += ' -w %s' % self.soft.params['cpus']
            cmd += ' -v'
            cmd += ' --sequence-db %s' % fp
            cmd += ' --models-dir %s/data/models' % folder
            cmd += ' -m %s all' % model
            res = '%s/macsyfinder.log' % out_dir
            self.out.append(out_dir)
            if not isfile(res):
                self.cmds.setdefault(self.sam, []).append(cmd)

    def prep_macsyfinder(self):
        folder = self.databases.macsyfinder['path']
        self.soft.io['I']['d'].add(folder)
        if self.pool in self.pools:
            for group in self.pools[self.pool]:
                o_dir, fp = get_out_dir(self.dir, self.inputs, self.pool, group)
                self.macsyfinder_cmd(fp, o_dir, folder)
        else:
            o_dir, fp = get_out_dir(self.dir, self.inputs, self.sam)
            fp_edit = '%s_edit.fasta' % fp.replace('.fasta', '')
            cmd = 'header_space_replace.py -i %s -o %s --n' % (fp, fp_edit)
            self.cmds.setdefault(self.sam, []).append(cmd)
            self.macsyfinder_cmd(fp_edit, o_dir, folder)

    def diamond_annot_cmd(self, o_dir, fp, tmp_dir):
        self.soft.io['I']['f'].add(fp)
        self.soft.dirs.add(o_dir)
        params = self.soft.params
        for k in params['k']:
            out = '%s/%s_k%s_%s.tsv' % (o_dir, self.sam, k, params['db_type'])
            cmd = 'diamond blastp'
            cmd += ' -d %s' % self.databases.diamond[params['db_type']]
            cmd += ' -q %s' % fp
            cmd += ' -o %s' % out
            cmd += ' -k 1'
            cmd += ' -p %s' % params['cpus']
            cmd += ' --id 80'
            cmd += ' -t %s' % tmp_dir
            self.cmds.setdefault(self.sam, []).append(cmd)
            self.soft.io['I']['f'].add(fp)
            self.soft.io['O']['f'].add(out)

    def prep_diamond_annot(self):
        tmp_dir = '$TMPDIR/diamond_annot_%s' % self.sam
        self.cmds[self.sam] = ['mkdir -p %s' % tmp_dir]
        if self.pool in self.pools:
            for group in self.pools[self.pool]:
                o_dir, fp = get_out_dir(self.dir, self.inputs, self.pool, group)
                self.diamond_annot_cmd(o_dir, fp, tmp_dir)
        else:
            o_dir, fp = get_out_dir(self.dir, self.inputs, self.sam)
            self.diamond_annot_cmd(o_dir, fp, tmp_dir)

    def prep_ioncom(self):
        # -------------------------------------------------------
        # --- WILL NEED TO BE REVISED AFTER I-TASSER DOWNLOAD ---
        # -------------------------------------------------------
        itasser_libs = self.databases.ioncom['itasser']
        itasser_soft = self.soft.params['itasser']
        ioncom_dir = '%s/IonCom_standalone' % self.databases.ioncom['path']
        self.soft.io['I']['d'].extend([ioncom_dir, itasser_soft, itasser_libs])
        if self.pool in self.pools:
            self.out[self.pool] = []
            for group in self.pools[self.pool]:
                o_dir, fp = get_out_dir(self.dir, self.inputs, self.pool, group)
                self.soft.io['I']['f'].add(fp)
                self.soft.io['O']['d'].add(o_dir)
                self.out[self.pool].append(o_dir)
                self.soft.dirs.add(o_dir)

                cmd = 'prepare_ioncom_inputs.py'
                cmd += ' -d %s -p %s -n 1' % (ioncom_dir, fp)
                cmd += '%s/run_IonCom.pl' % ioncom_dir
                self.cmds.setdefault(self.pool, []).append(cmd)

                output = '%s/output' % ioncom_dir
                cmd = 'for i in %s/*rehtml\n' % output
                cmd += 'do\n'
                # cmd += '    mkdir -p "%s/$(basename $(dirname "$i"))\n' % (
                #     output, out_dir)
                cmd += '    mv %s/*rehtml %s/.\n' % (output, o_dir)
                cmd += 'done\n'
                cmd += 'for i in %s/*/*/prob_*.txt\n' % output
                cmd += 'do\n'
                cmd += '    mkdir -p "$i"\n'
                cmd += 'done'

    def prep_integron_finder(self):
        self.out = {}
        if self.pool in self.pools:
            self.out[self.pool] = []
            for group in self.pools[self.pool]:
                o_dir, fp = get_out_dir(self.dir, self.inputs, self.pool, group)
                fp_out = fp.replace('.fasta', '_len2000.fasta')
                self.soft.io['O']['d'].add(o_dir)
                self.out[self.pool].append(o_dir)
                self.soft.dirs.add(o_dir)
                hmms_fp = write_hmms(self.dir, self.databases.hmms_dias)
                self.soft.io['I']['f'].update([fp, hmms_fp])
                cmd = 'filter_on_length.py -i %s -o %s -t 2000\n' % (fp, fp_out)
                cmd += 'integron_finder'
                cmd += ' --local-max'
                cmd += ' --outdir %s' % o_dir
                cmd += ' --cpu %s' % self.soft.params['cpus']
                if self.databases.hmms_dias:
                    cmd += ' --func-annot'
                    cmd += ' --path-func-annot %s' % hmms_fp
                cmd += ' --promoter-attI --pdf --gbk --mute -v'
                cmd += ' %s' % fp_out
                self.cmds.setdefault(self.pool, []).append(cmd)

    def prep_metawrap(self):
        self.out = {}
        binners = ['maxbin2', 'metabat2', 'concoct']
        for group in self.pools[self.pool]:
            tmp_dir = '$TMPDIR/mtwrp_%s_%s' % (self.pool, group)
            out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
            self.soft.dirs.add(out_dir)
            bin_tools = [basename(x).split('_bins')[0] for x
                         in glob.glob('%s/*' % out_dir) if x.endswith('_bins')]
            contigs_fp = self.inputs[self.pool][group][1]
            # if not isfile(contigs_fp):
            #     continue
            fastq_fps = [fastq for sam in self.pools[self.pool][group]
                         for fastq in self.config.fastq[sam]]
            bins_dirs = []
            if len(bin_tools) != 3:
                cmd = 'metawrap binning'
                cmd += ' -o %s' % out_dir
                cmd += ' -a %s' % contigs_fp
                for bin_tool in self.soft.params.get('binners', binners):
                    if bin_tool in bin_tools:
                        continue
                    cmd += ' --%s' % bin_tool
                    bins_dirs.append('%s/%s_bins' % (out_dir, bin_tool))
                cmd += ' -t 4'
                cmd += ' -l 1500'
                cmd += ' --universal'
                cmd += ' %s' % ' '.join(fastq_fps)
                if not isdir('%s/work_files' % out_dir) or self.config.force:
                    self.cmds.setdefault(group, []).append(cmd)
            bins_dirs = bins_dirs + ['work_files']
            self.out[group] = bins_dirs
            self.soft.io['I']['d'].add(tmp_dir)
            self.soft.io['I']['f'].update(([contigs_fp] + fastq_fps))
            self.soft.io['O']['d'].update(bins_dirs)

    def prep_metawrap_refine(self):
        self.out = {}
        c = self.softs['metawrap'].params['min_completion']
        x = self.softs['metawrap'].params['min_contamination']
        for group in self.pools[self.pool]:
            out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
            bin_folders = self.inputs[self.pool][group][:-1]
            self.soft.io['I']['d'].update(bin_folders)
            self.soft.dirs.add(out_dir)

            has_bins = 0
            cmd = 'metawrap bin_refinement'
            cmd += ' -o %s' % out_dir
            for fdx, folder in enumerate(bin_folders):
                if len(glob.glob('%s/bin.*.fa' % folder)) or self.config.force:
                    has_bins += 1
                cmd += ' -%s %s' % (['A', 'B', 'C'][fdx], folder)
            cmd += ' -t %s' % self.soft.params['cpus']
            cmd += ' -c %s -x %s' % (c, x)

            # is_file = isfile('%s/metawrap_%s_%s_bins.stats' % (out_dir, c, x))
            # if not is_file and has_bins == len(bin_folders):
            if has_bins == len(bin_folders):
                self.cmds.setdefault(self.pool, []).append(cmd)
            out = '%s/metawrap_%s_%s' % (out_dir, c, x)
            self.out[group] = ['%s.stats' % out, '%s_bins' % out]
            self.soft.io['O']['d'].update(['%s.stats' % out, '%s_bins' % out])

    def prep_metawrap_reassemble(self):
        self.out = {}
        modes = ['permissive', 'strict']
        c = self.softs['metawrap'].params['min_completion']
        x = self.softs['metawrap'].params['min_contamination']
        for group in self.pools[self.pool]:
            self.out[group] = {}
            metawrap_ref = self.inputs[self.pool][group][-1]
            self.soft.io['I']['d'].add(metawrap_ref)
            for sam in self.pools[self.pool][group]:
                out_dir = '%s/%s/%s/%s' % (self.dir, self.pool, group, sam)
                self.soft.io['I']['f'].update(self.config.fastq[sam])
                self.soft.io['O']['d'].add(out_dir)
                self.soft.dirs.add(out_dir)
                cmd = 'metawrap reassemble_bins'
                cmd += ' -o %s' % out_dir
                cmd += ' -1 %s' % self.config.fastq[sam][0]
                cmd += ' -2 %s' % self.config.fastq[sam][1]
                cmd += ' -t %s' % self.soft.params['cpus']
                cmd += ' -m %s' % self.soft.params['mem_num']
                cmd += ' -c %s -x %s' % (c, x)
                cmd += ' -b %s' % metawrap_ref
                cmd += ' --skip-checkm --parallel\n'
                for m in self.softs['metawrap'].params.get('reassembly', modes):
                    mode_dir = '%s/reassembled_bins_%s' % (out_dir, m)
                    if not isdir(mode_dir) or self.config.force:
                        cmd += 'mkdir %s\n' % mode_dir
                        cmd += 'mv %s/reassembled_bins/*%s.fa %s/.\n' % (
                            out_dir, m, mode_dir)
                        self.cmds.setdefault(self.pool, []).append(cmd)
                    self.out[group].setdefault(sam, []).append(mode_dir)
                    self.soft.io['I']['d'].add(mode_dir)

    def prep_metawrap_blobology(self):
        self.out = {}
        modes = self.softs['metawrap'].params.get('blobology', ['coassembly'])
        for group in self.pools[self.pool]:
            self.out[group] = {}
            metawrap_ref = self.inputs[self.pool][group][-1]
            contigs_fp = self.softs['spades'].outputs[self.pool][group][1]
            fastqs = {sam: self.config.fastq[sam] for sam
                      in self.pools[self.pool][group]}
            self.soft.io['I']['f'].add(contigs_fp)
            self.soft.io['I']['d'].add(metawrap_ref)
            for mode in modes:
                self.out[group][mode] = {}
                out_dir = '%s/%s/%s/%s' % (self.dir, self.pool, group, mode)
                sams_fastqs = get_sams_fastqs(mode, fastqs)
                for sam, fastq in sams_fastqs.items():
                    if sam:
                        out_dir += '/%s' % sam
                    self.soft.io['O']['d'].add(out_dir)
                    self.soft.io['I']['f'].update(fastq)
                    self.soft.dirs.add(out_dir)
                    is_file = isfile('%s/contigs.binned.blobplot' % out_dir)
                    if not is_file or self.config.force:
                        cmd = 'metawrap blobology'
                        cmd += ' -o %s' % out_dir
                        cmd += ' -a %s' % contigs_fp
                        cmd += ' -t %s' % self.soft.params['cpus']
                        cmd += ' --bins %s' % metawrap_ref
                        cmd += ' %s' % ' '.join(fastq)
                        self.cmds.setdefault(group, []).append(cmd)
                        self.out[group][mode][sam] = out_dir

    def prep_metawrap_quantify(self):
        self.out = {}
        modes = self.softs['metawrap'].params.get('blobology', ['coassembly'])
        for group in self.pools[self.pool]:
            self.out[group] = {}
            metawrap_ref = self.inputs[self.pool][group][-1]
            contigs_fp = self.softs['spades'].outputs[self.pool][group][1]
            fastqs = {sam: self.config.fastq[sam] for sam
                      in self.pools[self.pool][group]}
            self.soft.io['I']['f'].add(contigs_fp)
            self.soft.io['I']['d'].add(metawrap_ref)
            for mode in modes:
                self.out[group][mode] = {}
                out_dir = '%s/%s/%s/%s' % (self.dir, self.pool, group, mode)
                sams_fastqs = get_sams_fastqs(mode, fastqs)
                for sam, fastq in sams_fastqs.items():
                    if sam:
                        out_dir += '/%s' % sam
                    self.soft.io['O']['d'].add(out_dir)
                    self.soft.io['I']['f'].update(fastq)
                    self.soft.dirs.add(out_dir)
                    cmd = 'metawrap quant_bins'
                    cmd += ' -b %s' % metawrap_ref
                    cmd += ' -o %s' % out_dir
                    cmd += ' -a %s' % contigs_fp
                    cmd += ' -t %s' % self.soft.params['cpus']
                    cmd += ' %s' % ' '.join(fastq)
                    self.cmds.setdefault(group, []).append(cmd)
                    self.out[group][mode][sam] = out_dir

    def metawrap_classify_annotate(self, command):
        self.out = {}
        prev = self.soft.prev
        for group in self.pools[self.pool]:
            self.out[group] = {}
            if prev == 'metawrap_refine':
                bin_dirs = {'': [self.inputs[self.pool][group][1]]}
            elif prev == 'metawrap_reassemble':
                bin_dirs = self.inputs[self.pool][group]
            else:
                raise ValueError('No metawrap classify_bins after %s' % prev)
            self.soft.io['I']['d'].update(bin_dirs)
            for sam, bin_dir in bin_dirs.items():
                for bin in bin_dir:
                    odir = '%s/%s/%s' % (self.dir, self.pool, group)
                    if 'permissive' in bin:
                        odir += '/permissive'
                    elif 'strict' in bin:
                        odir += '/strict'
                    if sam:
                        odir += '/%s' % sam
                    self.soft.io['O']['d'].add(odir)
                    self.out[group].setdefault(sam, []).append(odir)
                    cmd = 'metawrap %s' % command
                    cmd += ' -b %s' % bin
                    cmd += ' -o %s' % odir
                    cmd += ' -t %s' % self.soft.params['cpus']
                    if command == 'annotate_bins':
                        out = glob.glob('%s/bin_funct_annotations/*.tab' % odir)
                        if not out or self.config.force:
                            self.cmds.setdefault(group, []).append(cmd)
                    elif command == 'classify_bins':
                        out = '%s/bin_taxonomy.tab' % odir
                        if not isfile(out) or self.config.force:
                            self.cmds.setdefault(group, []).append(cmd)

    def prep_metawrap_classify(self):
        self.metawrap_classify_annotate('classify_bins')

    def prep_metawrap_annotate(self):
        self.metawrap_classify_annotate('annotate_bins')

    def prep_drep(self):
        self.out = {}
        genomes = get_drep_bins(self.soft.prev, self.pools, self.inputs)
        for pool, group_paths in genomes.items():
            self.out[pool] = {}
            for stringency, sam_paths in group_paths.items():
                drep_dir = '%s/%s' % (self.dir, pool)
                if stringency:
                    drep_dir += '/%s' % stringency
                drep_in, paths = get_drep_inputs(drep_dir, sam_paths)
                self.soft.io['I']['f'].update(([drep_in] + list(paths)))
                for algorithm in ['fastANI', 'ANIn']:
                    drep_out = '%s/%s' % (drep_dir, algorithm)
                    self.soft.dirs.add(drep_out)
                    log = '%s/log' % drep_out
                    figure = '%s/figure' % drep_out
                    data_table = '%s/data_tables' % drep_out
                    dereplicated_genomes = '%s/dereplicated_genomes' % drep_out
                    outputs = [data_table, log, figure, dereplicated_genomes]
                    self.out[pool][algorithm] = outputs
                    if len(glob.glob('%s/*.fa' % dereplicated_genomes)):
                        continue
                    if isdir(drep_out):
                        self.soft.io['I']['d'].add(drep_out)
                    cmd = 'dRep dereplicate'
                    cmd += ' %s' % drep_out
                    cmd += ' --S_algorithm %s' % algorithm
                    cmd += ' --ignoreGenomeQuality'
                    if len(list(paths)) > 5000 and algorithm == 'fastANI':
                        cmd += ' --multiround_primary_clustering'
                        cmd += ' --greedy_secondary_clustering'
                        cmd += ' --run_tertiary_clustering'
                    cmd += ' -pa 0.9 -sa 0.98 -nc 0.3 -cm larger'
                    cmd += ' -p %s' % self.soft.params['cpus']
                    cmd += ' -g %s' % drep_in
                    self.soft.io['O']['d'].update(outputs)
                    self.cmds.setdefault(pool, []).append(cmd)

    def prep_metamarker(self):
        io, cmd, outputs, dirs = metamarker(self.dir, self.inputs)
        self.cmds[''] = cmd
        self.out = outputs
        self.soft.dirs.update(dirs)
        self.soft.io['I']['f'].update(io['I']['f'])
        self.soft.io['I']['d'].update(io['I']['d'])
        self.soft.io['O']['f'].update(io['O']['f'])
        self.soft.io['O']['d'].update(io['O']['d'])

    def prep_multiqc(self):
        self.cmds = {}

    def prep_atropos(self):
        self.cmds = {}

    def prep_kneaddata(self):
        self.cmds = {}

    def prep_human_filtering(self):
        self.cmds = {}

    def prep_mocat(self):
        self.cmds = {}

    def prep_mapDamage(self):
        self.cmds = {}

    def prep_seqtk(self):
        self.cmds = {}

    def prep_metaclade(self):
        self.cmds = {}

    def prep_ezTree(self):
        self.cmds = {}

    def prep_instrain(self):
        self.cmds = {}

    def prep_quast(self):
        self.cmds = {}

    def prep_map__cazy_spades(self):
        self.cmds = {}

    def prep_map__cazy_macsyfinder(self):
        self.cmds = {}

    def prep_map__spades_bins(self):
        self.cmds = {}

    def prep_map__drep(self):
        self.cmds = {}

    def prep_mag_data(self):
        self.cmds = {}

    def prep_yamb(self):
        self.cmds = {}

    def prep_summarize_yamb(self):
        self.cmds = {}

    def prep_mycc(self):
        self.cmds = {}

    def prep_gtdbtk(self):
        self.cmds = {}

    def prep_graspx(self):
        self.cmds = {}

    def prep_checkm(self):
        self.cmds = {}

    def prep_antismash(self):
        self.cmds = {}

    def prep_deepvirfinder(self):
        self.cmds = {}

    def prep_deepvirfinder_extractVirs(self):
        self.cmds = {}

    def prep_metaxa2(self):
        self.cmds = {}

    def prep_closed_ref_qiime1(self):
        self.cmds = {}

    def prep_kraken2(self):
        self.cmds = {}

    def prep_prokka(self):
        self.cmds = {}

    def prep_phyloflash(self):
        self.cmds = {}

    def prep_wish(self):
        self.cmds = {}

    def prep_thmm(self):
        self.cmds = {}

    def prep_cazy(self):
        self.cmds = {}

    def prep_itasser(self):
        self.cmds = {}

    def get_scratch(self):
        if self.soft.name in [
                'spades', 'metamarker', 'diamond_annot', 'ioncom', 'humann2',
                'cazy', 'human_filtering', 'woltka', 'mycc']:
            self.softs[self.soft.name].scratch = 1

    def register_command(self):
        self.softs[self.analysis].cmds = dict(self.cmds)
        self.cmds = {}

    # def add_transfer_localscratch(self):
    #     pass
    #     # if self.scratch:
    #     #     self.add_transfer_localscratch()

    # def post_processes(self):
    #     if self.name.startswith('map__spades') or self.name in [
    #         'count_reads', 'humann2', 'midas', 'cazy', 'prodigal',
    #         'count_reads_grep', 'diamond_custom', 'hmmer_custom'
    #     ]:
    #         out_merge_samples = merge_samples(soft, output_paths)
    #         self.outputs['%s_merge' % self.name] = out_merge_samples
    #     if self.name in ['cazy']:
    #         out_cazy = cazy_post_processing(soft, soft_prev, output_paths)
    #         output_paths['%s_data' % self.name] = out_cazy
    #         cazy_sequence_analysis(soft, soft_prev, output_paths)
    #     if self.name in ['midas']:
    #         out_postprocess_samples = postprocess(soft, output_paths)
    #     if self.name in ['spades']:
    #         outs = contigs_per_sample(soft, output_paths)
    #         output_paths['%s_per_sample' % self.name] = outs
    #         outs = longest_contigs_per_sample(soft, output_paths)
    #         output_paths['%s_longest' % self.name] = outs
    #         prep_quast(soft, output_paths)
    #     if self.name in [
    #             'metawrap_analysis_reassemble_bins', 'metawrap_ref', 'yamb']:
    #         rename_bins(soft, output_paths)
    #     if self.name in ['drep']:
    #         faindex_contig_to_bin_drep(soft, output_paths)
    #         outs = longest_contigs_per_drep_bin(soft, output_paths)
    #         output_paths['%s_longest' % self.name] = outs
    #         outs = contigs_per_drep_bin(soft, output_paths)
    #         output_paths['%s_bins_contigs' % self.name] = outs
    #         outs = genes_per_mag(soft, output_paths)
    #         output_paths['%s_genes_per_mag' % self.name] = outs
    #     if self.name in ['yamb']:
    #         outs = yamb_across(soft, output_paths)
    #         output_paths['%s_maps' % self.name] = outs
    #     if self.name in ['read_mapping']:
    #         collect_reads_per_contig(soft, soft_prev, output_paths)
    #
    #     return all_sh, output_paths, mem_n_u, procs, cnd, pooling_groups, IOs

