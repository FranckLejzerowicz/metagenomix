# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from importlib import import_module, util
from metagenomix._io_utils import mkdr, get_edit_fastq_cmd


class Commands(object):

    def __init__(self, config, databases, workflow):
        self.config = config
        self.databases = databases
        self.softs = workflow.softs
        self.analysis = ''
        self.dir = ''
        self.cmds = {}
        self.method = None
        self.soft = None
        self.out = []
        self.holistics = ['multiqc', 'simka', 'simka_pcoa', 'simka_procrustes',
                          'qiita_wol', 'ezTree', 'metamarker', 'woltka', 'drep',
                          'strainphlan', 'mag_data']
        print()
        print('**************** CONFIG *****************')
        for i, j in self.config.__dict__.items():
            print(i, '\t:\t', j)
        print('**************** CONFIG *****************')
        print()

        print('**************** pipeline ***************')
        for idx, soft in self.softs.items():
            print(idx, soft.prev, soft.name, soft.params)
        print('*****************************************')

    def collect(self):
        for softs in self.config.pipeline:
            self.analysis = softs[-1]
            self.soft = self.softs[self.analysis]
            self.get_method()
            self.get_inputs()
            self.get_dir()
            self.generic_command()

    def get_method(self):
        """Get the command-preparing method from this class (for the tools that
        are easy to deal with), or from auxillary modules located in the `tools`
        submodules path."""
        func = 'prep_%s' % self.analysis
        if util.find_spec(self.analysis, 'metagenomix.tools'):
            mod = import_module('metagenomix.tools.%s' % self.analysis)
            self.method = getattr(mod, func)
        elif func in dir(self) and callable(getattr(self, func)):
            self.method = getattr(self, func)
        else:
            raise ValueError('No method for software %s' % self.analysis)

    def get_inputs(self):
        """Update the `inputs` attribute of the software object."""
        if not self.soft.prev or self.soft.name == 'map__drep':
            inputs = self.config.fastq
        else:
            inputs = self.softs[self.soft.prev].outputs
        self.softs[self.soft.name].inputs = inputs

    def get_dir(self):
        self.dir = '%s/%s/after_%s' % (
            self.config.output, self.soft.name, self.soft.prev)

    def generic_command(self):
        if self.soft.name in self.holistics:
            self.prep_job()
        elif self.soft.name == 'pooling':
            for group in self.config.pooling_groups:
                for pool, pool_pd in self.config.meta.groupby(group):
                    self.prep_pooling(pool, pool_pd.index.tolist(), group)
        else:
            for sam, fqs in sorted(self.softs[self.soft.name].inputs):
                self.prep_job(sam)

    def prep_job(self, sam=None):
        self.method(self.soft, self.softs, self.config.pooling_groups, sam)
        if 'simka' in self.soft.name or self.soft.name in ['drep', 'mag_data']:
            self.softs[self.soft.name].outputs[self.soft.name] = self.out
        elif sam and self.out:
            self.softs[self.soft.name].outputs[sam] = self.out

        if self.scratch:
            tool.add_transfer_localscratch()
        self.io = tool.io

    # def get_n_procs(self) -> None:
    #     """
    #     Get the number of nodes and CPUs.
    #     """
    #     self.nodes = 1
    #     if self.procs > 12:
    #         self.nodes = (((self.procs - 1) // 12) + 1)
    #         p = self.procs // self.nodes
    #         if self.nodes * p < self.procs:
    #             p += 1
    #         self.procs = p

    # def get_sh(self, sam: str = None) -> str:
    #     """Get the path to the job file for one sample
    #     (or one all-samples analysis).
    #
    #     Parameters
    #     ----------
    #     sam : str
    #         Current sample name (could be None).
    #
    #     Returns
    #     -------
    #     run_sh : str
    #         Path to the job file for one sample (or one all-samples analysis).
    #     """
    #     run_dir = '%s/%s/jobs/staged' % (self.config.dir, self.name)
    #     mkdr(run_dir)
    #     run_sh = '%s/staged_%s_after_%s' % (run_dir, self.name, self.prev)
    #     if sam:
    #         run_sh += '_%s' % sam.replace(' ', '_')
    #     run_sh += '.sh'
    #     return run_sh
    #
    # def get_outdir(self) -> str:
    #     """Get the path to the output folder for the current analysis.
    #
    #     Returns
    #     -------
    #     out_dir : str
    #         Path to the output folder for the current analysis.
    #     """
    #     out_dir = '%s/%s/after_%s' % (self.config.dir, self.name, self.prev)
    #     mkdr(out_dir)
    #     return out_dir

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

    def edit_fastqs(self):
        self.analysis = 'edit_fastq'
        for sam, reads in self.config.fastq.items():
            cmd1 = get_edit_fastq_cmd(reads, 1)
            cmd2 = get_edit_fastq_cmd(reads, 2)
            self.cmds[sam] = [cmd1, cmd2]

    def prep_count_reads(self, soft, softs, sam, pools):
        pass

    def prep_fastqc(self):
        self.cmds = {}

    def prep_multiqc(self):
        self.cmds = {}

    def prep_atropos(self):
        self.cmds = {}

    def prep_kneaddata(self):
        self.cmds = {}

    def prep_cutadapt(self):
        self.cmds = {}

    def prep_human_filtering(self):
        self.cmds = {}

    def prep_mocat(self):
        self.cmds = {}

    def prep_shogun(self):
        self.cmds = {}

    def prep_woltka(self):
        self.cmds = {}

    def prep_pooling(self, pool, sams, group):
        self.softs[self.soft.name].outputs = {group: {}}
        out_dir = self.dir + '/' + group
        mkdr(out_dir)

        merges = ['']
        if self.soft.prev == 'flash':
            merges = ['extendedFrags', 'notCombined_1', 'notCombined_2']

        for merge in merges:
            paths = [x for s in sams for x in self.soft.inputs[s] if merge in x]
            print(paths)
            if len(sams) > 1:
                fas = '%s/%s.%s.fasta' % (out_dir, pool, merge)
                fas = fas.replace(' ', '_').replace('..', '.')
                for pdx, path in enumerate(paths):
                    if pdx:
                        pooling_cmd = 'cat %s >> %s' % (path, fas)
                    else:
                        pooling_cmd = '\n\ncat %s > %s' % (path, fas)
                    self.softs[self.soft.name].cmds.append(pooling_cmd)
                self.softs[self.soft.name].io['I']['f'].extend(paths)
                self.softs[self.soft.name].io['O']['f'].append(fas)
                self.softs[self.soft.name].io['O']['d'].append(out_dir)
            else:
                if len(paths) > 1:
                    raise ValueError('Error in pooling group...')
                fas = paths[0]
            self.softs[self.soft.name][group].setdefault(pool, []).append(fas)

    def prep_flash(self):
        self.cmds = {}

    def prep_count_reads_grep(self, soft, softs, sam, pools):
        for inputs in soft.inputs:
            if input_path.endswith('fasta'):
                break
            elif input_path.endswith('fastq'):
                break
            elif input_path.endswith('fastq.gz'):
                break
        input_path_rad = dirname(input_path)
        outp_dir_sam_out = '%s/pipeline_reads_counts_%s.tsv' % (
        input_path_rad, sam)
        #     if 1:
        if not isfile(outp_dir_sam_out):
            for idx, input_path in enumerate(input_paths):
                if input_path.endswith('fasta'):
                    cur_count = 'cur_count_%s' % idx
                    cur_cmd = "%s=`wc -l %s | cut -d ' ' -f 1 | awk '{x=$1/2; print x}'`" % (
                    cur_count, input_path)
                    # cur_cmd = '%s=`grep -c ">" %s`' % (cur_count, input_path)
                elif input_path.endswith('fastq'):
                    with open(input_paths) as f:
                        for line in f:
                            break
                    #                 to_grep = ':'.join(line.split(':')[:3])
                    cur_count = 'cur_count_%s' % idx
                    # cur_cmd = '%s=`grep -c "%s" %s`' % (cur_count, to_grep, input_path)
                    cur_cmd = "%s=`wc -l %s | cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`" % (
                    cur_count, input_path)
                elif input_path.endswith('fastq.gz'):
                    with gzip.open(input_path, 'rb') as f:
                        for line in f:
                            break
                            #                 to_grep = ':'.join(line.decode().split(':')[:3])
                    cur_count = 'cur_count_%s' % idx
                    cur_cmd = "%s=`zcat %s | wc -l | cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`" % (
                    cur_count, input_path)
                    # cur_cmd = '%s=`zcat %s | grep -c "%s"`' % (cur_count, input_path, to_grep)
                count_reads_grep_cmds.append(cur_cmd)
                if idx:
                    cmd = 'echo "%s,%s,%s,%s,${%s}" >> %s' % (
                        soft_prev, sam, input_path,
                        splitext(basename(input_path))[0],
                        cur_count, outp_dir_sam_out
                    )
                else:
                    cmd = 'echo "%s,%s,%s,%s,${%s}" > %s' % (
                        soft_prev, sam, input_path,
                        splitext(basename(input_path))[0],
                        cur_count, outp_dir_sam_out
                    )
                count_reads_grep_cmds.append(cmd)

            IO['I']['f'].extend(input_paths)
            IO['O']['f'].append(outp_dir_sam_out)

        # add the outputs
        count_reads_grep_outs = outp_dir_sam_out
        self.cmds = {}

    def prep_mapDamage(self):
        self.cmds = {}

    def prep_seqtk(self):
        self.cmds = {}

    def prep_simka(self):
        self.cmds = {}

    def prep_simka_pcoa(self):
        self.cmds = {}

    def prep_simka_procrustes(self):
        self.cmds = {}

    def prep_plass(self):
        self.cmds = {}

    def prep_metamarker(self):
        self.cmds = {}

    def prep_metaclade(self):
        self.cmds = {}

    def prep_ezTree(self):
        self.cmds = {}

    def prep_metaphlan(self):
        self.cmds = {}

    def prep_phylophlan(self):
        self.cmds = {}

    def prep_strainphlan(self):
        self.cmds = {}

    def prep_instrain(self):
        self.cmds = {}

    def prep_spades(self):
        self.cmds = {}

    def prep_quast(self):
        self.cmds = {}

    def prep_drep(self):
        self.cmds = {}

    def prep_read_mapping(self):
        self.cmds = {}

    def prep_map__cazy_spades(self):
        self.cmds = {}

    def prep_map__cazy_macsyfinder(self):
        self.cmds = {}

    def prep_map__spades_bins(self):
        self.cmds = {}

    def prep_map__drep(self):
        self.cmds = {}

    def prep_map__spades_prodigal(self):
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

    def prep_metawrap_bin(self):
        self.cmds = {}

    def prep_metawrap_ref(self):
        self.cmds = {}

    def prep_metawrap_analysis_blobology(self):
        self.cmds = {}

    def prep_metawrap_analysis_reassemble_bins(self):
        self.cmds = {}

    def prep_metawrap_analysis_quant_bins(self):
        self.cmds = {}

    def prep_metawrap_analysis_classify_bins(self):
        self.cmds = {}

    def prep_metawrap_analysis_annotate_bins(self):
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

    def prep_prodigal(self):
        self.cmds = {}

    def prep_diamond_annot(self):
        self.cmds = {}

    def prep_macsyfinder(self):
        self.cmds = {}

    def prep_diamond_custom(self):
        self.cmds = {}

    def prep_hmmer_custom(self):
        self.cmds = {}

    def prep_thmm(self):
        self.cmds = {}

    def prep_humann2(self):
        self.cmds = {}

    def prep_cazy(self):
        self.cmds = {}

    def prep_integron_finder(self):
        self.cmds = {}

    def prep_itasser(self):
        self.cmds = {}

    def prep_ioncom(self):
        self.cmds = {}

    def prep_midas(self):
        print('LOCAL')
        self.cmds = {}

    def get_scratch(self):
        if self.name in [
                'spades', 'metamarker', 'diamond_annot', 'ioncom', 'humann2',
                'cazy', 'human_filtering', 'woltka', 'mycc']:
            self.scratch = 1
