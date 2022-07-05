# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import sys
from skbio.io import read

from os.path import basename, dirname, isdir, isfile, splitext
from metagenomix._io_utils import (
    get_out_dir, write_hmms, io_update, reads_lines)

from metagenomix._cmds import caller


def prodigal_cmd(self, contigs_fp: str, out_dir: str, group: str = '') -> tuple:
    """Create command lines for Prodigal.

    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prodigal
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
    contigs_fp : str
        Contigs to predict
    out_dir : str
        Output directory
    group
        Pooling group
    """
    gbk = '%s/gene.coords.gbk' % out_dir
    proteins = '%s/protein.translations.fasta' % out_dir
    genes = '%s/potential.starts.fasta' % out_dir
    cmd = 'prodigal'
    cmd += ' -i %s' % contigs_fp
    cmd += ' -o %s' % gbk
    cmd += ' -a %s' % proteins
    cmd += ' -s %s' % genes
    cmd += ' -p %s' % self.soft.params['procedure']
    outputs = [gbk, proteins, genes]
    return cmd, outputs


def prodigal(self):
    """Create command lines for Prodigal.

    self : Commands class instance
        .dir : str
            Path to pipeline output folder for prodigal
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
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            contigs = self.inputs[self.pool][group][1]
            out_dir = '%s/%s/%s' % (self.dir, self.pool, group)
            cmd, outputs = prodigal_cmd(self, contigs, out_dir, group)
            self.outputs['outs'][group] = outputs
            if self.config.force or not isfile(outputs[0]):
                self.outputs['cmds'].setdefault(group, []).append(cmd)
                io_update(self, i_f=contigs, o_f=outputs, key=group)
    else:
        contigs = self.inputs[self.sam]
        out_dir = '%s/%s' % (self.dir, self.sam)
        cmd, outputs = prodigal_cmd(self, contigs, out_dir)
        self.outputs['outs'].extend(outputs)
        self.outputs['cmds'].append(cmd)
        io_update(self, i_f=contigs, o_f=outputs)


def macsyfinder_cmd(self, fp, o_dir, folder, group=None) -> tuple:
    cmd = ''
    outs = []
    io_update(self, i_f=fp, i_d=folder, key=group)
    self.outputs['dirs'].append(o_dir)
    for model in self.soft.params['models']:
        out_dir = '%s/%s' % (o_dir, model)
        outs.append(out_dir)
        self.outputs['dirs'].append(out_dir)
        io_update(self, o_d=out_dir, key=group)
        res = '%s/macsyfinder.log' % out_dir
        if self.config.force or not isfile(res):
            cmd += 'macsyfinder'
            cmd += ' --db-type %s' % self.soft.params['db_type']
            cmd += ' --replicon-topology %s' % self.soft.params['topology']
            cmd += ' --e-value-search %s' % self.soft.params['evalue']
            cmd += ' --coverage-profile %s' % self.soft.params['coverage']
            cmd += ' --out-dir %s' % out_dir
            cmd += ' --res-search-suffix _hmm.tsv'
            cmd += ' --res-extract-suffix _out.tsv'
            cmd += ' --worker %s' % self.soft.params['cpus']
            cmd += ' --sequence-db %s' % fp
            cmd += ' --models-dir %s/data/models' % folder
            cmd += ' --models %s all' % model
            cmd += ' --verbosity\n'
    return cmd, outs


def macsyfinder(self) -> None:
    """Create command lines for MacSyFinder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for macsyfinder
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
    folder = self.databases.paths['macsyfinder']
    if self.pool in self.pools:
        self.outputs['outs'] = {}
        for group in self.pools[self.pool]:
            o_dir, fp = get_out_dir(self, self.pool, group)
            cmd, outs = macsyfinder_cmd(self, fp, o_dir, folder, group)
            self.outputs['outs'][group] = outs
            if cmd:
                self.outputs['cmds'].setdefault(group, []).append(cmd)
    else:
        o_dir, fp = get_out_dir(self, self.sam)
        fp_edit = '%s_edit.fasta' % fp.replace('.fasta', '')
        cmd_edit = 'header_space_replace.py -i %s -o %s --n\n' % (fp, fp_edit)
        cmd, outs = macsyfinder_cmd(self, fp_edit, o_dir, folder)
        self.outputs['outs'].extend(outs)
        if cmd:
            cmd = cmd_edit + cmd
            self.outputs['cmds'].append(cmd)


def prepare_ioncom_inputs(self, ion_com: str, fp: str) -> str:
    """Create the preformatting command lines before running IonCom.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for IonCom
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
    ion_com : str
        Path to the IonCom standalone
    fp : str
        Path to the output file
    """
    cmd = 'prepare_ioncom_inputs.py'
    cmd += ' -d %s -p %s -n 1' % (ion_com, fp)
    cmd += '%s/run_IonCom.pl' % ion_com
    return cmd


def ioncom(self) -> None:
    """Create command lines for IonCom.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for IonCom
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
    # -------------------------------------------------------
    # --- WILL NEED TO BE REVISED AFTER I-TASSER DOWNLOAD ---
    # -------------------------------------------------------
    i_tasser_libs = self.databases.ioncom['itasser']
    i_tasser = self.soft.params['itasser']
    ion_com = '%s/IonCom_standalone' % self.databases.paths['ioncom']

    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            o_dir, fp = get_out_dir(self, self.pool, group)

            cmd = prepare_ioncom_inputs(self, ion_com, fp)
            self.outputs['cmds'].setdefault(group, []).append(cmd)

            output = '%s/output' % ion_com
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
            self.outputs['cmds'].setdefault(group, []).append(cmd)

            self.outputs['dir'].append(o_dir)
            self.outputs['outs'].setdefault(self.pool, []).append(o_dir)
            io_update(self, i_f=fp, i_d=[ion_com, i_tasser, i_tasser_libs],
                      o_d=o_dir, key=group)


def integron_finder(self):
    """Create command lines for integron_finder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder for integron_finder
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
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            o_dir, fp = get_out_dir(self, self.pool, group)
            fp_out = fp.replace('.fasta',
                                '_len%s.fasta' % self.soft.params['min_length'])
            hmms_fp = write_hmms(self)

            cmd = 'filter_on_length.py'
            cmd += ' -i %s' % fp
            cmd += ' -o %s' % fp_out
            cmd += ' -t %s\n' % self.soft.params['min_length']

            cmd += 'integron_finder'
            if self.soft.params['local_max']:
                cmd += ' --local-max'
            if self.soft.params['promoter_attI']:
                cmd += ' --promoter-attI'
            if self.soft.params['mute']:
                cmd += ' --mute'
            if self.soft.params['pdf']:
                cmd += ' --pdf'
            if self.soft.params['gbk']:
                cmd += ' --gbk'
            if self.soft.params['union_integrases']:
                cmd += ' --union-integrases'

            cmd += ' --verbose'
            cmd += ' --outdir %s' % o_dir
            cmd += ' --cpu %s' % self.soft.params['cpus']
            if self.databases.hmms_dias:
                cmd += ' --func-annot'
                cmd += ' --path-func-annot %s' % hmms_fp
                io_update(self, i_f=hmms_fp, key=group)
            if self.soft.params['prot_file']:
                cmd += ' --prot-file %s' % self.soft.params['prot_file']
            if self.soft.params['attc_model']:
                cmd += ' --attc-model %s' % self.soft.params['attc_model']
                cmd += ' --evalue-attc %s' % self.soft.params['evalue_attc']
                cmd += ' --max-attc-size %s' % self.soft.params['max_attc_size']
                cmd += ' --min-attc-size %s' % self.soft.params['min_attc_size']
            if self.soft.params['topology_file']:
                cmd += ' --topology-file %s' % self.soft.params['topology_file']
            else:
                cmd += ' --%s' % self.soft.params['topology']

            cmd += ' %s' % fp_out

            self.outputs['dirs'].append(o_dir)
            self.outputs['outs'].setdefault(self.pool, []).append(o_dir)
            self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=fp, o_d=o_dir, key=group)













def custom_cmd(self, input_file, out_dir, dia_hmm):
    outs = {}
    for target, gene_hmm_dia in self.databases.hmms_dias.items():
        if target == 'cazy':
            continue
        for gene, (hmm, dia) in gene_hmm_dia.items():
            sam_dir = '%s/%s' % (out_dir, gene)
            self.outputs['dirs'].append(sam_dir)
            out = '%s/%s.tsv' % (sam_dir, self.sam)
            if not isfile(out):
                if dia_hmm == 'hmm':
                    stdout = '%s_hmmer.out' % splitext(out)[0]
                    io_update(self, i_f=hmm, o_f=[out, stdout], key=target)
                    cmd = 'hmmsearch'
                    cmd += ' --cut_tc'
                    cmd += ' --tblout %s' % out
                    cmd += ' -o %s' % stdout
                    cmd += ' --cpu %s' % self.soft.params['cpus']
                    cmd += ' %s %s' % (hmm, input_file)
                else:
                    io_update(self, i_f=dia, o_f=out, key=target)
                    tmp = '%s/%s_tmp' % (sam_dir, self.sam)
                    cmd = 'mkdir -p %s\n' % tmp
                    cmd += 'diamond blastp'
                    cmd += ' -d %s' % dia
                    cmd += ' -q %s' % input_file
                    cmd += ' -o %s' % out
                    cmd += ' -k 1 --id 80'
                    cmd += ' -p %s' % self.soft.params['cpus']
                    cmd += ' -t %s' % tmp
                self.outputs['cmds'].setdefault(target, []).append(cmd)
            else:
                io_update(self, i_f=out, key=target)
            out2 = out.replace('.tsv', '_contigs.tsv')
            if not isfile(out2):
                io_update(self, o_f=out2, key=target)
                cmd = 'extract_custom_searched_contigs.py'
                cmd += ' -i %s' % out
                cmd += ' -o %s' % out2
                self.outputs['cmds'].setdefault(target, []).append(cmd)
            outs.setdefault(target, []).extend([out, out2])
    return outs


def prep_custom(self, dia_hmm):
    self.outputs['outs'] = {}
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            self.sam = group
            o_dir, fp = get_out_dir(self, self.pool, group)
            outs = self.custom_cmd(fp, o_dir, dia_hmm)
            io_update(self, i_f=fp, key=group)
            self.outputs['outs'][group] = outs
    else:
        o_dir, fp = get_out_dir(self, self.sam)
        outs = self.custom_cmd(fp, o_dir, dia_hmm)
        io_update(self, i_f=fp)
        self.outputs['outs'] = outs


def prep_diamond_custom(self):
    self.prep_custom('dia')


def prep_hmmer_custom(self):
    self.prep_custom('hmm')




















def diamond_annot_cmd(self, o_dir, fp, tmp_dir, group=None) -> tuple:
    cmd, outs = '', []
    io_update(self, i_f=fp, key=group)
    self.soft.dirs.add(o_dir)
    k = self.soft.params['k']
    for db, dmnds in self.soft.params['databases'].items():
        for dmnd in dmnds:
            base_db = basename(dmnd).rstrip('.dmnd')
            out = '%s/%s_k%s_%s.tsv' % (o_dir, self.sam, k, base_db)
            outs.append(out)
            cmd += 'mkdir -p %s\n' % tmp_dir
            cmd += 'diamond blastp'
            cmd += ' -d %s' % dmnd
            cmd += ' -q %s' % fp
            cmd += ' -o %s' % out
            cmd += ' -k %s' % k
            cmd += ' -p %s' % self.soft.params['cpus']
            cmd += ' --id %s' % self.soft.params['identity']
            cmd += ' -t %s\n' % tmp_dir
            io_update(self, o_f=out, key=group)

    return cmd, outs


def diamond(self):
    tmp_dir = '$TMPDIR/diamond_annot_%s' % self.sam
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            o_dir, fp = get_out_dir(self, self.pool, group)
            cmd, outs = diamond_annot_cmd(self, o_dir, fp, tmp_dir, group)
            self.outputs['outs'].setdefault(self.pool, []).extend(outs)
            if cmd:
                self.outputs['cmds'].setdefault(group, []).append(cmd)

    else:
        o_dir, fp = get_out_dir(self, self.sam)
        cmd, outs = diamond_annot_cmd(self, o_dir, fp, tmp_dir)
        self.outputs['outs'].extend(outs)
        if cmd:
            self.outputs['cmds'].append(cmd)


def search(self) -> None:
    """Create command lines for custom searches based on DIAMOND or HMMER.

    Parameters
    ----------
    self : Commands class instance
        Contains all the attributes needed for binning on the current sample
    """
    # This function splits the name of the software and calls as function (
    # from this module) the part of the software name that follows the first
    # underscore, e.g. software "search_diamond" would call `diamond()`
    caller(self, __name__)


def get_antismash_cmd(self, fa, base, out):
    cmd = '\nantismash'
    for boolean in ['rre', 'asf', 'cassis', 'pfam2go', 'tigrfam', 'fullhmmer',
                    'cc-mibig', 'cb_general', 'smcog_trees', 'clusterhmmer',
                    'cb_subclusters', 'cb_knownclusters']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --tta-threshold %s' % self.soft.params['--tta-threshold']
    cmd += ' --genefinding-tool %s' % self.soft.params['genefinding_tool']
    cmd += ' --html-description after_%s' % self.soft.prev
    cmd += ' --taxon %s' % self.soft.params['taxon']
    cmd += ' --cpus %s' % self.soft.params['cpus']
    cmd += ' --output-basename %s' % base
    cmd += ' --html-title "%s: %s"' % base
    cmd += ' --output-dir %s' % out
    cmd += ' %s\n' % fa
    return cmd


def antismash(self):
    if self.soft.prev != 'drep':
        sys.exit('[%s] Not possible (yet) after %s (only drep)' % (
            self.soft.name, self.soft.prev))
    for pool in self.pools:
        for algo, folder in self.inputs[pool].items():
            out_dir = '%s/%s/%s' % (self.dir, pool, algo)
            io_update(self, o_d=out_dir, key=pool)
            for fa in glob.glob('%s/*.fa' % folder):
                io_update(self, i_f=fa, key=pool)
                base = fa.split('/')[-1].replace('.fa', '')
                out = '%s/%s' % (out_dir, base)
                self.outputs['dirs'].append(out)
                if self.config.force or not glob.glob('%s/*' % out):
                    cmd = get_antismash_cmd(self, fa, base, out)
                    self.outputs['cmds'].setdefault(
                        (pool, algo), []).append(cmd)


def get_prokka_cmd(self, contigs, out_dir, prefix, config, cols):
    cmd = 'prokka'
    cmd += ' --mincontiglen %s' % self.soft.params['mincontiglen']
    cmd += ' --cpus %s' % self.soft.params['cpus']
    for boolean in ['metagenome', 'notrna', 'norrna', 'annotation']:
        if self.soft.params[boolean]:
            cmd += ' --%s' % boolean.replace('_', '-')
    cmd += ' --force'
    for col in cols:
        if col and config[col]:
            cmd += ' --%s %s' % (col, config[col])
    cmd += ' --prefix %s' % prefix
    cmd += ' --outdir %s' % out_dir
    cmd += ' %s' % contigs
    return cmd


def get_prokka_config(self, cols):
    configs = []
    cols += ['proteins']
    if 'taxa' in self.soft.params:
        with open(self.soft.params['taxa']) as f:
            for ldx, line in enumerate(f):
                ls = line.strip('\n').split('\t')
                if not ldx:
                    continue
                configs.append({cols[idx]: val for idx, val in enumerate(ls)})
    return configs


def prokka(self):
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            out = '%s/%s/%s' % (self.dir, self.pool, group)
            contigs = self.inputs[self.pool][group][1]
            if self.soft.params['config']:
                cols = ['genus', 'species', 'strain', 'plasmid']
                configs = get_prokka_config(self, cols)
            else:
                cols = ['']
                configs = {'': 'no_config'}

            cmd = ''
            for config in configs:
                prefix = '_'.join([config[x] for x in cols if config[x]])
                if config['proteins']:
                    prefix += '_%s' % splitext(basename(config['proteins']))[0]
                file_out = '%s/%s.out' % (out, prefix)
                if self.config.force or not isfile(file_out):
                    cmd += get_prokka_cmd(self, contigs, out,
                                          prefix, config, cols)

            self.outputs['dirs'].append(out)
            self.outputs['outs'].setdefault(group, []).append(out)
            if cmd:
                self.outputs['cmds'].setdefault(group, []).append(cmd)
            io_update(self, i_f=configs, o_d=out, key=group)


def get_barrnap_cmd(self, out, fasta):
    cmd = 'barrnap'
    cmd += ' --kingdom %s' % self.soft.params['kingdom']
    cmd += ' --threads %s' % self.soft.params['cpus']
    cmd += ' --reject %s' % self.soft.params['reject']
    cmd += ' --lencutoff %s' % self.soft.params['lencutoff']
    cmd += ' --evalue %s' % self.soft.params['evalue']
    if self.soft.params['incseq']:
        cmd += ' --incseq'
    cmd += ' --outseq %s' % out
    cmd += ' %s' % fasta
    return cmd


def get_ccmap_cmd(self, out, fasta):
    out_dir = dirname(out)
    cmd = ''
    if self.config.force:
        cmd += 'if [ -d "%s" ]; then rm -rf %s; fi\n' % (out_dir, out_dir)
    cmd += 'ccfind %s %s' % (fasta, out_dir)
    cmd += ' --terminal-fragment-size %s' % self.soft.params[
        'terminal_fragment_size']
    cmd += ' --min-percent-identity %s' % self.soft.params[
        'min_percent_identity']
    cmd += ' --min-aligned-length %s' % self.soft.params[
        'min_aligned_length']
    cmd += ' --ncpus %s' % self.soft.params['cpus']
    if self.soft.params['preserve_tmpdir']:
        cmd += ' --preserve-tmpdir'
    return cmd


def generic_on_fasta(self, tool):
    func = 'get_%s_cmd' % self.soft.name
    if self.pool in self.pools:
        for group in self.pools[self.pool]:
            fastas_d = {}
            if self.soft.prev == 'spades':
                fastas_d = {'': [self.inputs[self.pool][group][1]]}
            elif self.soft.prev == 'drep':
                fastas_d = self.inputs[self.pool][group]
            for _, fastas in fastas_d.items():
                for fasta in fastas:
                    out = '%s/%s/%s/%s.fas' % (self.dir, self.pool, group,
                                               splitext(basename(fasta))[0])
                    self.outputs['outs'].setdefault(
                        (self.pool, group), []).append(out)
                    self.outputs['dirs'].append(dirname(out))
                    io_update(self, i_f=fasta, o_d=out, key=group)
                    if tool == 'barrnap':
                        condition = isfile(out)
                    elif tool == 'ccmap':
                        condition = isdir(dirname(out))
                    if self.config.force or not condition:
                        cmd = globals()[func](self, out, fasta)
                        self.outputs['cmds'].setdefault(group, []).append(cmd)
    else:
        for fasta in self.inputs[self.sam]:
            out = '%s/%s/%s.fas' % (self.dir, self.sam,
                                    splitext(basename(fasta))[0])
            self.outputs['dirs'].append(dirname(out))
            self.outputs['outs'].append(out)
            io_update(self, i_f=fasta, o_d=out)
            if tool == 'barrnap':
                condition = isfile(out)
            elif tool == 'ccmap':
                condition = isdir(dirname(out))
            if self.config.force or not condition:
                cmd = globals()[func](self, out, fasta)
                self.outputs['cmds'].append(cmd)


def barrnap(self):
    generic_on_fasta(self, self.soft.name)


def ccmap(self):
    generic_on_fasta(self, self.soft.name)


def write_dbcan_subset(self, taxa: list, folder: str, name: str) -> str:
    """Write the fasta file suset to the target features and
    the command to make it a diamond database.

    Parameters
    ----------
    taxa : list
        Taxa from the user file.
    folder : str
        dbCAN-Seq folder.
    name : str
        Current subset name.

    Returns
    -------
    cmd : str
        Command to make the diamond db from the subsets fasta file.
    """
    path = self.databases['dbcan']
    cmd = ""
    fas_fp = '%s/%s.fa' % (folder, name)
    dia_fp = '%s.dmnd' % splitext(fas_fp)[0]
    io_update(self, i_f=fas_fp, o_f=dia_fp)
    if not isfile(dia_fp):
        with open(fas_fp, 'w') as o:
            taxon_found = False
            for taxon in taxa:
                meta_taxon_pd = self.databases.dbcan_meta.loc[
                    self.databases.dbcan_meta.genome_name.str.contains(taxon),:]
                if not meta_taxon_pd.shape[0]:
                    continue
                gcf_dir = '%s/dbCAN-seq/CAZyme_seq_list' % path
                if not isdir(gcf_dir):
                    os.makedirs(gcf_dir)
                for gcf in set(meta_taxon_pd.index):
                    gcf_fas = '%s/%s.fasta' % (gcf_dir, gcf)
                    if isfile(gcf_fas):
                        for e in read(gcf_fas, 'fasta'):
                            o.write('>%s\n%s\n' % (e.metadata['id'], e))
                        taxon_found = True
        if taxon_found:
            cmd = "diamond makedb --in %s -d %s\n" % (fas_fp, dia_fp)
    self.outputs['cmds'].append(cmd)


def set_dbcan_taxa(self) -> None:
    for name, fp in self.soft.params['taxa'].items():
        taxa = list(reads_lines(fp))
        if not taxa:
            continue
        folder = '%s/subsets' % self.dir
        if not isdir(folder):
            os.makedirs(folder)
        self.databases.cazys[name] = folder
        write_dbcan_subset(self, taxa, folder, name)


def cazy(self):
    set_dbcan_taxa(self)
