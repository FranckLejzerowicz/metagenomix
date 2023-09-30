# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import glob

import pkg_resources

import pandas as pd
from skbio.tree import TreeNode
from os.path import basename, exists, isdir, isfile, splitext

from metagenomix._io_utils import (
    mkdr, wget_pfam, get_hmm_dat, get_hmms_dias_cmd)

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')


class ReferenceDatabases(object):

    def __init__(self, config) -> None:
        self.config = config
        self.commands = {}
        self.cmds = {}
        self.paths = {}
        self.build = ''
        self.builds = {}
        self.content = {}
        self.db = ''
        self.path = ''
        self.length = 0
        self.messages = {}
        self.formats = []
        self.formatted = False
        self.format = ''
        self.fdir = ''
        self.hmms_pd = pd.DataFrame()
        self.pfams = {}
        self.dmnd = {}
        self.hmms = {}
        self.dbcan_meta = pd.DataFrame()
        self.cazys = {}
        self.gtdb = {}
        self.wol = {}

    def run(self) -> None:
        self.get_formats()
        self.check_path()
        self.check_no_default()
        self.check_params_dbs()
        if len(self.config.databases):
            print('  * Found databases config "%s"' % self.config.databases_yml)
            self.get_length()
            self.print_formats()
            self.set_databases()
        else:
            print('  * No database passed to option `-d`')

    def get_formats(self):
        """
        Set all possible database formats, which also corresponds to the names
        of folders where a database is built.

        For example, format "blastn" is for the path of a database for which
        there is a "blastn" sub-folder containing the indices for BLASTn:

            If the database yaml file has a database called:
                my_db: /the/absolute/path/to/my_db_folder

            Then there must exist one of these folder:
                /the/absolute/path/to/my_db_folder/blastn
                /the/absolute/path/to/my_db_folder/databases/blastn
            which would contain the BLASTn database index files

        For example, format "bowtie2" is for the path of a database for which
        there is a "bowtie2" sub-folder containing the .bt2 indices:

            Then there must exist one of these folder:
                /the/absolute/path/to/my_db_folder/bowtie2
                /the/absolute/path/to/my_db_folder/databases/bowtie2
            which would contain the .bt2 files
        """
        self.formats = ['bbmap', 'blastn', 'bowtie2', 'bracken', 'burst',
                        'centrifuge', 'diamond', 'hmmer', 'kraken2',
                        'minimap2', 'qiime2', 'utree', 'fasta']

    def check_path(self):
        for db, path in sorted(self.config.databases.items()):
            if not path:
                sys.exit('[databases] No value for "%s"' % db)
            if not isinstance(path, str):
                sys.exit('[databases] "%s" value not a path string' % db)
            if path[0] != '/':
                sys.exit('[databases] Path to "%s" not absolute' % db)

    def check_no_default(self):
        if 'default' in self.config.databases:
            print('[databases] Database name "default" not allowed (ignored)')
            del self.config.databases['default']

    def check_params_dbs(self):
        config_dbs = set(self.config.databases)
        tool_dbs = ['metaxa2']
        params_dbs = set([z for x, y in self.config.params_dbs.items()
                          for z in y if x not in tool_dbs])
        missing_dbs = params_dbs.difference(config_dbs)
        missing_dbs.discard('default')
        if missing_dbs:
            print('[databases] %s databases not in "%s" (ignored):' % (
                    len(missing_dbs), self.config.databases_yml))
            for ddx, db in enumerate(sorted(missing_dbs)):
                print('[databases]  - %s\t: %s' % (ddx, db))

    def get_length(self):
        for db, path in sorted(self.config.databases.items()):
            db_length = len(db) + 5
            if db_length > self.length:
                self.length = db_length

    def print_formats(self):
        print(' ' * self.length, end='\t')
        for fmt in self.formats:
            print('%s ' % fmt, end='')
        print()

    def set_databases(self) -> None:
        for db, path in sorted(self.config.databases.items()):
            self.db = db
            self.path = path.rstrip('/')
            self.messages = {}
            self.set_paths()

    def set_paths(self):
        if self.config.dev:
            self.set_database()
        else:
            # if the path of a database is not found, it will be ignored
            if exists(self.path):
                self.set_database()
            else:
                print('  - %s: "%s" not found (ignored)' % (self.db, self.path))

    def set_database(self):
        self.print_db()
        if hasattr(self, "set_%s" % self.db):
            # specific treatment for some databases
            getattr(self, "set_%s" % self.db)()
            print()
        else:
            self.set_format()
        if self.config.verbose:
            for fmt, mess in self.messages.get(self.db, {}).items():
                print('[databases] "%s": %s" (%s)' % (self.path, mess, fmt))
        self.set_dmnd()
        self.set_hmms()
        self.set_path()

    # ------------------------------------------------------
    # ------------------------------------------------------
    def set_dmnd(self):
        pass

    def set_hmms(self):
        pass
    # ------------------------------------------------------
    # ------------------------------------------------------

    def print_db(self):
        gaps = self.length - (len(self.db) + 5)
        print('    + %s%s' % (self.db, (' ' * gaps)), end='\t')

    def set_path(self) -> None:
        self.paths[self.db] = self.path

    def check_format(self):
        self.formatted = False
        if isdir(self.fdir):
            self.formatted = True
            if hasattr(self, "check_format_%s" % self.format):
                getattr(self, "check_format_%s" % self.format)()

    def set_format(self) -> None:
        self.builds[self.db] = {}
        for fmt in self.formats:
            self.format = fmt
            for sub_folder in ['', 'databases/']:
                self.fdir = self.path + '/%s' % sub_folder + self.format
                self.check_format()
                if self.config.dev or self.formatted:
                    self.register_build()
                    print('+%s' % (' ' * len(self.format)), end='')
                    break
            else:
                print('-%s' % (' ' * len(self.format)), end='')
        print()

    def register_build(self):
        if self.format == 'bowtie2':
            self.builds[self.db][self.format] = str(self.build)
        else:
            self.builds[self.db][self.format] = str(self.fdir)

    def register_command(self) -> None:
        self.commands.setdefault(self.db, {}).update(dict(self.cmds))
        self.cmds = {}

    def get_content(self, exp_files, exp_dirs):
        self.content = dict((x, 'file') for x in exp_files)
        self.content.update(dict((x, 'folder') for x in exp_dirs))

    def check_content(self):
        missing = {}
        for content, typ in self.content.items():
            x = '%s/%s' % (self.config.databases[self.db], content)
            if not exists(x):
                missing[content] = typ
        return missing

    def show_exit(self, missing):
        if missing:
            prints = sorted([' - %s (%s)' % x for x in missing.items()])
            sys.exit('\n"%s": missing content in "%s"\n%s' % (
                self.db, self.config.databases[self.db], '\n'.join(prints)))

    def check_tool_db(self, exp_files=None, exp_dirs=None):
        if exp_files is None:
            exp_files = []
        if exp_dirs is None:
            exp_dirs = []
        self.get_content(exp_files, exp_dirs)
        missing = self.check_content()
        self.show_exit(missing)

    def set_metaxa2(self):
        dbs = ['SSU', 'LSU', 'SSU_SILVA128', 'SSU_SILVA123.1',
               'SSU_Typestrains', 'ATP9-NAD9', 'COI', 'cpn60', 'EF1_alpha',
               'ITS2', 'matK', 'rbcL', 'rpb1', 'rpb2', 'trnH', 'trnL']
        for db in dbs:
            if isdir('%s/%s' % (self.fdir, db)):
                exp_files = ['%s/blast.cutoffs.txt' % db, '%s/blast.nhr' % db,
                             '%s/blast.nin' % db, '%s/blast.nsd' % db,
                             '%s/blast.nhr' % db, '%s/blast.nsi' % db,
                             '%s/blast.nsq' % db, '%s/blast.taxonomy.txt' % db]
                self.check_tool_db(exp_files)

    def set_midas(self):
        exp_files = ['exclude.txt', 'genome_info.txt', 'genome_taxonomy.txt',
                     'README.txt', 'species_info.txt', 'species_tree.newick']
        exp_dirs = ['marker_genes', 'ontologies', 'pan_genomes', 'rep_genomes']
        self.check_tool_db(exp_files, exp_dirs)

    def set_platon(self):
        exp_files = [
            'conjugation.h3f', 'conjugation.h3i', 'conjugation.h3m',
            'conjugation.h3p', 'inc-types.fasta', 'mobilization.h3f',
            'mobilization.h3i', 'mobilization.h3m', 'mobilization.h3p',
            'mps.dmnd', 'mps.tsv', 'ncbifam-amr.h3f', 'ncbifam-amr.h3i',
            'ncbifam-amr.h3m', 'ncbifam-amr.h3p', 'ncbifam-amr.tsv',
            'orit.nhr', 'orit.nin', 'orit.nsq', 'refseq-plasmids.nhr',
            'refseq-plasmids.nin', 'refseq-plasmids.nsq', 'refseq-plasmids.tsv',
            'replication.h3f', 'replication.h3i', 'replication.h3m',
            'replication.h3p', 'rRNA.i1f', 'rRNA.i1i', 'rRNA.i1m', 'rRNA.i1p']
        self.check_tool_db(exp_files)

    def set_itasser(self):
        exp_files = ['download_lib.pl', 'README.txt']
        exp_dirs = [
            'abs', 'bin', 'blast', 'COACH', 'COFACTOR', 'common', 'data',
            'example', 'file2html', 'I-TASSERmod', 'PSSpred', 'ResQ', 'src']
        self.check_tool_db(exp_files, exp_dirs)

    def set_ioncom(self):
        exp_files = ['readme.txt', 'run_IonCom.pl', 'runninglist']
        exp_dirs = ['bin', 'model', 'output']
        self.check_tool_db(exp_files, exp_dirs)

    def set_macsyfinder(self):
        exp_files = ['models/readme.txt']
        exp_dirs = ['models/CasFinder', 'models/CONJScan_plasmids',
                    'models/TFF-SF', 'models/TFFscan', 'models/TXSScan']
        self.check_tool_db(exp_files, exp_dirs)

    def set_squeezemeta(self):
        exp_files = [
            'db/arc.all.faa', 'db/arc.hmm', 'db/arc.scg.faa',
            'db/arc.scg.lookup', 'db/bac.all.faa', 'db/bac.hmm',
            'db/bac.scg.faa', 'db/bac.scg.lookup', 'db/bacar_marker.hmm',
            'db/DB_BUILD_DATE', 'db/eggnog.dmnd', 'db/euk.hmm',
            'db/keggdb.dmnd', 'db/marker.hmm', 'db/mito.hmm', 'db/nr.dmnd',
            'db/Pfam-A.hmm', 'db/ReadMe', 'db/selected_marker_sets.tsv',
            'db/silva.nr_v132.align', 'db/silva.nr_v132.tax',
            'db/taxon_marker_sets.tsv']
        exp_dirs = [
            'db/distributions', 'db/genome_tree', 'db/hmms', 'db/hmms_ssu',
            'db/img', 'db/LCA_tax', 'db/pfam', 'db/test_data']
        self.check_tool_db(exp_files, exp_dirs)

    def set_humann(self):
        exp_files = [
            'mpa_v30_CHOCOPhlAn_201901.1.bt2',
            'mpa_v30_CHOCOPhlAn_201901.2.bt2',
            'mpa_v30_CHOCOPhlAn_201901.3.bt2',
            'mpa_v30_CHOCOPhlAn_201901.4.bt2', 'mpa_v30_CHOCOPhlAn_201901.fna',
            'mpa_v30_CHOCOPhlAn_201901_marker_info.txt',
            'mpa_v30_CHOCOPhlAn_201901.md5', 'mpa_v30_CHOCOPhlAn_201901.pkl',
            'mpa_v30_CHOCOPhlAn_201901.rev.1.bt2',
            'mpa_v30_CHOCOPhlAn_201901.rev.2.bt2']
        self.check_tool_db(exp_files)

    def set_metaphlan(self):
        exp_files = [
            'mpa_v30_CHOCOPhlAn_201901.1.bt2',
            'mpa_v30_CHOCOPhlAn_201901.2.bt2',
            'mpa_v30_CHOCOPhlAn_201901.3.bt2',
            'mpa_v30_CHOCOPhlAn_201901.4.bt2', 'mpa_v30_CHOCOPhlAn_201901.fna',
            'mpa_v30_CHOCOPhlAn_201901_marker_info.txt',
            'mpa_v30_CHOCOPhlAn_201901.md5', 'mpa_v30_CHOCOPhlAn_201901.pkl',
            'mpa_v30_CHOCOPhlAn_201901.rev.1.bt2',
            'mpa_v30_CHOCOPhlAn_201901.rev.2.bt2']
        self.check_tool_db(exp_files)

    def set_midas2_gtdb(self):
        exp_files = ['genomes.tsv', 'md5sum.json', 'metadata.tsv', 'README.txt']
        exp_dirs = ['markers', 'markers_models', 'pangenomes']
        self.check_tool_db(exp_files, exp_dirs)

    def set_checkm(self):
        exp_files = ['selected_marker_sets.tsv', 'taxon_marker_sets.tsv']
        exp_dirs = ['distributions', 'genome_tree', 'hmms', 'hmms_ssu',
                    'img', 'pfam', 'test_data']
        self.check_tool_db(exp_files, exp_dirs)

    def set_gtdbtk(self):
        exp_files = [
            'fastani/create_genome_paths.sh', 'fastani/genome_paths.tsv',
            'masks/gtdb_r207_ar53.mask', 'masks/gtdb_r207_bac120.mask',
            'metadata/metadata.txt', 'mrca_red/gtdbtk_r207_ar53.tsv',
            'mrca_red/gtdbtk_r207_bac120.tsv', 'msa/gtdb_r207_ar53.faa',
            'msa/gtdb_r207_bac120.faa', 'radii/gtdb_radii.tsv',
            'taxonomy/ar53_taxonomy_r207_reps.tsv',
            'taxonomy/bac120_taxonomy_r207_reps.tsv',
            'taxonomy/gtdb_taxonomy.tsv']
        exp_dirs = [
            'fastani/database', 'markers/pfam', 'markers/tigrfam',
            'pplacer/gtdb_r207_ar53.refpkg', 'pplacer/gtdb_r207_bac120.refpkg',
            'split/backbone', 'split/class_level']
        self.check_tool_db(exp_files, exp_dirs)

    def set_genomad(self):
        exp_files = [
            'genomad_db', 'genomad_db.dbtype', 'genomad_db.index',
            'genomad_db.lookup', 'genomad_db.source', 'genomad_db_h',
            'genomad_db_h.dbtype', 'genomad_db_h.index', 'genomad_db_mapping',
            'genomad_db_taxonomy', 'genomad_integrase_db',
            'genomad_integrase_db.dbtype', 'genomad_integrase_db.index',
            'genomad_integrase_db.lookup', 'genomad_integrase_db.source',
            'genomad_integrase_db_h', 'genomad_integrase_db_h.dbtype',
            'genomad_integrase_db_h.index', 'genomad_marker_metadata.tsv',
            'genomad_mini_db', 'genomad_mini_db.dbtype',
            'genomad_mini_db.index', 'genomad_mini_db.lookup',
            'genomad_mini_db.source', 'genomad_mini_db_h',
            'genomad_mini_db_h.dbtype', 'genomad_mini_db_h.index',
            'genomad_mini_db_mapping', 'genomad_mini_db_taxonomy',
            'mini_set_ids', 'names.dmp', 'nodes.dmp',
            'plasmid_hallmark_annotation.txt',  'version.txt',
            'virus_hallmark_annotation.txt']
        self.check_tool_db(exp_files)

    def set_uniref(self):
        exp_files = ['uniref90.fasta', 'uniref90.xml']
        self.check_tool_db(exp_files)

    def set_wol(self) -> None:
        self.set_format()
        wol_dir = self.config.databases[self.db]
        # WOL is a fasta file used to build indices
        fna_dir = '%s/genomes' % wol_dir
        genomes = glob.glob('%s/fna/*' % fna_dir)
        if not genomes:
            metadata_fp = '%s/wol/metadata.tsv' % wol_dir
            batch_down_fp = '%s/wol/batch_down.sh' % RESOURCES
            make_down_list_py = '%s/wol/make_down_list.py' % RESOURCES
            cmd = 'mkdir -p %s\n' % fna_dir
            cmd += 'cd %s\n' % fna_dir
            cmd += '%s %s > download.list\n' % (make_down_list_py, metadata_fp)
            cmd += 'bash %s download.list' % batch_down_fp
            self.cmds['download_wol_genomes'] = [cmd]
        else:
            self.wol['fna'] = genomes

        lineages_fp = '%s/wol/lineages.txt' % RESOURCES
        lineages_pd = pd.read_csv(lineages_fp, header=None, sep='\t')
        lineages_pd.columns = ['Genome ID', 'taxonomy']
        lineages_pd['species'] = [x.split(';')[-1].strip().replace(
            ' ', '_').replace( '[', '').replace(
            ']', '') for x in lineages_pd['taxonomy']]
        self.wol['taxonomy'] = lineages_pd

        tree_fp = '%s/wol/tree.nwk' % RESOURCES
        self.wol['tree'] = TreeNode.read(tree_fp, format='newick')

        sizes_fp = '%s/wol/genome_sizes.txt' % RESOURCES
        sizes = pd.read_csv(sizes_fp, header=None, sep='\t')
        sizes.columns = ['gid', 'length']
        self.wol['sizes'] = sizes
        self.register_command()

    def set_dbcan(self) -> None:
        # dbcan has hmms and fasta file that can be turned to diamond indices
        self.get_dbcan_hmms()
        self.register_command()
        m = '%s/dbCAN-seq/metadata.txt' % self.config.databases['dbcan']
        if self.config.dev:
            self.dbcan_meta = pd.DataFrame([
                ['GCF_000007365.1', 'Buchnera aphidicola Sg'],
                ['GCF_000007725.1', 'Buchnera aphidicola Bp'],
                ['GCF_000021065.1', 'Buchnera aphidicola Tuc7'],
                ['GCF_000021085.1', 'Buchnera aphidicola 5A'],
                ['GCF_000090965.1', 'Buchnera aphidicola Cc'],
                ['GCF_001280225.1', 'Buchnera aphidicola SBA'],
                ['GCF_000218545.1', 'Cellulomonas gilvus ATCC']],
                columns=['gen', 'genome_name'],
            )
            self.dbcan_meta.set_index('gen', inplace=True)
        elif isfile(m):
            self.dbcan_meta = pd.read_csv(m, header=0, index_col=0, sep='\t')

    def get_dbcan_hmms(self) -> None:
        """Get all the .hmm files from the dbCAN database."""
        self.hmms[self.db] = {}
        self.dmnd[self.db] = {}
        for root, _, files in os.walk(self.config.databases['dbcan']):
            for fil in files:
                if fil.endswith('.hmm'):
                    # get the full HMM file path
                    hmm = root + '/' + fil
                    # get the extension-striped and basename
                    rad = splitext(hmm)[0]
                    base = basename(fil)
                    # collect the fasta and dmnd files for this HMM'ed sequences
                    fas, dia = '%s.fas' % rad, '%s.dmnd' % rad
                    self.hmms[self.db][base] = hmm
                    self.dmnd[self.db][base] = dia
                    if not isfile(dia):
                        cmd = 'diamond makedb --in %s -d %s\n' % (fas, dia)
                        self.cmds.setdefault(base, []).append(cmd)

    def pfam_models(self):
        pfam_terms_dir = '%s/pfam_terms' % RESOURCES
        self.pfams['res'] = pfam_terms_dir
        mkdr(pfam_terms_dir)
        if self.config.show_pfams or self.config.purge_pfams:
            term_dirs = sorted(glob.glob('%s/*' % pfam_terms_dir))
            if self.config.show_pfams:
                if term_dirs:
                    print(' > already extracted Pfam models', end='')
                    if self.config.purge_pfams:
                        print(' and removed:')
                    else:
                        print(':')
                    for term_dir in term_dirs:
                        hmm_fps = glob.glob('%s/*.hmm' % term_dir)
                        print('   - %s (%s HMMs)' % (term_dir, len(hmm_fps)))
                        if self.config.purge_pfams:
                            for hmm_fp in hmm_fps:
                                os.remove(hmm_fp)
                            os.rmdir(term_dir)
                else:
                    print(' > No Pfam models previously extracted')

    def set_pfam(self) -> None:
        cmd = wget_pfam(self.path)
        if cmd:
            self.cmds[''] = [cmd]
        else:
            self.hmms_pd = get_hmm_dat(self.path)
        self.register_command()
        self.pfam_models()

    # def get_pfam(
    #         self,
    #         descriptions=[],
    #         accessions=[],
    #         interpro=[]
    # ):
    #     """
    #
    #     Parameters
    #     ----------
    #     descriptions : list
    #         List of description terms
    #     accessions : list
    #         List of Pfam accessions
    #     interpro : list
    #         List of InterPro tsv export files
    #     """
    #     if descriptions
    #
    #
    #     hmm = '%s/Pfam-A.hmm' % self.paths['pfam']
    #     for term in params.get('terms', []):
    #         term_pd = self.hmms_pd.loc[
    #             self.hmms_pd['DE'].str.lower().str.contains(term.lower())
    #         ].copy()
    #         if not term_pd.shape[0]:
    #             continue
    #         print(term_pd)
    #         print(term_pddsa)
    #         hmms_dias, cmd = get_hmms_dias_cmd(
    #             hmm, term_pd, term, self.pfams['res'])
    #         terms_hmms_dias[term] = hmms_dias
    #         self.cmds[term] = [cmd]
    #     # collect the "pfams per term" for the parameters of current softwares
    #     # params['terms'] = terms
    #     # add these "pfams per term" to a full, database-level collection
    #     self.pfams.update(accessions)
    #     self.register_command()

    def check_format_bowtie2(self):
        bt2s = glob.glob('%s/*.bt2*' % self.fdir)
        if len(bt2s) < 6:
            self.messages[self.format] = 'Less than six "bt2*" files'
            self.formatted = False
        bt2_4s = glob.glob('%s/*.4.bt2*' % self.fdir)
        if bt2_4s:
            for bt2_4 in bt2_4s:
                self.build = bt2_4.rsplit('.4.bt2', 1)[0]
                for e in ['1.bt2', '2.bt2', '3.bt2', 'rev.1.bt2', 'rev.2.bt2']:
                    f = '%s.%s' % (self.build, e)
                    if not isfile(f) and not isfile('%sl' % f):
                        self.messages[self.format] = 'Missing "bt2*"'
                        self.formatted = False
        else:
            self.messages[self.format] = 'No "1.bt2*" file'
            self.formatted = False

    def check_format_bbmap(self):
        files = {'genome': 'info.txt', 'index': '*.block'}
        for fdx, (folder, fil) in enumerate(files.items()):
            fold = '%s/*/ref/%s/1' % (self.fdir, folder)
            if not glob.glob(fold):
                self.messages[self.format] = 'no "%s" folder' % folder
                self.formatted = False
            fils = glob.glob('%s/%s' % (fold, fil))
            if not fils:
                self.messages[self.format] = 'incomplete "%s"' % fold
                self.formatted = False

    def check_format_blastn(self):
        nhrs = glob.glob('%s/*.nhr' % self.fdir)
        if not nhrs:
            self.messages[self.format] = 'Not a single ".nhr" file'
            self.formatted = False
        for nhr in nhrs:
            for ext in ['nin', 'nsq']:
                if not isfile('%s.%s' % (splitext(nhr)[0], ext)):
                    self.messages[self.format] = 'No ".%s" file' % ext
                    self.formatted = False

    def check_format_bracken(self):
        kraken = '%s/database*mers.kraken' % self.fdir
        if not glob.glob(kraken):
            self.messages[self.format] = 'No "%s" file' % kraken
            self.formatted = False
        kmer_d = '%s/database*mers.kmer_distrib' % self.fdir
        if not glob.glob(kmer_d):
            self.messages[self.format] = 'No "%s" file' % kmer_d
            self.formatted = False
        for kraken in glob.glob('%s/database*mers.kraken' % self.fdir):
            k = '%s.kmer_distrib' % splitext(kraken)[0]
            if not isfile(k):
                self.messages[self.format] = 'No "%s"' % basename(k)
                self.formatted = False

    def check_format_burst(self):
        acxs = glob.glob('%s/*.acx' % self.fdir)
        if not acxs:
            self.messages[self.format] = 'No ".acx" file'
            self.formatted = False
        for acx in acxs:
            if not isfile('%s.edx' % splitext(acx)[0]):
                b = basename(splitext(acx)[0])
                self.messages[self.format] = 'No "%s.edx" file' % b
                self.formatted = False

    def check_format_centrifuge(self):
        cfs = glob.glob('%s/*.cf' % self.fdir)
        if len(cfs) < 4:
            self.messages[self.format] = 'Less than four "cf" files'
            self.formatted = False
        cf_1s = glob.glob('%s/*.1.cf' % self.fdir)
        if cf_1s:
            for cf_1 in cf_1s:
                cf_rad = cf_1.rsplit('.1.cf', 1)[0]
                for ext in ['.2.cf', '.3.cf', '.4.cf']:
                    if not isfile('%s.%s' % (cf_rad, ext)):
                        self.messages[self.format] = 'Missing "cf"'
                        self.formatted = False
        else:
            self.messages[self.db] = 'No ".1.cf" file'
            self.formatted = False

    def check_format_kraken2(self):
        h_k2d = '%s/*/hash.k2d' % self.fdir
        o_k2d = '%s/*/opts.k2d' % self.fdir
        t_k2d = '%s/*/taxo.k2d' % self.fdir
        if not glob.glob(h_k2d):
            self.messages[self.format] = 'No "hash.k2d" file'
            self.formatted = False
        elif not glob.glob(o_k2d):
            self.messages[self.format] = 'No "opts.k2d" file'
            self.formatted = False
        elif not glob.glob(t_k2d):
            self.messages[self.format] = 'No "taxo.k2d" file'
            self.formatted = False

    def check_format_minimap2(self):
        mmi = '%s/*.mmi' % self.fdir
        if not glob.glob(mmi):
            self.messages[self.format] = 'No ".mmi" file'
            self.formatted = False

    def check_format_qiime2(self):
        qza = '%s/*.qza' % self.fdir
        if not glob.glob(qza):
            self.messages[self.format] = 'No ".qza" file'
            self.formatted = False

    def check_format_utree(self):
        if not glob.glob('%s/*.ctr' % self.fdir):
            self.messages[self.format] = 'No ".ctr" file'
            self.formatted = False

    def check_format_fasta(self):
        fastas = []
        for ext in ['fa', 'fas', 'fasta']:
            fastas.extend(glob.glob('%s/*.%s' % (self.fdir, ext)))
        if not fastas:
            self.messages[self.format] = 'No fasta file'
            self.formatted = False
