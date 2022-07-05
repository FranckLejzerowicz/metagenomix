# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import yaml
import pkg_resources

import pandas as pd
from skbio.tree import TreeNode
from os.path import basename, isdir, isfile, splitext

from metagenomix._io_utils import (mkdr, get_pfam_wget_cmd, get_hmm_dat,
                                   get_pfams_cmd)

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')


class ReferenceDatabases(object):

    def __init__(self, config) -> None:
        self.config = config
        self.commands = {}
        self.cmds = {}
        self.paths = {}
        self.builds = {}
        self.database = ''
        self.pfams = {'terms': {}}
        self.hmms_dias = {}
        self.hmms_pd = pd.DataFrame()
        self.dbcan_meta = pd.DataFrame()
        self.cazys = {}
        self.gtdb = {}
        self.wol = {}
        self.formats = [
            'bbmap', 'blastn', 'bowtie2', 'burst', 'centrifuge', 'diamond',
            'hmmer', 'kraken2', 'minimap2', 'qiime2', 'utree']
        self.valid_databases = set()

    def run(self) -> None:
        self.get_formats()
        self.validate_databases()
        self.set_databases()

    def get_formats(self):
        """Override init (default) database formats with those from yaml file"""
        if 'formats' in self.config.databases:
            self.formats = self.config.databases['formats']

    def validate_databases(self) -> None:
        """Show and validate the presence of databases.
        Whether there are valid databases in the user-defined paths.
        """
        if len(self.config.databases):
            x = '%s\n' % ('#' * (11 + len(self.config.databases_yml)))
            print('\n%sDatabases: %s\n%s' % (x, self.config.databases_yml, x))
            for db, path in sorted(self.config.databases.items()):
                if db == 'formats':  # "formats" is not a database
                    continue
                if not path or not isinstance(path, str):  # must have a "path"
                    print('  - %s: path (char. string) missing (ignored)' % db)
                    continue
                if not self.config.dev:
                    if not isdir(path):  # not-found database will be ignored
                        print("  - %s: can't find %s (ignored)" % (db, path))
                    else:  # print database input if it found to exist
                        self.valid_databases.add(db)
                        print('  + %s' % yaml.dump({db: path}), end='')
                else:  # print database input if it found to exist in dev mode
                    self.valid_databases.add(db)
                    print('  + %s' % yaml.dump({db: path}), end='')
        else:
            print('No database passed to option `-d`')

    def set_databases(self) -> None:
        for database in sorted(self.valid_databases):
            self.database = database
            if hasattr(self, "set_%s" % database):
                # this is true for the following databases:
                #  - wol
                #  - mar
                #  - pfam
                #  - dbcan
                #  - ioncom
                #  - shogun
                getattr(self, "set_%s" % database)()
            else:
                self.set_path()

    def set_path(self) -> None:
        self.paths[self.database] = self.config.databases[self.database]
        self.builds[self.database] = self.get_db_formats()

    def register_command(self) -> None:
        self.commands.setdefault(self.database, {}).update(dict(self.cmds))
        self.cmds = {}

    # Below are the database-specific functions, potentially to make builds

    def get_db_formats(self):
        path = self.paths[self.database]
        formats = {}
        for db_format in self.formats:
            for subfolder in ['', 'databases/']:
                db_format_dir = path + '/%s' % subfolder + db_format
                if not self.config.dev and not isdir(db_format_dir):
                    continue
                formats[db_format] = db_format_dir
                break
        return formats

    def set_wol(self) -> None:
        wol_dir = self.config.databases[self.database]
        self.paths[self.database] = wol_dir
        self.builds[self.database] = self.get_db_formats()
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

    def set_pfam(self) -> None:
        pfams_dir = self.config.databases[self.database]
        self.paths[self.database] = pfams_dir
        self.pfams['dir'] = pfams_dir
        cmd = get_pfam_wget_cmd(pfams_dir)
        if cmd:
            self.cmds[''] = [cmd]
        else:
            self.hmms_pd = get_hmm_dat(pfams_dir)
        self.register_command()
        self.pfam_models()

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

    def get_pfam_terms(self, params):
        self.database = 'pfam'
        terms = {}
        hmm = '%s/Pfam-A.hmm' % self.pfams['dir']
        for term in params['terms']:
            term_pd = self.hmms_pd.loc[
                self.hmms_pd['DE'].str.lower().str.contains(term.lower())
            ].copy()
            if not term_pd.shape[0]:
                print('[hmmer] No Pfam model for term "%s" (ignored)' % term)
                continue
            pfams, cmd = get_pfams_cmd(hmm, term_pd, term, self.pfams['res'])
            terms[term] = pfams
            self.cmds[term] = [cmd]
        # collect the "pfams per term" for the parameters of the current tools
        params['terms'] = terms
        # add these "pfams per term" to a full, database-level collection
        self.pfams['terms'].update(terms)
        self.register_command()

    def get_dbcan_hmms(self) -> None:
        """Get all the .hmm files from the dbCAN database."""
        for root, _, files in os.walk(self.config.databases['dbcan']):
            for fil in files:
                if fil.endswith('.hmm'):
                    hmm = root + '/' + fil
                    rad = splitext(hmm)[0]
                    fas, dia = '%s.fas' % rad, '%s.dmnd' % rad
                    cmd = 'diamond makedb --in %s -d %s\n' % (fas, dia)
                    if 'cazy' not in self.hmms_dias:
                        self.hmms_dias['cazy'] = {}
                    self.hmms_dias['cazy'][basename(rad)] = [hmm, dia]
                    self.cmds.setdefault(basename(rad), []).append(cmd)

    def set_dbcan(self) -> None:
        self.get_dbcan_hmms()
        m = '%s/dbCAN-seq/metadata.txt' % self.config.databases['dbcan']
        if isfile(m):
            self.dbcan_meta = pd.read_csv(m, header=0, index_col=0, sep='\t')

    def set_shogun(self) -> None:
        shogun_config = self.config.databases['shogun']
        self.paths[self.database] = shogun_config
        for k, v in shogun_config.items():
            if k == 'rep82':
                if not isdir(v) or not glob.glob('%s/*' % v):
                    cmd = 'mkdir -p %s\ncd %s\n' % (v, v)
                    cmd += 'wget -i https://raw.githubusercontent.com/knights' \
                           '-lab/SHOGUN/master/docs/shogun_db_links.txt\n'
                    self.cmds['download_shogun'] = [cmd]


