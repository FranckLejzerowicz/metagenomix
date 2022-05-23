# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import glob
import subprocess
import pkg_resources

import pandas as pd
from skbio.io import read
from skbio.tree import TreeNode
from os.path import basename, dirname, isdir, isfile, splitext

from metagenomix._io_utils import (mkdr, get_pfam_file, get_hmm_dat,
                                   get_pfams_cmd, reads_lines)

RESOURCES = pkg_resources.resource_filename('metagenomix', 'resources')


class ReferenceDatabases(object):

    def __init__(self, config) -> None:
        self.config = config
        self.databases_commands = {}
        self.tmps = []
        self.cmds = {}
        self.database = ''
        self.hmms_dias = {}
        self.hmms_pd = pd.DataFrame()
        self.wol = {}
        self.dbcan_meta = pd.DataFrame()
        self.krakens = {}
        self.cazys = {}
        self.midas = {}
        self.metaphlan = {}
        self.macsyfinder = {}
        self.humann = {}
        self.shogun = {}
        self.diamond = {}
        self.ioncom = {}

    def init(self) -> None:
        if self.show_loaded_databases():
            self.set_databases()

    def show_loaded_databases(self) -> bool:
        """Show and validate the presence of databases.

        Returns
        -------
        set_databases : bool
            Whether there are valid databases in the user-defined paths.
        """
        set_databases = False
        if len(self.config.databases):
            print('Databases in', self.config.databases_yml)
            invalid_keys = []
            for database, database_d in sorted(self.config.databases.items()):
                # if 'path' not in database_d or not isdir(database_d['path']):
                if 'path' not in database_d:
                    print('[Database] "%s" folder not found' % database)
                    invalid_keys.append(database)
                    continue
                else:
                    if isinstance(database_d['path'], dict):
                        print('- %s:' % database)
                        for k, v in database_d['path'].items():
                            print('  - %s: %s' % (k, v))
                    else:
                        print('- %s: %s' % (database, database_d['path']))
            if invalid_keys:
                for invalid_key in invalid_keys:
                    del self.config.databases[invalid_key]
            if len(self.config.databases):
                set_databases = True
        else:
            print('No database passed using option `-d`')
        return set_databases

    def set_databases(self) -> None:
        for database, keys in self.config.databases.items():
            db_method = getattr(self, "set_%s" % database)
            db_method()

    def set_squeezemeta(self) -> None:
        self.database = 'squeezemeta'

    def set_mar(self) -> None:
        self.database = 'mar'
        path = self.config.databases['mar']['path']
        mol_type_per_db = {
            'BLAST': [['nucleotides', 'fna'], ['proteins', 'faa']],
            'Genomes': [['genomic', 'fa'], ['protein', 'faa']]}
        for db_type, types in mol_type_per_db.items():
            # Build the database from concatenation of MArDB and MarRef
            for (typ, ext) in types:
                name = '%s_%s' % (db_type,typ)
                out_fp = '%s/%s/%s.%s' % (path, name, name, ext)
                if isfile(out_fp):
                    continue
                out_dir = dirname(out_fp)
                mkdr(out_dir)
                cmd = ''
                for idx, mar in enumerate(['MarRef', 'MarDB']):
                    if db_type == 'Genomes':
                        to_glob = '%s/%s/%s/*/*%s.%s' % (
                            path, mar, db_type, typ, ext)
                    elif db_type == 'BLAST':
                        to_glob = '%s/%s/%s/%s/*%s_V3.%s' % (
                            path, mar, db_type, typ, typ, ext)
                    for path in glob.glob(to_glob):
                        if idx:
                            cmd += 'cat %s >> %s\n' % (path, out_fp)
                        else:
                            cmd += 'cat %s > %s\n' % (path, out_fp)
                self.cmds.setdefault(name, []).append(cmd)
        self.register_command()

    def set_pfam(self) -> None:
        self.database = 'pfam'
        pfam_dir = self.config.databases['pfam']['path']
        pfam_job = '%s/jobs' % pfam_dir
        mkdr(pfam_job)
        fas_dir = '%s/fastas' % pfam_dir
        mkdr(fas_dir)
        hmm = '%s/Pfam-A.hmm' % pfam_dir
        fas = '%s/Pfam-A.fasta' % pfam_dir
        dat = '%s/Pfam-A.hmm.dat' % pfam_dir
        tsv = '%s/Pfam-A.hmm.dat.tsv' % pfam_dir
        if get_pfam_file(hmm) and get_pfam_file(dat) and get_pfam_file(fas):
            self.hmms_pd = get_hmm_dat(dat, tsv)
            if self.config.databases['pfam']['terms']:
                for t in self.config.databases['pfam']['terms']:
                    term_pd = self.hmms_pd.loc[
                        self.hmms_pd['DE'].str.lower().str.contains(t)].copy()
                    if not term_pd.shape[0]:
                        continue
                    pfams, cmd = get_pfams_cmd(hmm, term_pd, pfam_dir, t)
                    self.hmms_dias[t] = pfams
                    self.cmds[t] = cmd
        self.register_command()

    def get_dbcan_hmms(self) -> None:
        """Get all the .hmm files from the dbCAN database."""
        for root, _, files in os.walk(self.config.databases['dbcan']['path']):
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
        path = self.config.databases['dbcan']['path']
        cmd = ""
        fas_fp = '%s/%s.fa' % (folder, name)
        dia_fp = '%s.dmnd' % splitext(fas_fp)[0]
        if not isfile(dia_fp):
            with open(fas_fp, 'w') as o:
                taxon_found = False
                for taxon in taxa:
                    meta_taxon_pd = self.dbcan_meta.loc[
                        self.dbcan_meta.genome_name.str.contains(taxon), :]
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
        return cmd

    def set_dbcan_taxa(self) -> None:
        for name, fp in self.config.databases['dbcan'].get('taxa', {}).items():
            taxa = list(reads_lines(fp))
            if not taxa:
                continue
            folder = '%s/subset' % self.config.databases['dbcan']['path']
            mkdr(folder)
            self.cazys[name] = folder
            self.write_dbcan_subset(taxa, folder, name)

    def set_dbcan(self) -> None:
        self.database = 'dbcan'
        self.get_dbcan_hmms()
        m = '%s/dbCAN-seq/metadata.txt' % self.config.databases['dbcan']['path']
        if isfile(m):
            self.dbcan_meta = pd.read_csv(m, header=0, index_col=0, sep='\t')
            self.set_dbcan_taxa()

    def set_wol(self) -> None:
        self.database = 'wol'
        wol_dir = self.config.databases['wol']['path']
        self.wol['path'] = wol_dir

        fna_dir = '%s/fna' % wol_dir
        mkdr(fna_dir)
        if not glob.glob('%s/*' % fna_dir):
            metadata_fp = '%s/wol/metadata.tsv' % RESOURCES
            batch_down_fp = '%s/wol/batch_down.sh' % RESOURCES
            make_down_list_py = '%s/wol/make_down_list.py' % RESOURCES
            cmd = 'cd %s\n' % fna_dir
            cmd += '%s %s > download.list\n' % (make_down_list_py, metadata_fp)
            cmd += 'bash %s download.list' % batch_down_fp
            self.cmds['download_wol_genomes'] = cmd
        else:
            self.wol['fna'] = glob.glob('%s/*' % fna_dir)

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

    def set_midas(self) -> None:
        midas_db_path = self.config.databases['midas']['path']
        self.midas['path'] = midas_db_path
        self.config.midas_foci['all'] = (midas_db_path, '')

    def set_humann(self) -> None:
        self.humann['path'] = self.config.databases['humann']['path']

    def set_metaphlan(self) -> None:
        self.metaphlan['path'] = self.config.databases['metaphlan']['path']

    def set_macsyfinder(self) -> None:
        self.macsyfinder['path'] = self.config.databases['macsyfinder']['path']

    def set_diamond(self) -> None:
        paths = self.config.databases['diamond']['path']
        self.diamond['uniref'] = paths['uniref']
        self.diamond['genome'] = paths['genome']

    def set_ioncom(self) -> None:
        self.ioncom['path'] = self.config.databases['ioncom']['path']
        self.ioncom['itasser'] = self.config.databases['ioncom']['itasser']

    def set_shogun(self) -> None:
        shogun_config = self.config.databases['shogun']['path']
        for k, v in shogun_config.items():
            self.shogun[k] = v
            if k == 'rep82':
                if not isdir(v) or not glob.glob('%s/*' % v):
                    cmd = 'mkdir -p %s\ncd %s\n' % (v, v)
                    cmd += 'wget -i https://raw.githubusercontent.com/knights-' \
                           'lab/SHOGUN/master/docs/shogun_db_links.txt\n'
                    self.cmds['download_shogun'] = cmd

    def register_command(self) -> None:
        self.databases_commands[self.database] = dict(self.cmds)
        self.cmds = {}


