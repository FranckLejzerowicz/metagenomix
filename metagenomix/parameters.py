# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import numpy as np
from os.path import dirname, isdir, isfile

from metagenomix._io_utils import read_yaml
from metagenomix.tools.alignment import *


def check_int(param: str, value, name):
    """Verifies that the value of the current parameter is an integer."""
    if not str(value).isdigit():
        raise sys.exit('Param "%s" for "%s" must be integer' % (param, name))


def check_mems(param, value, name):
    """Verifies that the requested "memory" is in kb, mb or gb."""
    mems = ['kb', 'mb', 'gb']
    if value not in mems:
        raise sys.exit('Param "%s" for "%s" must be of %s' % (
            param, name, str(mems)))


def check_env(config, value, name):
    """Verifies that the conda environment or module to load do exist."""
    if name not in config.modules and value not in config.conda_envs:
        sys.exit('"%s" not a module or conda env' % value)


def check_path(config, value, name):
    """Verifies that a path exists (usually for databases)."""
    if not config.dev and not isfile(value) and not isdir(value):
        sys.exit('path "%s" for "%s" do not exist' % (value, name))


def check_scratch(value, name):
    """Verifies that the scratch location is specified correctly."""
    if value not in ['scratch', 'userscratch'] and not isinstance(value, int):
        sys.exit('"%s" for "%s" must be "scratch", "userscratch" or an int' % (
            value, name))


def show_valid_params(param, values, name):
    """Verifies that the parameter value is authorized."""
    m = '[%s] user parameter "%s" must be among:\n' % (name, param)
    for value in values:
        m += '  - %s\n' % value
    sys.exit(m)


def check_nums(params, defaults, vals, dtype, tool, mi=None, ma=None):
    for param in vals:
        if param not in params:
            params[param] = defaults[param]
        else:
            if dtype == int:
                if not isinstance(params[param], dtype):
                    sys.exit('[%s] Param "%s" not of %s' % (tool, param, dtype))
            else:
                try:
                    float(params[param])
                except ValueError:
                    sys.exit('[%s] Param "%s" not a number' % (tool, param))
        if mi and ma and not mi <= dtype(params[param]) <= ma:
            sys.exit('[%s] Param "%s" not in [%s-%s]' % (tool, param, mi, ma))


def check_databases(name, params, databases):
    """Verifies that the databases are set in the user parameters.
        - If only one database exists for a given list of databases,
    then the non-existing databases are ignored and the tool will use
    only this one database.
        - If no database is found to exist at all, then metagenomix
    will exit with a helpful message."""
    if 'databases' not in params:
        sys.exit('[%s] "databases" must be a parameter' % name)
    dbs_existing = []
    dbs_missing = []
    for db in params['databases']:
        if db in databases.paths or db == 'default':
            dbs_existing.append(db)
        else:
            dbs_missing.append(db)
    if not dbs_existing:
        sys.exit('[%s] No databases: %s' % (name, '; '.join(dbs_missing)))
    elif dbs_missing:
        print('[%s] Missing databases: "%s"' % (name, '", "'.join(dbs_missing)))
    return dbs_existing


def check_default(params, defaults, name, let_go: list = [], multi: list = []):
    """Verifies that the parameters given by the used for the current tool
    and that are to be checked (i.e., not in the `let_go` list) are valid.
    This mean that they must be one (or sevreral) of the values matching the
    passed defaults. If a default parameter is not given by the user,
    then the parameter is created and set to the first as value of the
    passed defaults.
    Note: a list of parameters for which the value has to be a list (such as
          "aligners" or "models") will have the default value set to a list.
    """
    for param, values in defaults.items():
        if param in let_go:  # skip params checked specifically per software
            continue
        if param not in params:
            if param in ['aligners', 'blobology', 'reassembly', 'models']:
                params[param] = [values[0]]
            else:
                params[param] = values[0]
        else:
            if isinstance(params[param], list):
                if set(sorted(params[param])).difference(values):
                    show_valid_params(param, values, name)
            elif params[param] not in values:
                show_valid_params(param, values, name)
        if param in multi:
            if not isinstance(params[param], list):
                sys.exit('[%s] Param "%s" must be a list' % (name, param))
        elif isinstance(params[param], (list, dict)):
            sys.exit('[%s] Param "%s" must not be a list/dict' % (name, param))


class Parameters(object):
    """
    This class shall contains the checks for all tools specifically,
    but since these rely on the systematic elaboration of `defaults` and
    `dtypes` data structures that are useful later, it is best to consign
    these inside class instances.
    """

    def __init__(self, soft):
        self.name = soft.name
        self.defaults = {}
        self.dtypes = {}
        self.let_go = {}
        self.params = {}

    def check_default(self):
        """Verifies that the parameters given by the used for the current tool
        and that are to be checked (i.e., not in the `let_go` list) are valid.
        This mean that they must be one (or sevreral) of the values matching the
        passed defaults. If a default parameter is not given by the user,
        then the parameter is created and set to the first as value of the
        passed defaults.
        Note: a list of parameters for which the value has to be a list (such as
              "aligners" or "models") will have the default value set to a list.
        """
        for param, values in self.defaults.items():
            if param in self.let_go:  # skip params checked per software
                continue
            if param not in self.params:
                if param in ['aligners', 'blobology', 'reassembly',
                             'models']:
                    self.params[param] = [values[0]]
                else:
                    self.params[param] = values[0]
            else:
                if isinstance(self.params[param], list):
                    if set(sorted(self.params[param])).difference(values):
                        self.show_valid_params(param, values)
                elif self.params[param] not in values:
                    self.show_valid_params(param, values)

    def show_valid_params(self, param, values):
        """Verifies that the parameter value is authorized."""
        m = '[%s] user parameter "%s" must be among:\n' % (self.name, param)
        for value in values:
            m += '  - %s\n' % value
        sys.exit(m)

    def check_nums(self, vals, dtype, tool, mi=None, ma=None):
        for param in vals:
            if param not in self.params:
                self.params[param] = self.defaults[param]
            else:
                if dtype == int:
                    if not isinstance(self.params[param], dtype):
                        sys.exit('[%s] Param "%s" not of %s' % (
                            tool, param, dtype))
                else:
                    try:
                        float(self.params[param])
                    except ValueError:
                        sys.exit('[%s] Param "%s" not a number' % (
                            tool, param))
            if mi and ma and not mi <= dtype(self.params[param]) <= ma:
                sys.exit('[%s] Param "%s" not in [%s-%s]' % (
                    tool, param, mi, ma))

    def check_databases(self, databases):
        """Verifies that the databases are set in the user parameters.
            - If only one database exists for a given list of databases,
        then the non-existing databases are ignored and the tool will use
        only this one database.
            - If no database is found to exist at all, then metagenomix
        will exit with a helpful message."""
        if 'databases' not in self.params:
            sys.exit('[%s] "databases" must be a parameter' % self.name)
        dbs_existing = []
        dbs_missing = []
        for db in self.params['databases']:
            if db in databases.paths or db == 'default':
                dbs_existing.append(db)
            else:
                dbs_missing.append(db)
        placeholders = (self.name, '; '.join(dbs_missing))
        if not dbs_existing:
            sys.exit('[%s] No databases: %s' % placeholders)
        elif dbs_missing:
            print('[%s] Missing databases: "%s"' % placeholders)
        return dbs_existing


# ============================================= #
#  Below are the function that are called from  #
#  class method Workflow.set_params()           #
#       (see module metagenomix.pipeline)       #
# ============================================= #


def get_diamond_hmmer_databases(self, tool, params):
    """Collect the paths of the passed databases that have .dmnd files within
    a diamond folder. This function replaces that value of the 'databases'
    key, from a list of database names to a dict with as keys those names
    that are valid, i.e., for which the value can be the list of .dmnd files."""
    if tool == 'diamond':
        ext = '.dmnd'
    elif tool == 'hmmer':
        ext = '.hmm'
    else:
        sys.exit('Searching only possibly using diamond or hmmsearch')
    valid_dbs = {}
    dbs_existing = check_databases(tool, params, self.databases)
    for db in dbs_existing:
        if self.databases.builds[db].get(tool):
            files = '%s/*%s' % (self.databases.builds[db][tool], ext)
            if self.config.dev:
                valid_dbs[db] = [files]
            else:
                file_paths = glob.glob(files)
                if file_paths:
                    valid_dbs[db] = file_paths
    params['databases'] = valid_dbs


def get_hmmer_databases(self, params):
    if 'terms' in params:
        terms = params['terms']
        if not self.databases.hmms_pd.shape[0]:
            print('[search] Empty Pfam data file (hmmer step ignored)')
        else:
            terms_hmms_dias = self.databases.get_pfam_terms(params)
            not_found_terms = set(terms).difference(set(terms_hmms_dias))
            if not_found_terms:
                print('[search] No HMM for hmmer "terms"')
                for term in not_found_terms:
                    print('[search]    - %s' % term)
                if len(not_found_terms) == len(set(terms)):
                    print('[search]    (i.e. all terms, hmmer step ignored)')
    else:
        print('[search] Params "hmmer:terms" missing (ignored)')


def check_search(self, params, soft):
    tool = soft.name.rsplit('_')[1]
    defaults = {}
    defaults.update(search_tool(params, soft, tool))
    get_diamond_hmmer_databases(self, tool, params)
    defaults['databases'] = '<list of databases>'
    return defaults


def check_ccmap(self, params, soft):
    defaults = {
        'preserve_tmpdir': [False, True], 'terminal_fragment_size': 500,
        'min_percent_identity': 94, 'min_aligned_length': 50}
    check_nums(params, defaults,
               ['terminal_fragment_size', 'min_aligned_length'], int, soft.name)
    check_nums(params, defaults, ['min_percent_identity'],
               int, soft.name, 0, 100)
    let_go = [
        'terminal_fragment_size', 'min_aligned_length', 'min_percent_identity']
    check_default(params, defaults, soft.name, let_go)
    return defaults


def check_barrnap(self, params, soft):
    defaults = {
        'evalue': 1e-06, 'reject': 0.25, 'lencutoff': 0.8,
        'incseq': [False, True], 'kingdom': ['bac', 'mito', 'arc', 'euk']}
    floats = ['evalue', 'reject', 'lencutoff']
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, floats)
    return defaults


def check_integron_finder(self, params, soft):
    defaults = {
        'evalue_attc': 1, 'min_length': 1500,
        'min_attc_size': 40, 'max_attc_size': 200,
        'pdf': [True, False], 'gbk': [True, False], 'mute': [True, False],
        'local_max': [True, False], 'promoter_attI': [True, False],
        'union_integrases': [False, True], 'topology': ['linear', 'circ'],
    }
    check_nums(params, defaults, ['min_length', 'min_attc_size',
                                  'max_attc_size'], int, soft.name)
    check_nums(params, defaults, ['evalue_attc'], float, soft.name, 0, 100)
    let_go = ['min_length', 'min_attc_size', 'max_attc_size', 'evalue_attc']
    check_default(params, defaults, soft.name, let_go)
    defaults['prot_file'] = '<Path to proteins file used for annotations>'
    defaults['attc_model'] = '<Path to the attc model (Covariance Matrix)>'
    defaults['topology_file'] = '<Path to file with each replicon topology>'
    return defaults


def check_prokka_config(params):
    t = params['config']
    cols = ['genus', 'species', 'strain', 'plasmid', 'proteins']
    with open(t) as f:
        for ldx, line in enumerate(f):
            if ldx == 10:  # only check 10 rows
                break
            ls = line.strip('\n').split('\t')
            if not ldx:
                if ls != cols:
                    sys.exit('[prokka] Wrong headers in "config" param file '
                             '("%s"): must be "%s"' % (t, '\\t'.join(cols)))
            if len(ls) != 5:
                sys.exit('[prokka] "config" param file ("%s") do not have 5 '
                         'tab-separated fields (could be empty):\t%s' % (t, ls))


def check_prokka(self, params, soft):
    defaults = {
        'kingdom': ['Archaea', 'Bacteria', 'Mitochondria', 'Viruses'],
        'evalue': 1e-09, 'coverage': 80, 'mincontiglen': 100,
        'gcode': [11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        'notrna': [False, True], 'norrna': [False, True],
        'metagenome': [True, False]
    }
    check_nums(params, defaults, ['evalue'], float, soft.name, 0, 100)
    check_nums(params, defaults, ['coverage'], int, soft.name, 0, 100)
    check_nums(params, defaults, ['mincontiglen'], int, soft.name)
    check_default(params, defaults, soft.name, ['evalue', 'coverage',
                                                'mincontiglen'])
    if 'config' in params:
        if not isfile(params['config']):
            sys.exit('[prokka] Param "config" must be an existing file')
        check_prokka_config(params)
    else:
        params['config'] = None
    defaults['config'] = '<Path to config file for prokka>'
    return defaults


def check_antismash(self, params, soft):
    defaults = {
        'tta-threshold': 0.65,
        'rre': [True, False],
        'asf': [True, False],
        'cassis': [True, False],
        'pfam2go': [True, False],
        'tigrfam': [True, False],
        'cc-mibig': [True, False],
        'fullhmmer': [True, False],
        'cb_general': [True, False],
        'smcog_trees': [True, False],
        'clusterhmmer': [True, False],
        'cb_subclusters': [True, False],
        'cb_knownclusters': [True, False],
        'taxon': ['bacteria', 'fungi'],
        'genefinding_tool': [
            'prodigal-m', 'prodigal', 'glimmerhmm', 'none', 'error']
    }
    check_nums(params, defaults, ['tta-threshold'], float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, ['tta-threshold'])
    return defaults


def check_quast(self, params, soft):
    defaults = {
        'min_contig': 1000,
        'min_alignment': 65,
        'min_identity': 95.0,
        'circos': [False, True],
        'glimmer': [False, True],
        'rna_finding': [False, True],
        'conserved_genes_finding': [False, True],
        'space_efficient': [False, True],
        'ambiguity_usage': ['all', 'none', 'one'],
        'fast': [True, False],  # activates all below
        'no_check': [True, False],
        'no_plots': [True, False],
        'no_html': [True, False],
        'no_icarus': [True, False],
        'no_snps': [True, False],
        'no_gc': [True, False],
        'no_sv': [True, False],
        'no_read_stats': [True, False],
    }
    ints = ['min_contig', 'min_alignment']
    floats = ['min_identity']
    check_nums(params, defaults, ints, int, soft.name)
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    if 'label' in params:
        label = params['label']
        if label not in set(self.config.meta):
            sys.exit('[quast] Params "label" must be a valid metadata column')
    check_default(params, defaults, soft.name, (ints + floats + ['label']))
    defaults['label'] = '<an existing metadata column>'
    return defaults


def check_filtering(self, params, soft):
    db = 'databases'
    if db not in params or not isinstance(params[db], dict):
        sys.exit('[filtering] Param "%s" not a name:path bowtie2 db dict' % db)
    else:
        for (d, index) in params[db].items():
            if not self.config.dev and not glob.glob('%s.*' % index):
                sys.exit('[filtering] Param "%s" no bowtie2 file %s*' % (db, d))
    return {'databases': '<list of databases>'}


def check_atropos(self, params, soft):
    defaults = {
        'q': 15,
        'overlap': 3,
        'max_reads': 0,
        'indel_cost': 1,
        'error_rate': 0.1,
        'nextseq_trim': 5,
        'minimum_length': 50,
        'quality_cutoff': 15,
        'pair_filter': ['any', 'both'],
        'aligner': ['adapter', 'insert'],
        'stats': ['none', 'pre', 'post', 'both'],
        'report_formats': ['txt', 'json', 'yaml', 'pickle'],
        'a': ['GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 'GGGGGGGGGG'],
        'A': ['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', 'GGGGGGGGGG']
    }
    ints = ['q', 'overlap', 'max_reads', 'indel_cost', 'nextseq_trim',
            'minimum_length', 'quality_cutoff']
    floats = ['error_rate']
    check_nums(params, defaults, ints, int, soft.name)
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, (floats + ints + ['a', 'A']))
    for aA in ['a', 'A']:
        if aA not in params:
            params[aA] = []
        if not isinstance(params[aA], list):
            sys.exit('[atropos] Param "%s" must be a list' % aA)
    return defaults


def check_kneaddata(self, params, soft):
    defaults = {'trimmomatic': None, 'bowtie2': None, 'databases': [],
                'purge': [True, False]}
    tools = ['trimmomatic', 'bowtie2']
    for tool in tools:
        if tool not in params:
            params[tool] = defaults[tool]
        elif not self.config.dev and not glob.glob('%s*' % params[tool]):
            print('[kneaddata] Param "%s": no path found (use $PATH)' % tool)
    if 'databases' not in params or not isinstance(params['databases'], list):
        sys.exit('[kneaddata] Params "databases" must a list of existing paths')
    check_default(params, defaults, soft.name, (tools + ['databases']))
    defaults['databases'] = '<list of databases>'
    defaults['bowtie2'] = '<path to bowtie2 (default to $PATH)>'
    defaults['trimmomatic'] = '<path to trimmomatic (default to $PATH)>'
    return defaults


def del_search_tool_keys(params):
    for name in ['diamond', 'hmmer']:
        if name in params:
            del params[name]


def expand_search_params(params, defaults, name):
    if name not in params:
        for k, v in defaults.items():
            if isinstance(v, list):
                params[k] = v[0]
            else:
                params[k] = v
    else:
        for k, v in params[name].items():
            params[k] = v
    del_search_tool_keys(params)


def search_diamond(params, soft) -> dict:
    defaults = {
        'mode': [False, 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive',
                 'very-sensitive', 'ultra-sensitive'],
        'strand': ['both', 'minus', 'plus'],
        'masking': ['tantan', 'none', 'seg'],
        'evalue': 0.001,
        'top': 0,
        'max_hsps': 1,
        'max_target_seqs': 1,
        'id': 80,
        'query_cover': 80,
        'subject_cover': 0
    }
    expand_search_params(params, defaults, 'diamond')
    floats, ints_ = ['evalue'], ['max_target_seqs', 'top', 'max_hsps']
    ints = ['id', 'query_cover', 'subject_cover']
    check_nums(params, defaults, ints_, int, soft.name)
    check_nums(params, defaults, ints, int, soft.name, 0, 100)
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, (ints_ + ints + floats))
    return defaults


def search_hmmer(params, soft) -> dict:
    defaults = {
        'cut_ga': [False, True],
        'cut_nc': [False, True],
        'cut_tc': [False, True],
        'noali': [True, False],
        'nobias': [False, True],
        'max': [False, True],
        'textw': 120,
        'domZ': 10,
        'domE': 1,
        'Z': 10,
        'E': 1,
        'F1': 0.02,
        'F2': 1e-3,
        'F3': 1e-5,
    }
    expand_search_params(params, defaults, 'hmmer')
    e_vals, z_vals = ['E', 'domE'], ['Z', 'domZ', 'textw']
    floats = ['F1', 'F2', 'F3']
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    check_nums(params, defaults, e_vals, float, soft.name, 0, 100)
    check_nums(params, defaults, z_vals, int, soft.name, 0, 1000)
    check_default(params, defaults, soft.name,
                  (e_vals + z_vals + floats + ['terms']))
    defaults['terms'] = '<List of strings to search in Pfam data (for hmmer)>'
    return defaults


def search_tool(params, soft, tool) -> dict:
    if tool == 'diamond':
        defaults = search_diamond(params, soft)
    else:
        defaults = search_hmmer(params, soft)
    return defaults


def check_metaxa2(self, params, soft) -> dict:
    defaults = {'databases': ['default'],
                'align': ['none', 'all', 'uncertain'],
                'mode': ['metagenome', 'genome'],
                'megablast': [True, False],
                'graphical': [True, False],
                'reltax': [True, False],
                'plus': [True, False],
                'T': '0,60,70,75,85,90,97',
                'E': 1, 'S': 12, 'N': 2, 'M': 5, 'R': 75,
                'r': 0.8, 'd': 0.7, 'l': 50}
    if 'T' not in params:
        params['T'] = defaults['T']
    elif len([x for x in str(params['T']).split(',') if x.isdigit()]) != 7:
        sys.exit('[metaxa2] Param "T" must be 7 tab-separated integers [0-100]')
    check_databases('metaxa2', params, self.databases)
    check_nums(params, defaults, ['r', 'd', 'E'], float, soft.name, 0, 1)
    check_nums(params, defaults, ['R'], int, soft.name, 0, 100)
    check_nums(params, defaults, ['S', 'N', 'M', 'l'], int, soft.name)
    check_default(params, defaults, soft.name, ['databases', 'r', 'd', 'E', 'R',
                                                'S', 'N', 'M', 'l', 'T'])
    defaults['databases'].append('<list of databases>')
    return defaults


def check_count(self, params, soft):
    defaults = {'cat': ['zcat', 'gzcat']}
    check_default(params, defaults, soft.name)
    return defaults


def check_kraken2(self, params, soft):
    defaults = {'databases': ['default'], 'confidence': 0.5}
    check_databases('bowtie2', params, self.databases)
    check_nums(params, defaults, ['confidence'], float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, ['confidence', 'databases'])
    defaults['databases'].append('<list of databases>')
    return defaults


def check_shogun(self, params, soft):
    defaults = {'aligners': ['bowtie2', 'burst', 'utree']}
    check_default(params, defaults, soft.name, [], ['aligners'])
    if 1:
        valid_dbs = {}
        dbs_existing = check_databases('shogun', params, self.databases)
        for db in dbs_existing:
            path = self.databases.paths[db]

            yamls = []
            for folder in ['', 'databases/']:
                yaml = '%s/%sshogun/metadata.yaml' % (path, folder)
                yamls.append(yaml)
                if isfile(yaml):
                    break
            else:
                if not self.config.dev:
                    raise IOError(
                        '[shogun] file must exist: %s' % ' or '.join(yamls))
            metadata = read_yaml(yaml)
            for aligner in list(params['aligners']):
                if aligner in metadata:
                    if not self.config.dev:
                        ali = metadata[aligner]
                        if ali[0] == '/':
                            formatted = ali
                        else:
                            formatted = '%s/%s' % (dirname(yaml), ali)
                        if glob.glob('%s*' % formatted):
                            valid_dbs.setdefault(db, []).append(aligner)
                        else:
                            print('[shogun] No formatted db "%s"' % formatted)
                    else:
                        valid_dbs.setdefault(db, []).append(aligner)
        params['databases'] = valid_dbs
        if not params['databases']:
            print('[shogun] No database formatted for shogun: will be skipped')
    defaults['databases'] = '<list of databases>'
    return defaults


def check_bowtie2(self, params, soft):
    defaults = {
        # default: 'paired'
        'pairing': ['paired', 'concat', 'single'],
        'discordant': [True, False],
        'k': '16', 'np': '1',
        'mp': '1,1', 'rdg': '0,1', 'rfg': '0,1',
        'score-min': 'L,0,-0.05'
    }
    let_go = []
    check_bowtie_k_np(soft, params, defaults, let_go)
    check_bowtie_mp_rdg_rfg(soft, params, defaults, let_go)
    check_bowtie_score_min(soft, params, defaults, let_go)
    check_default(params, defaults, soft.name, let_go)
    dbs_existing = check_databases('bowtie2', params, self.databases)
    valid_dbs = {}
    for db in dbs_existing:
        if 'bowtie2' in self.databases.builds[db]:
            bt2_path = '%s/*.*.bt2*' % self.databases.builds[db]['bowtie2']
            if not self.config.dev:
                bt2_paths = glob.glob(bt2_path)
                if bt2_paths:
                    valid_dbs[db] = bt2_paths[0].rsplit('.', 2)[0]
            else:
                valid_dbs[db] = bt2_path.rsplit('.', 2)[0]
    params['databases'] = valid_dbs
    defaults['databases'] = '<list of databases>'
    return defaults


def check_spades(self, params, soft):
    defaults = {'k': ['33', '55', '77', '99', '127'],
                'bio': [False, True],
                'meta': [True, False],
                'plasmid': [False, True],
                'only_assembler': [True, False]}
    if 'k' not in params:
        params['k'] = defaults['k']
    else:
        kerrors = [x for x in params['k'] if not str(x).isdigit()]
        if len(kerrors):
            sys.exit('[spades] "k" must be integers (%s)' % ','.join(kerrors))
    check_default(params, defaults, soft.name, ['k'])
    return defaults


def check_viralverify(self, params, soft):
    if 'path' not in params:
        sys.exit("[viralverify] Please provide path to software's 'bin' folder")
    if not self.config.dev and not isdir(params['path']):
        sys.exit("[viralverify] Please provide path to software's 'bin' folder")
    defaults = {'thr': 7, 'p': [False, True], 'db': [False, True]}
    check_nums(params, defaults, ['thr'], int, 'viralverify')
    check_default(params, defaults, soft.name, ['thr'])
    defaults['path'] = '<Path to the software "bin" folder>'
    return defaults


def simka_error(key: str, kmer_reads: dict) -> None:
    """Throw an error is the user-defined simka params are not right.

    Parameters
    ----------
    key : str
        "" or ""
    kmer_reads : dict
        Params for sampling of the kmers of reads size in Simka
    """
    for p in ['start', 'end', 'size']:
        if p not in kmer_reads or not isinstance(kmer_reads[p], int):
            sys.exit('Please give an integer for simka subparam "%s" (param '
                     '"%s") in your input .yml params file' % (p, key))


def check_simka(self, params, soft):
    defaults = {
        'simkaMin': [True, False],
        'kmer': np.linspace(15, 80, 6),
        'log_reads': np.logspace(3, 7, 3),
        'nb_kmers': 50000,
        'min_read_size': 100,
    }
    if 'kmer' not in params:
        params['kmer'] = defaults['kmer']
    else:
        kmer = params['kmer']
        simka_error('kmer', kmer)
        params['kmer'] = np.linspace(kmer['start'], kmer['end'], kmer['size'])
    if 'log_reads' not in params:
        params['log_reads'] = defaults['log_reads']
    else:
        n = params['log_reads']
        simka_error('log_reads', n)
        params['log_reads'] = np.logspace(n['start'], n['end'], n['size'])
    check_nums(params, defaults, ['nb_kmers', 'min_read_size'], int, soft.name)
    check_default(params, defaults, soft.name,
                  ['kmer', 'log_reads', 'nb_kmers', 'min_read_size'])
    params['kmer'] = [int(x) for x in params['kmer']]
    params['log_reads'] = [int(x) for x in params['log_reads']]
    defaults['path'] = '<Path to the SimKa program folder>'
    defaults['kmer'] = [str(x) for x in np.linspace(15, 80, 6)]
    defaults['log_reads'] = [str(x) for x in np.logspace(3, 7, 3)]
    return defaults


def check_metamarker(self, params, soft):
    defaults = {'identity': 0.9}
    check_nums(params, defaults, ['identity'], float, 'metamarker', 0, 1)
    if 'groups' not in params:
        sys.exit('[metamarker] Param "groups" not found (metadata variable(s))')
    for group in params['groups']:
        if group not in set(self.config.meta.columns):
            sys.exit('[metamarker] metadata variable "%s" not found' % group)
    defaults['groups'] = '<List of metadata columns>'
    return defaults


def check_metawrap(self, params, soft):
    defaults = {
        'binners': ['maxbin2', 'metabat2', 'concoct'],
        'min_completion': 25, 'min_contamination': 5, 'min_contig_length': 1000,
        'min_completion_reassembly': 25, 'min_contamination_reassembly': 5,
        'reassembly': ['permissive', 'strict'],
        'blobology': ['coassembly', 'sample']
    }
    if 'binners' not in params:
        params['binners'] = defaults['binners']
    mins = ['min_completion', 'min_completion_reassembly',
            'min_contamination', 'min_contamination_reassembly']
    check_nums(params, defaults, mins, int, 'metawrap:binning', 0, 100)
    check_nums(params, defaults, ['min_contig_length'], int, 'metawrap:binning')
    mins.append('min_contig_length')
    check_default(params, defaults, soft.name, mins,
                  ['binners', 'reassembly', 'blobology'])
    return defaults


def check_drep(self, params, soft):
    defaults = {
        'n_PRESET': ['normal', 'tight'],
        'coverage_method': ['larger', 'total'],
        'S_algorithm': ['fastANI', 'ANIn', 'ANImf', 'gANI', 'goANI'],
        'clusterAlg': ['average', 'weighted', 'single', 'median',
                       'ward', 'centroid', 'complete'],
        'SkipMash': [False, True],
        'SkipSecondary': [False, True],
        'run_tertiary_clustering': [False, True],
        'greedy_secondary_clustering': [False, True],
        'multiround_primary_clustering': [False, True],
        'primary_chunksize': 5000,
        'MASH_sketch': 1000,
        'P_ani': 0.9,
        'S_ani': 0.95,
        'cov_thresh': 0.1,
        'warn_dist': 0.25,
        'warn_sim': 0.98,
        'warn_aln': 0.25
    }
    if 'anicalculator' not in params:
        sys.exit('[checkm] Param "anicalculator" needed: path to binary folder')
    if 'S_algorithm' not in params:
        params['S_algorithm'] = ['fastANI', 'ANIn']
    ints = ['MASH_sketch', 'primary_chunksize']
    check_nums(params, defaults, ints, int, soft.name)
    flts = ['P_ani', 'S_ani', 'cov_thresh', 'warn_dist', 'warn_sim', 'warn_aln']
    check_nums(params, defaults, flts, float, soft.name, 0, 1)
    check_default(params, defaults, soft.name, (ints + flts), ['S_algorithm'])
    defaults['anicalculator'] = '<path to folder with "anicalculator" binary>'
    return defaults


def check_checkm(self, params, soft):
    defaults = {
        'min_seq_len': 1500,
        'all_reads': [True, False],
        'min_align': 0.98,
        'max_edit_dist': 0.02,
        'aai_strain': 0.9,
        'length': 0.7,
        'min_qc': 15,
        'e_value': 1e-10,
        'coverage': [False, True],
        'reduced_tree': [False, True],
        'ali': [False, True],
        'nt': [False, True],
        'genes': [False, True],
        'tab_table': [True, False],
        'individual_markers': [False, True],
        'skip_adj_correction': [False, True],
        'skip_pseudogene_correction': [False, True],
        'ignore_thresholds': [False, True],
    }
    if 'data' not in params:
        sys.exit('[checkm] Param "data" needed: path for "checkm data setRoot"')
    ints = ['min_seq_len', 'min_qc']
    check_nums(params, defaults, ints, int, soft.name)
    floats = ['min_align', 'max_edit_dist', 'aai_strain', 'length']
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    floats_ = ['e_value']
    check_nums(params, defaults, floats_, float, soft.name)
    check_default(params, defaults, soft.name, (ints + floats + floats_))
    defaults['data'] = '<path to the CheckM reference data>'
    return defaults


def check_gtdbtk(self, params, soft):
    defaults = {}
    check_nums(params, defaults, [''], int, soft.name)
    check_nums(params, defaults, [''], float, soft.name)
    check_default(params, defaults, soft.name)
    return defaults


def check_prodigal(self, params, soft):
    defaults = {'procedure': ['meta', 'single']}
    check_default(params, defaults, soft.name)
    return defaults


def check_fastp(self, params, soft):
    defaults = {'average_qual': 20,
                'trim_poly_g': [True, False]}
    check_nums(params, defaults, ['average_qual'], int, 'fastp')
    check_default(params, defaults, soft.name, ['average_qual'])
    return defaults


def check_humann(self, params, soft):
    defaults = dict({
        'evalue': 0.1,
        'identity': 90,
        'query_cov': 80,
        'subject_cov': 80,
        'prescreen_threshold': 0.001,
        'uniref': ['uniref90', 'uniref50'],
        'skip_translated': [False, True],
        'nucleotide_db': '/path/to/full_chocophlan.v0.1.1/chocophlan',
        'protein_db': '/path/to/uniref90',
    })
    let_go = ['profiles', 'nucleotide_db', 'protein_db']
    param_dtype = [('evalue', float), ('prescreen_threshold', float),
                   ('identity', int), ('query_cov', int), ('subject_cov', int)]
    for (param, dtype) in param_dtype:
        let_go.append(param)
        if param not in params:
            params[param] = defaults[param]
        else:
            if not isinstance(params[param], dtype):
                sys.exit('[humann] Param "%s" not an int' % param)
        if not 0 <= float(params[param]) <= 100:
            sys.exit('[humann] Param "%s" not in [0-100]' % param)

    if not self.config.dev:
        if 'nucleotide_db' not in params or not isdir(params['nucleotide_db']):
            sys.exit('[humann] Param "nucleotide_db" must be an existing path')
        if 'protein_db' not in params or not isdir(params['protein_db']):
            sys.exit('[humann] Param "protein_db" must be an existing path')
    if 'profiles' in params:
        if not isinstance(params['profiles'], dict):
            sys.exit('[humann] Param "profiles" must be a dict structure')
        for key, value in params['profiles'].items():
            if self.config.dev:
                continue
            if not isinstance(value, str):
                sys.exit('[humann] Param "profiles::%s" must be a string' % key)
            elif not isfile(value):
                sys.exit('[humann] Param "profiles::%s" does not exist' % key)
    check_default(params, defaults, soft.name, let_go)
    defaults['profiles'] = '<dict of taxonomic profiles name: file>'
    return defaults


def check_midas(self, params, soft):
    defaults = {
        'focus': {'all': ''},
        's': ['very-sensitive', 'sensitive', 'very-fast', 'fast'],
        'm': ['local', 'global'],
        'species_cov': 3.0,
        'species_topn': 10,
        'word_size': 28,
        'aln_cov': 0.75,
        'mapid': 94.0,
        'readq': 20,
        'trim': 0,
        'n': 0
    }
    if 'tracking' not in params:
        params['tracking'] = []
    else:
        if not isinstance(params['tracking'], list):
            sys.exit('[midas] Param "tracking" must be a list structure')
        meta_vars = set(self.config.meta.columns)
        params['tracking'] = [x for x in params['tracking'] if x in meta_vars]
        v = [x for x in params['tracking'] if x not in meta_vars]
        if len(v):
            print('[midas] Param "tracking" no metadata var "%s" (ignored)' % v)
        if not params['tracking']:
            print('[midas] Param "tracking" no metadata variable (ignored)')

    if 'focus' not in params:
        params['focus'] = defaults['focus']
    else:
        if not isinstance(params['focus'], dict):
            sys.exit('[midas] Param "focus" must be a dict structure')
        for k, v in params['focus'].items():
            if not isinstance(v, str):
                sys.exit('[midas] Param "focus::%s" must be a char. string' % k)
            if not self.config.dev and not isfile(v):
                sys.exit('[midas] Param "focus::%s::%s" not a file' % (k, v))
    ints = ['species_topn', 'readq', 'trim', 'n', 'word_size']
    check_nums(params, defaults, ints, int, soft.name)
    floats = ['mapid', 'aln_cov']
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    floats_ = ['species_cov']
    check_nums(params, defaults, floats_, float, soft.name)
    check_default(params, defaults, soft.name,
                  (ints + floats + floats_ + ['tracking', 'focus']))
    defaults['focus'] = '<dict of key:value pair some_name: /path/to/spc.txt'
    defaults['tracking'] = '<list of metadata columns>'
    return defaults


def check_macsyfinder(self, params, soft):
    defaults = {
        'db_type': (['unordered', 'ordered_replicon', 'gembase'], str),
        'replicon_topology': (['linear', 'circular'], str),
        'models': (['TXSS', 'TFF-SF', 'CAS'], list),
        'evalue': (0.1, float),
        'coverage': (0.5, float)
    }
    check_nums(params, defaults, ['evalue', 'coverage'], float, soft.name, 0, 1)
    let_go = ['evalue', 'coverage']
    check_default(params, defaults, soft.name, let_go)
    return defaults


def check_coconet(self, params, soft):
    defaults = {
        'quiet': [False, True],
        'no_rc': [False, True],
        'silent': [False, True],
        'continue': [False, True],
        'patience': [False, True],
        'recruit_small_contigs': [False, True],
        'algorithm': ['leiden', 'spectral'],
        'features': ['coverage' 'composition'],
        'tlen_range': [],
        'debug': 20,
        'flag': 3596,
        'min_ctg_len': 2048,
        'min_prevalence': 2,
        'min_mapping_quality': 30,
        'min_aln_coverage': 50,
        'min_dtr_size': 10,
        'fragment_step': 128,
        'test_ratio': 0.1,
        'n_train': 4000000,
        'n_test': 10000,
        'learning_rate': 0.001,
        'batch_size': 256,
        'test_batch': 400,
        'load_batch': 100,
        'compo_neurons': [64, 32],
        'cover_neurons': [64, 32],
        'cover_filters': 16,
        'cover_kernel': 4,
        'cover_stride': 2,
        'merge_neurons': 32,
        'kmer': 4,
        'wsize': 64,
        'wstep': 32,
        'n_frags': 30,
        'max_neighbors': 250,
        'vote_threshold': 0,
        'theta': 0.8,
        'gamma1': 0.3,
        'gamma2': 0.4,
        'n_clusters': 0,
        'fragment_length': -1,
    }
    for p in ['tlen_range', 'compo_neurons', 'cover_neurons']:
        if p not in params:
            params[p] = defaults[p]
        else:
            ints = [int(x) for x in params[p] if x.isdigit()]
            if not isinstance(params[p], list) or len(ints) != 2:
                sys.exit('[coconet] Param "%s" must be a 2-integers list' % p)
            if ints[0] > ints[1]:
                sys.exit('[coconet] Param "%s": integers must be sorted' % p)
    ints = ['debug', 'flag', 'min_ctg_len', 'min_prevalence',
            'min_mapping_quality', 'min_aln_coverage', 'min_dtr_size',
            'fragment_step', 'n_train', 'n_test', 'batch_size', 'test_batch',
            'load_batch', 'cover_filters', 'cover_kernel', 'cover_stride',
            'merge_neurons', 'kmer', 'wsize', 'wstep', 'n_frags',
            'max_neighbors', 'vote_threshold', 'n_clusters', 'fragment_length']
    check_nums(params, defaults, ints, int, soft.name)
    floats = ['test_ratio', 'learning_rate', 'theta', 'gamma1', 'gamma2']
    check_nums(params, defaults, floats, float, soft.name, 0, 1)
    let_go = ['tlen_range', 'compo_neurons', 'cover_neurons']
    check_default(params, defaults, soft.name, (ints + floats + let_go))
    return defaults


def check_tiara(self, params, soft):
    defaults = {
        'min_len': 3000, 'first_stage_kmer': 6, 'second_stage_kmer': 7,
        'to_fasta': ['all', 'mit', 'pla', 'bac', 'arc', 'euk', 'unk', 'pro'],
        'prob_cutoff': [0.65, 0.65],
        'probabilities': [True, False],
        'verbose': [True, False],
        'gzip': [False, True]
    }
    fas = 'to_fasta'
    if fas not in params:
        params[fas] = ['all']
    else:
        if not isinstance(params[fas], list):
            sys.exit('[tiara] Param "%s" not of %s' % (fas, defaults[fas]))
        to_fasta = [x for x in params[fas] if x in defaults[fas]]
        if not to_fasta:
            sys.exit('[tiara] Param "%s" not in %s' % (fas, defaults[fas]))
        if 'all' in to_fasta:
            params['to_fasta'] = ['all']
    prob = 'prob_cutoff'
    if prob not in params:
        params[prob] = defaults[prob]
    else:
        if not isinstance(params[prob], list) or len(params[prob]) > 2:
            sys.exit('[tiara] Param "%s" not a max 2-floats [0-1] list' % prob)
        for c in params[prob]:
            try:
                float(c)
            except ValueError:
                sys.exit('[tiara] Param "%s": %s not [0-1] float' % (prob, c))
    params[prob] = map(str, params[prob])
    ints = ['min_len', 'first_stage_kmer', 'second_stage_kmer']
    check_nums(params, defaults, ints, int, soft.name)
    check_default(params, defaults, soft.name, (ints + [fas, prob]))
    return defaults


def check_bracken(self, params, soft):
    defaults = {
        'read_len': 100,
        'level': ['S', 'D', 'P', 'C', 'O', 'F', 'G', 'S1'],
        'threshold': 0
    }
    check_nums(params, defaults, ['read_len', 'threshold'], int, soft.name)
    check_default(params, defaults, soft.name, ['read_len', 'threshold'])
    return defaults


def check_plasmidfinder(self, params, soft):
    defaults = {
        'extented_output': [True, False],
        'methodPath': ['kma', 'blastn'],
        'threshold': 70,
        'mincov': 75,
    }
    if 'kma' not in params:
        sys.exit('[plasmidfinder] Param "kma" missing (kma binary)')
    elif not self.config.dev and not isfile(params['kma']):
        sys.exit('[plasmidfinder] "kma" binary not found: "%s"' % params['kma'])

    if 'path' not in params:
        sys.exit('[plasmidfinder] Param "path" missing (plasmidfinder.py path)')
    elif not self.config.dev and not isfile(params['path']):
        sys.exit('[plasmidfinder] "%s" not found' % params['path'])

    if 'databasePath' not in params:
        sys.exit('[plasmidfinder] Param "databasePath" missing (databases dir)')
    elif not self.config.dev and not isdir(params['db_path']):
        sys.exit('[plasmidfinder] no "databasePath": %s' % params['db_path'])

    if 'databases' in params:
        if not isinstance(params['databases'], str):
            dbs = params['databases']
            sys.exit('[plasmidfinder] Param "databases" not a string: %s' % dbs)

    ints = ['threshold', 'mincov']
    check_nums(params, defaults, ints, int, soft.name, 0, 100)
    check_default(params, defaults, soft.name, ints)
    defaults['kma'] = '<Path to the "kma" binary>'
    defaults['path'] = '<Path to the "plasmidfinder.py" script>'
    defaults['databasePath'] = '<Path databases folder (with a "config" file)>'
    defaults['databases'] = "<Comma-separated databases (first field of the " \
                            "databasePath's 'config' file)>"
    return defaults


def check_plasforest(self, params, soft):
    defaults = {
        'size_of_batch': 0,
        'b': [True, False],
        'f': [True, False],
        'r': [False, True]
    }
    check_nums(params, defaults, ['size_of_batch'], int, soft.name)
    check_default(params, defaults, soft.name, ['size_of_batch'])
    return defaults


# def check_TOOL(self, params, soft):
#     defaults = {}


# def check_TOOL(self, params, soft):
#     defaults = {}


# def check_TOOL(self, params, soft):
#     defaults = {}


# def check_TOOL(self, params, soft):
#     defaults = {}


# def check_TOOL(self, params, soft):
#     defaults = {}


# def check_TOOL(self, params, soft):
#     defaults = {}
