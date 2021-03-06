# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob
import numpy as np
from os.path import dirname, isdir, isfile

from metagenomix._io_utils import read_yaml


def tech_params(self, tech: str):
    params = {}
    for param, values in self.soft.params.items():
        if isinstance(values, dict) and tech in set(values):
            params[param] = values[tech]
        else:
            params[param] = values
    return params


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


def check_num_sys_exit(tool: str, err: str, tech: str = None):
    if tech:
        sys.exit('[%s] %s (%s)' % (tool, err, tech))
    else:
        sys.exit('[%s] %s' % (tool, err))


def check_nums_generic(param, val, dtype, tool, mi=None, ma=None, tech=None):
    if dtype == int:
        if not isinstance(val, dtype):
            err = 'Param "%s" not of %s' % (param, dtype)
            check_num_sys_exit(tool, err, tech)
    else:
        try:
            float(val)
        except ValueError:
            err = 'Param "%s" not a number' % param
            check_num_sys_exit(tool, err, tech)
    if mi and ma and not mi <= dtype(val) <= ma:
        err = 'Param "%s" not in [%s-%s]' % (param, mi, ma)
        check_num_sys_exit(tool, err, tech)


def check_nums(self, params, defaults, vals, dtype, tool, mi=None, ma=None):
    for param in vals:
        if param not in params:
            params[param] = defaults[param]
        else:
            val = params[param]
            if isinstance(val, dict) and set(val).issubset(self.config.techs):
                for tech, v in val.items():
                    check_nums_generic(param, v, dtype, tool, mi, ma, tech)
            else:
                check_nums_generic(param, val, dtype, tool, mi, ma)


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


def show_valid_params(param, values, name, tech=None):
    """Verifies that the parameter value is authorized."""
    m = '[%s] user parameter "%s" must be among' % (name, param)
    if tech:
        m += ' (%s):' % tech
    m += '\n'
    for value in values:
        m += '  - %s\n' % value
    sys.exit(m)


def check_default_generic(param, vals, values, tool, tech=None):
    if isinstance(vals, list):
        if set(sorted(vals)).difference(values):
            show_valid_params(param, values, tool, tech)
    elif vals not in values:
        show_valid_params(param, values, tool, tech)


def check_type(vals, param, multi, tool, tech=None):
    if param in multi:
        if not isinstance(vals, list):
            check_type_sys_exit(tool, param, False, tech)
    elif isinstance(vals, (list, dict)):
        check_type_sys_exit(tool, param, True, tech)


def check_default_type(self, params, param, multi, tool):
    vals = params[param]
    if isinstance(vals, dict) and set(vals).issubset(self.config.techs):
        for tech, v in vals.items():
            check_type(v, param, multi, tool, tech)
    else:
        check_type(vals, param, multi, tool)


def check_type_sys_exit(tool, param, no_dict, tech=None):
    if no_dict:
        m = '[%s] Param "%s" must not be a list/dict' % (tool, param)
    else:
        m = '[%s] Param "%s" must be a list' % (tool, param)
    if tech:
        m += ' (%s)' % tech
    sys.exit(m)


def check_default(self, params, defaults, tool, let_go: list = [],
                  multi: list = []):
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
            vals = params[param]
            if isinstance(vals, dict) and set(vals).issubset(self.config.techs):
                for tech, v in vals.items():
                    check_default_generic(param, v, values, tool, tech)
            else:
                check_default_generic(param, vals, values, tool)
        check_default_type(self, params, param, multi, tool)


def check_binary(self, t, params, defaults, opt):
    if [x for x in self.config.modules.get(t, []) if x.lower().startswith(t)]:
        return None
    if opt == 'path':
        message = 'Param "path" for path to binaries folder missing'
        defaults['path'] = '<Path to folder containing the %s binary>' % t
        isfile_or_isdir = isdir
    if opt in ['binary', 'trimmomatic', 'bowtie2']:
        message = 'Param "binary" for path to binary/executable'
        defaults['binary'] = '<Path to the %s binary>' % t
        isfile_or_isdir = isfile
    if opt not in params:
        sys.exit('[%s] %s' % (t, message))
    if not self.config.dev and not isfile_or_isdir(params[opt]):
        sys.exit('[%s] Please provide valid path to param "%s"' % (opt, t))


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
    defaults.update(search_tool(self, params, soft, tool))
    get_diamond_hmmer_databases(self, tool, params)
    defaults['databases'] = '<list of databases>'
    return defaults


def check_ccmap(self, params, soft):
    defaults = {
        'preserve_tmpdir': [False, True],
        'terminal_fragment_size': 500,
        'min_percent_identity': 94,
        'min_aligned_length': 50}
    check_nums(self, params, defaults,
               ['terminal_fragment_size', 'min_aligned_length'], int, soft.name)
    check_nums(self, params, defaults, ['min_percent_identity'],
               int, soft.name, 0, 100)
    let_go = [
        'terminal_fragment_size', 'min_aligned_length', 'min_percent_identity']
    check_default(self, params, defaults, soft.name, let_go)
    return defaults


def check_barrnap(self, params, soft):
    defaults = {
        'evalue': 1e-06,
        'reject': 0.25,
        'lencutoff': 0.8,
        'incseq': [False, True],
        'kingdom': ['bac', 'mito', 'arc', 'euk']
    }
    floats = ['evalue', 'reject', 'lencutoff']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, floats)
    return defaults


def check_integron_finder(self, params, soft):
    defaults = {
        'evalue_attc': 1,
        'min_length': 1500,
        'min_attc_size': 40,
        'max_attc_size': 200,
        'pdf': [True, False],
        'gbk': [True, False],
        'mute': [True, False],
        'local_max': [True, False],
        'promoter_attI': [True, False],
        'union_integrases': [False, True],
        'topology': ['linear', 'circ'],
    }
    check_nums(self, params, defaults, ['min_length', 'min_attc_size',
                                        'max_attc_size'], int, soft.name)
    check_nums(self, params, defaults, ['evalue_attc'],
               float, soft.name, 0, 100)
    let_go = ['min_length', 'min_attc_size', 'max_attc_size', 'evalue_attc']
    check_default(self, params, defaults, soft.name, let_go)
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
        'evalue': 1e-09,
        'coverage': 80,
        'mincontiglen': 100,
        'gcode': [11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        'notrna': [False, True],
        'norrna': [False, True],
        'metagenome': [True, False]
    }
    check_nums(self, params, defaults, ['evalue'], float, soft.name, 0, 100)
    check_nums(self, params, defaults, ['coverage'], int, soft.name, 0, 100)
    check_nums(self, params, defaults, ['mincontiglen'], int, soft.name)
    check_default(self, params, defaults, soft.name,
                  ['evalue', 'coverage', 'mincontiglen'])
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
    check_nums(self, params, defaults, ['tta-threshold'],
               float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, ['tta-threshold'])
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
    check_nums(self, params, defaults, ints, int, soft.name)
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    if 'label' in params:
        label = params['label']
        if label not in set(self.config.meta):
            sys.exit('[quast] Params "label" must be a valid metadata column')
    check_default(self, params, defaults, soft.name,
                  (ints + floats + ['label']))
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
    check_nums(self, params, defaults, ints, int, soft.name)
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name,
                  (floats + ints + ['a', 'A']))
    for aA in ['a', 'A']:
        if aA not in params:
            params[aA] = []
        if not isinstance(params[aA], list):
            sys.exit('[atropos] Param "%s" must be a list' % aA)
    return defaults


def check_kneaddata(self, params, soft):
    defaults = {
        'trimmomatic': None,
        'bowtie2': None,
        'purge': [True, False]}
    tools = ['trimmomatic', 'bowtie2']
    for tool in tools:
        if tool not in params:
            params[tool] = defaults[tool]
        else:
            check_binary(self, 'kneaddata', params, defaults, tool)
    if 'databases' not in params or not isinstance(params['databases'], list):
        sys.exit('[kneaddata] Params "databases" must a list of existing paths')
    check_default(self, params, defaults, soft.name, (tools + ['databases']))
    defaults['databases'] = '<list of databases>'
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


def search_diamond(self, params, soft) -> dict:
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
    check_nums(self, params, defaults, ints_, int, soft.name)
    check_nums(self, params, defaults, ints, int, soft.name, 0, 100)
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints_ + ints + floats))
    return defaults


def search_hmmer(self, params, soft) -> dict:
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
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_nums(self, params, defaults, e_vals, float, soft.name, 0, 100)
    check_nums(self, params, defaults, z_vals, int, soft.name, 0, 1000)
    check_default(self, params, defaults, soft.name,
                  (e_vals + z_vals + floats + ['terms']))
    defaults['terms'] = '<List of strings to search in Pfam data (for hmmer)>'
    return defaults


def search_tool(self, params, soft, tool) -> dict:
    if tool == 'diamond':
        defaults = search_diamond(self, params, soft)
    else:
        defaults = search_hmmer(self, params, soft)
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
    check_nums(self, params, defaults, ['r', 'd', 'E'], float, soft.name, 0, 1)
    check_nums(self, params, defaults, ['R'], int, soft.name, 0, 100)
    check_nums(self, params, defaults, ['S', 'N', 'M', 'l'], int, soft.name)
    check_default(self, params, defaults, soft.name,
                  ['databases', 'r', 'd', 'E', 'R', 'S', 'N', 'M', 'l', 'T'])
    defaults['databases'].append('<list of databases>')
    return defaults


def check_count(self, params, soft):
    defaults = {'cat': ['zcat', 'gzcat']}
    check_default(self, params, defaults, soft.name)
    return defaults


def check_kraken2(self, params, soft):
    defaults = {
        'databases': ['default'],
        'confidence': 0.5
    }
    check_databases('kraken2', params, self.databases)
    check_nums(self, params, defaults, ['confidence'], float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name,
                  ['confidence', 'databases'])
    defaults['databases'].append('<list of databases>')
    return defaults


def check_shogun(self, params, soft):
    defaults = {'aligners': ['bowtie2', 'burst', 'utree']}
    check_default(self, params, defaults, soft.name, [], ['aligners'])
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
        'pairing': ['paired', 'concat', 'single'],
        'discordant': [True, False],
        'k': '16', 'np': '1',
        'mp': '1,1', 'rdg': '0,1', 'rfg': '0,1',
        'score-min': 'L,0,-0.05'
    }
    for opt in ['k', 'np']:
        if opt in params and not str(params[opt]).isdigit():
            sys.exit('[bowtie2] "%s" option invalid' % opt)
        else:
            params[opt] = defaults[opt]

    for opt in ['mp', 'rdg', 'rfg']:
        if opt in params:
            if len([x.isdigit() for x in str(params[opt]).split(',')]) != 2:
                sys.exit('[bowtie2] "%s" option invalid' % opt)
        else:
            params[opt] = defaults[opt]

    if 'score-min' in params:
        s = params['score-min'].split(',')
        if len(s) != 3 or s[0] not in 'CLSG':
            sys.exit('[bowtie2] "score-min" option invalid')
        else:
            for r in [1, 2]:
                try:
                    float(s[r])
                except ValueError:
                    sys.exit('[bowtie2] "score-min" option invalid')
    else:
        params['score-min'] = defaults['score-min']

    let_go = ['k', 'np', 'mp', 'rdg', 'rfg', 'score-min']
    check_default(self, params, defaults, soft.name, let_go)
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
    defaults = {
        'hybrid': ['illumina', 'nanopore', 'pacbio'],
        'k': ['33', '55', '77', '99', '127'],
        'checkpoints': ['all', 'last'],
        'cov_cutoff': ['off', 'auto'],
        'meta': [True, False],
        'careful': [False, True],
        'continue': [True, False],
        'disable_rr': [False, True],
        'only_assembler': [True, False],
        'disable_gzip_output': [False, True],
        'only_error_correction': [False, True],
    }
    if 'hybrid' not in params:
        params['hybrid'] = defaults['hybrid']
    if 'k' not in params:
        params['k'] = defaults['k']
    else:
        kerrors = [x for x in params['k'] if not str(x).isdigit()]
        if len(kerrors):
            sys.exit('[spades] "k" must be integers (%s)' % ','.join(kerrors))
    check_default(self, params, defaults, soft.name, ['k'], ['hybrid'])
    return defaults


def check_viralverify(self, params, soft):
    defaults = {'thr': 7, 'p': [False, True], 'db': [False, True]}
    check_nums(self, params, defaults, ['thr'], int, 'viralverify')
    check_default(self, params, defaults, soft.name, ['thr'])
    check_binary(self, soft.name, params, defaults, 'path')
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
    check_nums(self, params, defaults, ['nb_kmers', 'min_read_size'],
               int, soft.name)
    check_default(self, params, defaults, soft.name,
                  ['kmer', 'log_reads', 'nb_kmers', 'min_read_size'])
    params['kmer'] = [int(x) for x in params['kmer']]
    params['log_reads'] = [int(x) for x in params['log_reads']]
    defaults['kmer'] = [str(x) for x in np.linspace(15, 80, 6)]
    defaults['log_reads'] = [str(x) for x in np.logspace(3, 7, 3)]
    check_binary(self, 'kneaddata', params, defaults, 'path')
    return defaults


def check_metamarker(self, params, soft):
    defaults = {'identity': 0.9}
    check_nums(self, params, defaults, ['identity'], float, 'metamarker', 0, 1)
    if 'groups' not in params:
        sys.exit('[metamarker] Param "groups" not found (metadata variable(s))')
    for group in params['groups']:
        if group not in set(self.config.meta.columns):
            sys.exit('[metamarker] metadata variable "%s" not found' % group)
    defaults['groups'] = '<List of metadata columns>'
    return defaults


def check_metawrap(self, params, soft):
    defaults = {
        'min_completion': 25,
        'min_contamination': 5,
        'min_contig_length': 1000,
        'min_completion_reassembly': 25,
        'min_contamination_reassembly': 5,
        'binners': ['maxbin2', 'metabat2', 'concoct'],
        'reassembly': ['permissive', 'strict'],
        'blobology': ['coassembly', 'sample']
    }
    if 'binners' not in params:
        params['binners'] = defaults['binners']
    mins = ['min_completion', 'min_completion_reassembly',
            'min_contamination', 'min_contamination_reassembly']
    check_nums(self, params, defaults, mins, int, 'metawrap:binning', 0, 100)
    check_nums(self, params, defaults, ['min_contig_length'],
               int, 'metawrap:binning')
    mins.append('min_contig_length')
    check_default(self, params, defaults, soft.name, mins,
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
    check_nums(self, params, defaults, ints, int, soft.name)
    flts = ['P_ani', 'S_ani', 'cov_thresh', 'warn_dist', 'warn_sim', 'warn_aln']
    check_nums(self, params, defaults, flts, float, soft.name, 0, 1)
    check_default(
        self, params, defaults, soft.name, (ints + flts), ['S_algorithm'])
    check_binary(self, soft.name, params, defaults, 'path')
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
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['min_align', 'max_edit_dist', 'aai_strain', 'length']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    floats_ = ['e_value']
    check_nums(self, params, defaults, floats_, float, soft.name)
    check_default(self, params, defaults, soft.name, (ints + floats + floats_))
    defaults['data'] = '<path to the CheckM reference data>'
    return defaults


def check_gtdbtk(self, params, soft):
    defaults = {}
    check_nums(self, params, defaults, [''], int, soft.name)
    check_nums(self, params, defaults, [''], float, soft.name)
    check_default(self, params, defaults, soft.name)
    return defaults


def check_prodigal(self, params, soft):
    defaults = {
        'f': ['gbk', 'gff', 'sco'],
        'p': ['meta', 'single'],
        'c': [False, True],
        'm': [False, True],
        'n': [False, True],
        'q': [False, True],
    }
    check_default(self, params, defaults, soft.name)
    return defaults


def check_fastp(self, params, soft):
    defaults = {
        'split': 0,
        'split_prefix_digits': 4.,
        'split_by_lines': [False, True],
        'overrepresentation_analysis': [False, True],
        'overrepresentation_sampling': 20,
        'correction': [False, True],
        'overlap_len_require': 30,
        'overlap_diff_limit': 5,
        'complexity_threshold': 30,
        'filter_by_index_threshold': 0,
        'low_complexity_filter': [False, True],
        'disable_length_filtering': [False, True],
        'length_required': 15,
        'length_limit': 0,
        'phred64': [False, True],
        'compression': 4,
        'reads_to_process': 0,
        'dont_overwrite': [False, True],
        'verbose': [False, True],
        'disable_adapter_trimming': [False, True],
        'detect_adapter_for_pe': [True, False],
        'trim_front1': 0,
        'trim_tail1': 0,
        'max_len1': 0,
        'trim_front2': 0,
        'trim_tail2': 0,
        'max_len2': 0,
        'dedup': [False, True],
        'dup_calc_accuracy': 0,
        'dont_eval_duplication': [False, True],
        'trim_poly_g': [True, False],
        'poly_g_min_len': 10,
        'disable_trim_poly_g': [False, True],
        'trim_poly_x': [False, True],
        'poly_x_min_len': 10,
        'cut_window_size': 4,
        'cut_mean_quality': 20,
        'cut_front': [False, True],
        'cut_front_window_size': 4,
        'cut_front_mean_quality': 20,
        'cut_tail': [False, True],
        'cut_tail_window_size': 4,
        'cut_tail_mean_quality': 20,
        'cut_right': [False, True],
        'cut_right_window_size': 4,
        'cut_right_mean_quality': 20,
        'disable_quality_filtering': [False, True],
        'qualified_quality_phred': 15,
        'unqualified_percent_limit': 40,
        'n_base_limit': 5,
        'average_qual': 20,
    }
    int1 = [
        'overlap_len_require', 'overlap_diff_limit','filter_by_index_threshold',
        'length_required', 'length_limit', 'average_qual', 'n_base_limit',
        'qualified_quality_phred', 'cut_window_size', 'cut_front_window_size',
        'cut_tail_window_size', 'cut_right_window_size', 'poly_g_min_len',
        'poly_x_min_len', 'trim_front1', 'trim_tail1', 'max_len1',
        'trim_front2', 'trim_tail2', 'max_len2', 'reads_to_process']
    for adapter in ['adapter_sequence', 'adapter_sequence_r2', 'adapter_fasta']:
        if adapter not in params:
            params[adapter] = None
    check_nums(self, params, defaults, int1, int, 'fastp')
    int2 = ['split']
    check_nums(self, params, defaults, int2, int, 'fastp', 2, 999)
    int3 = ['split_prefix_digits']
    check_nums(self, params, defaults, int3, int, 'fastp', 1, 10)
    int4 = ['overrepresentation_sampling']
    check_nums(self, params, defaults, int4, int, 'fastp', 1, 10000)
    int5 = ['complexity_threshold', 'unqualified_percent_limit']
    check_nums(self, params, defaults, int5, int, 'fastp', 1, 100)
    int6 = ['cut_mean_quality', 'cut_front_mean_quality',
            'cut_tail_mean_quality', 'cut_right_mean_quality']
    check_nums(self, params, defaults, int6, int, 'fastp', 1, 36)
    int7 = ['dup_calc_accuracy']
    check_nums(self, params, defaults, int7, int, 'fastp', 1, 6)
    int8 = ['compression']
    check_nums(self, params, defaults, int8, int, 'fastp', 1, 9)
    let_go = int1 + int2 + int3 + int4 + int5 + int6 + int7 + int8
    check_default(self, params, defaults, soft.name, let_go)
    defaults['adapter_sequence'] = '<the adapter for read1>'
    defaults['adapter_sequence_r2'] = '<the adapter for read2>'
    defaults['adapter_fasta'] = '<Path to fasta file with sequences to trim>'
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
    check_default(self, params, defaults, soft.name, let_go)
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
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['mapid', 'aln_cov']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    floats_ = ['species_cov']
    check_nums(self, params, defaults, floats_, float, soft.name)
    check_default(self, params, defaults, soft.name,
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
    check_nums(self, params, defaults, ['evalue', 'coverage'],
               float, soft.name, 0, 1)
    let_go = ['evalue', 'coverage']
    check_default(self, params, defaults, soft.name, let_go)
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
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['test_ratio', 'learning_rate', 'theta', 'gamma1', 'gamma2']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    let_go = ['tlen_range', 'compo_neurons', 'cover_neurons']
    check_default(self, params, defaults, soft.name, (ints + floats + let_go))
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
    params[prob] = [str(x) for x in params[prob]]
    ints = ['min_len', 'first_stage_kmer', 'second_stage_kmer']
    check_nums(self, params, defaults, ints, int, soft.name)
    check_default(self, params, defaults, soft.name, (ints + [fas, prob]))
    return defaults


def check_bracken(self, params, soft):
    defaults = {
        'read_len': 100,
        'level': ['S', 'D', 'P', 'C', 'O', 'F', 'G', 'S1'],
        'threshold': 0
    }
    check_nums(self, params, defaults, ['read_len', 'threshold'],
               int, soft.name)
    check_default(self, params, defaults, soft.name, ['read_len', 'threshold'])
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

    db_path = 'databasePath'
    if db_path not in params:
        sys.exit('[plasmidfinder] Param "databasePath" missing (databases dir)')
    elif not self.config.dev and not isdir(params[db_path]):
        sys.exit('[plasmidfinder] no "databasePath" %s' % params[db_path])

    if 'databases' in params:
        if not isinstance(params['databases'], str):
            dbs = params['databases']
            sys.exit('[plasmidfinder] Param "databases" not a string: %s' % dbs)

    ints = ['threshold', 'mincov']
    check_nums(self, params, defaults, ints, int, soft.name, 0, 100)
    check_default(self, params, defaults, soft.name, ints)
    check_binary(self, soft.name, params, defaults, 'binary')
    defaults['kma'] = '<Path to the "kma" binary>'
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
    check_nums(self, params, defaults, ['size_of_batch'], int, soft.name)
    check_default(self, params, defaults, soft.name, ['size_of_batch'])
    check_binary(self, soft.name, params, defaults, 'binary')
    return defaults


def check_flye(self, params, soft):
    defaults = {
        'stop_after': [
            None, 'assembly', 'consensus', 'repeat', 'trestle', 'polishing'],
        'resume_from': [
            None, 'assembly', 'consensus', 'repeat', 'trestle', 'polishing'],
        'read_error': 0.,
        'iterations': 1,
        'min_overlap': 0,
        'asm_coverage': 0,
        'meta': [True, False],
        'resume': [False, True],
        'scaffold': [True, False],
        'polish_target': [False, True],
        'keep_haplotypes': [False, True],
        'long_reads': [
            'pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw',
            'nano-corr', 'nano-hq'],
    }
    ints = ['iterations', 'asm_coverage', 'min_overlap']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['read_error']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + floats))
    defaults['stop_after'] = '<Stage name that is completed and to stop after>'
    defaults['genome_size'] = '<estimated genome size (e.g., 5m or 2.6g)>'
    return defaults


def check_canu(self, params, soft):
    defaults = {
        'stage': [
            None, 'haplotype', 'correct', 'trim', 'assemble', 'trim-assemble'],
        'corMhapSensitivity': ['high', 'normal', 'low'],
        'processing': [None, 'corrected', 'trimmed'],
        'correctedErrorRate': 0.,
        'minOverlapLength': 500,
        'rawErrorRate': 0.,
        'maxInputCoverage': 10000,
        'minReadLength': 1000,
        'corOutCoverage': 10000,
        'corMinCoverage': 0,
        'redMemory': 32,
        'oeaMemory': 32,
        'batMemory': 200,
    }
    ints = ['minReadLength', 'minOverlapLength', 'redMemory', 'oeaMemory',
            'batMemory', 'maxInputCoverage']
    check_nums(self, params, defaults, ints, int, soft.name)
    ints_ = ['corMinCoverage', 'corOutCoverage']
    check_nums(self, params, defaults, ints_, int, soft.name, 0, 4)
    floats = ['correctedErrorRate', 'rawErrorRate']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + ints_ + floats))
    check_binary(self, soft.name, params, defaults, 'path')
    defaults['genome_size'] = '<estimated genome size (e.g., 5m or 2.6g)>'
    defaults['specifications'] = '<Path to assembly option specifications file>'
    return defaults


def check_unicycler(self, params, soft):
    defaults = {
        'hybrid': [True, False],
        'min_bridge_qual': 25.0,
        'start_gene_cov': 95.0,
        'start_gene_id': 90.0,
        'max_kmer_frac': 0.95,
        'min_kmer_frac': 0.2,
        'depth_filter': 0.25,
        'min_component_size': 1000,
        'min_dead_end_size': 1000,
        'min_polish_size': 10000,
        'min_anchor_seg_len': 0,
        'min_fasta_length': 100,
        'kmer_count': 10,
        'linear_seqs': 0,
        'low_score': 0,
        'kmers': 0,
        'vcf': [False, True],
        'no_pilon': [False, True],
        'no_rotate': [True, False],
        'no_correct': [False, True],
        'no_miniasm': [False, True],
        'largest_component': [False, True],
        'keep': ['2', '1', '3', '0'],
        'verbosity': ['2', '3', '1', '0'],
        'mode': ['conservative', 'normal', 'bold'],
        'bowtie2_build_path': ['bowtie2-build'],
        'bcftools_path': ['bcftools'],
        'makeblastdb_path': ['makeblastdb'],
        'samtools_path': ['samtools'],
        'spades_path': ['spades.py'],
        'tblastn_path': ['tblastn'],
        'bowtie2_path': ['bowtie2'],
        'racon_path': ['racon'],
        'pilon_path': ['pilon'],
        'java_path': ['java'],
        'scores': '3,-6,-5,-2',
    }
    if 'kmers' not in params:
        params['kmers'] = defaults['kmers']
    else:
        k_errors = [x for x in params['kmers'] if not str(x).isdigit()]
        if len(k_errors):
            sys.exit('[spades] "k" must be integers (%s)' % ','.join(k_errors))
    sc = 'scores'
    if sc not in params:
        params[sc] = defaults[sc]
    else:
        if not isinstance(params[sc], str) or params[sc].count(',') != 3:
            sys.exit('[unicycler] Param "scores" must be 4 ","-separated ints')
    ints = ['min_component_size', 'min_dead_end_size', 'min_polish_size',
            'min_anchor_seg_len', 'min_fasta_length', 'kmer_count',
            'linear_seqs', 'low_score', 'kmers']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['min_bridge_qual']
    check_nums(self, params, defaults, floats, float, soft.name)
    floats_ = ['start_gene_cov', 'start_gene_id']
    check_nums(self, params, defaults, floats_, float, soft.name, 0, 100)
    floats__ = ['max_kmer_frac', 'min_kmer_frac', 'depth_filter']
    check_nums(self, params, defaults, floats__, float, soft.name, 0, 1)
    paths = [x for x in defaults if '_path' in x]
    check_default(self, params, defaults, soft.name,
                  (ints + floats + floats_ + floats__ + paths))
    defaults['contamination'] = '<Path to fasta file of known contaminations>'
    defaults['existing_long_read_assembly'] = '<GFA long-read assembly>'
    for path in ['bowtie2_build_path', 'bowtie2_path', 'bcftools_path',
                 'samtools_path', 'makeblastdb_path', 'tblastn_path',
                 'spades_path', 'racon_path', 'pilon_path', 'java_path']:
        defaults[path] = '<Path to the "%s" executable>' % path
    defaults['start_genes'] = '<Path to fasta file of genes for start point ' \
                              'of rotated replicons>'
    return defaults


def check_ngmerge(self, params, soft):
    defaults = {
        'p': 0.10,
        'm': 20,
        'e': 50,
        'q': 33,
        'u': 40,
        'a': [False, True],
        'd': [False, True],
        's': [False, True],
        'i': [False, True],
        'g': [False, True]
    }
    ints = ['m', 'e', 'q', 'u']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['p']
    check_nums(self, params, defaults, floats, float, soft.name)
    check_default(self, params, defaults, soft.name, (ints + floats))
    check_binary(self, soft.name, params, defaults, 'path')
    return defaults


def check_flash(self, params, soft):
    defaults = {
        'min_overlap': 10,
        'max_overlap': 65,
        'mismatch': 0.25
    }
    ints = ['min_overlap', 'max_overlap']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['mismatch']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + floats))
    return defaults


def check_pear(self, params, soft):
    defaults = {
        'max_uncalled_base': 1.,
        'p_value': 0.01,
        'min_assembly_length': 50,
        'max_assembly_length': 0,
        'quality_threshold': 0,
        'min_trim_length': 1,
        'score_method': 2,
        'min_overlap': 10,
        'test_method': 1,
        'phred_base': 33,
        'cap': 40,
        'empirical_freqs': [True, False],
        'keep_original': [False, True],
        'stitch': [False, True],
        'nbase': [False, True],
    }
    ints = ['min_assembly_length', 'max_assembly_length', 'quality_threshold',
            'min_trim_length', 'score_method', 'min_overlap', 'test_method',
            'phred_base', 'cap']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['p_value', 'max_uncalled_base']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + floats))
    return defaults


def check_bbmerge(self, params, soft):
    defaults = {
        'interleaved': [False, True],
        'reads': -1,
        'nzo': [True, False],
        'showhiststats': [True, False],
        'ziplevel': 2,
        'ordered': [False, True],
        'mix': [False, True],
        'qtrim': [False, True],
        'qtrim2': [False, True],
        'trimq': 10,
        'minlength': 1,
        'maxlength': -1,
        'tbo': [False, True],
        'minavgquality': 0,
        'maxexpectederrors': 0,
        'forcetrimleft': 0,
        'forcetrimright': 0,
        'forcetrimright2': 0,
        'forcetrimmod': 5,
        'ooi': [False, True],
        'trimpolya': [True, False],
        'usejni': [False, True],
        'merge': [True, False],
        'ecco': [False, True],
        'trimnonoverlapping': [False, True],
        'useoverlap': [True, False],
        'mininsert': 35,
        'mininsert0': 35,
        'minoverlap': 12,
        'minoverlap0': 8,
        'minq': 9,
        'maxq': 41,
        'entropy': [True, False],
        'efilter': 6,
        'pfilter': 0.00004,
        'kfilter': 0,
        'ouq': [False, True],
        'owq': [True, False],
        'usequality': [True, False],
        'iupacton': [False, True],
        'ratiomode': [True, False],
        'maxratio': 0.09,
        'ratiomargin': 5.5,
        'ratiooffset': 0.55,
        'maxmismatches': 20,
        'ratiominoverlapreduction': 3,
        'minsecondratio': 0.1,
        'forcemerge': [False, True],
        'flatmode': [False, True],
        'margin': 2,
        'mismatches': 3,
        'requireratiomatch': [False, True],
        'trimonfailure': [True, False],
        'strictness': [None, 'strict', 'verystrict', 'ultrastrict', 'maxstrict',
                       'loose', 'veryloose', 'ultraloose', 'maxloose', 'fast'],
        'k': 31,
        'extend': 0,
        'extend2': 0,
        'iterations': 1,
        'rem': [False, True],
        'rsem': [False, True],
        'ecctadpole': [False, True],
        'reassemble': [False, True],
        'removedeadends': [False, True],
        'removebubbles': [False, True],
        'mindepthseed': 3,
        'mindepthextend': 2,
        'branchmult1': 20,
        'branchmult2': 3,
        'branchlower':3,
        'ibb': [False, True],
        'prealloc': 1.,
        'prefilter': 0,
        'filtermem': 0,
        'minprob': 0.5,
        'minapproxoverlap': 26,
        'eccbloom': [False, True],
        'testmerge': [False, True],
        'eoom': [True, False],
        'da': [False, True],
    }
    ints = ['reads', 'ziplevel', 'trimq', 'minlength', 'maxlength',
            'minavgquality', 'maxexpectederrors', 'forcetrimleft',
            'forcetrimright', 'forcetrimright2', 'forcetrimmod', 'mininsert',
            'mininsert0', 'minoverlap', 'minoverlap0', 'minq', 'maxq',
            'efilter', 'kfilter', 'maxmismatches', 'ratiominoverlapreduction',
            'margin', 'mismatches', 'k', 'extend', 'extend2', 'iterations',
            'mindepthseed', 'mindepthextend', 'branchmult1', 'branchmult2',
            'branchlower', 'prefilter', 'filtermem', 'minapproxoverlap']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['pfilter', 'maxratio', 'ratiomargin',
              'ratiooffset', 'minsecondratio']
    check_nums(self, params, defaults, floats, float, soft.name)
    floats_ = ['prealloc', 'minprob']
    check_nums(self, params, defaults, floats_, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + floats + floats_))
    return defaults


def check_necat(self, params, soft):
    defaults = {
        'min_read_length': 3000,
        'prep_output_coverage': 40,
        'ovlp_fast_options_n': 500,
        'ovlp_fast_options_z': 20,
        'ovlp_fast_options_b': 2000,
        'ovlp_fast_options_e': 0.5,
        'ovlp_fast_options_j': 0,
        'ovlp_fast_options_u': 1,
        'ovlp_fast_options_a': 1000,
        'ovlp_sensitive_options_n': 500,
        'ovlp_sensitive_options_z': 10,
        'ovlp_sensitive_options_e': 0.5,
        'ovlp_sensitive_options_j': 0,
        'ovlp_sensitive_options_u': 1,
        'ovlp_sensitive_options_a': 1000,
        'cns_fast_options_a': 2000,
        'cns_fast_options_x': 4,
        'cns_fast_options_y': 12,
        'cns_fast_options_l': 1000,
        'cns_fast_options_e': 0.5,
        'cns_fast_options_p': 0.8,
        'cns_fast_options_u': 0,
        'cns_sensitive_options_a': 2000,
        'cns_sensitive_options_x': 4,
        'cns_sensitive_options_y': 12,
        'cns_sensitive_options_l': 1000,
        'cns_sensitive_options_e': 0.5,
        'cns_sensitive_options_p': 0.8,
        'cns_sensitive_options_u': 0,
        'trim_ovlp_options_n': 100,
        'trim_ovlp_options_z': 10,
        'trim_ovlp_options_b': 2000,
        'trim_ovlp_options_e': 0.5,
        'trim_ovlp_options_j': 1,
        'trim_ovlp_options_u': 1,
        'trim_ovlp_options_a': 400,
        'asv_ovlp_options_n': 100,
        'asv_ovlp_options_z': 10,
        'asv_ovlp_options_b': 2000,
        'asv_ovlp_options_e': 0.5,
        'asv_ovlp_options_j': 1,
        'asv_ovlp_options_u': 0,
        'asv_ovlp_options_a': 400,
        'num_iter': 2,
        'cns_output_coverage': 30,
        'cleanup': 1,
        'polish_contigs': [True, False]
    }
    ints = ['min_read_length', 'prep_output_coverage', 'num_iter',
            'cns_output_coverage', 'cleanup']
    ints += [x for x in defaults if x[-1] in 'axylnzbjua']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = [x for x in defaults if x[-2:] in ['_e', '_p']]
    check_nums(self, params, defaults, floats, float, soft.name)
    check_default(self, params, defaults, soft.name, (ints + floats))
    check_binary(self, soft.name, params, defaults, 'path')
    defaults['genome_size'] = '<estimated genome size (in bases)>'
    return defaults


def check_pooling(self, params, soft):
    param = 'pool_single_and_merged'
    defaults = {param: [True, False]}
    if param in params:
        if set(params[param]).issubset(self.config.techs):
            sys.exit('[pooling] Param "%s" can not be set per tech' % param)
    check_default(self, params, defaults, soft.name)
    return defaults


def check_cutadapt(self, params, soft):
    defaults = {
        'times': 1,
        'overlap': 3,
        'error_rate': 0.1,
        'quality_base': 33,
        'trim_n': [False, True],
        'revcomp': [False, True],
        'zero_cap': [False, True],
        'no_indels': [False, True],
        'pair_adapters': [False, True],
        'discard_trimmed': [False, True],
        'discard_untrimmed': [False, True],
        'match_read_wildcards': [False, True],
        'no_match_adapter_wildcards': [False, True],
        'report': ['full', 'minimal'],
        'pair-filter': ['any', 'both', 'first'],
        'action': ['trim', 'retain', 'mask', 'lowercase', 'none'],
    }
    ints = ['times', 'overlap', 'quality_base']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['error_rate']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (ints + floats))
    return defaults


def check_megahit(self, params, soft):
    defaults = {
        'mem_flag': [1, 0],
        'bubble_level': [2, 1, 0],
        'prune_level': [2, 1, 0, 3],
        'presets': [None, 'meta-sensitive', 'meta-large'],
        'k_list': '21,29,39,59,79,99,119,141',
        'k_min': 21,
        'k_max': 141,
        'k_step': 12,
        'min_count': 2,
        'prune_depth': 2,
        'cleaning_rounds': 5,
        'min_contig_len': 200,
        'merge_level': '20,0.95',
        'disconnect_ratio': 0.1,
        'low_local_ratio': 0.2,
        'memory': 0.9,
        'no_mercy': [False, True],
        'no_local': [False, True],
        'continue': [False, True],
        'kmin_1pass': [False, True],
        'no_hw_accel': [False, True],
        'keep_tmp_files': [False, True],
    }
    p = 'k_list'
    let_go = [p]
    if p not in params:
        params[p] = defaults[p]
    else:
        if ',' not in params[p]:
            sys.exit('[megahit] Param "%s" must be ","-separated values' % p)
        vals = params[p].split(',')
        if [v for v in vals if not v.isdigit() or 15 < int(v) or int(v) > 255]:
            sys.exit('[megahit] Param "%s" invalid: %s' % (p, params[p]))
    p = 'merge_level'
    let_go.append(p)
    if p not in params:
        params[p] = defaults[p]
    else:
        if params[p].count(',') != 1:
            sys.exit('[megahit] Param "%s" must be 2 ","-separated values' % p)
        l, s = params[p].split(',')
        if not str(l).isdigit():
            sys.exit('[megahit] Param "%s" invalid: %s' % (p, params[p]))
        try:
            float(str(s))
        except ValueError:
            sys.exit('[megahit] Param "%s" invalid: %s' % (p, params[p]))
    ints = ['k_min', 'k_max', 'k_step', 'min_count', 'prune_depth',
            'cleaning_rounds', 'min_contig_len']
    check_nums(self, params, defaults, ints, int, soft.name)
    floats = ['memory', 'low_local_ratio', 'disconnect_ratio']
    check_nums(self, params, defaults, floats, float, soft.name, 0, 1)
    check_default(self, params, defaults, soft.name, (let_go + ints + floats))
    return defaults


# def check_circlator(self, params, soft):
#     defaults = {
#         '': [],
#         '':,
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults


# def check_circlator(self, params, soft):
#     defaults = {
#         '': [],
#         '':,
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
