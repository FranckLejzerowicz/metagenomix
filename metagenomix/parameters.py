# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import sys

import numpy as np
from os.path import isdir

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


def check_default(params, defaults, name, let_go: list = []):
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


def check_search(self, params, soft):
    tool = soft.name.rsplit('_')[1]
    defaults = {}
    if tool == 'diamond':
        defaults.update(search_diamond(params, soft))
        dbs_existing = check_databases(tool, params, self.databases)
        valid_dbs = {}
        for db in dbs_existing:
            path = self.databases.paths[db]
            dmnds = '%s/diamond/*.dmnd' % path
            if not self.config.dev:
                dmnd_paths = glob.glob(dmnds)
                if dmnd_paths:
                    valid_dbs[db] = dmnd_paths
            else:
                valid_dbs[db] = [dmnds]
        params['databases'] = valid_dbs
        print(valid_dbs)

    if tool == 'hmmer':
        defaults.update(search_hmmer(params, soft))
        if 'terms' not in params:
            print('[search] Params "%s:terms" missing (search ignored)' % tool)
        else:
            terms = str(params['terms'])
            self.databases.get_pfam_terms(params)
            if not params['terms']:
                print('[search] Params "%s:terms" has no HMM profile for %s' % (
                    tool, terms))
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
    defaults = {'databases': None}
    if 'databases' not in params or not isinstance(params['databases'], list):
        sys.exit('[filtering] Param "databases" must a list of bowtie2 dbs')
    else:
        for d in params['databases']:
            if not self.config.dev and not glob.glob('%s/*' % d):
                sys.exit('[filtering] Param "databases" no "%s" bowtie2 db' % d)
    defaults['databases'] = '<list of databases>'
    return defaults


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
    defaults = {'mode': [False, 'fast', 'mid-sensitive', 'sensitive',
                         'more-sensitive', 'very-sensitive', 'ultra-sensitive'],
                'k': 1, 'top': 2, 'evalue': 0.001,
                'identity': 80, 'query_cov': 80}
    expand_search_params(params, defaults, 'diamond')
    ints, evals, ktop = ['identity', 'query_cov'], ['evalue'], ['k', 'top']
    check_nums(params, defaults, ints, int, soft.name, 0, 100)
    check_nums(params, defaults, evals, float, soft.name, 0, 1)
    check_nums(params, defaults, ktop, int, soft.name)
    check_default(params, defaults, soft.name, (ints + evals + ktop))
    return defaults


def search_hmmer(params, soft) -> dict:
    defaults = {'max': [True, False], 'nobias': [True, False],
                'Z': 10, 'E': 1, 'domZ': 10, 'domE': 1}
    expand_search_params(params, defaults, 'hmmer')
    e_vals, z_vals = ['E', 'domE'], ['Z', 'domZ']
    check_nums(params, defaults, e_vals, float, soft.name, 0, 100)
    check_nums(params, defaults, z_vals, int, soft.name, 0, 1000)
    check_default(params, defaults, soft.name, (e_vals + z_vals + ['terms']))
    return defaults


def check_metaxa2(self, params, soft) -> dict:
    defaults = {'align': ['none', 'all', 'uncertain'],
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
    check_nums(params, defaults, ['r', 'd', 'E'], float, soft.name, 0, 1)
    check_nums(params, defaults, ['R'], int, soft.name, 0, 100)
    check_nums(params, defaults, ['S', 'N', 'M', 'l'], int, soft.name)
    check_default(params, defaults, soft.name, ['r', 'd', 'E', 'R',
                                                'S', 'N', 'M', 'l', 'T'])
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
    check_default(params, defaults, soft.name)
    if 1:
        valid_dbs = {}
        dbs_existing = check_databases('shogun', params, self.databases)
        for db in dbs_existing:
            path = self.databases.paths[db]
            metadata_yml = '%s/shogun/metadata.yaml' % path
            if not self.config.dev and not isfile(metadata_yml):
                raise IOError('[shogun] file "%s" must exist' % metadata_yml)
            metadata = read_yaml(metadata_yml)
            for aligner in list(params['aligners']):
                if aligner in metadata:
                    if not self.config.dev:
                        ali = metadata[aligner]
                        if ali[0] == '/':
                            formatted = ali
                        else:
                            formatted = '%s/shogun/%s' % (path, ali)
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
        path = self.databases.paths[db]
        bt2_path = '%s/bowtie2/*.*.bt2' % path
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
    defaults = {'simkaMin': [True, False],
                'kmer': np.linspace(15, 80, 6),
                'log_reads': np.logspace(3, 7, 3)}
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
    check_default(params, defaults, soft.name, ['kmer', 'log_reads'])
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
    check_default(params, defaults, soft.name, mins)
    return defaults


def check_flash(self, params, soft):
    defaults = {'min_overlap': 30, 'max_overlap': 100, 'mismatch': 0}
    check_nums(params, defaults, list(defaults.keys()), int, soft.name)
    mi, ma = params['min_overlap'], params['max_overlap']
    if mi > ma:
        sys.exit('[flash] "min_overlap" (%s) > "max_overlap" (%s)' % (mi, ma))
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
        'uniref': ['uniref90'],
        'skip_translated': [False, True],
        'nucleotide_db': '/path/to/full_chocophlan.v0.1.1/chocophlan',
        'protein_db': '/path/to/uniref90/uniref',
    })
    let_go = ['profiles']
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
    defaults = {'focus': {'all': [self.databases.paths['midas'], '']}}
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
            if not isinstance(v, list) or len(v) != 2:
                sys.exit('[midas] Param "focus::%s" must be 2-items list' % k)
            if not self.config.dev and not isfile(v[1]):
                sys.exit('[midas] Param "focus::%s::%s" not a file' % (k, v[1]))
    defaults = {'focus': '<dict of >'}
    return defaults


def check_macsyfinder(self, params, soft):
    defaults = {'db_type': (['unordered', 'ordered_replicon', 'gembase'], str),
                'replicon_topology': (['linear', 'circular'], str),
                'models': (['TXSS', 'TFF-SF', 'CAS'], list),
                'evalue': (0.1, float),
                'coverage': (0.5, float)}
    check_nums(params, defaults, ['evalue', 'coverage'], float, soft.name, 0, 1)
    let_go = ['evalue', 'coverage']
    check_default(params, defaults, soft.name, let_go)
    return defaults


# def check_TOOL(params, soft, databases, config):
#     defaults = {}
