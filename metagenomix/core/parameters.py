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


def tech_params(self, tech: str, soft: str = None):

    params = self.soft.params
    if soft in params and isinstance(params[soft], dict):
        tool_params = params.pop(soft)
        for param in list(tool_params):
            if param in params:
                del tool_params[param]
        params.update(tool_params)

    if tech not in ['illumina', 'pacbio', 'nanopore']:
        return params

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
    if not config.dev and value is not None:
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


def check_nums(self, params, defaults, vals, dtype, mi=None, ma=None):
    for param in vals:
        if param not in params:
            if defaults[param] == [None]:
                params[param] = None
            else:
                params[param] = defaults[param]
        else:
            val = params[param]
            if isinstance(val, dict) and set(val).issubset(self.config.techs):
                for tech, v in val.items():
                    check_nums_generic(param, v, dtype, self.name, mi, ma, tech)
            else:
                check_nums_generic(param, val, dtype, self.name, mi, ma)


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


def check_default_type(self, params, param, multi):
    vals = params[param]
    if isinstance(vals, dict) and set(vals).issubset(self.config.techs):
        for tech, v in vals.items():
            check_type(v, param, multi, self.name, tech)
    else:
        check_type(vals, param, multi, self.name)


def check_type_sys_exit(tool, param, no_dict, tech=None):
    if no_dict:
        m = '[%s] Param "%s" must not be a list/dict' % (tool, param)
    else:
        m = '[%s] Param "%s" must be a list' % (tool, param)
    if tech:
        m += ' (%s)' % tech
    sys.exit(m)


def add_skip_samples(
        params: dict,
        defaults: dict,
        tool: str
) -> None:
    if 'skip_samples' not in params:
        params['skip_samples'] = [None]
    elif not isinstance(params['skip_samples'], list):
        sys.exit('[%s] Param "skip_samples" not a (samples) list' % tool)
    defaults['skip_samples'] = '<list of samples / co-assembly groups to skip>'


def check_default(
        self,
        params,
        defaults,
        let_go: list = [],
        multi: list = []
):
    """Verifies that the parameters given by the used for the current tool
    and that are to be checked (i.e., not in the `let_go` list) are valid.
    This mean that they must be one (or several) of the values matching the
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
                    check_default_generic(param, v, values, self.name, tech)
            else:
                check_default_generic(param, vals, values, self.name)
        check_default_type(self, params, param, multi)
    add_skip_samples(params, defaults, self.name)


def check_binary(self, params, defaults, opt, n=None):
    if n is None:
        n = self.name
    if [x for x in self.config.modules.get(n, []) if x.lower().startswith(n)]:
        return None
    if opt == 'path':
        message = 'Param "path" for path to binaries folder missing'
        defaults['path'] = '<Path to folder containing the %s binary>' % n
        file_or_dir = isdir
    if opt in ['binary', 'trimmomatic', 'bowtie2', 'anicalculator']:
        message = 'Param "%s" for path to binary/executable' % opt
        defaults[opt] = '<Path to the %s binary>' % self.name
        file_or_dir = isfile
    if opt not in params:
        sys.exit('[%s] %s' % (n, message))
    if not (self.config.dev or self.config.workshop):
        if params[opt] != 'module' and not file_or_dir(params[opt].split()[0]):
            sys.exit('[%s] Please provide valid path to param "%s"' % (opt, n))


class Parameters(object):
    """
    This class shall contains the checks for all softwares specifically,
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
            if self.config.dev or self.config.workshop:
                valid_dbs[db] = [files]
            else:
                file_paths = glob.glob(files)
                if file_paths:
                    valid_dbs[db] = file_paths
    params['databases'] = valid_dbs

#
# def get_hmmer_databases(self, params):
#     empties = set()
#     acc = self.databases.get_pfam(descriptions=params.get('descriptions', []))
#     acc = self.databases.get_pfam(accessions=params.get('accessions', []))
#     acc = self.databases.get_pfam(interpro=params.get('interpro', []))
#
#         not_found_terms = set(terms).difference(set(terms_hmms_dias))
#         if not_found_terms:
#             print('[search] No HMM for hmmer "terms"')
#             for term in not_found_terms:
#                 print('[search]    - %s' % term)
#             if len(not_found_terms) == len(set(terms)):
#                 print('[search]    (i.e. all terms, hmmer step ignored)')
#     else:
#         empties.add('[search] Params "hmmer:descriptions" missing (ignored)')
#
#
# def check_search(self, params, soft):
#     tool = soft.name.rsplit('_')[1]
#     defaults = {}
#     defaults.update(search_tool(self, params, soft, tool))
#     get_diamond_hmmer_databases(self, tool, params)
#     # if tool == 'hmmer':
#     #     if self.databases.hmms_pd.shape[0]:
#     #         get_hmmer_databases(self, params)
#     #     else:
#             print('[search] Empty Pfam data file (hmmer step ignored)')
#     defaults['databases'] = '<list of databases>'
#     return defaults


def check_ccfind(self, params):
    defaults = {
        'preserve_tmpdir': [False, True],
        'terminal_fragment_size': 500,
        'min_percent_identity': 94,
        'min_aligned_length': 50
    }
    int1 = ['terminal_fragment_size', 'min_aligned_length']
    check_nums(self, params, defaults, int1, int)
    int2 = ['min_percent_identity']
    check_nums(self, params, defaults, int2, int, 0, 100)
    check_default(self, params, defaults, (int1 + int2))
    check_binary(self, params, defaults, 'path')
    check_binary(self, params, defaults, 'binary')
    defaults['path'] = '<path to the ssearch36 installation folder>'
    defaults['binary'] = '<path to the ccfind binary>'
    return defaults


def check_barrnap(self, params):
    defaults = {
        'evalue': 1e-06,
        'reject': 0.25,
        'lencutoff': 0.8,
        'incseq': [True, False],
        'kingdom': ['bac', 'mito', 'arc', 'euk']
    }
    floats = ['evalue', 'reject', 'lencutoff']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, floats)
    return defaults


def check_integronfinder(self, params):
    defaults = {
        'distance_threshold': 400,
        'evalue_attc': 1,
        'min_length': 1500,
        'min_attc_size': 40,
        'max_attc_size': 200,
        'calin_threshold': 2,
        'pdf': [True, False],
        'gbk': [True, False],
        'mute': [True, False],
        'gembase': [False, True],
        'local_max': [True, False],
        'promoter_attI': [True, False],
        'union_integrases': [False, True],
        'keep_palindromes': [False, True],
        'no_proteins': [False, True],
        'func_annot': [True, False],
        'circ': [False, True],
        'linear': [True, False],
        'keep_tmp': [True, False]
    }
    paths = [
        'annot_parser',
        'attc_model',
        'topology_file',
        'path_func_annot'
    ]
    nf = ['parallel_integron_finder.nf', 'nextflow.config']
    if 'path' in params and not (self.config.dev or self.config.workshop):
        if sum([isfile('%s/%s' % (params['path'], x)) for x in nf]) != 2:
            print('[integronfinder] Will run without using nextflow (slower)')

    for p in paths:
        if p in params:
            if not isfile(str(params[p])):
                sys.exit('[integronfinder] Path "%s" do not exist' % params[p])
        else:
            params[p] = None
    ints = ['distance_threshold', 'min_length', 'min_attc_size',
            'max_attc_size', 'calin_threshold']
    check_nums(self, params, defaults, ints, int)
    ints1 = ['evalue_attc']
    check_nums(self, params, defaults, ints1, float, 0, 100)
    let_go = ints + ints1 + ["path"]
    check_default(self, params, defaults, let_go)
    defaults['annot_parser'] = '<Parser to use to get info from protein file>'
    defaults['attc_model'] = '<Path to the attc model (Covariance Matrix)>'
    defaults['topology_file'] = '<Path to file with each replicon topology>'
    defaults['path_func_annot'] = '<Path containing all hmm bank paths>'
    defaults['path'] = '<Path to folder containing "%s">' % '", "'.join(nf)
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


def check_prokka(self, params):
    defaults = {
        'kingdom': ['Archaea', 'Bacteria', 'Mitochondria', 'Viruses'],
        'evalue': 1e-09,
        'coverage': 80,
        'mincontiglen': 100,
        'gcode': [11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        'notrna': [False, True],
        'norrna': [False, True],
        'metagenome': [True, False],
        'addgenes': [False, True],
        'addmrna': [False, True],
        'gffver': 3,
        'compliant': [False, True],
        'centre': [''],
        'gram': [''],
        'accver': 1,
        'rawproduct': [False, True],
        'cdsrnaolap': [False, True],
        'fast': [False, True],
        'noanno': [False, True],
        'rnammer': [False, True],
        'rfam': 0
    }
    floats = ['evalue']
    int1 = ['mincontiglen', 'gffver', 'accver', 'rfam']
    int2 = ['coverage']
    check_nums(self, params, defaults, floats, float, 0, 100)
    check_nums(self, params, defaults, int1, int)
    check_nums(self, params, defaults, int2, int, 0, 100)
    check_default(self, params, defaults, (floats + int1 + int2))
    if 'config' in params and params['config'] is not None:
        if not isfile(params['config']):
            sys.exit('[prokka] Param "config" must be an existing file')
        check_prokka_config(params)
    else:
        params['config'] = None
    defaults['config'] = '<Path to config file for prokka>'
    return defaults


def check_antismash(self, params):
    defaults = {
        'tta_threshold': 0.65,
        'rre': [False, True],
        'asf': [False, True],
        'cassis': [False, True],
        'pfam2go': [False, True],
        'tigrfam': [False, True],
        'cc_mibig': [False, True],
        'fullhmmer': [False, True],
        'cb_general': [False, True],
        'smcog_trees': [False, True],
        'clusterhmmer': [False, True],
        'cb_subclusters': [False, True],
        'cb_knownclusters': [False, True],
        'taxon': ['bacteria', 'fungi'],
        'genefinding_tool': [
            'prodigal-m', 'prodigal', 'glimmerhmm', 'none', 'error']
    }
    floats = ['tta_threshold']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, floats)
    return defaults


def check_quast(self, params):
    defaults = {
        'no_check': [False, True],
        'no_plots': [False, True],
        'no_html': [False, True],
        'no_icarus': [False, True],
        'no_snps': [False, True],
        'no_gc': [False, True],
        'no_sv': [False, True],
        'no_read_stats': [False, True],
        'fast': [False, True],
        'sv_bedpe': [None],
        'sam': [None],
        'bam': [None],
        'ref_sam': [None],
        'ref_bam': [None],
        'space_efficient': [False, True],
        'memory_efficient': [False, True],
        'plots_format': ['pdf', 'emf', 'eps', 'png', 'ps',
                         'raw', 'rgba', 'svg', 'svgz'],
        'report_all_metrics': [False, True],
        'est_insert_size': [None],
        'upper_bound_min_con': [None],
        'upper_bound_assembly': [False, True],
        'scaffold_gap_max_size': 10000,
        'unaligned_part_size': 500,
        'skip_unaligned_mis_contigs': [False, True],
        'fragmented': [False, True],
        'fragmented_max_indent': 200,
        'local_mis_size': 200,
        'extensive_mis_size': 1000,
        'strict_NA': [False, True],
        'unique_mapping': [False, True],
        'ambiguity_score': 0.99,
        'ambiguity_usage': ['one', 'none', 'all'],
        'min_identity': 95.0,
        'min_alignment': 65,
        'use_all_alignments': [False, True],
        'reuse_combined_alignments': [False, True],
        'x_for_Nx': 90,
        'contig_thresholds': '0,1000,5000,10000,25000,50000',
        'use_input_ref_order': [False, True],
        'blast_db': [None],
        'max_ref_number': 50,
        'operons': [None],
        'conserved_genes_finding': [True, False],
        'rna_finding': [True, False],
        'gene_thresholds': '0,300,1500,3000',
        'glimmer': [False, True],
        'gene_finding': [True, False],
        'circos': [True, False],
        'k_mer_size': 101,
        'k_mer_stats': [False, True],
        'large': [False, True],
        'fungus': [False, True],
        'eukaryote': [False, True],
        'split_scaffolds': [False, True],
        'min_contig': 500,
        'r': [None],
        'references_list': [None],
        'g': [None]
    }
    for param in ['contig_thresholds', 'gene_thresholds']:
        if param in params:
            if [x for x in str(params[param]).split(',') if not x.isdigit()]:
                sys.exit('[quast] "%s" invalid: "%s"' % (param, params[param]))
        else:
            params[param] = defaults[param]
    ints = ['min_contig', 'min_alignment', 'k_mer_size', 'max_ref_number',
            'extensive_mis_size', 'local_mis_size', 'fragmented_max_indent',
            'unaligned_part_size', 'scaffold_gap_max_size']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['x_for_Nx']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    floats = ['min_identity']
    check_nums(self, params, defaults, floats, float, 80, 100)
    floats2 = ['ambiguity_score']
    check_nums(self, params, defaults, floats2, float, 0.8, 1)
    if 'label' in params:
        if params['label'] not in set(self.config.meta):
            sys.exit('[quast] Params "label" must be a valid metadata column')
    let_go = ints + ints2 + floats + floats2 + [
        'label', 'contig_thresholds', 'gene_thresholds']
    check_default(self, params, defaults, let_go)
    defaults['label'] = '<an existing metadata column>'
    check_binary(self, params, defaults, 'path')

    if 'circos_path' not in params:
        params['circos_path'] = None
    else:
        db_dir = params['circos_path']
        if not (self.config.dev or self.config.workshop):
            if db_dir and not isdir(db_dir):
                sys.exit('[quast] Params "circos_path": path do not exist')

    defaults['path'] = '<path to the quast installation folder>'
    defaults['circos_path'] = '<path to the circos installation folder>'
    return defaults


def check_bowtie2(self, params, no_database=False):
    defaults = {
        'presets': ['sensitive', 'very-sensitive', 'very-fast', 'fast', None],
        'paired': [True, False],
        'fr': ['fr', 'rf', 'ff'],
        'i': 'S,1,1.15',
        'n_ceil': 'L,0,0.15',
        'rdg': '5,3',
        'rfg': '5,3',
        'score_min': [None, 'G,20,8', 'L,0,-0.05'],
        'skip': 0,
        'upto': 0,
        'trim5': 0,
        'trim3': 0,
        'trim_to': 0,
        'N': 0,
        'L': 22,
        'dpad': 15,
        'gbar': 4,
        'ma': -1000,
        'k': 0,
        'np': 1,
        'mp': 6,
        'D': 15,
        'R': 2,
        'minins': 0,
        'maxins': 500,
        'met': 240,
        'seed': 12345,
        'q': [True, False],
        'tab5': [False, True],
        'tab6': [False, True],
        'qseq': [False, True],
        'f': [False, True],
        'r': [False, True],
        'c': [False, True],
        'phred33': [True, False],
        'phred64': [False, True],
        'int_quals': [False, True],
        'ignore_quals': [False, True],
        'nofw': [False, True],
        'norc': [False, True],
        'no_1mm_upfront': [False, True],
        'end_to_end': [True, False],
        'local': [False, True],
        'all': [False, True],
        'no_mixed': [False, True],
        'no_discordant': [False, True],
        'dovetail': [False, True],
        'no_contain': [False, True],
        'no_overlap': [False, True],
        'align_paired_reads': [False, True],
        'preserve_tags': [False, True],
        't': [False, True],
        'un': [False, True],
        'al': [False, True],
        'un_conc': [False, True],
        'al_conc': [False, True],
        'quiet': [False, True],
        'met_file': [False, True],
        'met_stderr': [False, True],
        'no_unal': [True, False],
        'no_head': [False, True],
        'no_sq': [False, True],
        'omit_sec_seq': [False, True],
        'sam_no_qname_trunc': [False, True],
        'xeq': [False, True],
        'soft_clipped_unmapped_tlen': [False, True],
        'sam_append_comment': [False, True],
        'reorder': [False, True],
        'mm': [False, True],
        'qc_filter': [False, True],
        'non_deterministic': [False, True]
    }
    if 'mp' in params:
        if ',' in str(params['mp']):
            if len([x.isdigit() for x in str(params['mp']).split(',')]) != 2:
                sys.exit('[bowtie2] Dual-value "mp" option must be two INTs')
        elif not str(params['mp']).isdigit():
            sys.exit('[bowtie2] Single-value "mp" option must be an int')
    else:
        params['mp'] = defaults['mp']

    for param in ['rdg', 'rfg']:
        if param in params:
            if len([x.isdigit() for x in str(params[param]).split(',')]) != 2:
                sys.exit('[bowtie2] "%s" option invalid' % param)
        else:
            params[param] = defaults[param]

    for param in ['i', 'n_ceil', 'score_min']:
        if param in params:
            s = params[param].split(',')
            if len(s) != 3 or s[0] not in 'CLSG':
                sys.exit('[bowtie2] "%s" option invalid' % param)
            else:
                for r in [1, 2]:
                    try:
                        float(s[r])
                    except ValueError:
                        sys.exit('[bowtie2] "%s" option invalid' % param)
        else:
            if param == 'score_min':
                params[param] = None
            else:
                params[param] = defaults[param]

    ints = ['skip', 'upto', 'trim5', 'trim3', 'trim_to', 'dpad', 'gbar',
            'k', 'np', 'D', 'R', 'minins', 'maxins', 'met', 'seed', 'ma']
    check_nums(self, params, defaults, ints, int)
    check_nums(self, params, defaults, ['N'], int, 0, 1)
    check_nums(self, params, defaults, ['L'], int, 4, 31)

    let_go = ints + ['N', 'L', 'i', 'n_ceil', 'score_min', 'rdg', 'rfg', 'mp']
    check_default(self, params, defaults, let_go)

    if not no_database:
        dbs_existing = check_databases(self.name, params, self.databases)
        valid_dbs = {}
        for db in dbs_existing:
            if 'bowtie2' in self.databases.builds[db]:
                bt2_path = '%s*.*.bt2*' % self.databases.builds[db]['bowtie2']
                if self.config.dev or self.config.workshop:
                    valid_dbs[db] = bt2_path.rsplit('.', 2)[0]
                else:
                    bt2_paths = glob.glob(bt2_path)
                    if bt2_paths:
                        valid_dbs[db] = bt2_paths[0].rsplit('.', 2)[0]
        params['databases'] = valid_dbs
        defaults['databases'] = '<list of databases>'
    return defaults


def check_mapping(self, params):
    defaults = {
        'aligners': ['minimap2', 'bowtie2', 'bwa', 'bbmap'],
        'per_tech': [False, True],
        'all_reads': [False, True]
    }
    for aligner in defaults['aligners']:
        aligner_params = dict(params.get(aligner, {}))
        check_aligner = 'check_%s' % aligner
        if check_aligner in globals():
            func = globals()[check_aligner]
            aligner_defaults = func(self, aligner_params, True)
            defaults[aligner] = aligner_defaults
            params[aligner] = aligner_params
    check_default(self, params, defaults, defaults['aligners'], ['aligners'])
    return defaults


def check_filtering(self, params):
    defaults = {
        'aligner': ['bowtie2', 'minimap2', 'bwa', 'bbmap'],
    }
    if "databases" not in params or not isinstance(params["databases"], list):
        sys.exit('[filtering] Param "databases" not a list of databases')

    for aligner in defaults['aligner']:
        aligner_params = dict(params.get(aligner, {}))
        check_aligner = 'check_%s' % aligner
        if check_aligner in globals():
            func = globals()[check_aligner]
            aligner_defaults = func(self, aligner_params, True)
            defaults[aligner] = aligner_defaults
            params[aligner] = aligner_params

    check_default(self, params, defaults, defaults['aligner'])

    aligner = params['aligner']
    for db in params["databases"]:
        if self.config.dev or self.config.workshop:
            continue
        if aligner == 'bowtie2':
            if db not in self.databases.paths:
                sys.exit('[filtering] DB "%s" not in databases config' % db)
            elif db not in self.databases.builds:
                sys.exit('[filtering] DB "%s" not built for bowtie2' % db)
            elif not self.databases.builds[db].get('bowtie2'):
                sys.exit('[filtering] DB "%s" not built for bowtie2' % db)
            dbs = glob.glob('%s*.bt2*' % self.databases.builds[db]['bowtie2'])
            if len(dbs) != 6:
                sys.exit('[filtering] Param "databases" bowtie2 files '
                         'missing for database "%s"' % db)
    defaults['databases'] = '<list of databases>'
    return defaults


def check_atropos(self, params):
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
    for aA in ['a', 'A']:
        if aA not in params:
            params[aA] = []
        if not isinstance(params[aA], list):
            sys.exit('[atropos] Param "%s" must be a list' % aA)

    ints = ['q', 'overlap', 'max_reads', 'indel_cost', 'nextseq_trim',
            'minimum_length', 'quality_cutoff']
    floats = ['error_rate']
    check_nums(self, params, defaults, ints, int)
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (floats + ints + ['a', 'A']))
    return defaults


def check_kneaddata(self, params):
    defaults = {
        'trimmomatic': None,
        'bowtie2': None,
        'purge': [True, False]}
    tools = ['trimmomatic', 'bowtie2']
    for tool in tools:
        if tool not in params:
            params[tool] = defaults[tool]
        else:
            check_binary(self, params, defaults, tool, 'kneaddata')
    if 'databases' not in params or not isinstance(params['databases'], list):
        sys.exit('[kneaddata] Params "databases" must a list of existing paths')
    check_default(self, params, defaults, (tools + ['databases']))
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


def search_diamond(self, params) -> dict:
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
    check_nums(self, params, defaults, ints_, int)
    check_nums(self, params, defaults, ints, int, 0, 100)
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints_ + ints + floats))
    return defaults


def search_hmmer(self, params) -> dict:
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
    evals, zvals = ['E', 'domE'], ['Z', 'domZ', 'textw']
    floats = ['F1', 'F2', 'F3']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_nums(self, params, defaults, evals, float, 0, 100)
    check_nums(self, params, defaults, zvals, int, 0, 1000)
    check_default(self, params, defaults, (evals + zvals + floats + ['terms']))
    defaults['descriptions'] = '<List of strings to search Pfam descriptions>'
    defaults['accessions'] = '<List of strings to search Pfam accessions>'
    defaults['interpro'] = '<List of tsv files exported from InterPro>'
    return defaults


def search_tool(self, params, tool) -> dict:
    if tool == 'diamond':
        defaults = search_diamond(self, params)
    else:
        defaults = search_hmmer(self, params)
    return defaults


def check_metaxa2(self, params) -> dict:
    defaults = {
        'g': ['ssu', 'lsu', 'string'],
        # 'f': ['auto', 'fasta', 'fastq', 'paired-end', 'paired-fasta'],
        'mode': ['metagenome', 'genome', 'auto'],
        'p': [None],
        'date': [False, True],
        'plus': [False, True],
        'q': 20,
        'databases': ['default'],
        'quality_percent': 10,
        'quality_filter': [False, True],
        'quality_trim': [False, True],
        'ignore_paired_read': [True, False],
        'distance': 150,
        't': ['bacteria', 'b', 'archaea', 'a', 'eukaryota', 'e', 'all', 'A',
              'mitochondrial', 'm', 'chloroplast', 'c', 'other', 'o'],
        'E': 1,
        'S': 12,
        'N': 2,
        'M': 5,
        'R': 75,
        'T': '0,60,70,75,85,90,97',
        'H': 5,
        'selection_priority': ['score', 'domains', 'eval', 'sum'],
        'scoring_model': ['new', 'old'],
        'search_eval': [None],
        'search_score': 0,
        'blast_eval': 1e-5,
        'blast_score': [None],
        'blast_wordsize': 14,
        'blast_bin': [None],
        'allow_single_domain': '1e-10,0',
        'allow_reorder': [True, False],
        'complement': [True, False],
        'multi_thread': [False, True],
        'heuristics': [True, False],
        'megablast': [False, True],
        'ref_identity': 99,
        'summary': [True, False],
        'graphical': [True, False],
        'fasta': [True, False],
        'split_pairs': [False, True],
        'table': [False, True],
        'taxonomy': [True, False],
        'reltax': [False, True],
        'taxlevel': 0,
        'not_found': [False, True],
        'align': ['none', 'all', 'uncertain'],
        'truncate': [True, False],
        'guess_species': [False, True],
        'silent': [False, True],
        'graph_scale': 0,
        'save_raw': [False, True],
        'usearch': 0,
        'usearch_bin': [None],
        'ublast': [True, False],
        'r': 0,
        'l': 50,
        'd': 0,
        'm': 0,
        'n': 1,
        's': [False, True],
        'remove_na': [True, False],
        'lists': [True, False],
        'separate': [True, False],
        'unknown': [False, True],
        'reference': [None],
    }
    if 'T' not in params:
        params['T'] = defaults['T']
    elif len([x for x in str(params['T']).split(',') if x.isdigit()]) != 7:
        sys.exit('[metaxa2] Param "T" must be 7 tab-separated integers [0-100]')

    p = 'allow_single_domain'
    if p not in params:
        params[p] = defaults[p]
    else:
        if params[p].count(',') != 1:
            sys.exit('[metaxa2] Param "%s" must have two values' % p)
        e, s = params[p].split(',')
        try:
            float(s)
        except ValueError:
            sys.exit('[metaxa2] Param "%s": first value not float' % p)
        if not str(s).isdigit():
            sys.exit('[metaxa2] Param "%s": second value not integer' % p)

    # check_databases('metaxa2', params, self.databases)
    floats = ['E', 'usearch', 'blast_eval', 'search_eval']
    check_nums(self, params, defaults, floats, float, 0, 1)
    ints = ['S', 'N', 'M', 'q', 'distance', 'search_score', 'blast_wordsize',
            'taxlevel', 'graph_scale', 'l', 'm', 'n']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['R', 'H', 'quality_percent', 'ref_identity', 'r', 'd']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    check_default(self, params, defaults, (floats + ints + ints2 + [
        'databases', 'T', 'allow_single_domain', 'blast_bin']))
    check_binary(self, params, defaults, 'path')
    defaults['databases'].append('<list of databases>')
    return defaults


def check_kraken2(self, params):
    defaults = {
        'databases': ['default'],
        'confidence': 0.5
    }
    check_databases(self.name, params, self.databases)
    check_nums(self, params, defaults, ['confidence'], float, 0, 1)
    check_default(self, params, defaults, ['confidence', 'databases'])
    defaults['databases'].append('<list of databases>')
    return defaults


def check_shogun(self, params):
    defaults = {'aligners': ['bowtie2', 'burst', 'utree']}
    check_default(self, params, defaults, [], ['aligners'])
    if 1:
        valid_dbs = {}
        dbs_existing = check_databases(self.name, params, self.databases)
        for db in dbs_existing:
            path = self.databases.paths[db]
            yamls = []
            for folder in ['', 'databases/']:
                yaml = '%s/%sshogun/metadata.yaml' % (path, folder)
                yamls.append(yaml)
                if isfile(yaml):
                    break
            else:
                if not self.config.dev or self.config.workshop:
                    raise IOError(
                        '[shogun] file must exist: %s' % ' or '.join(yamls))
            metadata = read_yaml(yaml)
            for aligner in list(params['aligners']):
                if aligner in metadata:
                    if not self.config.dev or self.config.workshop:
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


def check_minimap2(self, params, no_database=False):
    defaults = {
        'paired': [True, False],
        'k': 15,
        'w': 10,
        'H': [False, True],
        'I': '10G',
        'idx_no_seq': [False, True],
        'd': [''],
        'f': 0.0002,
        'min_occ_floor': 0,
        'g': 10000,
        'r': 500,
        'n': 3,
        'm': 40,
        'D': [False, True],
        'P': [False, True],
        'dual': [True, False],
        'X': [False, True],
        'p': 0.8,
        'N': 5,
        'G': '200k',
        'F': 800,
        'M': 0.5,
        'hard_mask_level': [False, True],
        'max_chain_skip': 25,
        'max_chain_iter': 5000,
        'no_long_join': [False, True],
        'lj_min_ratio': 0.5,
        'splice': [False, True],
        'sr': [False, True],
        'split_prefix': [None],
        'frag': [False, True],
        'for_only': [False, True],
        'rev_only': [False, True],
        'heap_sort': [False, True],
        'no_pairing': [False, True],
        'A': 2,
        'B': 4,
        'O': '4,24',
        'E': '2,1',
        'C': 0,
        'z': '400,200',
        's': 80,
        'u': ['n', 'f', 'b'],
        'end_bonus': 0,
        'score_N': 1,
        'splice_flank': [True, False],
        'junc_bed': [''],
        'junc_bonus': 9,
        'end_seed_pen': 6,
        'no_end_flt': [False, True],
        'cap_sw_mem': '100m',
        'Q': [False, True],
        'L': [False, True],
        'R': [''],
        'y': [False, True],
        'c': [False, True],
        'cs': [None, 'short', 'long'],
        'MD': [False, True],
        'eqx': [False, True],
        'Y': [False, True],
        'seed': 11,
        '2': [False, True],
        'K': '500M',
        'secondary': [True, False],
        'max_qlen': [None],
        'paf_no_hit': [False, True],
        'sam_hit_only': [True, False],
        'x': [None, 'map-pb', 'map-ont', 'map-hifi', 'ava-pb', 'ava-ont',
              'asm5', 'asm10', 'asm20', 'splice', 'splice:hq', 'sr'],
        'no_kalloc': [False, True],
        'print_qname': [False, True],
        'print_seeds': [False, True],
        'a': [True, False]
    }

    for p in ['O', 'K', 'z']:
        if p not in params:
            params[p] = defaults[p]
        elif len([x.isdigit() for x in str(params[p]).split(',')]) != 2:
            sys.exit('[minimap2] "%s" option invalid (INT,INT)' % p)

    for p in ['I', 'E', 'G', 'cap_sw_mem']:
        if p not in params:
            params[p] = defaults[p]
        elif not str(params[p][:-1]).isdigit():
            sys.exit('[minimap2] "%s" option invalid (NUM)' % p)

    ints = ['max_qlen', 'w', 'g', 'F', 'r', 'n', 'm', 'N', 'A', 'B', 's',
            'max_chain_skip', 'max_chain_iter', 'C', 's', 'end_bonus',
            'score_N', 'junc_bonus', 'seed', 'min_occ_floor', 'end_seed_pen']
    check_nums(self, params, defaults, ints, int)
    check_nums(self, params, defaults, ['k'], int, 3, 28)

    floats = ['f', 'p', 'M', 'lj_min_ratio']
    check_nums(self, params, defaults, floats, float, 0, 1)

    check_default(self, params, defaults, (ints[1:]+floats+['k']))

    if not no_database:
        dbs_existing = check_databases(self.name, params, self.databases)
        valid_dbs = {}
        for db in dbs_existing:
            if 'minimap2' in self.databases.builds[db]:
                db_path = '%s/*.mmi' % self.databases.builds[db]['minimap2']
                if not self.config.dev or self.config.workshop:
                    bt2_paths = glob.glob(db_path)
                    if bt2_paths:
                        valid_dbs[db] = bt2_paths[0]
                else:
                    valid_dbs[db] = db_path
        params['databases'] = valid_dbs
        defaults['databases'] = '<list of databases>'
    return defaults


def check_bbmap(self, params, no_database=False):
    defaults = {
        'paired': [True, False],
    }
    if not no_database:
        dbs_existing = check_databases(self.name, params, self.databases)
        valid_dbs = {}
        for db in dbs_existing:
            if 'bbmap' in self.databases.builds[db]:
                db_path = '%s/*' % self.databases.builds[db]['bbmap']
                if not self.config.dev or self.config.workshop:
                    bt2_paths = glob.glob(db_path)
                    if bt2_paths:
                        valid_dbs[db] = bt2_paths[0]
                else:
                    valid_dbs[db] = db_path
        params['databases'] = valid_dbs
        defaults['databases'] = '<list of databases>'
    return defaults


def check_bwa(self, params, no_database=False):
    defaults = {
        'paired': [True, False],
    }
    if not no_database:
        dbs_existing = check_databases(self.name, params, self.databases)
        valid_dbs = {}
        for db in dbs_existing:
            if 'bwa' in self.databases.builds[db]:
                db_path = '%s/*' % self.databases.builds[db]['bwa']
                if not self.config.dev or self.config.workshop:
                    bt2_paths = glob.glob(db_path)
                    if bt2_paths:
                        valid_dbs[db] = bt2_paths[0]
                else:
                    valid_dbs[db] = db_path
        params['databases'] = valid_dbs
        defaults['databases'] = '<list of databases>'
    return defaults


def check_spades(self, params):
    defaults = {
        'hybrid': ['illumina', 'nanopore', 'pacbio'],
        'k': ['auto'],
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
    if 'hybrid' in params:
        invalids = [x for x in params['hybrid'] if x not in defaults['hybrid']]
        if invalids:
            sys.exit('[%s] Param "hybrid" has invalid terms: %s' % (
                self.name, ','.join(invalids)))
    else:
        params['hybrid'] = ['illumina']
    if 'k' not in params:
        params['k'] = defaults['k']
    else:
        if isinstance(params['k'], str) and params['k'] == 'auto':
            params['k'] = ['auto']
        elif isinstance(params['k'], list):
            kerr = [x for x in params['k'] if not str(x).isdigit()]
            if len(kerr):
                sys.exit('[%s] "k" must be integers (%s)' % (
                    self.name, ','.join(kerr)))
    check_default(self, params, defaults, ['k', 'hybrid'])
    return defaults


def check_viralverify(self, params):
    defaults = {'thr': 7, 'p': [False, True], 'db': [False, True]}
    check_nums(self, params, defaults, ['thr'], int, 'viralverify')
    check_default(self, params, defaults, ['thr'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<path to the ViralVerify installation folder>'
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


def check_simka(self, params):
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
    ints = ['nb_kmers', 'min_read_size']
    check_nums(self, params, defaults, ints, int)
    defs = ['kmer', 'log_reads'] + ints
    check_default(self, params, defaults, defs)
    params['kmer'] = [int(x) for x in params['kmer']]
    params['log_reads'] = [int(x) for x in params['log_reads']]
    defaults['kmer'] = [str(x) for x in np.linspace(15, 80, 6)]
    defaults['log_reads'] = [str(x) for x in np.logspace(3, 7, 3)]
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<path to the SimKa installation folder>'
    return defaults


def check_metamarker(self, params):
    defaults = {'identity': 0.9}
    check_nums(self, params, defaults, ['identity'], float, 0, 1)
    if 'groups' not in params:
        sys.exit('[metamarker] Param "groups" not found (metadata variable(s))')
    for group in params['groups']:
        if group not in set(self.config.meta.columns):
            sys.exit('[metamarker] metadata variable "%s" not found' % group)
    defaults['groups'] = '<List of metadata columns>'
    return defaults


def check_metawrap(self, params):
    defaults = {
        'min_completion': 25,
        'max_contamination': 5,
        'min_contig_length': 1000,
        'min_completion_reassembly': 25,
        'max_contamination_reassembly': 5,
        'binners': ['maxbin2', 'metabat2', 'concoct'],
        'reassembly': ['permissive', 'strict'],
        'blobology': ['coassembly', 'sample']
    }
    if 'binners' not in params:
        params['binners'] = defaults['binners']
    mins = ['min_completion', 'min_completion_reassembly',
            'max_contamination', 'max_contamination_reassembly']
    check_nums(self, params, defaults, mins, int, 0, 100)
    check_nums(self, params, defaults, ['min_contig_length'], int)
    mins.append('min_contig_length')
    check_default(
        self, params, defaults, mins, ['binners', 'reassembly', 'blobology'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<path to the metaWRAP installation folder>'
    return defaults


def check_drep(self, params):
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
    if 'S_algorithm' not in params:
        params['S_algorithm'] = ['fastANI', 'ANIn']
    ints = ['MASH_sketch', 'primary_chunksize']
    check_nums(self, params, defaults, ints, int)
    flts = ['P_ani', 'S_ani', 'cov_thresh', 'warn_dist', 'warn_sim', 'warn_aln']
    check_nums(self, params, defaults, flts, float, 0, 1)
    check_default(self, params, defaults, (ints + flts), ['S_algorithm'])
    check_binary(self, params, defaults, 'anicalculator')
    return defaults


def check_checkm(self, params):
    defaults = {
        'multi': 10,
        'unique': 10,
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
        'force_domain': [False, True],
        'no_refinement': [False, True],
    }
    ints = ['min_seq_len', 'min_qc', 'multi', 'unique']
    check_nums(self, params, defaults, ints, int)
    floats = ['min_align', 'max_edit_dist', 'aai_strain', 'length']
    check_nums(self, params, defaults, floats, float, 0, 1)
    floats_ = ['e_value']
    check_nums(self, params, defaults, floats_, float)
    check_default(self, params, defaults, (ints + floats + floats_))
    return defaults


def check_checkm2(self, params):
    defaults = {
        'lowmem': [False, True],
        'general': [False, True],
        'specific': [False, True],
        'allmodels': [True, False],
        'genes': [False, True],
        'force': [True, False],
        'dbg_cos': [False, True],
        'dbg_vectors': [False, True]
    }
    if 'database_path' in params:
        if not (self.config.dev or self.config.workshop):
            if not isdir(params['database_path']):
                sys.exit('[checkm2] Param "database_path" do not exist')
    else:
        params['database_path'] = None
    check_default(self, params, defaults, ['database_path'])
    defaults['database_path'] = '<path to an alternative reference database>'
    return defaults


def check_prodigal(self, params):
    defaults = {
        'f': ['gff', 'gbk', 'sco'],
        'p': ['meta', 'single'],
        'c': [False, True],
        'm': [False, True],
        'n': [False, True],
        'q': [False, True],
    }
    check_default(self, params, defaults)
    return defaults


def check_fastp(self, params):
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
        'overlap_len_require', 'overlap_diff_limit',
        'filter_by_index_threshold', 'length_required', 'length_limit',
        'average_qual', 'n_base_limit', 'qualified_quality_phred',
        'cut_window_size', 'cut_front_window_size', 'cut_tail_window_size',
        'cut_right_window_size', 'poly_g_min_len', 'poly_x_min_len',
        'trim_front1', 'trim_tail1', 'max_len1', 'trim_front2', 'trim_tail2',
        'max_len2', 'reads_to_process']
    for adapter in ['adapter_sequence', 'adapter_sequence_r2', 'adapter_fasta']:
        if adapter not in params:
            params[adapter] = None
    check_nums(self, params, defaults, int1, int)
    int2 = ['split']
    check_nums(self, params, defaults, int2, int, 2, 999)
    int3 = ['split_prefix_digits']
    check_nums(self, params, defaults, int3, int, 1, 10)
    int4 = ['overrepresentation_sampling']
    check_nums(self, params, defaults, int4, int, 1, 10000)
    int5 = ['complexity_threshold', 'unqualified_percent_limit']
    check_nums(self, params, defaults, int5, int, 1, 100)
    int6 = ['cut_mean_quality', 'cut_front_mean_quality',
            'cut_tail_mean_quality', 'cut_right_mean_quality']
    check_nums(self, params, defaults, int6, int, 1, 36)
    int7 = ['dup_calc_accuracy']
    check_nums(self, params, defaults, int7, int, 1, 6)
    int8 = ['compression']
    check_nums(self, params, defaults, int8, int, 1, 9)
    let_go = int1 + int2 + int3 + int4 + int5 + int6 + int7 + int8
    check_default(self, params, defaults, let_go)
    defaults['adapter_sequence'] = '<the adapter for read1>'
    defaults['adapter_sequence_r2'] = '<the adapter for read2>'
    defaults['adapter_fasta'] = '<Path to fasta file with sequences to trim>'
    return defaults


def check_humann(self, params):
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

    if not self.config.dev or self.config.workshop:
        if 'nucleotide_db' not in params or not isdir(params['nucleotide_db']):
            sys.exit('[humann] Param "nucleotide_db" must be an existing path')
        if 'protein_db' not in params or not isdir(params['protein_db']):
            sys.exit('[humann] Param "protein_db" must be an existing path')
    if 'profiles' in params:
        if not isinstance(params['profiles'], dict):
            sys.exit('[humann] Param "profiles" must be a dict structure')
        for key, value in params['profiles'].items():
            if self.config.dev or self.config.workshop:
                continue
            if not isinstance(value, str):
                sys.exit('[humann] Param "profiles::%s" must be a string' % key)
            elif not isfile(value):
                sys.exit('[humann] Param "profiles::%s" does not exist' % key)
    check_default(self, params, defaults, let_go)
    defaults['profiles'] = '<dict of taxonomic profiles name: file>'
    return defaults


def check_midas(self, params):
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
            if not (self.config.dev or self.config.workshop) and not isfile(v):
                sys.exit('[midas] Param "focus::%s::%s" not a file' % (k, v))
    ints = ['species_topn', 'readq', 'trim', 'n', 'word_size']
    check_nums(self, params, defaults, ints, int)
    floats = ['mapid', 'aln_cov']
    check_nums(self, params, defaults, floats, float, 0, 1)
    floats_ = ['species_cov']
    check_nums(self, params, defaults, floats_, float)
    defs = ['tracking', 'focus']
    check_default(self, params, defaults, (ints + floats + floats_ + defs))
    check_binary(self, params, defaults, 'path')
    defaults['focus'] = '<dict of key:value pair some_name: /path/to/spc.txt'
    defaults['tracking'] = '<list of metadata columns>'
    defaults['path'] = '<path to the MIDAS installation folder>'
    return defaults


def check_macsyfinder(self, params):
    defaults = {
        'models': [
            'CasFinder', 'CONJScan_plasmids', 'TFFscan', 'TFF-SF', 'TXSScan'],
        'e_value_search': 0.1,
        'i_evalue_sel': 0.001,
        'coverage_profile': 0.5,
        'mandatory_weight': 1.0,
        'accessory_weight': 0.5,
        'exchangeable_weight': 0.8,
        'redundancy_penalty': 1.5,
        'out_of_cluster': 0.7
    }
    if 'models' not in params:
        params['models'] = defaults['models']
    else:
        if not isinstance(params['models'], list):
            sys.exit('[models] Param "models" not a list')
        elif not set(params['models']).issubset(defaults['models']):
            sys.exit('[models] Param "models" not in: %s' % params['models'])
    flo1 = ['e_value_search', 'i_evalue_sel', 'coverage_profile',
            'mandatory_weight', 'accessory_weight', 'exchangeable_weight',
            'out_of_cluster']
    check_nums(self, params, defaults, flo1, float, 0, 1)
    flo2 = ['redundancy_penalty']
    check_nums(self, params, defaults, flo2, float)
    check_default(self, params, defaults, (flo1 + flo2 + ['models']))
    return defaults


def check_coconet(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['test_ratio', 'learning_rate', 'theta', 'gamma1', 'gamma2']
    check_nums(self, params, defaults, floats, float, 0, 1)
    let_go = ['tlen_range', 'compo_neurons', 'cover_neurons']
    check_default(self, params, defaults, (ints + floats + let_go))
    return defaults


def check_tiara(self, params):
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
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, (ints + [fas, prob]))
    return defaults


def check_bracken(self, params):
    defaults = {
        'read_len': 100,
        'level': ['S', 'D', 'P', 'C', 'O', 'F', 'G', 'S1'],
        'threshold': 0
    }
    check_nums(self, params, defaults, ['read_len', 'threshold'], int)
    check_default(self, params, defaults, ['read_len', 'threshold'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<path to the Bracken installation folder>'
    return defaults


def check_plasmidfinder(self, params):
    defaults = {
        'extented_output': [True, False],
        'threshold': 70,
        'mincov': 75,
    }
    if 'methodPath' not in params:
        sys.exit('[plasmidfinder] Param "methodPath" missing (kma or blastn '
                 'binary)')
    elif not (self.config.dev or self.config.workshop):
        if not isfile(params['methodPath']):
            sys.exit('[plasmidfinder] "methodPath" binary not found'
                     ': "%s"' % params['methodPath'])

    db_path = 'databasePath'
    if db_path not in params:
        sys.exit('[plasmidfinder] Param "databasePath" missing (databases dir)')
    elif not (self.config.dev or self.config.workshop):
        if not isdir(params[db_path]):
            sys.exit('[plasmidfinder] no "databasePath" %s' % params[db_path])

    if 'databases' in params:
        if not isinstance(params['databases'], str):
            dbs = params['databases']
            sys.exit('[plasmidfinder] Param "databases" not a string: %s' % dbs)

    ints = ['threshold', 'mincov']
    check_nums(self, params, defaults, ints, int, 0, 100)
    check_default(self, params, defaults, ints)
    check_binary(self, params, defaults, 'binary')
    defaults['methodPath'] = '<Path to the "kma" or "blastn" binary>'
    defaults['binary'] = '<Path to the "plasmidfinder.py" binary>'
    defaults['databasePath'] = '<Path databases folder (with a "config" file)>'
    defaults['databases'] = "<Comma-separated databases (first field of the " \
                            "databasePath's 'config' file)>"
    return defaults


def check_plasforest(self, params):
    defaults = {
        'size_of_batch': 0,
        'b': [True, False],
        'f': [True, False],
        'r': [False, True]
    }
    check_nums(self, params, defaults, ['size_of_batch'], int)
    check_default(self, params, defaults, ['size_of_batch'])
    check_binary(self, params, defaults, 'binary')
    defaults['binary'] = '<Path to the "PlasForest.py" binary>'
    return defaults


def check_flye(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['read_error']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    defaults['stop_after'] = '<Stage name that is completed and to stop after>'
    defaults['genome_size'] = '<estimated genome size (e.g., 5m or 2.6g)>'
    return defaults


def check_canu(self, params):
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
    check_nums(self, params, defaults, ints, int)
    ints_ = ['corMinCoverage', 'corOutCoverage']
    check_nums(self, params, defaults, ints_, int, 0, 4)
    floats = ['correctedErrorRate', 'rawErrorRate']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + ints_ + floats))
    check_binary(self, params, defaults, 'path')
    defaults['genome_size'] = '<estimated genome size (e.g., 5m or 2.6g)>'
    defaults['specifications'] = '<Path to assembly option specifications file>'
    return defaults


def check_unicycler(self, params):
    defaults = {
        'hybrid': [None, 'pacbio', 'nanopore'],
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
    check_nums(self, params, defaults, ints, int)
    floats = ['min_bridge_qual']
    check_nums(self, params, defaults, floats, float)
    floats_ = ['start_gene_cov', 'start_gene_id']
    check_nums(self, params, defaults, floats_, float, 0, 100)
    floats__ = ['max_kmer_frac', 'min_kmer_frac', 'depth_filter']
    check_nums(self, params, defaults, floats__, float, 0, 1)
    paths = [x for x in defaults if '_path' in x]
    check_default(self, params, defaults,
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


def check_ngmerge(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['p']
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<absolute path to the folder containing "NGmerge">'
    return defaults


def check_flash(self, params):
    defaults = {
        'min_overlap': 10,
        'max_overlap': 65,
        'mismatch': 0.25
    }
    ints = ['min_overlap', 'max_overlap']
    check_nums(self, params, defaults, ints, int)
    floats = ['mismatch']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'path')
    return defaults


def check_pear(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['p_value', 'max_uncalled_base']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    return defaults


def check_bbmerge(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['pfilter', 'maxratio', 'ratiomargin',
              'ratiooffset', 'minsecondratio']
    check_nums(self, params, defaults, floats, float)
    floats_ = ['prealloc', 'minprob']
    check_nums(self, params, defaults, floats_, float, 0, 1)
    check_default(self, params, defaults, (ints + floats + floats_))
    return defaults


def check_necat(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = [x for x in defaults if x[-2:] in ['_e', '_p']]
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'path')
    defaults['genome_size'] = '<estimated genome size (in bases)>'
    return defaults


def check_pooling(self, params):
    defaults = {'pool_single_and_merged': [True, False]}
    check_default(self, params, defaults)
    return defaults


def check_cutadapt(self, params):
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
    check_nums(self, params, defaults, ints, int)
    floats = ['error_rate']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    return defaults


def check_megahit(self, params):
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
        if [v for v in vals if not v.isdigit() or int(v) < 15 or int(v) > 255]:
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
    check_nums(self, params, defaults, ints, int)
    floats = ['memory', 'low_local_ratio', 'disconnect_ratio']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (let_go + ints + floats))
    return defaults


def check_plass(self, params):
    defaults = {
        'type': ['assemble', 'nuclassemble'],
        'add_self_matches': [False, True],
        'alph_size': 13,
        'spaced_kmer_mode': [False, True],
        'mask': [False, True],
        'mask_lower_case': [False, True],
        'k': 14,
        'split_memory_limit': 0,
        'wrapped_scoring': [False, True],
        'e': 0.,
        'c': 0.,
        'a': 0,
        'cov_mode': [3, 0, 1, 2, 4, 5],
        'min_seq_id': 0.9,
        'min_aln_len': 0,
        'seq_id_mode': [0, 1, 2],
        'kmer_per_seq': 60,
        'kmer_per_seq_scale': 'nucl:0.200,aa:0.000',
        'adjust_kmer_len': [False, True],
        'hash_shift': 67,
        'include_only_extendable': [True, False],
        'ignore_multi_kmer': [True, False],
        'num_iterations': 12,
        'rescore_mode': [3, 0, 1, 2, 4],
        'min_length': 45.,
        'max_length': 32734,
        'max_gaps': 2147483647,
        'contig-start-mode': [2, 1, 0],
        'contig_end_mode': [2, 1, 0],
        'orf_start_mode': [1, 0, 2],
        'forward_frames': '1,2,3',
        'reverse_frames': '1,2,3',
        'translation_table': list(range(1, 32)),
        'translate': [False, True],
        'use_all_table_starts': [False, True],
        'id_offset': 0,
        'protein_filter_threshold': 0.2,
        'filter_proteins': [True, False],
        'dbtype': [0, 1, 2],
        'shuffle': [True, False],
        'createdb_mode': [False, True],
        'db_load_mode': [0, 1, 2, 3],
        'compressed': [False, True],
        'v': [0, 1, 2, 3],
        'max_seq_len': 65535,
        'delete_tmp_inc': [False, True],
        'remove_tmp_files': [False, True],
        'filter_hits': [False, True],
        'sort_results': [False, True],
        'create_lookup': [False, True],
        'write_lookup': [True, False],
    }
    p = 'kmer_per_seq_scale'
    if p not in params:
        params[p] = defaults[p]
    else:
        if 'nucl:' not in params[p] or ',aa:' not in params[p]:
            sys.exit('[plass] Params "%s" invalid format (nucl:#,aa:#)' % p)
        vals = params[p].split(',aa:')
        try:
            float(vals[-1])
            float(vals[0].split('nucl:')[-1])
        except ValueError:
            sys.exit('[plass] Params "%s" invalid format (nucl:#,aa:#)' % p)

    for p in ['forward_frames', 'reverse_frames']:
        if p not in params:
            params[p] = defaults[p]
        else:
            if not isinstance(params[p], str):
                sys.exit('[plass] Prams "%s" invalid format (not a string)' % p)
            vals = params[p].split(',')
            if [v for v in vals if not v.isdigit() or int(v) not in [1, 2, 3]]:
                sys.exit('[plass] Prams "%s" not ","-separated "1","2","3"' % p)
    let_go = ['kmer_per_seq_scale', 'forward_frames', 'reverse_frames']
    ints = ['a', 'split_memory_limit', 'k', 'min_aln_len', 'kmer_per_seq',
            'hash_shift', 'min_length', 'max_length', 'max_gaps',
            'id_offset', 'max_seq_len']
    check_nums(self, params, defaults, ints, int)
    int1 = ['num_iterations']
    check_nums(self, params, defaults, int1, int, 1, 1000000000)
    int2 = ['alph_size']
    check_nums(self, params, defaults, int2, int, 2, 21)
    floats = ['e', 'c', 'min_seq_id', 'protein_filter_threshold']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(
        self, params, defaults, (ints + int1 + int2 + floats + let_go))
    return defaults


def check_pirate(self, params):
    defaults = {
        'steps': '50,60,70,80,90,95,98',
        'features':	['CDS', 'mRNA'],
        'nucl': [False, True],
        'align': [False, True],
        'rplots': [False, True],
        'diamond': [False, True],
        'para_off': [False, True],
        'hsp_prop': [False, True],
        'cd_core_off': [True, False],
        'z': [1, 2, 0],
        'perc': 98,
        'cd_low': 98,
        'cd_step': 0.5,
        'evalue': 1e-6,
        'flat': 1.5,
    }
    p = 'steps'
    if p not in params:
        params[p] = defaults[p]
    else:
        if not isinstance(params[p], str):
            sys.exit('[plass] Prams "%s" invalid format (not a string)' % p)
        vals = params[p].split(',')
        if [v for v in vals if not v.isdigit() or int(v) not in range(101)]:
            sys.exit('[plass] Prams "%s" not ","-separated integers' % p)
    int1 = ['perc', 'cd_low']
    check_nums(self, params, defaults, int1, int, 0, 100)
    flo1 = ['flat']
    check_nums(self, params, defaults, flo1, float)
    flo2 = ['cd_step', 'evalue']
    check_nums(self, params, defaults, flo2, float, 0, 1)
    check_default(self, params, defaults, (int1 + flo1 + flo2))
    return defaults


def check_deeparg(self, params):
    defaults = {
        'min_prob': 0.8,
        'arg_alignment_overlap': 0.8,
        'arg_alignment_evalue': 1e-10,
        'arg_alignment_identity': 50,
        'arg_num_alignments_per_entry': 1000,
        'model_version': ['v2'],
        'deeparg_identity': 80,
        'deeparg_probability': 0.8,
        'deeparg_evalue': 1e-10,
        'gene_coverage': 1,
        'bowtie_16s_identity': 0.8
    }
    ints = ['arg_num_alignments_per_entry']
    check_nums(self, params, defaults, ints, int)
    ints1 = ['arg_alignment_identity', 'deeparg_identity']
    check_nums(self, params, defaults, ints1, int, 0, 100)
    floats = ['min_prob', 'arg_alignment_overlap', 'arg_alignment_evalue',
              'deeparg_probability', 'deeparg_evalue', 'gene_coverage',
              'bowtie_16s_identity']
    check_nums(self, params, defaults, floats, float, 0, 1)
    let_go = ints + ints1 + floats
    check_default(self, params, defaults, let_go)
    defaults['model'] = '<Model to use (SS: for reads, LS: for genes)>'
    db = 'db_dir'
    if db not in params or (not (self.config.dev or self.config.workshop) and not isdir(params[db])):
        sys.exit('[%s] Params "%s" must be an existing path' % (self.name, db))
    defaults[db] = '<Path to the installed deepARG database folder>'

    return defaults


def check_lorikeet(self, params):
    defaults = {
        'mapper': [
            "minimap2-sr", "bwa-mem", "bwa-mem2", "ngmlr-ont", "ngmlr-pb",
            "minimap2-ont", "minimap2-pb", "minimap2-hifi", "minimap2-no-preset"
        ],
        'longread_mapper': [
            "minimap2-ont", "minimap2-pb", "minimap2-hifi",
            "minimap2-no-preset", "ngmlr-ont", "ngmlr-pb"],
        'minimap2_params': [""],
        'bwa_params': [""],
        'ngmlr_params': [""],
        'parallel_genomes': 4,
        'min_read_aligned_length': 0,
        'min_read_aligned_length_pair': 0,
        'contig_end_exclusion': 0,
        'ploidy': 1,
        'qual_by_depth_filter': 25,
        'qual_threshold': 150,
        'depth_per_sample_filter': 5,
        'min_long_read_size': 1500,
        'min_long_read_average_base_qual': 20,
        'min_base_quality': 10,
        'min_mapq': 20,
        'base_quality_score_threshold': 18,
        'max_input_depth': 200000,
        'min_contig_size': 2500,
        'min_sv_qual': 3,
        'phred_scaled_global_read_mismapping_rate': 45,
        'pair_hmm_gap_continuation_penalty': 10,
        'min_assembly_region_size': 50,
        'max_assembly_region_size': 300,
        'assembly_region_padding': 100,
        'min_dangling_branch_length': 4,
        'min_prune_factor': 2,
        'num_pruning_samples': 1,
        'max_unpruned_variants': 100,
        'max_prob_propagation_distance': 50,
        'max_mnp_distance': 0,
        'min_read_percent_identity': 0,
        'min_read_aligned_percent': 0,
        'min_read_percent_identity_pair': 0,
        'min_read_aligned_percent_pair': 0,
        'trim_min': 0.00,
        'trim_max': 1.00,
        'heterozygosity': 0.001,
        'heterozygosity_stdev': 0.01,
        'indel_heterozygosity': 0.000125,
        'standard_min_confidence_threshold_for_calling': 30.0,
        'active_probability_threshold': 0.002,
        'initial_error_rate_for_pruning': 0.001,
        'pruning_log_odds_threshold': 1.0,
        'pcr_indel_model': ["conservative"],
        'quiet': [True, False],
        'verbose': [False, True],
        'sharded': [False, True],
        'split_bams': [False, True],
        'calculate_fst': [False, True],
        'calculate_dnds': [False, True],
        'do_not_call_svs': [False, True],
        'include_secondary': [False, True],
        'proper_pairs_only': [False, True],
        'exclude_supplementary': [False, True],
        'dont_use_soft_clipped_bases': [False, True],
        'minimap2_reference_is_index': [False, True],
        'use_posteriors_to_calculate_qual': [False, True],
        'annotate_with_num_discovered_alleles': [False, True],
        'dont_increase_kmer_sizes_for_cycles': [False, True],
        'allow_non_unique_kmers_in_ref': [False, True],
        'recover_all_dangling_branches': [False, True],
        'do_not_run_physical_phasing': [False, True],
        'disable_optimizations': [False, True],
        'use_adaptive_pruning': [False, True],
        'discard_unmapped': [False, True],
        'disable_avx': [False, True],
        'force': [True, False],
        'kmer_sizes': "10 25",
        'limiting_interval': "",
        'features_vcf': [None],
        'reference': [None],
        'genome_fasta_directory': [None],
    }
    if 'kmer_sizes' not in params:
        params['kmer_sizes'] = defaults['kmer_sizes']
    elif [x for x in params['kmer_sizes'].split() if not x.isdigit()]:
        sys.exit("[lorikeet] Param 'kmer_sizes' should be INTs only")

    if 'limiting_interval' not in params:
        params['limiting_interval'] = defaults['limiting_interval']
    elif params['limiting_interval'].count('-') != 1:
        sys.exit("[lorikeet] Param 'limiting_interval' should contain 1 hyphen")
    elif [x for x in params['limiting_interval'].split('-') if not x.isdigit()]:
        sys.exit("[lorikeet] Param 'limiting_interval' should be INT-INT")

    ints = ['parallel_genomes', 'min_read_aligned_length',
            'min_read_aligned_length_pair', 'contig_end_exclusion', 'ploidy',
            'qual_by_depth_filter', 'qual_threshold', 'depth_per_sample_filter',
            'min_long_read_size', 'min_long_read_average_base_qual',
            'min_base_quality', 'min_mapq', 'base_quality_score_threshold',
            'max_input_depth', 'min_contig_size', 'min_sv_qual',
            'phred_scaled_global_read_mismapping_rate',
            'pair_hmm_gap_continuation_penalty', 'min_assembly_region_size',
            'max_assembly_region_size', 'assembly_region_padding',
            'min_dangling_branch_length', 'min_prune_factor',
            'num_pruning_samples', 'max_unpruned_variants',
            'max_prob_propagation_distance', 'max_mnp_distance']
    check_nums(self, params, defaults, ints, int)
    flo = ['heterozygosity', 'heterozygosity_stdev', 'indel_heterozygosity',
           'standard_min_confidence_threshold_for_calling',
           'initial_error_rate_for_pruning']
    check_nums(self, params, defaults, flo, float)
    flo1 = ['min_read_percent_identity', 'min_read_aligned_percent',
            'min_read_percent_identity_pair', 'min_read_aligned_percent_pair']
    check_nums(self, params, defaults, flo1, float, 0, 100)
    flo2 = ['trim_min', 'trim_max', 'active_probability_threshold',
            'pruning_log_odds_threshold']
    check_nums(self, params, defaults, flo2, float, 0, 1)
    check_default(self, params, defaults, (
            ints + flo + flo1 + flo2 + ['kmer_sizes', 'limiting_interval']))
    return defaults


def check_metamarc(self, params):
    defaults = {
        'coverage': 80,
        'dedup': [False, True],
        'evalue': 10,
        'kmer': 0,
        'level': ['1', '2', '3'],
        'multicorrect': [False, True]
    }
    if 'level' not in params:
        params['level'] = defaults['level']
    int1 = ['evalue', 'kmer']
    check_nums(self, params, defaults, int1, int)
    int2 = ['coverage']
    check_nums(self, params, defaults, int2, int, 0, 100)
    check_default(self, params, defaults, (int1 + int2), ['level'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the meta-marc software folder (from github)>'
    return defaults


def check_hifiadapterfilt(self, params):
    defaults = {'l': 44, 'm': 97}
    check_nums(self, params, defaults, ['l'], int)
    check_nums(self, params, defaults, ['m'], float, 0, 100)
    check_default(self, params, defaults, ['l', 'm'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the folder containing "pbadapterfilt.sh">'
    return defaults


def check_karga(self, params):
    defaults = {
        'k': 17,
        'r': ['y', 'n', 'yes', 'no'],
        'm': ['y', 'n', 'yes', 'no'],
        'i': 25000,
        's': 12345
    }
    if 'databases' not in params:
        sys.exit('[karga] Param "databases" missing (name: /path/to/db.fasta')
    check_nums(self, params, defaults, ['k', 'i', 's'], int)
    check_default(self, params, defaults, ['k', 'i', 's'])
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the folder containing "KARGA.class">'
    defaults['databases'] = '<Paths to the ARG/MGE fasta files (with ' \
                            'resistance annotation in header) per name (dict)>'
    return defaults


def check_kargva(self, params):
    defaults = {
        'k': 9,
        'm': ['y', 'n', 'yes', 'no'],
        'i': 125000,
    }
    check_nums(self, params, defaults, ['k', 'i'], int)
    check_default(self, params, defaults, ['k', 'i'])
    check_binary(self, params, defaults, 'path')
    if 'databases' not in params:
        params['databases'] = {
            'default': '%s/kargva_db_v5.fasta' % params['path']}
    defaults['path'] = '<Path to the folder containing "KARGA.class">'
    defaults['databases'] = '<Paths to the ARG/MGE fasta files (with ' \
                            'resistance mutations in header) per name (dict)>'
    return defaults


def check_abricate(self, params):
    defaults = {
        'setupdb': [False, True],
        'noheader': [False, True],
        'csv': [False, True],
        'nopath': [False, True],
        'minid': 80,
        'mincov': 80,
        'databases': ['vfdb', 'card', 'ecoh', 'resfinder', 'ncbi', 'megares',
                      'argannot', 'plasmidfinder', 'ecoli_vf']
    }
    if 'databases' not in params:
        params['databases'] = defaults['databases']
    ints = ['minid', 'mincov']
    check_nums(self, params, defaults, ints, int, 0, 100)
    check_default(self, params, defaults, ints, ['databases'])
    return defaults


def check_amrplusplus2(self, params):
    defaults = {
        'leading': 10,
        'trailing': 3,
        'slidingwindow': 4,
        'minlen': 36,
        'threshold': 80,
        'min': 5,
        'max': 100,
        'skip': 5,
        'samples': 1
    }
    paths = {
        'adapters': "data/adapters/nextera.fa",
        'fqc_adapters': "data/adapters/nextera.tab",
        'host_index': "",
        'host': "data/host/chr21.fasta.gz",
        'kraken_db': "",
        'amr_index': "",
        'amr': "data/amr/megares_database_v1.02.fasta",
        'annotation': "data/amr/megares_annotations_v1.02.csv",
        'snp_annotation': "data/amr/snp_location_metadata.csv",
        'snp_confirmation': "bin/snp_confirmation.py"
    }
    for param, filename in paths.items():
        if param not in params:
            params[param] = filename
    ints0 = ['leading', 'trailing', 'minlen', 'min', 'max', 'skip', 'samples']
    check_nums(self, params, defaults, ints0, int)
    ints1 = ['slidingwindow']
    check_nums(self, params, defaults, ints1, int, 4, 15)
    ints2 = ['threshold']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    check_default(self, params, defaults,
                  (ints0 + ints1 + ints2 + list(paths.keys())))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to folder containing "main_AmrPlusPlus_v2.nf">'
    defaults['adapters'] = '<Path to adapter sequences file>'
    defaults['fqc_adapters'] = '<Path to tab delimited adapter sequences file>'
    defaults['host_index'] = '<Path to host genome index files>'
    defaults['host'] = '<Path to host genome file>'
    defaults['kraken_db'] = '<Path to Kraken database>'
    defaults['amr_index'] = '<Path to amr index files>'
    defaults['amr'] = '<Path to antimicrobial resistance (MEGARes) database>'
    defaults['annotation'] = '<Path to amr annotation file>'
    defaults['snp_annotation'] = '<Path to SNP metadata file>'
    defaults['snp_confirmation'] = '<Path to "snp_confirmation.py" script>'
    return defaults


def check_metaphlan(self, params):
    defaults = {
        'tax_lev': ['k', 'p', 'c', 'o', 'f', 'g', 's', 'a'],
        'bt2_ps': ['very-sensitive', 'sensitive', 'sensitive-local',
                   'very-sensitive-local'],
        'stat': ['tavg_g', 'avg_g', 'avg_l', 'tavg_l', 'wavg_g', 'wavg_l',
                 'med'],
        't': ['rel_ab', 'rel_ab_w_read_stats', 'clade_profiles',
              'marker_ab_table', 'marker_pres_table',
              'clade_specific_strain_tracker'],
        'perc_nonzero': 0.33,
        'stat_q': 0.2,
        'pres_th': 0.,
        'min_ab': 0.,
        'min_mapq_val': 5,
        'min_cu_len': 2000,
        'read_min_len': 70,
        'min_alignment_len': 0,
        'ignore_eukaryotes': [False, True],
        'ignore_bacteria': [False, True],
        'ignore_archaea': [False, True],
        'add_viruses': [False, True],
        'avoid_disqm': [False, True]
    }
    if 't' not in params:
        params['t'] = ['rel_ab', 'rel_ab_w_read_stats', 'clade_profiles',
                       'marker_ab_table', 'marker_pres_table',
                       'clade_specific_strain_tracker']
    if 'tax_lev' not in params:
        params['tax_lev'] = ['a']
    if 'ignore_markers' in params and not isfile(params['ignore_markers']):
        sys.exit('[metaphlan] Param "ignore_markers" must be a file path')
    ints = ['min_mapq_val', 'min_cu_len', 'min_alignment_len', 'read_min_len']
    check_nums(self, params, defaults, ints, int)
    floats = ['stat_q', 'perc_nonzero', 'pres_th', 'min_ab']
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats), ['tax_lev', 't'])
    defaults['ignore_markers'] = '<File containing a list of markers to ignore>'
    return defaults


def check_mapdamage2(self, params):
    defaults = {
        'downsample': [None],
        'downsample_seed': 1234,
        'merge_reference_sequences': [False, True],
        'length': 70,
        'around': 10,
        'min_basequal': 0,
        'fasta': [False, True],
        'plot_only': [False, True],
        'quiet': [False, True],
        'verbose': [False, True],
        'mapdamage_modules': [False, True],
        'ymax': 0.3,
        'readplot': 25,
        'refplot': 10,
        'rand': 30,
        'burn': 10000,
        'adjust': 10,
        'iter': 50000,
        'forward': [False, True],
        'reverse': [False, True],
        'var_disp': [False, True],
        'jukes_cantor': [False, True],
        'diff_hangs': [False, True],
        'fix_nicks': [False, True],
        'use_raw_nick_freq': [False, True],
        'single_stranded': [False, True],
        'theme_bw': [False, True],
        'seq_length': 12,
        'stats_only': [False, True],
        'no_stats': [True, False],
        'check_R_packages': [False, True],
        'rescale': [False, True],
        'rescale_only': [False, True],
        'all_reads': [False, True],
        'rescale_length_5p': 12,
        'rescale_length_3p': 12
    }
    if 'downsample' not in params:
        params['downsample'] = None
    else:
        try:
            downsample = float(params['downsample'])
            if downsample > 1:
                params['downsample'] = np.ceil(downsample)
            elif downsample > 0:
                params['downsample'] = downsample
            else:
                sys.exit('[mapdamage2] Param "downsample" must be >0')
        except:
            sys.exit('[mapdamage2] Param "downsample" must be numeric')
    ints = ['downsample_seed', 'length', 'around', 'min_basequal', 'readplot',
            'refplot', 'rand', 'burn', 'adjust', 'iter', 'seq_length',
            'rescale_length_5p', 'rescale_length_3p']
    check_nums(self, params, defaults, ints, int)
    floats = ['ymax']
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats + ['downsample']))
    return defaults


def check_amrfinderplus(self, params):
    defaults = {
        'database': [None],
        'annotation_format': [None, 'bakta', 'genbank', 'microscope', 'patric',
                              'pgap', 'prokka', 'pseudomonasdb', 'rast'],
        'pgap': [False, True],
        'ident_min': -1,
        'coverage_min': 0.5,
        'organism': [None, 'Acinetobacter_baumannii', 'Burkholderia_cepacia',
                     'Burkholderia_pseudomallei', 'Campylobacter',
                     'Clostridioides_difficile', 'Enterococcus_faecalis',
                     'Enterococcus_faecium', 'Escherichia', 'Klebsiella',
                     'Neisseria', 'Pseudomonas_aeruginosa', 'Salmonella',
                     'Staphylococcus_aureus', 'Staphylococcus_pseudintermedius',
                     'Streptococcus_agalactiae', 'Streptococcus_pneumoniae',
                     'Streptococcus_pyogenes', 'Vibrio_cholerae'],
        'translation_table': 11,
        'plus': [False, True],
        'report_common': [False, True],
        'blast_bin': [None],
        'report_all_equal': [False, True],
        'nucleotide_flank5_size': 0,
        'quiet': [False, True],
        'gpipe_org': [False, True],
        'parm': [None, 'nosame', 'noblast', 'skip_hmm_check', 'bed'],
    }
    if 'database' not in params:
        params['database'] = '$AMRFINDER_DB'
    if 'blast_bin' not in params:
        params['blast_bin'] = '$BLAST_BIN'
    if 'parm' not in params:
        params['parm'] = [None]
    int1 = ['translation_table']
    check_nums(self, params, defaults, int1, int, 1, 33)
    int2 = ['nucleotide_flank5_size']
    check_nums(self, params, defaults, int2, int)
    floats = ['ident_min', 'coverage_min']
    check_nums(self, params, defaults, floats, float, -1, 1)
    check_default(self, params, defaults,
                  (int1 + int2 + floats + ['database', 'blast_bin']), ['parm'])
    return defaults


def check_squeezemeta(self, params):
    defaults = {
        'm': ['sequential', 'coassembly', 'merged'],
        'restart': [False, True],
        'step': [None],
        'force_overwrite': [False, True],
        'cleaning': [False, True],
        'cleaning_options': [None],
        'a': ['megahit', 'spades', 'rnaspades', 'canu', 'flye'],
        'assembly_options': [None],
        'contiglen': 200,
        'extassembly': [None],
        'singletons': [False, True],
        'contigid': [None],
        'norename': [False, True],
        'nocog': [False, True],
        'nokegg': [False, True],
        'nopfam': [False, True],
        'euk': [False, True],
        'consensus': 50,
        'extdb': [None],
        'doublepass': [False, True],
        'map': ['bowtie', 'bwa', 'minimap2-ont', 'minimap2-pb', 'minimap2-sr'],
        'mapping_options': [None],
        'nobins': [False, True],
        'binners': ['maxbin', 'metabat', 'concoct'],
        'taxbinmode': ['s', 'c', 's+c', 'c+s'],
        't': 12,
        'block_size': [None],
        'canumem': 32,
        'lowmem': [False, True],
        'minion': [False, True],
        'empty': [False, True]
    }
    if 'binners' not in params:
        params['binners'] = defaults['binners']
    int1 = ['step']
    check_nums(self, params, defaults, int1, int, 1, 21)
    int2 = ['contiglen', 'consensus', 't', 'canumem', 'block_size']
    check_nums(self, params, defaults, int2, int)
    check_default(self, params, defaults, (int1 + int2), ['binners'])
    return defaults


def check_binspreader(self, params):
    defaults = {
        'l': 0,
        'e': 1e-5,
        'n': 5000,
        'la': 0.6,
        'metaalpha': 0.6,
        'bin_weight': 0.1,
        'm': [False, True],
        'Smax': [False, True],
        'Rcorr': [False, True],
        'cami': [False, True],
        'zero_bin': [False, True],
        'tall_multi': [False, True],
        'bin_dist': [False, True],
        'sparse_propagation': [False, True],
        'no_unbinned_bin': [False, True],
        'length_threshold': [None],
        'distance_bound': [None],
        'dataset': [None],
        'reads': [False, True]
    }
    ints = ['l', 'n']
    check_nums(self, params, defaults, ints, int)
    floats = ['e', 'la', 'metaalpha', 'bin_weight']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to BinSPreader installation folder>'
    return defaults


def check_midas2(self, params):
    defaults = {
        'force': [False, True],
        'word_size': 28,
        'marker_reads': 2,
        'marker_covered': 2,
        'prebuilt_bowtie2_indexes': [None],
        'prebuilt_bowtie2_species': [None],
        'species_list': {'all_species': ''},
        'select_by': 'median_marker_coverage,unique_fraction_covered',
        'select_threshold': '2.0,0.5',
        'aln_speed': ['very-fast', 'fast', 'sensitive', 'very-sensitive'],
        'aln_mode': ['local', 'global'],
        'aln_interleaved': [False, True],
        'fragment_length': 5000,
        'min_copy': 0.35,
        'cluster_pid': ['95', '75', '80', '85', '90', '99'],
        'max_reads': [None],
        'aln_mapid': [None],
        'aln_readq': 20,
        'aln_cov': 0.75,
        'min_cov': 1.0,
        'read_depth': 2,
        'chunk_size': [None],
        'aln_mapq': [None],
        'aln_baseq': 30,
        'aln_trim': 0,
        'paired_only': [False, True],
        'snp_maf': 0.1,
        'ignore_ambiguous': [False, True],
        'analysis_ready': [False, True],
        'genome_depth': [None],
        'genome_coverage': 0.4,
        'sample_counts': [None],
        'site_depth': 2,
        'site_ratio': 3.0,
        'site_prev': 0.9,
        'snv_type': ['common', 'rare'],
        'snp_pooled_method': ['prevalence', 'abundance'],
        'snp_type':  ['bi', 'tri', 'quad', 'any', 'mono'],
        'locus_type': ['any', 'CDS', 'IGR'],
        'advanced': [False, True],
        'robust_chunk': [False, True],
    }
    if 'species_list' not in params:
        params['species_list'] = defaults['species_list']
    else:
        if not isinstance(params['species_list'], dict):
            sys.exit('[%s] Param "species_list" not name:path dict' % self.name)

    if 'snp_type' not in params:
        params['snp_type'] = ['bi', 'tri', 'quad']

    if 'select_by' not in params:
        params['select_by'] = defaults['select_by']
    else:
        allowed = {'marker_read_counts', 'median_marker_coverage',
                   'marker_coverage', 'marker_relative_abundance',
                   'unique_fraction_covered'}
        for param in params['select_by'].split(','):
            if param not in allowed:
                sys.exit('[%s] Param "select_by" error ("%s" not in %s)' %
                         (self.name, param, str(allowed)))

    if 'select_threshold' not in params:
        params['select_threshold'] = defaults['select_threshold']
    else:
        for p in params['select_threshold'].split(','):
            try:
                float(p)
            except ValueError:
                sys.exit('[midas2] Param "select_threshold" error'
                         '("%s" not FLOAT,FLOAT,...)' % p)

    ints = ['word_size', 'marker_reads', 'marker_covered', 'aln_readq',
            'read_depth', 'aln_baseq', 'aln_trim', 'sample_counts',
            'fragment_length', 'site_depth', 'max_reads', 'chunk_size',
            'aln_mapq']
    check_nums(self, params, defaults, ints, int)
    floats = ['aln_cov', 'genome_coverage', 'site_prev']
    check_nums(self, params, defaults, floats, float, 0, 1)
    floats2 = ['aln_mapid', 'min_cov', 'min_copy', 'genome_depth', 'site_ratio']
    check_nums(self, params, defaults, floats2, float)
    floats3 = ['snp_maf']
    check_nums(self, params, defaults, floats3, float, 0, 0.5)
    check_default(self, params, defaults, (
        ints[:-3] + floats + floats2[1:] + floats3 + [
            'select_by', 'select_threshold', 'species_list']), ['snp_type'])
    defaults['databases'] = '<Path the folder containing the MIDAS2 databases>'
    defaults['genome_depth'] = '<1 for merge_genes; 5 for merge_snps>'
    defaults['sample_counts'] = '<1 for merge_genes; 2 for merge_snps>'
    return defaults


def check_woltka(self, params):
    defaults = {
        'taxa': [
            'phylum', 'order', 'class', 'family', 'genus', 'species', 'none'],
        'classifications': ['', 'go', 'eggnog', 'metacyc', 'kegg'],
        'go': ['process', 'function', 'component'],
        'demux': [False, True],
        'trim_sub': [None],
        'nodes_dmp': [None],
        'newick': [None],
        'lineage': [None],
        'columns': [None],
        'map_as_rank': [False, True],
        'names_dmp': [None],
        'rank': [None],
        'uniq': [False, True],
        'major': [None],
        'above': [False, True],
        'subok': [False, True],
        'coords': [None],
        'overlap': 80,
        'stratify': [None],
        'sizes': [None],
        'frac': [False, True],
        'scale': [None],
        'digits': [None],
        'to_biom': [None],
        'unassigned': [False, True],
        'name_as_id': [True, False],
        'add_rank': [True, False],
        'add_lineage': [False, True],
        'outmap': [None],
        'zipmap': ['gz', 'none', 'bz2', 'xz'],
        'outcov': [None],
        'chunk': [None],
        'cache': 1024,
        'no_exe': [False, True],
        'min_count': 0,
        'min_percent': 0,
        'divide': [False, True],
        'field': 2,
        'threshold': [None]
    }
    if 'classifications' not in params:
        params['classifications'] = ['go', 'eggnog', 'metacyc', 'kegg']
    if 'taxa' not in params:
        params['taxa'] = defaults['taxa']
    if 'go' not in params:
        params['go'] = defaults['go']

    ints = ['cache', 'min_count', 'field']
    check_nums(self, params, defaults, ints, int)
    floats = ['overlap', 'min_percent']
    check_nums(self, params, defaults, floats, float, 0, 100)
    check_default(self, params, defaults,
                  (ints + floats), ['classifications', 'taxa', 'go'])
    defaults['to_biom'] = [True, False]
    return defaults


def check_mobsuite(self, params):
    defaults = {
        'multi': [True, False],
        'debug': [False, True],
        'filter_db': [None],
        'genome_filter_db_prefix': [None],
        'prefix': [None],
        'mash_genome_neighbor_threshold': 0.002,
        'max_contig_size': 450000,
        'max_plasmid_size': 450000,
        'min_rep_evalue': 1e-05,
        'min_mob_evalue': 1e-05,
        'min_con_evalue': 1e-05,
        'min_rpp_evalue': 1e-05,
        'min_length': [None],
        'min_rep_ident': 80,
        'min_mob_ident': 80,
        'min_con_ident': 80,
        'min_rpp_ident': 80,
        'min_rep_cov': 80,
        'min_mob_cov': 80,
        'min_con_cov': [None],
        'min_rpp_cov': 80,
        'min_overlap': 10,
        'unicycler_contigs': [False, True],
        'run_overhang': [False, True],
        'keep_tmp': [False, True],
        'plasmid_db': [None],
        'plasmid_mash_db': [None],
        'plasmid_meta': [None],
        'plasmid_db_type': ['blastn'],
        'plasmid_replicons': [None],
        'repetitive_mask': [None],
        'plasmid_mob': [None],
        'plasmid_mpf': [None],
        'plasmid_orit': [None],
        'database_directory': [None],
        'mode': [None, 'Build', 'Update'],
        'mob_typer_file': [None],
        'taxonomy': [None],
        'ref_cluster_file': [None],
        'ref_fasta_file': [None],
        'primary_cluster_dist': 0.06,
        'secondary_cluster_dist': 0.025,
        'cluster': [False, True],
        'new_plasmids': [None]
    }
    if 'db_dir' not in params:
        params['db_dir'] = None
    else:
        db_dir = params['db_dir']
        if not (self.config.dev or self.config.workshop):
            if db_dir and isdir(db_dir):
                sys.exit('[mobsuite] Params "db_dir": path do not exist')

    if 'min_length' not in params:
        params['min_length'] = None
    else:
        min_length = params['min_length']
        if not (self.config.dev or self.config.workshop):
            if min_length and not str(min_length).isdigit():
                sys.exit("[mobsuite] 'min_length' not of <class 'int'>")

    if 'min_con_cov' not in params:
        params['min_con_cov'] = None
    else:
        m_con_cov = params['min_con_cov']
        if not (self.config.dev or self.config.workshop) and m_con_cov:
            if not str(m_con_cov).isdigit() or not 0 <= int(m_con_cov) <= 100:
                sys.exit("[mobsuite] 'min_con_cov' not of <class 'int'>")

    ints = ['max_contig_size', 'max_plasmid_size']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['min_rep_ident', 'min_mob_ident', 'min_con_ident', 'min_rpp_ident',
             'min_rep_cov', 'min_mob_cov', 'min_rpp_cov', 'min_overlap']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    floats = ['mash_genome_neighbor_threshold', 'min_rep_evalue',
              'min_mob_evalue', 'min_con_evalue', 'min_rpp_evalue',
              'primary_cluster_dist', 'secondary_cluster_dist']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (
            ints + ints2 + floats + ['db_dir', 'min_length', 'min_con_cov']))
    defaults['db_dir'] = '<Path to the MOB-suite database installed locally>'
    return defaults


def check_yamb(self, params):
    defaults = {'l': 0, 'm': 0}
    check_nums(self, params, defaults, ['l', 'm'], int)
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the YAMB folder>'
    return defaults


def check_salmon(self, params):
    defaults = {
        'useAlignments': [True, False],
        'kmerLen': 31,
        'gencode': [False, True],
        'features': [False, True],
        'keepDuplicates': [False, True],
        'keepFixedFasta': [False, True],
        'filterSize': -1,
        'sparse': [False, True],
        'decoys': [None],
        'ont': [False, True],
        'no_clip': [False, True],
        'type': ['puff'],
        'seqBias': [False, True],
        'gcBias': [False, True],
        'posBias': [False, True],
        'incompatPrior': 0,
        'geneMap': [None],
        'auxTargetFile': [None],
        'meta': [False, True],
        'discardOrphans': [False, True],
        'discardOrphansQuasi': [False, True],
        'validateMappings': [False, True],
        'consensusSlack': 0.35,
        'preMergeChainSubThresh': 0.75,
        'postMergeChainSubThresh': 0.9,
        'orphanChainSubThresh': 0.95,
        'scoreExp': 1,
        'numErrorBins': 6,
        'noErrorModel': [False, True],
        'sampleOut': [False, True],
        'sampleUnaligned': [False, True],
        'mappingCacheMemoryLimit': 2000000,
        'minScoreFraction': 0.65,
        'mismatchSeedSkip': 3,
        'disableChainingHeuristic': [False, True],
        'decoyThreshold': 1,
        'ma': 2,
        'mp': -4,
        'go': 6,
        'ge': 2,
        'bandwidth': 15,
        'allowDovetail': [False, True],
        'recoverOrphans': [False, True],
        'mimicBT2': [False, True],
        'mimicStrictBT2': [False, True],
        'softclip': [False, True],
        'softclipOverhangs': [False, True],
        'fullLengthAlignment': [False, True],
        'hardFilter': [False, True],
        'minAlnProb': 1e-05,
        'writeMappings': [None],
        'writeQualities': [False, True],
        'hitFilterPolicy': ['BEFORE', 'AFTER', 'BOTH', 'NONE'],
        'alternativeInitMode': [False, True],
        'auxDir': [None],
        'skipQuant': [False, True],
        'dumpEq': [False, True],
        'dumpEqWeights': [False, True],
        'minAssignedFrags': 10,
        'reduceGCMemory': [False, True],
        'biasSpeedSamp': 5,
        'fldMax': 1000,
        'fldMean': 250,
        'fldSD': 25,
        'forgettingFactor': 0.65,
        'initUniform': [False, True],
        'maxOccsPerHit': 1000,
        'maxReadOcc': 200,
        'maxRecoverReadOcc': 2500,
        'noLengthCorrection': [False, True],
        'noEffectiveLengthCorrection': [False, True],
        'noSingleFragProb': [False, True],
        'noFragLengthDist': [False, True],
        'noBiasLengthThreshold': [False, True],
        'numBiasSamples': 2000000,
        'numAuxModelSamples': 5000000,
        'numPreAuxModelSamples': 5000,
        'useEM': [False, True],
        'useVBOpt': [True, False],
        'rangeFactorizationBins': 4,
        'numGibbsSamples': 0,
        'noGammaDraw': [False, True],
        'numBootstraps': 0,
        'bootstrapReproject': [False, True],
        'thinningFactor': 16,
        'perTranscriptPrior': [False, True],
        'perNucleotidePrior': [False, True],
        'sigDigits': 3,
        'vbPrior': 0.01,
        'writeOrphanLinks': [False, True],
        'writeUnmappedNames': [False, True],
        'column': ['tpm', 'len', 'elen', 'numreads'],
        'genes': [False, True],
        'missing': [None]
    }
    ints = ['sigDigits', 'thinningFactor', 'numBootstraps', 'numGibbsSamples',
            'rangeFactorizationBins', 'numBiasSamples', 'numAuxModelSamples',
            'numPreAuxModelSamples', 'maxOccsPerHit', 'maxReadOcc',
            'maxRecoverReadOcc', 'minAssignedFrags', 'reduceGCMemory',
            'biasSpeedSamp', 'fldMax', 'fldMean', 'fldSD', 'minAssignedFrags',
            'kmerLen', 'filterSize', 'scoreExp', 'numErrorBins',
            'mappingCacheMemoryLimit', 'mismatchSeedSkip', 'ma', 'mp', 'go',
            'ge', 'bandwidth']
    check_nums(self, params, defaults, ints, int)
    floats = ['vbPrior', 'minAlnProb', 'decoyThreshold', 'consensusSlack',
              'preMergeChainSubThresh', 'postMergeChainSubThresh',
              'orphanChainSubThresh', 'minScoreFraction', 'incompatPrior']
    check_nums(self, params, defaults, floats, float, 0, 1)
    floats2 = ['forgettingFactor']
    check_nums(self, params, defaults, floats2, float, 0.5, 1)
    check_default(self, params, defaults, (ints + floats + floats2))
    return defaults


def check_kallisto(self, params):
    defaults = {
        'kmer_size': 31,
        'make_unique': [False, True],
        'bias': [False, True],
        'bootstrap_samples': 0,
        'seed': 42,
        'plaintext': [False, True],
        'fusion': [False, True],
        'single': [False, True],
        'single_overhang': [False, True],
        'fr_stranded': [False, True],
        'rf_stranded': [False, True],
        'ec_file': [None],
        'fragment_file': [None],
        'fragment_length': [None],
        'genemap': [None],
        'sd': [None],
        'pseudobam': [False, True],
        'genomebam': [False, True],
        'paired': [False, True],
        'unstranded': [False, True],
        'num': [False, True],
        'bam': [False, True],
        'gtf': [None],
        'chromosomes': [None],
        'batch': [None],
        'technology': [None],
    }
    ints = ['kmer_size', 'bootstrap_samples', 'seed']
    check_nums(self, params, defaults, ints, int)
    floats = []
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats))
    return defaults


def check_metaclade2(self, params):
    defaults = {
        'arch': [False, True],
        'user_cfg': [None],
        'remove_temp': [False, True],
        'evalue_cutoff': 1e-3,
        'evalue_cutconf': 1e-10,
        'overlappingAA': 30,
        'overlappingMaxDomain': 50,
        'sge': [False, True],
    }
    ints = ['overlappingAA', 'overlappingMaxDomain']
    check_nums(self, params, defaults, ints, int)
    floats = ['evalue_cutoff', 'evalue_cutconf']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the metaclade2 installation folder>'
    return defaults


def check_mycc(self, params):
    defaults = {
        't': 1000,
        'lt': 0.7,
        'ct': 0,
        'meta': ['single', 'meta'],
        'p': 20,
        'pm': 500,
        'mask': [False, True],
        'keep': [False, True],
    }
    ints = ['pm', 'ct', 't']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['p']
    check_nums(self, params, defaults, ints2, int, 5, 50)
    floats = ['lt']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + ints2 + floats))
    return defaults


def check_argsoap(self, params):
    defaults = {
        'f': ['fq'],
        'q': [False, True],
        'z': [False, True],
        'x': 1e-10,
        'y': 3,
        'v': 0.45,
        'l': 25,
        'e': 1e-7,
        'd': 80,
        'p': 'stage2output',
        'r': [False, True]
    }
    ints = ['y', 'l']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['d']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    floats = ['x', 'v', 'e']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + ints2 + floats))
    return defaults


def check_deeptmhmm(self, params):
    defaults = {}
    return defaults


def check_gtdbtk(self, params):
    defaults = {
        'prot_model': ['WAG', 'JTT', 'LG'],
        'rnd_seed': 42,
        'min_af': 0.65,
        'min_perc_aa': 10,
        'cols_per_gene': 42,
        'min_consensus': 25,
        'max_consensus': 95,
        'min_perc_taxa': 50,
        'genes': [None],
        'batchfile': [None],
        'taxa_filter': [None],
        'outgroup_taxon': [None],
        'custom_taxonomy_file': [None],
        'gtdbtk_classification_file': [None],
        'force': [False, True],
        'gamma': [False, True],
        'archaea': [False, True],
        'bacteria': [True, False],
        'full_tree': [False, True],
        'no_support': [False, True],
        'skip_gtdb_refs': [False, True],
        'keep_intermediates': [False, True],
        'custom_msa_filters': [False, True],
        'write_single_copy_genes': [False, True],
    }
    if params.get('outgroup_taxon'):
        if not isinstance(params.get['outgroup_taxon'], dict):
            sys.exit('[gtdbtk] Param "outgroup_taxon" must be a dict')
        if not params['outgroup_taxon'].get('bacteria'):
            params['outgroup_taxon']['bacteria'] = 'p__Cyanobacteria'
        if not params['outgroup_taxon'].get('archaea'):
            params['outgroup_taxon']['archaea'] = 'p__Undinarchaeota'
    else:
        params['outgroup_taxon'] = {}
        if params.get('bacteria', True):
            params['outgroup_taxon']['bacteria'] = 'p__Cyanobacteria'
        if params.get('archaea', False):
            params['outgroup_taxon']['archaea'] = 'p__Undinarchaeota'

    ints = ['rnd_seed', 'cols_per_gene']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['min_perc_aa', 'min_consensus', 'max_consensus', 'min_perc_taxa']
    check_nums(self, params, defaults, ints2, int, 0, 100)
    floats = ['min_af']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults,
                  (ints + ints2 + floats + ['outgroup_taxon']))
    defaults['outgroup_taxon'] = {'bacteria': 'p__Cyanobacteria',
                                  'archaea': 'p__Undinarchaeota'}
    return defaults


# def check_itassermtd(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
# def check_itasser(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
# def check_graspx(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
# def check_hmmgraspx(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
def check_rundbcan(self, params):
    defaults = {
        'inputType': ['protein', 'prok', 'meta'],
        'AuxillaryFile': [None],
        'Tools': ['all', 'diamond', 'hmmer', 'eCAMI'],
        'dbCANFile': [None],
        'dia_eval': 1e-102,
        'hmm_eval': 1e-15,
        'stp_eval': 1e-4,
        'tf_eval': 1e-4,
        'eCAMI_kmer_db': [None],
        'eCAMI_k_mer': 8,
        'eCAMI_jobs': 8,
        'eCAMI_important_k_mer_number': 5,
        'eCAMI_beta': 2,
        'hmm_cov': 0.35,
        'tf_cov': 0.35,
        'stp_cov': 0.3,
        'out_pre': [None],
        'db_dir': [None],
        'cgc_dis': 2,
        'use_signalP': [False, True],
        'signalP_path': [None],
        'gram': ['all', 'n', 'p']
    }
    ints = ['eCAMI_k_mer', 'eCAMI_jobs', 'eCAMI_important_k_mer_number',
            'eCAMI_beta']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['cgc_dis']
    check_nums(self, params, defaults, ints2, int, 2, 10)
    floats = ['hmm_cov', 'tf_cov', 'stp_cov', 'dia_eval',
              'hmm_eval', 'stp_eval', 'tf_eval']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + ints2 + floats))
    return defaults


def check_ioncom(self, params):
    defaults = {}
    return defaults


# def check_anvio(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
# def check_mocat2(self, params, soft):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int, soft.name)
#     floats = []
#     check_nums(self, params, defaults, floats, float, soft.name)
#     check_default(self, params, defaults, soft.name, (ints + floats))
#     defaults[''] = '<>'
#     return defaults
#
#
def check_phyloflash(self, params):
    defaults = {
        'id': 70,
        'taxlevel': 4,
        'clusterid': 97,
        'maxinsert': 1200,
        'readlength': 100,
        'amplimit': 500000,
        'ref_minlength': 800,
        'emirge': [None],
        'trusted': [None],
        'sortmerna': [None],
        'readlimit': [None],
        'sc': [False, True],
        'log': [False, True],
        'zip': [False, True],
        'html': [True, False],
        'crlf': [False, True],
        'poscov': [False, True],
        'treemap': [False, True],
        'keeptmp': [False, True],
        'everything': [False, True],
        'skip_spades': [False, True],
        'decimalcomma': [False, True],
        'almosteverything': [False, True]
    }
    ints = ['readlength', 'amplimit', 'taxlevel', 'ref_minlength']
    check_nums(self, params, defaults, ints, int)
    ints2 = ['id']
    check_nums(self, params, defaults, ints2, int, 50, 98)
    ints3 = ['clusterid']
    check_nums(self, params, defaults, ints3, int, 50, 100)
    ints4 = ['maxinsert']
    check_nums(self, params, defaults, ints4, int, 0, 1200)
    check_default(self, params, defaults, (ints + ints2 + ints3 + ints4))
    return defaults


def check_instrain(self, params):
    defaults = {
        'max_insert_relative': 3,
        'min_scaffold_reads': 1,
        'min_genome_coverag': 0,
        'rarefied_coverage': 50,
        'window_length': 10000,
        'min_read_ani': 0.95,
        'min_insert': 50,
        'min_freq': 0.05,
        'min_mapq': -1,
        'min_cov': 5,
        'fdr': 1e-06,
        'min_snp': 20,
        'use_full_fasta_header': [False, True],
        'force_compress': [False, True],
        'pairing_filter': ['paired_only', 'non_discordant', 'all_reads'],
        'priority_reads': [None],
        'detailed_mapping_info': [False, True],
        'gene_file': [None],
        'stb': [None],
        'mm_level': [False, True],
        'skip_mm_profiling': [False, True],
        'database_mode': [False, True],
        'store_everything': [False, True],
        'scaffolds_to_profile': [None],
        'skip_genome_wide': [False, True],
        'skip_plot_generation': [False, True],
        'breadth': 0.5,
        'scaffolds': [None],
        'genome': [None],
        'store_coverage_overlap': [False, True],
        'store_mismatch_locations': [False, True],
        'include_self_comparisons': [False, True],
        'group_length': 10000000,
        'ani_threshold': 0.99999,
        'coverage_treshold': 0.1,
        'clusterAlg': ['average', 'median', 'weighted', 'complete', 'ward',
                       'centroid', 'single'],
        'bams': [None],
        'skip_popANI': [False, True],
        'breadth_cutoff': 0.5,
        'stringent_breadth_cutoff': 0.0
    }
    ints = ['max_insert_relative', 'min_scaffold_reads', 'group_length',
            'min_genome_coverag', 'rarefied_coverage', 'window_length',
            'min_insert', 'min_mapq', 'min_cov', 'min_snp', 'ani_threshold']
    check_nums(self, params, defaults, ints, int)
    floats = ['min_read_ani', 'min_freq', 'fdr', 'breadth', 'ani_threshold',
              'coverage_treshold', 'breadth_cutoff', 'stringent_breadth_cutoff']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    return defaults


def check_closedref(self, params):
    defaults = {}
    return defaults


def check_threecac(self, params):
    defaults = {}
    return defaults


def check_plasclass(self, params):
    defaults = {}
    return defaults


def check_deepvirfinder(self, params):
    defaults = {'len': 1}
    if 'model_dir' not in params:
        sys.exit('[deepvirfinder] Params "model_dir" must exist')
    ints = ['len']
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, ints)
    defaults['model_dir'] = '<Path to the DeepVirFinder model directory>'
    return defaults


def check_wish(self, params):
    defaults = {
        'k': 8,
        'a': 16,
        'b': 1,
        'c': ['build', 'predict'],
        'p': [False, True],
        'z': [False, True],
        'n': [None]
    }
    ints = ['k', 'a', 'b']
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, ints)
    return defaults


def check_eggnogmapper(self, params):
    defaults = {
        'x': [False, True],
        'D': [False, True],
        'F': [False, True],
        'P': [False, True],
        'M': [False, True],
        'H': [False, True],
        'y': [False, True],
        'f': [False, True],
        's': [False, True],
        'q': [False, True],
        'resume': [False, True],
        'override': [False, True],
        'translate': [False, True],
        'outfmt_short': [False, True],
        'dmnd_ignore_warnings': [False, True],
        'usemem': [False, True],
        'report_no_hits': [False, True],
        'cut_ga': [False, True],
        'no_annot': [False, True],
        'dbmem': [False, True],
        'report_orthologs': [False, True],
        'md5': [False, True],
        'no_file_comments': [False, True],
        'excel': [False, True],
        'overlap_tol': 0.0,
        'evalue': 0.001,
        'start_sens': 3,
        'sens_steps': 3,
        'final_sens': 7,
        'port': 51700,
        'end_port': 53200,
        'num_servers': 1,
        'num_workers': 1,
        'hmm_maxhits': 1,
        'hmm_maxseqlen': 5000,
        'Z': 40000000,
        'seed_ortholog_evalue': 0.001,
        'mp_start_method': ['spawn', 'fork', 'forkserver'],
        'itype': ['proteins', 'CDS', 'genome', 'metagenome'],
        'genepred': ['search', 'prodigal'],
        'allow_overlaps': ['none', 'strand', 'diff_frame', 'all'],
        'm': ['diamond', 'mmseqs', 'hmmer', 'no_search', 'cache', 'novel_fams'],
        'sensmode': [
            'sensitive', 'default', 'fast', 'mid-sensitive', 'more-sensitive'],
        'dmnd_algo': ['auto', '0', '1', 'ctg'],
        'dmnd_iterate': ['yes', 'no'],
        'qtype': ['seq', 'hmm'],
        'dbtype': ['hmmdb', 'seqdb'],
        'clean_overlaps': [
            'none', 'all', 'clans', 'hmmsearch_all', 'hmmsearch_clans'],
        'tax_scope': [
            'auto', 'auto_broad', 'all_narrow', 'archaea', 'bacteria',
            'bacteria_broad', 'eukaryota', 'eukaryota_broad',
            'prokaryota_broad', 'none'],
        'tax_scope_mode': [
            'inner_narrowest', 'broadest', 'inner_broadest', 'narrowest'],
        'target_orthologs': [
            'all', 'one2one', 'many2one', 'one2many', 'many2many'],
        'go_evidence': ['non-electronic', 'experimental', 'all'],
        'pfam_realign': ['none', 'realign', 'denovo'],
        'decorate_gff': ['no', 'yes'],
        'decorate_gff_ID_field': ['ID'],
        'taxa': [None],
        'dbname': [None],
        'taxids': [None],
        'annotate_hits_table': [None],
        'cache': [None],
        'trans_table': [None],
        'training_genome': [None],
        'training_file': [None],
        'pident': [None],
        'query_cover': [None],
        'subject_cover': [None],
        'score': [None],
        'dmnd_db': [None],
        'matrix': [None, 'BLOSUM62', 'BLOSUM90', 'BLOSUM80', 'BLOSUM50',
                   'BLOSUM45', 'PAM250', 'PAM70', 'PAM30'],
        'dmnd_frameshift': [None],
        'gapopen': [None],
        'gapextend': [None],
        'block_size': [None],
        'index_chunks': [None],
        'mmseqs_db': [None],
        'mmseqs_sub_mat': [None],
        'database': [None],
        'servers_list': [None],
        'seed_ortholog_score': [None],
        'target_taxa': [None],
        'excluded_taxa': [None],
    }
    if 'data_dir' not in params:
        sys.exit('[eggnogmapper] Params "data_dir" must exist')
    elif not isdir(params['data_dir']):
        if not (self.config.dev or self.config.workshop):
            sys.exit('[eggnogmapper] Params "data_dir" must be database folder')

    ints = ['start_sens', 'sens_steps', 'final_sens', 'port', 'end_port',
            'num_servers', 'num_workers', 'hmm_maxhits', 'hmm_maxseqlen', 'Z']
    check_nums(self, params, defaults, ints, int)
    floats = ['overlap_tol', 'evalue', 'seed_ortholog_evalue']
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats + ['data_dir']))
    defaults['data_dir'] = '<Path to the "eggnog-mapper-data" database folder>'
    return defaults


def check_ngless(self, params):
    defaults = {
        'script': [None],
        'validate_only': [False, True],
        'print_last': [False, True],
        'strict_threads': [False, True],
        'create_report': [False, True],
        'keep_temporary_files': [False, True],
        'no_header': [False, True],
        'subsample': [False, True],
        'experimental_features': [False, True],
        'export_json': [None],
        'export_cwl': [None],
        'check_deprecation': [False, True],
        'search_dir': [None],
        'search_path': [None],
        'index_path': [None],
        'color': ['auto', 'force', 'no'],
        'trace': [False, True],
    }
    check_default(self, params, defaults)
    return defaults


def check_motus(self, params):
    defaults = {
        'B': [False, True],
        'v': ['3', '1', '2', '4'],
        'I': [None],
        'M': [None],
        'e': [False, True],
        'c': [False, True],
        'p': [False, True],
        'u': [False, True],
        'q': [False, True],
        'C': [None, 'precision', 'recall', 'parenthesis'],
        'A': [False, True],
        'k': ['mOTU', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'],
        'g': 3,
        'l': 75,
        'y': ['insert.scaled_counts', 'base.coverage', 'insert.raw_counts'],
        'b': [False, True],
        'fb': 80.0,
        'fd': 5.0,
        'fm': 2,
        'fp': 0.50,
        'fc': 5.0,
        'K': [False, True]
    }
    ints = ['g']
    check_nums(self, params, defaults, ints, int, 1, 10)
    ints2 = ['l']
    check_nums(self, params, defaults, ints, int, 0, 100)
    floats = ['fb', 'fd', 'fm']
    check_nums(self, params, defaults, ints, int, 1, 10)
    floats2 = ['fp']
    check_nums(self, params, defaults, floats2, float, 0, 1)
    floats3 = ['fc']
    check_nums(self, params, defaults, floats3, float)
    check_default(
        self, params, defaults, (ints + ints2 + floats + floats2 + floats3))
    return defaults


def check_metacompare(self, params):
    defaults = {}
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the metacompare installation folder>'
    return defaults


def check_oritfinder(self, params):
    defaults = {}
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the oritfinder installation folder>'
    return defaults


def check_diting(self, params):
    defaults = {'noclean': [False, True]}
    if 'samples' not in params:
        params['samples'] = None
    else:
        if not isinstance(params['samples'], list):
            if params['samples'] not in ['all', None]:
                sys.exit('[diting] Param "samples" not list, "all" or "null"')
    check_default(self, params, defaults)
    defaults['samples'] = '<list of samples/co-assembly groups, "all" or null>'
    return defaults


def check_abritamr(self, params):
    defaults = {
        'identity': 0.9,
        'sop': ['general', 'plus'],
        'species': [None, 'Burkholderia_cepacia', 'Acinetobacter_baumannii',
                    'Streptococcus_pyogenes','Streptococcus_agalactiae',
                    'Streptococcus_pneumoniae', 'Enterococcus_faecium',
                    'Pseudomonas_aeruginosa', 'Staphylococcus_pseudintermedius',
                    'Clostridioides_difficile', 'Klebsiella', 'Neisseria',
                    'Campylobacter', 'Salmonella', 'Escherichia',
                    'Staphylococcus_aureus', 'Burkholderia_pseudomallei',
                    'Enterococcus_faecalis']
    }
    if 'samples' not in params:
        params['samples'] = None
    else:
        if not isinstance(params['samples'], list):
            if params['samples'] not in ['all', None]:
                sys.exit('[abritamr] Param "samples" not list, "all" or "null"')
    check_nums(self, params, defaults, ['identity'], float, 0, 1)
    check_default(self, params, defaults, ['identity'])
    defaults['samples'] = '<list of samples/co-assembly groups, "all" or null>'
    return defaults


def check_keggcharter(self, params):
    defaults = {
        'input': [None],
        'metabolic_maps': [
            '00680', '00020', '00030', '00040', '00051', '00052', '00053',
            '00061', '00062', '00071', '00073', '00100', '00120', '00121',
            '00130', '00140', '00190', '00195', '00196', '00220', '00230',
            '00232', '00240', '00250', '00253', '00254', '00260', '00261',
            '00270', '00280', '00281', '00290', '00300', '00310', '00311',
            '00330', '00331', '00332', '00333', '00340', '00350', '00360',
            '00361', '00362', '00363', '00364', '00365', '00380', '00400',
            '00401', '00402', '00403', '00404', '00405', '00410', '00430',
            '00440', '00450', '00460', '00470', '00480', '00500', '00510',
            '00511', '00512', '00513', '00514', '00515', '00520', '00521',
            '00522', '00523', '00524', '00525', '00531', '00532', '00533',
            '00534', '00540', '00541', '00542', '00543', '00550', '00552',
            '00561', '00562', '00563', '00564', '00565', '00571', '00572',
            '00590', '00591', '00592', '00600', '00601', '00603', '00604',
            '00620', '00621', '00622', '00623', '00624', '00625', '00626',
            '00627', '00630', '00633', '00640', '00642', '00643', '00650',
            '00660', '00670', '00010', '00710', '00720', '00730', '00740',
            '00750', '00760', '00770', '00780', '00785', '00790', '00791',
            '00830', '00860', '00900', '00901', '00902', '00903', '00904',
            '00905', '00906', '00908', '00909', '00910', '00920', '00930',
            '00940', '00941', '00942', '00943', '00944', '00945', '00946',
            '00950', '00960', '00965', '00966', '00970', '00980', '00981',
            '00982', '00983', '00984', '00996', '00997', '00998', '00999',
            '01010', '01040', '01051', '01052', '01053', '01054', '01055',
            '01056', '01057', '01059', '01060', '01061', '01062', '01063',
            '01064', '01065', '01066', '01070', '01100', '01110', '01120',
            '01200', '01210', '01212', '01220', '01230', '01232', '01240',
            '01250', '01501', '01502', '01503', '01521', '01522', '01523',
            '01524', '02010', '02020', '02024', '02025', '02026', '02030',
            '02040', '02060', '03008', '03010', '03013', '03015', '03018',
            '03020', '03022', '03030', '03040', '03050', '03060', '03070',
            '03230', '03240', '03250', '03260', '03261', '03262', '03263',
            '03264', '03265', '03266', '03267', '03268', '03320', '03410',
            '03420', '03430', '03440', '03450', '03460', '04010', '04011',
            '04012', '04013', '04014', '04015', '04016', '04020', '04022',
            '04024', '04060', '04061', '04062', '04064', '04066', '04068',
            '04070', '04071', '04072', '04075', '04080', '04110', '04111',
            '04112', '04113', '04114', '04115', '04120', '04122', '04130',
            '04136', '04137', '04138', '04139', '04140', '04141', '04142',
            '04144', '04145', '04146', '04150', '04151', '04152', '04210',
            '04211', '04212', '04213', '04214', '04215', '04216', '04217',
            '04218', '04260', '04261', '04270', '04310', '04320', '04330',
            '04340', '04341', '04350', '04360', '04361', '04370', '04371',
            '04380', '04390', '04391', '04392', '04510', '04512', '04514',
            '04520', '04530', '04540', '04550', '04610', '04611', '04612',
            '04613', '04614', '04620', '04621', '04622', '04623', '04624',
            '04625', '04626', '04630', '04640', '04650', '04657', '04658',
            '04659', '04660', '04662', '04664', '04666', '04668', '04670',
            '04672', '04710', '04711', '04712', '04713', '04714', '04720',
            '04721', '04722', '04723', '04724', '04725', '04726', '04727',
            '04728', '04730', '04740', '04742', '04744', '04745', '04750',
            '04810', '04910', '04911', '04912', '04913', '04914', '04915',
            '04916', '04917', '04918', '04919', '04920', '04921', '04922',
            '04923', '04924', '04925', '04926', '04927', '04928', '04929',
            '04930', '04931', '04932', '04933', '04934', '04935', '04936',
            '04940', '04950', '04960', '04961', '04962', '04964', '04966',
            '04970', '04971', '04972', '04973', '04974', '04975', '04976',
            '04977', '04978', '04979', '05010', '05012', '05014', '05016',
            '05017', '05020', '05022', '05030', '05031', '05032', '05033',
            '05034', '05100', '05110', '05111', '05120', '05130', '05131',
            '05132', '05133', '05134', '05135', '05140', '05142', '05143',
            '05144', '05145', '05146', '05150', '05152', '05160', '05161',
            '05162', '05163', '05164', '05165', '05166', '05167', '05168',
            '05169', '05170', '05171', '05200', '05202', '05203', '05204',
            '05205', '05206', '05207', '05208', '05210', '05211', '05212',
            '05213', '05214', '05215', '05216', '05217', '05218', '05219',
            '05220', '05221', '05222', '05223', '05224', '05225', '05226',
            '05230', '05231', '05235', '05310', '05320', '05321', '05322',
            '05323', '05330', '05332', '05340', '05410', '05412', '05414',
            '05415', '05416', '05417', '05418', '07011', '07012', '07013',
            '07014', '07015', '07016', '07017', '07018', '07019', '07020',
            '07021', '07023', '07024', '07025', '07026', '07027', '07028',
            '07029', '07030', '07031', '07032', '07033', '07034', '07035',
            '07036', '07037', '07038', '07039', '07040', '07041', '07042',
            '07043', '07044', '07045', '07046', '07047', '07048', '07049',
            '07050', '07051', '07052', '07053', '07054', '07055', '07056',
            '07057', '07110', '07112', '07114', '07117', '07211', '07212',
            '07213', '07214', '07215', '07216', '07217', '07218', '07219',
            '07220', '07221', '07222', '07223', '07224', '07225', '07226',
            '07227', '07228', '07229', '07230', '07231', '07232', '07233',
            '07234', '07235'],
        'genomic_columns': [None],
        'transcriptomic_columns': [None],
        'taxa_list': [None],
        'number_of_taxa': 10,
        'kegg_column': ['kegg'],
        'ko_column': ['ko'],
        'ec_column': ['ec'],
        'input_quantification': [False, True],
        'input_taxonomy': [None],
        'taxa_column': ['genus', None],
        'resume': [False, True],
        'step': 150,
        'taxonomic_levels': [
            'phylum', 'order', 'class', 'family', 'genus', 'species', 'none'],
    }
    if 'taxonomic_levels' not in params:
        params['taxonomic_levels'] = [
            'phylum', 'order', 'class', 'family', 'genus']

    if 'metabolic_maps' not in params:
        params['metabolic_maps'] = defaults['metabolic_maps'][:1]
    else:
        if not isinstance(params['metabolic_maps'], list):
            sys.exit('[keggcharter] Params "metabolic_maps" not a list')
        maps = list(map(str, params['metabolic_maps']))
        defs = defaults['metabolic_maps']
        diffs = list(set(maps).difference(set(defs)))
        if diffs:
            sys.exit('[keggcharter] maps not allowed: %s' % '", "'.join(diffs))

    ints = ['number_of_taxa', 'step']
    check_nums(self, params, defaults, ints, int)
    check_default(
        self, params, defaults, ints, ['metabolic_maps', 'taxonomic_levels'])
    return defaults


def check_metamic(self, params):
    defaults = {
        'mlen': 5000,
        'st': 0.8,
        'slen': 1000,
        'nb': [None],
        'rb': [None],
        'at': [None]
    }
    ints = ['mlen', 'slen']
    check_nums(self, params, defaults, ints, int)
    floats = ['st']
    check_nums(self, params, defaults, floats, float)
    check_default(self, params, defaults, (ints + floats))
    return defaults


def check_trf(self, params):
    defaults = {
        'match': 2,
        'mismatch': 7,
        'delta': 7,
        'match_probability': 80,
        'indel_probability': 10,
        'min_score': 50,
        'max_period': 500,
        'm': [True, False],
        'f': [True, False],
        'd': [True, False],
        'h': [True, False],
        'r': [False, True],
        'l': 10
    }
    ints = ['match', 'mismatch', 'delta', 'min_score', 'max_period', 'l']
    check_nums(self, params, defaults, ints, int)
    probas = ['match_probability', 'indel_probability']
    check_nums(self, params, defaults, probas, int, 0, 100)
    check_default(self, params, defaults, (ints + probas))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the trf binary folder>'
    return defaults


def check_kmerssr(self, params):
    defaults = {
        'a': 'A,C,G,T',
        'p': '4,5,6,7,8',
        'l': 100,
        'L': 500000000,
        'n': 16,
        'N': 10000,
        'r': 2,
        'R': 10000,
        'Q': 1000,
        'A': [False, True],
        'e': [False, True],
        'd': [False, True],
        's': [None]
    }
    if 'a' not in params:
        params['a'] = defaults['a']
    else:
        if not isinstance(params['a'], str) or ',' not in params['a']:
            sys.exit('[kmerssr] Param "a" not a string or has no comma')
        not_acgt = [x for x in params['a'].split(',') if x not in 'ACGT']
        if not_acgt:
            sys.exit('[kmerssr] Param "a" must be of "A", "C", "G", or "T"')

    if 'p' not in params:
        params['p'] = defaults['p']
    else:
        if not isinstance(params['p'], str) or ',' not in params['p']:
            sys.exit('[kmerssr] Param "p" not a string or has no comma')
        not_digits = [x for x in params['p'].split(',') if not x.isdigit()]
        if not_digits:
            sys.exit('[kmerssr] Param "p" must be digits')

    ints = ['l', 'L', 'n', 'N', 'r', 'R', 'Q']
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, (ints + ['a', 'p']))
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to the Kmer-SSR binary folder>'
    return defaults


def check_divissr(self, params):
    defaults = {
        'min_motif_size': 1,
        'max_motif_size': 6,
        'min_length': 12,
        'filter_reads': [False, True],
        'analyse': [True, False],
        'compound': [True, False],
        'annotate': [None],
        'comp_dist': 0,
        'anno_format': ['GFF', 'GTF'],
        'gene_key': 'gene',
        'up_promoter': 1000,
        'down_promoter': 1000
    }
    ints = ['min_motif_size', 'max_motif_size', 'min_length',
            'comp_dist', 'up_promoter', 'down_promoter']
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, ints)
    return defaults


def check_staramr(self, params):
    defaults = {
        'pointfinder_organism': [
            None, 'helicobacter_pylori', 'campylobacter', 'salmonella',
            'enterococcus_faecium', 'escherichia_coli', 'enterococcus_faecalis'
        ],
        'plasmidfinder_database_type': [
            None, 'gram_positive', 'enterobacteriaceae', 'enterobacteriales'],
        'mlst_scheme': [None],
        'exclude_genes_file': [None],
        'genome_size_lower_bound': 4000000,
        'genome_size_upper_bound': 6000000,
        'minimum_N50_value': 10000,
        'minimum_contig_length': 300,
        'unacceptable_number_contigs': 1000,
        'pid_threshold': 98.0,
        'percent_length_overlap_resfinder': 60.0,
        'percent_length_overlap_pointfinder': 95.0,
        'percent_length_overlap_plasmidfinder': 60.0,
        'ignore_invalid_files': [False, True],
        'no_exclude_gene': [False, True],
        'exclude_negatives': [False, True],
        'exclude_resistance_phenotypes': [False, True],
        'report_all_blast': [False, True]
    }
    if 'exclude_genes_file' in params and params['exclude_genes_file']:
        if not isfile(params['exclude_genes_file']):
            sys.exit('[staramr] Param "exclude_genes_file" is not a file')
    else:
        params['exclude_genes_file'] = None
    ints = ['genome_size_lower_bound', 'genome_size_upper_bound',
            'minimum_N50_value', 'minimum_contig_length',
            'unacceptable_number_contigs']
    check_nums(self, params, defaults, ints, int)
    floats = ['pid_threshold', 'percent_length_overlap_resfinder',
              'percent_length_overlap_pointfinder',
              'percent_length_overlap_plasmidfinder']
    check_nums(self, params, defaults, floats, float)
    check_default(
        self, params, defaults, (ints + floats + ['exclude_genes_file']))
    return defaults


def check_platon(self, params):
    defaults = {
        'mode': ['accuracy', 'sensitivity', 'specificity'],
        'characterize': [False, True],
        'meta': [True, False],
    }
    check_default(self, params, defaults)
    return defaults


def check_genomad(self, params):
    defaults = {
        'presets': [None, 'conservative', 'relaxed'],
        'cleanup': [False, True],
        'restart': [False, True],
        'disable_find_proviruses': [False, True],
        'disable_nn_classification': [False, True],
        'enable_score_calibration': [False, True],
        'conservative_taxonomy': [False, True],
        'sensitivity': 4.2,
        'splits': 0,
        'skip_integrase_identification': [False, True],
        'skip_trna_identification': [False, True],
        'composition': ['auto', 'metagenome', 'virome'],
        'force_auto': [False, True],
        'single_window': [False, True],
        'batch_size': 128,
        'min_score': 0.7,
        'max_fdr': 0.1,
        'min_plasmid_marker_enrichment': 0.0,
        'min_virus_marker_enrichment': 0.0,
        'min_plasmid_hallmarks': 0,
        'min_plasmid_hallmarks_short_seqs': 1,
        'min_virus_hallmarks': 0,
        'min_virus_hallmarks_short_seqs': 1,
        'max_uscg': 4,
        'evalue': 0.001,
        'use_minimal_db': [False, True],
        'crf_threshold': 0.4,
        'marker_threshold': 12.0,
        'marker_threshold_integrase': 8.0,
        'marker_threshold_edge': 8.0,
        'max_integrase_distance': 10000,
        'max_trna_distance': 5000,
    }
    ints = [
        'splits', 'min_plasmid_hallmarks', 'min_plasmid_hallmarks_short_seqs',
        'min_virus_hallmarks', 'min_virus_hallmarks_short_seqs', 'max_uscg',
        'max_integrase_distance', 'max_trna_distance', 'batch_size'
    ]
    check_nums(self, params, defaults, ints, int)
    floats1 = ['sensitivity', 'min_score', 'max_fdr', 'marker_threshold',
               'marker_threshold_integrase', 'marker_threshold_edge',
               'min_plasmid_marker_enrichment', 'min_virus_marker_enrichment']
    check_nums(self, params, defaults, floats1, float)
    floats2 = ['evalue', 'crf_threshold']
    check_nums(self, params, defaults, floats2, float, 0, 1)
    check_default(self, params, defaults, (ints + floats1 + floats2))
    return defaults


def check_hamronization(self, params):
    defaults = {
        'format': ['tsv', 'json'],
        'analysis_software_version': {
            'deeparg': '0.20',
            'staramr': '0.9.1',
            'abricate': '1.0.0',
            'amrfinderplus': '3.10.16'
        },
        'reference_database_version': {
            'deeparg': '3.2.5',
            'staramr': '2.0.2',
            'abricate': '3.2.5',
            'amrfinderplus': '2022-08-09.1'
        },
        'input_file_name': {
            'deeparg': 'nucleotide.sequences.fasta'
        },
        'reference_database_name': {}
    }
    tools = ['deeparg', 'abricate', 'abritamr', 'amrplusplus', 'ariba',
             'csstar', 'fargene', 'groot', 'kmerresistance', 'resfams',
             'resfinder', 'mykrobe', 'kmerresistance', 'pointfinder', 'rgi',
             'srax', 'srst2', 'staramr', 'tbprofiler']
    versions = ['analysis_software_version', 'reference_database_version']
    for version in versions:
        if version not in params:
            params[version] = defaults[version]
        else:
            if not isinstance(params[version], dict):
                sys.exit('[%s] Param "%s" not a dict (version per tool)' % (
                    self.name, version))
            for tool, v in params[version].items():
                if tool not in tools:
                    sys.exit('[%s] Tool "%s" not recognized (param "%s")' % (
                        self.name, tool, version))
                if '.' not in v:
                    sys.exit('[%s] Tool "%s" version invalid (no "."): "%s"' % (
                        self.name, tool, v))

    if 'input_file_name' not in params:
        params['input_file_name'] = defaults['input_file_name']

    exclude = versions + ['input_file_name', 'reference_database_name']
    check_default(self, params, defaults, exclude)
    return defaults


def check_metadmg(self, params):
    defaults = {
        'min_similarity_score': [None],
        'max_similarity_score': [None],
        'min_edit_dist': [None],
        'max_edit_dist': [None],
        'min_mapping_quality': 0,
        'custom_database': [False, True],
        'lca_rank': [None, 'family', 'genus', 'species'],
        'max_position': 15,
        'min_reads': 0,
        'weight_type': 1,
        'parallel_samples': 1,
        'cores_per_sample': 1,
        'damage_mode': ['lca', 'local', 'global'],
        'forward_only': [False, True],
        'bayesian': [False, True],
        'long_name': [False, True],
        'overwrite': [False, True],
        'all_reads': [False, True]
    }
    ints = ['min_edit_dist', 'max_edit_dist', 'min_mapping_quality',
            'max_position', 'min_reads', 'weight_type', 'parallel_samples',
            'cores_per_sample']
    check_nums(self, params, defaults, ints, int)
    floats = ['min_similarity_score', 'max_similarity_score']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults, (ints + floats))
    check_binary(self, params, defaults, 'binary')
    defaults['taxnames'] = '<Path to the (NCBI) names-mdmg.dmp>'
    defaults['taxnodes'] = '<Path to the (NCBI) nodes-mdmg.dmp>'
    defaults['acc2tax'] = '<Path to the (NCBI) acc2tax.gz>'
    defaults['binary'] = '<Path (+ option) needed to metaDMG-cpp binary>'
    return defaults


def check_flanker(self, params):
    defaults = {
        'database': ['ncbi', 'argannot', 'card', 'ecoh', 'ecoli_vf', 'megares',
                     'plasmidfinder', 'resfinder', 'vfdb'],
        'closest_match': [False, True],
        'flank': ['both', 'upstream', 'downstream'],
        'mode': ['default', 'mm', 'sm'],
        'circ': [False, True],
        'include_gene': [False, True],
        'verbose': 0,
        'window': 1000,
        'window_stop': 0,
        'window_step': 0,
        'cluster': [False, True],
        'threshold': 0.001,
        'kmer_length': 21,
        'sketch_size': 1000,
    }
    if 'gene' not in params:
        params['gene'] = []
    if 'list_of_genes' not in params:
        params['list_of_genes'] = ''
    int1 = ['verbose']
    check_nums(self, params, defaults, int1, int, 0, 2)
    int2 = ['window', 'window_stop', 'window_step',
            'kmer_length', 'sketch_size']
    check_nums(self, params, defaults, int2, int)
    floats = ['threshold']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults,
                  (int1 + int2 + floats + ['gene', 'list_of_genes']))
    defaults['gene'] = '<List of genes to target in the abricate db>'
    defaults['list_of_genes'] = '<File of genes to target in the abricate db>'
    return defaults


def check_plasx(self, params):
    defaults = {
        'splits': 1,
        'anvio_annot': [False, True],
        'save_tmp': [False, True],
        'databases': ['cog14', 'pfam32']
    }
    ints = ['splits']
    if 'databases' not in params:
        params['databases'] = defaults['databases']
    check_databases(self.name, params, self.databases)
    check_nums(self, params, defaults, ints, int)
    check_default(self, params, defaults, (ints + ['databases']))
    return defaults


def check_mobmess(self, params):
    defaults = {'min_similarity': 0.9, 'min_coverage': 0.9, 'min_length': 2000}
    if 'contigs_intersection' not in params:
        params['contigs_intersection'] = []
    if 'circularity_tool' not in params:
        params['circularity_tool'] = None
    if 'contigs_names' not in params:
        params['contigs_names'] = None
    else:
        if not isfile(params['contigs_names']):
            sys.exit('[mobmess] Param "contigs_names" is not a file')

    ints = ['min_length']
    check_nums(self, params, defaults, ints, int)
    floats = ['min_similarity', 'min_coverage']
    check_nums(self, params, defaults, floats, float, 0, 1)
    check_default(self, params, defaults,
                  (floats + ints + ['contigs_intersection', 'contigs_names']))
    defaults['contigs_intersection'] = '<Plasmid-detection tool for which to ' \
                                       'take the contigs intersection>'
    defaults['circularity_tool'] = '<Tool that yielded contig circularity info>'
    defaults['contigs_names'] = '<Path to a file containing contigs to use>'
    return defaults


def check_coverm(self, params):
    defaults = {
        'genome_fasta_extension': ['fna'],
        'separator': [None],
        'single_genome': [False, True],
        'use_full_contig_names': [False, True],
        'dereplicate': [False, True],
        'min_completeness': [None],
        'max_contamination': [None],
        'dereplication_ani': 95,
        'dereplication_aligned_fraction': 15,
        'dereplication_fragment_length': 3000,
        'dereplication_quality_formula': [
            'Parks2020_reduced', 'completeness-4contamination',
            'completeness-5contamination', 'dRep'],
        'dereplication_prethreshold_ani': 90,
        'dereplication_precluster_method': ['skani', 'dashing', 'finch'],
        'dereplication_cluster_method': ['skani', 'fastani'],
        'exclude_genomes_from_deshard': [False, True],
        'sharded': [False, True],
        'mapper': ['minimap2-sr', 'bwa-mem', 'bwa-mem2', 'minimap2-ont',
                   'minimap2-pb', 'minimap2-hifi', 'minimap2-no-preset'],
        'minimap2_params': [None],
        'minimap2_reference_is_index': [False, True],
        'bwa_params': [None],
        'strobealign_params': [None],
        'strobealign_use_index': [False, True],
        'min_read_aligned_length': 0,
        'min_read_percent_identity': 0,
        'min_read_aligned_percent': 0,
        'min_read_aligned_length_pair': 0,
        'min_read_percent_identity_pair': 0,
        'min_read_aligned_percent_pair': 0,
        'proper_pairs_only': [False, True],
        'exclude_supplementary': [False, True],
        'include_secondary': [False, True],
        'methods': ['relative_abundance', 'mean', 'trimmed_mean', 'tpm',
                    'coverage_histogram', 'covered_bases', 'variance',
                    'length', 'count', 'metabat', 'reads_per_base', 'rpkm'],
        'min_covered_fraction': 0,
        'contig_end_exclusion': 75,
        'trim_min': 5,
        'trim_max': 95,
        'output_format': ['sparse', 'dense'],
        'no_zeros': [False, True],
        'discard_unmapped': [False, True]
    }
    ints1 = ['min_covered_fraction', 'contig_end_exclusion', 'trim_min',
             'trim_max', 'min_read_percent_identity',
             'min_read_aligned_percent', 'min_read_percent_identity_pair',
             'min_read_aligned_percent_pair', 'dereplication_prethreshold_ani',
             'dereplication_aligned_fraction', 'dereplication_ani']
    check_nums(self, params, defaults, ints1, int, 0, 100)
    ints2 = ['min_read_aligned_length', 'min_read_aligned_length_pair',
             'dereplication_fragment_length']
    check_nums(self, params, defaults, ints2, int)
    check_default(self, params, defaults, (ints1 + ints2))
    defaults['genome_fasta_files'] = '<Path(s) to FASTA files of each genome>'
    defaults['genome_fasta_directory'] = '<Directory of genome FASTA files>'
    defaults['genome_fasta_list'] = '<FASTA file paths (one per line)>'
    defaults['reference'] = '<FASTA file of contigs>'
    defaults['genome_definition'] = '<list of "genome_name<tab>contig" lines>'
    defaults['checkm2_quality_report'] = '<CheckM version 2 quality_report.tsv>'
    defaults['checkm_tab_table'] = "<CheckM tab table (checkm's --tab_table)>"
    defaults['genome_info'] = '<dRep style genome info table defining quality>'
    return defaults


def check_carveme(self, params):
    defaults = {
        'dna': [False, True],
        'egg': [False, True],
        'refseq': [False, True],
        'diamond_args': [None],
        'recursive': [False, True],
        'universe': ['bacteria'],
        'universe_file': [None],
        'cobra': [False, True],
        'fbc2': [False, True],
        'ensemble': [False, True],
        'gapfill': [None],
        'init': [None],
        'mediadb': [None],
        'verbose': [False, True],
        'debug': [False, True],
        'soft': [False, True],
        'hard': [False, True],
        'reference': [None],
        'solver': ['cplex', 'gurobi']
    }
    check_default(self, params, defaults)
    check_binary(self, params, defaults, 'path')
    defaults['path'] = '<Path to folder with `carve` and `merge_community`>'
    return defaults


# def check_skani(self, params):
#     defaults = {
#     }
#     ints = []
#     check_nums(self, params, defaults, ints, int)
#     floats = []
#     check_nums(self, params, defaults, floats, float)
#     check_default(self, params, defaults, (ints + floats))
#     defaults[''] = '<>'
#     return defaults


# def check_SOFTNAME(self, params, soft):
#     # list here all the parameters and their default values
#     defaults = {
#         '<PARAM_1>': 10,
#         '<PARAM_2>': 0.6,
#         '<PARAM_3>': [True, False],
#         '<PARAM_4>': 'something,requiring,manual,check',
#         ...
#     }
#
#     # list the names of the parameters that have unbounded integers as values
#     ints = [<UNBOUNDED INTEGER-VALUE PARAMETER NAMES (from "defaults" DICT)>]
#     check_nums(self, params, defaults, ints, int, soft.name)
#
#     # list the names of the parameters that have unbounded floats as values
#     floats = [<UNBOUNDED FLOAT-VALUE PARAMETER NAMES (from "defaults" DICT)>]
#     check_nums(self, params, defaults, floats, float, soft.name)
#
#     # list the names of the parameters that have bounded integers as values
#     ints2 = [<BOUNDED INTEGER-VALUE PARAMETER NAMES (from "defaults" DICT)>]
#     # same as above but replace <m> and <M> by min and max integer values
#     check_nums(self, params, defaults, ints2, int, soft.name, <m>, <M>)
#
#     # list the names of the parameters that have bounded floats as values
#     floats2 = [<BOUNDED FLOAT-VALUE PARAMETER NAMES (from "defaults" DICT)>]
#     # same as above but replace <m> and <M> by min and max float values
#     check_nums(self, params, defaults, floats2, float, soft.name, <m>, <M>)
#
#     # make manual checks (see examples above - it will depend on the tool...)
#     manu = [<MANUALLY-CHECKED PARAMETER NAMES (from "defaults" DICT)>]
#
#     # make a list with all the above-checked parameters
#     let_go = ints + floats + ints2 + floats2 + manu
#
#     # add the "let_go" list at the end
#     check_default(self, params, defaults, soft.name, let_go)
#     # if the parameter value can be a list, get these params and add it too:
#     # multi = [<PARAMETERS THAT CAN BE A LIST>]
#     # check_default(self, params, defaults, soft.name, let_go, multi)
#
#     # Finally add to "defaults" dict those params that need not to be checked
#     defaults['<PARAM_NAME>'] = '<SOME USEFUL DESCRIPTION (for --show-params)>'
#     return defaults
