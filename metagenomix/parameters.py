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
from os.path import isdir, isfile

from metagenomix._io_utils import read_yaml
from metagenomix.tools.alignment import *


def check_ints(param, value, soft):
    if not str(value).isdigit():
        raise IOError('Param "%s" for "%s" must be integer' % (
            param, soft.name))


def check_mems(param, value, soft):
    mems = ['kb', 'mb', 'gb']
    if value not in mems:
        raise IOError('Param "%s" for "%s" must be of %s' % (
            param, soft.name, str(mems)))


def check_env(config, value, soft):
    if soft.name not in config.modules and value not in config.conda_envs:
        raise EnvironmentError('"%s" not a module or conda env' % value)


def check_path(value):
    if not isfile(value) and not isdir(value):
        raise IOError('"%s" do not exist' % value)


def show_valid_params(param, values, soft):
    m = '[%s] user parameter "%s" must be among:\n' % (soft.name, param)
    for value in values:
        m += '  - %s' % value
    sys.exit(m)


def check_databases(name, user_params, databases):
    if 'databases' not in user_params:
        raise IOError('[%s] "databases" must be a parameter' % name)
    dbs_existing = []
    dbs_missing = []
    for db in user_params['databases']:
        if db in databases.paths:
            dbs_existing.append(db)
        else:
            dbs_missing.append(db)
    if not dbs_existing:
        raise IOError('[%s] No databases: %s' % (name, '; '.join(dbs_missing)))
    elif dbs_missing:
        print('[%s] Missing databases: "%s"' % (name, '", "'.join(dbs_missing)))
    return dbs_existing


def check_generic(defaults, user_params, soft, let_go: list = []):
    for param, values in defaults.items():
        if param in let_go:  # skip params checked specifically per software
            continue
        if param not in user_params:
            soft.params[param] = values[0]
        else:
            if isinstance(user_params[param], list):
                if set(sorted(user_params[param])).difference(values):
                    show_valid_params(param, values, soft)
            else:
                if user_params[param] not in values:
                    show_valid_params(param, values, soft)


# ============================================= #
#  Below are the function that are called from  #
#  class method Workflow.set_user_params()      #
#       (see module metagenomix.pipeline)       #
# ============================================= #


def check_shogun(user_params, soft, databases, config):
    defaults = {
        # default: 'paired'
        'pairing': ['paired', 'concat', 'single'],
        # default: ['bowtie2']
        'aligners': ['bowtie2', 'burst', 'utree']
    }
    check_generic(defaults, user_params, soft)
    if 1:
        valid_dbs = {}
        dbs_existing = check_databases('shogun', user_params, databases)
        for db in dbs_existing:
            path = databases.paths[db]
            metadata_yml = '%s/shogun/metadata.yaml' % path
            if not config.dev and not isfile(metadata_yml):
                raise IOError('[shogun] file "%s" must exist' % metadata_yml)
            metadata = read_yaml(metadata_yml)
            for aligner in list(user_params['aligners']):
                if aligner in metadata:
                    if not config.dev:
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
        user_params['databases'] = valid_dbs
        if not user_params['databases']:
            print('[shogun] No database formatted for shogun: will be skipped')


def check_bowtie2(user_params, soft, databases, config):
    defaults = {
        # default: 'paired'
        'pairing': ['paired', 'concat', 'single'],
        'discordant': [True, False],
        'k': '16', 'np': '1',
        'mp': '1,1', 'rdg': '0,1', 'rfg': '0,1',
        'score-min': 'L,0,-0.05'
    }
    let_go = []
    check_bowtie_k_np(user_params, defaults, let_go)
    check_bowtie_mp_rdg_rfg(user_params, defaults, let_go)
    check_bowtie_score_min(user_params, defaults, let_go)
    check_generic(defaults, user_params, soft, let_go)
    dbs_existing = check_databases('bowtie2', user_params, databases)
    valid_dbs = {}
    for db in dbs_existing:
        path = databases.paths[db]
        bt2_path = '%s/bowtie2/*.*.bt2' % path
        if not config.dev:
            bt2_paths = glob.glob(bt2_path)
            if bt2_paths:
                valid_dbs[db] = bt2_paths[0].rsplit('.', 2)[0]
        else:
            valid_dbs[db] = bt2_path.rsplit('.', 2)[0]
    user_params['databases'] = valid_dbs


def check_kraken2(user_params, soft, databases, config):
    defaults = {'databases': (['default'] + sorted(databases.paths))}
    check_generic(defaults, user_params, soft)


def check_spades(user_params, soft, databases, config):
    defaults = {'k': ['33', '55', '77', '99', '127'],
                'meta': [True, False], 'only_assembler': [True, False]}
    if 'k' not in user_params:
        user_params['k'] = defaults['k']
    else:
        kerrors = [x for x in user_params['k'] if not str(x).isdigit()]
        if len(kerrors):
            sys.exit('[spades] "k" must be integers (%s)' % ','.join(kerrors))
    check_generic(defaults, user_params, soft, ['k'])


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


def check_simka(user_params, soft, databases, config):
    defaults = {'simkaMin': [True, False],
                'kmer': np.linspace(15, 80, 6),
                'log_reads': np.logspace(3, 7, 3)}
    if 'kmer' not in user_params:
        user_params['kmer'] = defaults['kmer']
    else:
        kmer = user_params['kmer']
        simka_error('kmer', kmer)
        user_params['kmer'] = np.linspace(
            kmer['start'], kmer['end'], kmer['size'])

    if 'log_reads' not in user_params:
        user_params['log_reads'] = defaults['log_reads']
    else:
        reads = user_params['log_reads']
        simka_error('log_reads', reads)
        user_params['log_reads'] = np.logspace(
            reads['start'], reads['end'], reads['size'])

    check_generic(defaults, user_params, soft, ['kmer', 'log_reads'])


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}


# def check_TOOL(user_params, soft, databases, config):
#     default = {}
