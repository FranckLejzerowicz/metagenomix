# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import glob
from os.path import isdir, isfile

from metagenomix._io_utils import read_yaml


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


def check_generic(defaults, user_params, soft):
    for param, values in defaults.items():
        if param not in user_params:
            soft.params[param] = values[0]
        else:
            if isinstance(user_params[param], list):
                if set(sorted(user_params[param])).difference(values):
                    show_valid_params(param, values, soft)
            else:
                if user_params[param] not in values:
                    show_valid_params(param, values, soft)


def check_databases(name, user_params, databases):
    if 'databases' not in user_params:
        raise IOError('[%s] "databases" must be a parameter' % name)
    for db in user_params['databases']:
        if db not in databases.paths:
            raise IOError('[%s] database "%s" must have a path' % (name, db))


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
        if not config.dev:
            check_databases('shogun', user_params, databases)
        for db in user_params['databases']:
            path = databases.paths[db]
            metadata_yml = '%s/shogun/metadata.yaml' % path
            if not config.dev and not isfile(metadata_yml):
                raise IOError('[shogun] database "%s" must have a path' % db)
            metadata = read_yaml(metadata_yml)
            final_aligners = []
            for aligner in list(user_params['aligners']):
                if aligner in metadata:
                    if not config.dev:
                        ali = metadata[aligner]
                        if ali[0] == '/':
                            aligner_db = ali
                        else:
                            aligner_db = '%s/shogun/%s' % (path, ali)
                        if glob.glob('%s*' % aligner_db):
                            final_aligners.append(aligner)
                    else:
                        final_aligners.append(aligner)
            user_params['aligners'] = final_aligners


def check_bowtie2(user_params, soft, databases, config):
    defaults = {
        # default: 'paired'
        'pairing': ['paired', 'concat', 'single'],
        'header': ['yes', 'no'],
    }
    check_generic(defaults, user_params, soft)


def check_kraken2(user_params, soft, databases, config):
    defaults = {'databases': (['default'] + sorted(databases.paths))}
    check_generic(defaults, user_params, soft)
