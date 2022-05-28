# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import isdir, isfile


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
        raise EnvironmentError('"%s" not a conda env' % value)


def check_path(value):
    if not isfile(value) and not isdir(value):
        raise IOError('"%s" do not exist' % value)


def check_soft_params():
    pass


def set_user_params(config, soft):
    ints = ['time', 'procs', 'mem_num', 'chunks']
    soft_params = config.user_params.get(soft.name, {})
    for param, value in soft_params.items():
        soft.params[param] = value
        if param in ints:
            check_ints(param, value, soft)
        elif param == 'mem_dim':
            check_mems(param, value, soft)
        elif param == 'env':
            check_env(config, value, soft)
        elif param == 'path':
            check_path(value)
        else:
            check_soft_params()
        soft.params[param] = value
