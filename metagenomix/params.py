# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import isdir, isfile


def set_user_params(config, soft):
    mem_dim = ['kb', 'mb', 'gb']
    integer_params = ['time', 'procs', 'mem_num', 'chunks']
    soft_params = config.user_params.get(soft.name, {})
    for param, value in soft_params.items():
        print(param, value)
        soft.params[param] = value
        if param in integer_params:
            if not str(value).isdigit():
                raise IOError('Param "%s" for "%s" must be integer' % (
                    param, soft.name))
        elif param == 'mem_dim':
            if value not in mem_dim:
                raise IOError('Param "%s" for "%s" must be of %s' % (
                    param, soft.name, str(mem_dim)))
        elif param == 'env':
            if value in config.conda_envs:
                raise EnvironmentError('"%s" not a conda env' % value)
        elif param == 'path':
            if not isfile(value) and not isdir(value):
                raise IOError('"%s" do not exist' % value)
        soft.params[param] = value
