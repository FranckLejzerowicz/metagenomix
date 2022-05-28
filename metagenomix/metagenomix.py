# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

from metagenomix.pipeline import Workflow
from metagenomix.config import AnalysesConfig
from metagenomix.databases import ReferenceDatabases
from metagenomix.commands import Commands
from metagenomix.jobs import CreateScripts


def metagenomix(**kwargs):
    """Main metagenomics workflow.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """

    # Collect general config and initialization checks
    config = AnalysesConfig(**kwargs)
    config.init()
    # print()
    # print('**************** CONFIG *****************')
    # for i, j in config.__dict__.items():
    #     print(i, '\t:\t', j)
    # print('**************** CONFIG *****************')
    # print()

    databases = ReferenceDatabases(config)
    databases.init()
    # print()
    # print('**************** database ***************')
    # for db in databases.databases_commands:
    #     for sam, cmds in databases.databases_commands[db].items():
    #         print('\n'.join(cmds))
    # print('*****************************************')

    workflow = Workflow(config)
    workflow.init()
    # print('**************** pipeline ***************')
    # for idx, soft in workflow.softs.items():
    #     print(idx, soft.prev, soft.name, soft.params)
    # print('*****************************************')
    # print()

    commands = Commands(config, databases, workflow)
    commands.collect()
    workflow.make_dirs()
    # print()
    # print('**************** analyses ***************')
    # for soft in commands.analyses_commands:
    #     for sam, cmds in commands.analyses_commands[soft].items():
    #         print('\n'.join(cmds))
    # print('*****************************************')

    scripting = CreateScripts(config)
    scripting.database_cmds(databases)
    scripting.software_cmds(commands)
