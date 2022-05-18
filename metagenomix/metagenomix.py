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
    workflow = Workflow(config)
    workflow.init()
    # print('**************** pipeline ***************')
    # for idx, soft in workflow.softs.items():
    #     print(idx, soft.prev, soft.name, soft.params)
    # print('*****************************************')
    # print()
    cmds = Commands(config, databases, workflow)
    cmds.collect()
    workflow.make_dirs()

    if 1:
        print('**************** commands ***************')
        for softs in config.pipeline:
            print()
            print(softs[-1])
            # print(workflow.softs[softs[-1]].outputs)
            for sam, cmds in workflow.softs[softs[-1]].cmds.items():
                print(sam)
                print('\n'.join(cmds))
        print('*****************************************')

    # scripting = CreateScripts(config)
    # scripting.write_scripts(databases.databases_commands)
    # scripting.write_scripts(analysis.analyses_commands)
