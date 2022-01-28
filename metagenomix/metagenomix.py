# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

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

    databases = ReferenceDatabases(config)
    databases.init()

    workflow = Workflow(config)
    workflow.init()

    cmds = Commands(config, databases, workflow)
    cmds.collect()

    # scripting = CreateScripts(config)
    # scripting.write_scripts(databases.databases_commands)
    # scripting.write_scripts(analysis.analyses_commands)
