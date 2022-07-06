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
    """Run metagenomix to account for the passed commanbd line a.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """

    # Parse the command line arguments and validate computing environment
    config = AnalysesConfig(**kwargs)
    config.run()

    # Validate and get commands to set up the passed databases if need be
    databases = ReferenceDatabases(config)
    databases.run()

    # Read and validate the workflow of tools to run as a pipeline
    workflow = Workflow(config, databases)
    workflow.run()

    # Collect the command line and  the workflow of tools to run as a pipeline
    commands = Commands(config, databases, workflow)
    commands.run()

    workflow.make_dirs()

    # Make .sh and scheduler (.slm or .pbs) scripts to
    scripting = CreateScripts(config)
    scripting.database_cmds(databases)  # build the databases
    scripting.software_cmds(commands)   # run the analysis pipeline
    scripting.display()  # show the scripts to run
