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
    print('* Reading configurations')
    config.run()

    # Validate and get commands to set up the passed databases if need be
    databases = ReferenceDatabases(config)
    print('* Checking databases')
    databases.run()

    # Read and validate the workflow of tools to run as a pipeline
    workflow = Workflow(config, databases)
    print('* Reading pipeline')
    workflow.visit()
    print('* Setting up graph and output paths')
    workflow.setup()
    print('* Checking parameters')
    workflow.parametrize()

    # Collect the command line and  the workflow of tools to run as a pipeline
    commands = Commands(config, databases, workflow)
    print('* Collecting command lines')
    commands.run()
    print('* Creating output folders')
    commands.make_dirs()

    # Make .sh and scheduler (.slm or .pbs) scripts to
    scripting = CreateScripts(config, workflow, databases, commands)
    # print('* Writing database formatting commands')
    # scripting.database_cmds()  # build the databases
    print('* Writing pipeline command lines')
    scripting.software_cmds()   # run the analysis pipeline
    if len(scripting.run['database']) or len(scripting.run['software']):
        print('< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >')
        scripting.display()  # show the scripts to run
