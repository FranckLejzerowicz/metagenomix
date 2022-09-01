# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix.core.pipeline import Workflow
from metagenomix.core.config import AnalysesConfig
from metagenomix.core.databases import ReferenceDatabases
from metagenomix.core.commands import Commands


def metagenomix(**kwargs) -> tuple:
    """Run metagenomix to account for the passed command-line arguments.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults

    Returns
    -------
    instances : tuple
        Class instances for pipeline configuration, databases, workflow,
        and commands
    """

    # Parse the command line arguments and validate computing environment
    config = AnalysesConfig(**kwargs)
    print('* Reading configurations')
    config.run()

    # Validate and get commands to set up the passed databases if need be
    databases = ReferenceDatabases(config)
    print('* Checking databases')
    databases.run()

    # Read and validate the workflow of softwares to run as a pipeline
    workflow = Workflow(config, databases)
    print('* Reading pipeline')
    workflow.visit()
    print('* Setting up graph and output paths')
    workflow.setup()
    print('* Checking parameters')
    workflow.parametrize()

    # Collect command line and workflow of softwares to run as a pipeline
    commands = Commands(config, databases, workflow)
    print('* Collecting command lines')
    commands.collect()

    instances = (config, databases, workflow, commands)
    return instances
