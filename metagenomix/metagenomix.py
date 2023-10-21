# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys

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
    logging = kwargs['logging']
    # Parse the command line arguments and validate computing environment
    config = AnalysesConfig(**kwargs)
    logging.info('* Reading configurations')
    config.run(logging)

    # Validate and get commands to set up the passed databases if need be
    databases = ReferenceDatabases(config)
    logging.info('* Checking databases')
    databases.run(logging)

    # Read and validate the workflow of softwares to run as a pipeline
    workflow = Workflow(config, databases)
    logging.info('* Checking pipeline')
    workflow.visit()
    logging.info('* Setting up pipeline graph')
    workflow.setup()
    logging.info('* Checking parameters')
    workflow.parametrize()
    logging.info('* Preparing workflow')
    workflow.prepare()

    if config.show_params:
        workflow.show_params(logging)
        sys.exit('Do not use "-s"/"--show-params" to proceed with create')

    # Collect command line and workflow of softwares to run as a pipeline
    commands = Commands(config, databases, workflow)
    logging.info('* Collecting command lines')
    commands.collect(logging)

    instances = (config, databases, workflow, commands)
    return instances
