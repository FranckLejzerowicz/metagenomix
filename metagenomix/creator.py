# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import logging
from metagenomix.metagenomix import metagenomix
from metagenomix.core.jobs import Created


def set_log():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    # Using StreamHandler writing to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    # Add the two Handlers
    logger.addHandler(ch)


def creator(**kwargs):
    """Run metagenomix to account for the passed command-line arguments.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    set_log()
    kwargs['logging'] = logging
    logging.info('\n>>> `metagenomix create` started >>>\n')
    # Collect all command and init the script creating instance
    kwargs['command'] = 'create'
    creating = Created(*metagenomix(**kwargs))
    # print('* Writing database formatting commands')
    # creating.database_cmds()  # build the databases
    logging.info('* Creating output folders')
    creating.make_dirs()
    logging.info('* Writing pipeline command lines')
    creating.software_cmds()   # run the analysis pipeline
    if len(creating.run['database']) or len(creating.run['software']):
        m = '\n< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >'
        logging.info(m)
        creating.display()  # show the scripts to run
    creating.write_links()
    creating.write_logs()
    logging.info('\n<<< `metagenomix create` completed <<<\n')
