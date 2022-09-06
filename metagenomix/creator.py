# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix.metagenomix import metagenomix
from metagenomix.core.jobs import Created


def creator(**kwargs):
    """Run metagenomix to account for the passed command-line arguments.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix creator ===\n')
    # Collect all command and init the script creating instance
    creating = Created(*metagenomix(**kwargs))
    # print('* Writing database formatting commands')
    # scripting.database_cmds()  # build the databases
    print('* Creating output folders')
    creating.make_dirs()
    print('* Writing pipeline command lines')
    creating.software_cmds()   # run the analysis pipeline
    if len(creating.run['database']) or len(creating.run['software']):
        print('< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >')
        creating.display()  # show the scripts to run
    creating.write_logs()
    print('\nCompleted.')
