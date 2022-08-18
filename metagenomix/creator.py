# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix.metagenomix import metagenomix
from metagenomix.jobs import CreateScripts


def creator(**kwargs):
    """Run metagenomix to account for the passed command-line arguments.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    # Make .sh and scheduler (.slm or .pbs) scripts to
    scripting = CreateScripts(*metagenomix(**kwargs))
    # print('* Writing database formatting commands')
    # scripting.database_cmds()  # build the databases
    print('* Creating output folders')
    scripting.make_dirs()
    print('* Writing pipeline command lines')
    scripting.software_cmds()   # run the analysis pipeline
    if len(scripting.run['database']) or len(scripting.run['software']):
        print('< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY >')
        scripting.display()  # show the scripts to run
