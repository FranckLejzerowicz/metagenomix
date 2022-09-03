# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix.metagenomix import metagenomix
from metagenomix.core.jobs import CreateScripts


def merger(**kwargs):
    """Combine the per-sample outputs into feature tables.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix exporter ===\n')
    # Collect all command and init the script creating instance
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
    scripting.get_hash()
    scripting.versioning()
    print('\nCompleted.')
