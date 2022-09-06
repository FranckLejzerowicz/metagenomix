# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from metagenomix.metagenomix import metagenomix
from metagenomix.core.jobs import Created


def merger(**kwargs):
    """Combine the per-sample outputs into feature tables.

    Parameters
    ----------
    kwargs : dict
        All arguments passed in command line, including defaults
    """
    print('\n === metagenomix exporter ===\n')
    # Collect all command and init the script creating instance
    merging = Created(*metagenomix(**kwargs))
    # print('* Writing database formatting commands')
    # scripting.database_cmds()  # build the databases
    print('* Writing pipeline command lines')
    merging.software_cmds()   # run the analysis pipeline
