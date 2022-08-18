# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix import __version__

from metagenomix.scripts.modules import create
from metagenomix.scripts.modules import check
from metagenomix.scripts.modules import manage


@click.group(help="Metagenomix command line manager")
@click.version_option(__version__, prog_name="metagenomix")
def standalone_metagenomix():
    pass


standalone_metagenomix.add_command(create.create)
standalone_metagenomix.add_command(check.check)
standalone_metagenomix.add_command(manage.manage)


if __name__ == '__main__':
    standalone_metagenomix()
