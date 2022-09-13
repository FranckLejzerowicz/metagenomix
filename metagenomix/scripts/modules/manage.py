# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.manager import manager
from metagenomix import __version__


@click.command()
@click.option("-i", "--folder", required=True,
              help='Path to pipeline output folder (`-o` for "create" module)')
@click.option("-o", "--summary",
              help='Output summary filename (will be in "<`-i`>/_managemnent")')
@click.option("-p", "--pipeline",
              help="Path to the file containing the softwares to run in order")
@click.option("-s", "--software", multiple=True,
              help="Software(s) to manage (or all in `-i/-p`)")
@click.option("--jobs/--no-jobs", default=True,
              help="[Management task] Enable job output management")
@click.option("--remove/--no-remove", default=False,
              help="[Management task] Enable completed output removal")
@click.option("--rename/--no-rename", default=False,
              help="[Management task] Enable output renaming")
@click.option("--store/--no-store", default=False,
              help="[Management task] Enable output storage")
@click.option("--confirm/--no-confirm", default=True, show_default=True,
              help="Whether to ask for confirmation before applying task")
@click.version_option(__version__, prog_name="metagenomix")
def manage(
        folder,
        summary,
        pipeline,
        software,
        jobs,
        remove,
        rename,
        store,
        confirm
):
    """Edit the contents of your pipeline output folder."""
    manager(
        dir=folder,
        out=summary,
        pipeline=pipeline,
        softwares=software,
        jobs=jobs,
        remove=remove,
        rename=rename,
        store=store,
        confirm=confirm
    )
