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
@click.option("-o", "--output", default=None, show_default=True,
              help='Output folder (will be in "<`-i`>/_management")')
@click.option("-p", "--pipeline", default=None, show_default=True,
              help="Path to the file containing the softwares to run in order")
@click.option("-s", "--software", multiple=True,
              help="Software(s) to manage (or all in `-i/-p`)")
@click.option("-d", "--disk", default='',
              help="Path to a storage disk for `--store` task")
@click.option("-x", "--chunks", type=int, show_default=False, default=None,
              help="Number of scripts for each `--store` task")
@click.option("--jobs/--no-jobs", show_default=True,
              help="[Task] Enable job output management")
@click.option("--remove/--no-remove", default=False, show_default=True,
              help="[Task] Enable output removal")
@click.option("--rename/--no-rename", default=False, show_default=True,
              help="[Task] Enable output renaming")
@click.option("--store/--no-store", default=False, show_default=True,
              help="[Task] Enable output storage")
@click.option("--confirm/--no-confirm", default=True, show_default=True,
              help="Whether to ask for confirmation before applying task")
@click.version_option(__version__, prog_name="metagenomix")
def manage(
        folder,
        output,
        pipeline,
        software,
        disk,
        chunks,
        jobs,
        remove,
        rename,
        store,
        confirm
):
    """Deal with the contents of your pipeline output folder."""
    manager(
        dir=folder,
        out=output,
        pipeline=pipeline,
        softwares=software,
        disk=disk,
        chunks=chunks,
        jobs=jobs,
        remove=remove,
        rename=rename,
        store=store,
        confirm=confirm
    )
