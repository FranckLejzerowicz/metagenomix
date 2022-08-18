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
@click.option("-i", "--folder", required=True, help="Path to pipeline folder")
@click.option("-o", "--summary", help="Summary of changes")
@click.option("-p", "--pipeline", help="Softwares making the pipeline")
@click.option("-t", "--tools", help="Software(s) to manage (or all in `-p`)")
@click.option("--remove/--no-remove", default=True, help="Enable removals")
@click.option("--jobs/--no-jobs", default=True, help="Enable jobs clearing")
@click.option("--rename/--no-rename", default=True, help="Enable renaming")
@click.version_option(__version__, prog_name="metagenomix")
def manage(
        folder,
        summary,
        pipeline,
        tools,
        remove,
        jobs,
        rename,
):

    manager(
        dir=folder,
        out=summary,
        pipeline=pipeline,
        tools=tools,
        remove=remove,
        jobs=jobs,
        rename=rename,
    )
