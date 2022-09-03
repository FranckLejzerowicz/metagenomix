# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.merger import merger
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
@click.option("--remove/--no-remove", default=True,
              help="Enable assistance to manage folder(s)/file(s) removals")
@click.option("--jobs/--no-jobs", default=True,
              help="Enable assistance to manage job stdour/stderr files")
@click.option("--rename/--no-rename", default=True,
              help="Enable assistance to manage folder(s)/file(s) renaming")
@click.option("--rename/--no-rename", default=True,
              help="Enable assistance to manage folder(s)/file(s) renaming")
@click.option("--verbose/--no-verbose", default=True, show_default=True,
              help="Whether to show input/outputs and other details")
@click.version_option(__version__, prog_name="metagenomix")
def merge(
        folder,
        summary,
        pipeline,
        software,
        remove,
        jobs,
        rename,
        verbose
):
    """Combine the per-sample outputs into feature tables."""
    merger(
        folder=folder,
        out=summary,
        pipeline=pipeline,
        softwares=software,
        remove=remove,
        jobs=jobs,
        rename=rename,
        verbose=verbose
    )
