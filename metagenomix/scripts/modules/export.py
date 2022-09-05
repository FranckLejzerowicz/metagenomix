# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.exporter import exporter
from metagenomix import __version__


@click.command()
@click.option("-i", "--folder", required=True,
              help='Path to pipeline output folder (`-o` for "create" module)')
@click.option("-o", "--output", required=True, help="Output archive file")
@click.option("-p", "--pipeline",
              help="Path to the file containing the softwares to run in order")
@click.option("-s", "--software", multiple=True,
              help="Software(s) to manage (or all in `-i/-p`)")
@click.option("-e", "--extension", multiple=True, default=None,
              show_default=True, help="Extension to select files")
@click.option("-r", "--regex", multiple=True, default=None,
              show_default=False, help="Regex for file names to select")
@click.option("-l", "--location", default='USERWORK',
              help='If not creating the tar locally, create it there')
@click.option("--local/--no-local", default=False,
              help="Creates the tar locally, and not in the location of `-l`")
@click.option("-a", "--account", show_default=False, default=None,
              help="User account for your HPC (in use for Slurm)")
@click.option("--jobs/--no-jobs", default=True, show_default=True,
              help="Whether to make the archive in a job")
@click.option("--torque/--no-torque", default=False, show_default=True,
              help="Whether to prepare Torque jobs instead of Slurm")
@click.option("--verbose/--no-verbose", default=True, show_default=True,
              help="Whether to show input/outputs and other details")
@click.version_option(__version__, prog_name="metagenomix")
def export(
        folder,
        output,
        pipeline,
        software,
        extension,
        regex,
        location,
        local,
        account,
        jobs,
        torque,
        verbose
):
    """Prepare an archive for specific pipeline outputs."""
    exporter(
        folder=folder,
        out=output,
        pipeline=pipeline,
        softwares=software,
        exts=extension,
        regex=regex,
        location=location,
        local=local,
        account=account,
        jobs=jobs,
        torque=torque,
        verbose=verbose
    )
