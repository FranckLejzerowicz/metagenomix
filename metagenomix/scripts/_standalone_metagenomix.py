# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.metagenomix import metagenomix
from metagenomix import __version__


@click.command()
@click.option(
    "-i", "--i-fastq-dir", required=True, multiple=True,
    help="Path to fastq files folder.")
@click.option(
    "-m", "--m-metadata", required=True, help="Path to the metadata file.")
@click.option(
    "-o", "--o-out-dir", required=True, help="Path to pipeline output folder.")
@click.option(
    "-n", "--p-project-name", required=True, help="Name for your project.")
@click.option(
    "-d", "--p-databases", show_default=True, help="Databases (yaml file).")
@click.option(
    "-p", "--p-pipeline", show_default=True, help="Tools to pipeline.")
@click.option(
    "-u", "--p-run-params", show_default=True, help="Server run parameters.")
@click.option(
    "-c", "--p-co-assembly", show_default=True,
    help="Metadata column(s) defining the co-assembly groups (yaml file).")
@click.option(
    "-s", "--p-strains", show_default=True, help="Species for strain-level"
                                                 " analyses (yaml file).")
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists).")
@click.option(
    "--jobs/--no-jobs", default=True, show_default=True,
    help="Whether to prepare Torque jobs from scripts.")
@click.option(
    "--torque/--no-torque", default=False, show_default=True,
    help="Whether to prepare Torque jobs instead of Slurm.")
@click.option(
    "-l", "--localscratch", type=int, show_default=False, default=None,
    help="Use localscratch with the provided memory amount (in GB)")
@click.option(
    "--scratch/--no-scratch", default=False, show_default=True,
    help="Use the scratch folder to move files and compute")
@click.option(
    "--userscratch/--no-userscratch", default=False, show_default=True,
    help="Use the userscratch folder to move files and compute")
@click.option(
    "--verbose/--no-verbose", default=False, show_default=True,
    help="Whether to show expected input and outputs and other tools' details.")
@click.version_option(__version__, prog_name="metagenomix")


def standalone_metagenomix(
        m_metadata,
        i_fastq_dir,
        o_out_dir,
        p_project_name,
        p_co_assembly,
        p_databases,
        p_pipeline,
        p_run_params,
        p_strains,
        force,
        jobs,
        torque,
        localscratch,
        scratch,
        userscratch,
        verbose,
):

    metagenomix(
        meta_fp=m_metadata,
        fastq_dirs=i_fastq_dir,
        output_dir=o_out_dir,
        project=p_project_name,
        coassembly_yml=p_co_assembly,
        databases_yml=p_databases,
        pipeline_tsv=p_pipeline,
        user_params_yml=p_run_params,
        strains_yml=p_strains,
        force=force,
        jobs=jobs,
        torque=torque,
        localscratch=localscratch,
        scratch=scratch,
        userscratch=userscratch,
        verbose=verbose,
    )


if __name__ == "__main__":
    standalone_metagenomix()
