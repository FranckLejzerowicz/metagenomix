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
    "-M", "--p-modules", show_default=True, help="modules to use per software"
                                                 " analyses (yaml file).")
@click.option(
    "-a", "--p-account", show_default=False, default=None,
    help="User account on HPC")
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
    "--show-params/--no-show-params", default=False, show_default=False,
    help="Show all possible parameters for all tools of your pipeline.")
@click.option(
    "--show-pfams/--no-show-pfams", default=False, show_default=False,
    help="Show terms for which Pfam HMM models were already extracted before.")
@click.option(
    "--purge-pfams", "--no-purge-pfams", default=None, show_default=False,
    help="Remove terms for Pfam HMM models that were already extracted before.")
@click.option(
    "--verbose/--no-verbose", default=False, show_default=True,
    help="Whether to show expected input and outputs and other tools' details.")
@click.option(
    "--dev/--no-dev", default=False, show_default=True,
    help="For development...")
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
        p_modules,
        p_account,
        force,
        jobs,
        torque,
        localscratch,
        scratch,
        userscratch,
        show_params,
        show_pfams,
        purge_pfams,
        verbose,
        dev
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
        modules_yml=p_modules,
        force=force,
        jobs=jobs,
        torque=torque,
        account=p_account,
        localscratch=localscratch,
        scratch=scratch,
        userscratch=userscratch,
        show_params=show_params,
        show_pfams=show_pfams,
        purge_pfams=purge_pfams,
        verbose=verbose,
        dev=dev
    )


if __name__ == "__main__":
    standalone_metagenomix()
