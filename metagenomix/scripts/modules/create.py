# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.creator import creator
from metagenomix import __version__


@click.command()
@click.option(
    "-i", "--fastq-dir-illumina", multiple=True,
    help="Path to short Illumina reads fastq files folder(s)")
@click.option(
    "-j", "--fastq-dir-pacbio", multiple=True,
    help="Path to long PacBio reads fastq files folder(s)")
@click.option(
    "-k", "--fastq-dir-nanopore", multiple=True,
    help="Path to long MinION Nanopore reads fastq files folder")
@click.option(
    "-o", "--output-dir", required=True, help="Path to pipeline output folder")
@click.option(
    "-m", "--metadata", required=True, help="Path to the metadata file")
@click.option(
    "-p", "--pipeline", required=True, show_default=True,
    help="Path to the file containing the softwares to run in order")
@click.option(
    "-d", "--databases", show_default=True, help="Databases (yaml file)")
@click.option(
    "-u", "--user-params", show_default=True,
    help="Parameters for the softwares of the pipeline (yaml file)")
@click.option(
    "-c", "--co-assembly", show_default=True,
    help="Metadata column(s) defining the co-assembly groups (yaml file)")
@click.option(
    "-t", "--strains", show_default=True,
    help="Species for strain-level analyses (yaml file)")
@click.option(
    "-n", "--project-name", required=True, help="Name for your project")
@click.option(
    "-a", "--account", show_default=False, default=None,
    help="User account for your HPC (in use for Slurm)")
@click.option(
    "-M", "--modules", show_default=True,
    help="modules to use per software analyses (yaml file)")
@click.option(
    "-x", "--chunks", type=int, show_default=False, default=None,
    help="Number of jobs to split the commands into for each tool")
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists)")
@click.option(
    "--jobs/--no-jobs", default=True, show_default=True,
    help="Whether to prepare Torque jobs from scripts")
@click.option(
    "--torque/--no-torque", default=False, show_default=True,
    help="Whether to prepare Torque jobs instead of Slurm")
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
    "--move-back/--no-move-back", default=True, show_default=True,
    help="Do not move back from scrach (makes sense only for --userscratch)")
@click.option(
    "--show-params/--no-show-params", default=False, show_default=False,
    help="Show all possible parameters for each software")
@click.option(
    "--show-status/--no-show-status", default=False, show_default=False,
    help="Show status (needed inputs, done/to do outputs) for each software")
@click.option(
    "--show-pfams/--no-show-pfams", default=False, show_default=False,
    help="Show terms for which Pfam HMM models were already extracted before")
@click.option(
    "--purge-pfams/--no-purge-pfams", default=False, show_default=False,
    help="Remove terms for Pfam HMM models that were already extracted before")
@click.option(
    "--verbose/--no-verbose", default=False, show_default=True,
    help="Whether to show input/outputs and other details")
@click.option(
    "--dev/--no-dev", default=False, show_default=True,
    help="For development...")
@click.version_option(__version__, prog_name="metagenomix")

def create(
        fastq_dir_illumina,
        fastq_dir_pacbio,
        fastq_dir_nanopore,
        output_dir,
        metadata,
        pipeline,
        databases,
        user_params,
        co_assembly,
        strains,
        project_name,
        modules,
        account,
        chunks,
        force,
        jobs,
        torque,
        localscratch,
        scratch,
        userscratch,
        move_back,
        show_params,
        show_status,
        show_pfams,
        purge_pfams,
        verbose,
        dev
):
    """Write jobs for your pipeline configuration."""
    creator(
        illumina_dirs=fastq_dir_illumina,
        pacbio_dirs=fastq_dir_pacbio,
        nanopore_dirs=fastq_dir_nanopore,
        output_dir=output_dir,
        meta_fp=metadata,
        pipeline_tsv=pipeline,
        databases_yml=databases,
        user_params_yml=user_params,
        coassembly_yml=co_assembly,
        strains_yml=strains,
        modules_yml=modules,
        project=project_name,
        force=force,
        jobs=jobs,
        torque=torque,
        account=account,
        chunks=chunks,
        localscratch=localscratch,
        scratch=scratch,
        userscratch=userscratch,
        move_back=move_back,
        show_params=show_params,
        show_status=show_status,
        show_pfams=show_pfams,
        purge_pfams=purge_pfams,
        verbose=verbose,
        dev=dev
    )
