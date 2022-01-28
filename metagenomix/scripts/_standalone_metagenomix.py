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
    "-k", "--p-scratch", show_default=True, help='Abs path to scratch folder.')
@click.option("-a2", "--p-a2", show_default=True)
@click.option("-a3", "--p-a3", show_default=True)
@click.option("-a4", "--p-a4", show_default=True)
@click.option("-a5", "--p-a5", show_default=True)
@click.option("-a6", "--p-a6", show_default=True)
@click.option("-a7", "--p-a7", show_default=True)
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Force the re-writing of scripts for all commands"
         "(default is to not re-run if output file exists).")
@click.option(
    "--jobs/--no-jobs", default=True, show_default=True,
    help="Whether to prepare Torque jobs from scripts.")
@click.option(
    "--slurm/--no-slurm", default=False, show_default=True,
    help="Whether to prepare Slurm and not Torque jobs.")
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
        p_scratch,
        p_a2,
        p_a3,
        p_a4,
        p_a5,
        p_a6,
        p_a7,
        force,
        jobs,
        slurm,
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
        scratch=p_scratch,
        a2=p_a2,
        a3=p_a3,
        a4=p_a4,
        a5=p_a5,
        a6=p_a6,
        a7=p_a7,
        force=force,
        jobs=jobs,
        slurm=slurm,
        verbose=verbose,
    )


if __name__ == "__main__":
    standalone_metagenomix()
