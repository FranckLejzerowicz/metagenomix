# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from metagenomix.monitoring import monitoring
from metagenomix import __version__


@click.command()
@click.option(
    "-i", "--i-fastq-illumina", multiple=True,
    help="Path to short Illumina reads fastq files folder(s).")
@click.option(
    "-j", "--i-fastq-pacbio", multiple=True,
    help="Path to long PacBio reads fastq files folder(s).")
@click.option(
    "-k", "--i-fastq-nanopore", multiple=True,
    help="Path to long MinION Nanopore reads fastq files folder.")
@click.option(
    "-m", "--m-metadata", required=True, help="Path to the metadata file.")
@click.option(
    "-o", "--o-out-dir", required=True, help="Path to pipeline output folder.")
@click.option(
    "-s", "--o-summary", required=True, help="Path to summary output file.")
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
    "-t", "--p-strains", show_default=True,
    help="Species for strain-level analyses (yaml file).")
@click.option(
    "--show-params/--no-show-params", default=False, show_default=False,
    help="Show all possible parameters for all tools of your pipeline")
@click.option(
    "--force/--no-force", default=False, show_default=True,
    help="Check as if the re-writing of scripts for all commands was planned")
@click.option(
    "--verbose/--no-verbose", default=False, show_default=True,
    help="Whether to show input/outputs and other details")
@click.option(
    "--dev/--no-dev", default=False, show_default=True,
    help="For development...")
@click.version_option(__version__, prog_name="metagenomix")
def monitor(
        m_metadata,
        i_fastq_illumina,
        i_fastq_pacbio,
        i_fastq_nanopore,
        o_out_dir,
        o_summary,
        p_co_assembly,
        p_databases,
        p_pipeline,
        p_run_params,
        p_strains,
        show_params,
        force,
        verbose,
        dev
):
    """Check IO/job status of your pipeline configuration."""
    monitoring(
        meta_fp=m_metadata,
        illumina_dirs=i_fastq_illumina,
        pacbio_dirs=i_fastq_pacbio,
        nanopore_dirs=i_fastq_nanopore,
        output_dir=o_out_dir,
        summary_fp=o_summary,
        coassembly_yml=p_co_assembly,
        databases_yml=p_databases,
        pipeline_tsv=p_pipeline,
        user_params_yml=p_run_params,
        strains_yml=p_strains,
        show_params=show_params,
        force=force,
        verbose=verbose,
        dev=dev
    )
