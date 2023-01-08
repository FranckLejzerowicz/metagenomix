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
    "-i", "--fastq-illumina", multiple=True,
    help="Path to short Illumina reads fastq files folder(s)")
@click.option(
    "-j", "--fastq-pacbio", multiple=True,
    help="Path to long PacBio reads fastq files folder(s)")
@click.option(
    "-k", "--fastq-nanopore", multiple=True,
    help="Path to long MinION Nanopore reads fastq files folder")
@click.option(
    "-o", "--out-dir", required=True, help="Path to pipeline output folder")
@click.option(
    "-r", "--result", help="Path to current results output file")
@click.option(
    "-s", "--software", multiple=True, help="Software(s) to manage (or 'all')")
@click.option(
    "-m", "--metadata", required=True, help="Path to the metadata file")
@click.option(
    "-d", "--databases", show_default=True, help="Databases (yaml file)")
@click.option(
    "-p", "--pipeline", show_default=True, help="Softwares to pipeline")
@click.option(
    "-u", "--user-params", show_default=True, help="Server run parameters")
@click.option(
    "-c", "--co-assembly", show_default=True,
    help="Metadata column(s) defining the co-assembly groups (yaml file)")
@click.option(
    "-t", "--strains", show_default=True,
    help="Species for strain-level analyses (yaml file)")
@click.option(
    "--show-params/--no-show-params", default=False, show_default=False,
    help="Show all possible parameters for all softwares of your pipeline")
@click.option(
    "--jobs/--no-jobs", default=False, show_default=True,
    help="Monitor the mem/cpu usage of previous jobs")
@click.option(
    "--dev/--no-dev", default=False, show_default=True,
    help="For development...")
@click.version_option(__version__, prog_name="metagenomix")
def monitor(
        metadata,
        fastq_illumina,
        fastq_pacbio,
        fastq_nanopore,
        out_dir,
        result,
        software,
        co_assembly,
        databases,
        pipeline,
        user_params,
        strains,
        show_params,
        jobs,
        dev
):
    """Check IO/job status of your pipeline configuration."""
    monitoring(
        meta_fp=metadata,
        illumina_dirs=fastq_illumina,
        pacbio_dirs=fastq_pacbio,
        nanopore_dirs=fastq_nanopore,
        output_dir=out_dir,
        result_fp=result,
        softwares=software,
        coassembly_yml=co_assembly,
        databases_yml=databases,
        pipeline_tsv=pipeline,
        user_params_yml=user_params,
        strains_yml=strains,
        show_params=show_params,
        jobs=jobs,
        dev=dev
    )
