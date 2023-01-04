:construction:

# Monitoring

This explains how to use the `metagenomix monitor` command to follow the 
progress of a 
[created](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#creating)
pipeline. 

## Usage

```
Usage: metagenomix monitor [OPTIONS]

  Check IO/job status of your pipeline configuration.

Options:
  -i, --fastq-illumina TEXT       Path to short Illumina reads fastq files
                                  folder(s)
  -j, --fastq-pacbio TEXT         Path to long PacBio reads fastq files
                                  folder(s)
  -k, --fastq-nanopore TEXT       Path to long MinION Nanopore reads fastq
                                  files folder
  -o, --out-dir TEXT              Path to pipeline output folder  [required]
  -s, --summary TEXT              Path to summary output file
  -m, --metadata TEXT             Path to the metadata file  [required]
  -d, --databases TEXT            Databases (yaml file)
  -p, --pipeline TEXT             Softwares to pipeline
  -u, --user-params TEXT          Server run parameters
  -c, --co-assembly TEXT          Metadata column(s) defining the co-assembly
                                  groups (yaml file)
  -t, --strains TEXT              Species for strain-level analyses (yaml
                                  file)
  --show-params / --no-show-params
                                  Show all possible parameters for all
                                  softwares of your pipeline
  --force / --no-force            Check as if the re-writing of scripts for
                                  all commands was planned  [default: no-
                                  force]
  --verbose / --no-verbose        Whether to show input/outputs and other
                                  details  [default: no-verbose]
  --dev / --no-dev                For development...  [default: no-dev]
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```