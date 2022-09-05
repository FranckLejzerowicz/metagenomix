# Exporting

This explains how to use the `metagenomix export` command to wrap up the 
current outputs of all or specific softwares of a shotgun metagenomics pipeline.

## Input

To export a set of files, the **mandatory** `-i` argument must be a path to a 
[pipeline output folder](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#output)
(i.e., the `-o` argument to `metagenomix create`).

### Softwares

By default, all the files in all the softwares present in the input folder 
will be exported. It is possible to restrict the export to the union of 
softwares names that can be provided using either/or:
* a [pipeline configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md)
  file, using the `-p` (or  `--pipeline`) argument.
* a list a softwares, using the `-s` (or  `--software`) arguments, possibly 
  multiple times and using regular expressions. 

For example, the following is a valid `metagenomix export` command:
```
metagenomix export \
    -i /path/to/pipeline/output \
    -o /path/to/exports.txt \
    -s checkm2 \
    -s metawrap_*
```
where one export file will be created for all the outputs of `checkm2`, as 
well as `metawrap_`, `metawrap_`, 

### Softwares


## Output

Exporting essentially consists of copying and archiving target files while
conserving the structure of the folders from the input location.

### Archived export

By default - if not extension is used to target files - then all the
[job output](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#job-outputs)
files will also be included.

### Job script

Exporting many files can take a long time, and therefore, it is key

### Terminal display

```
 === metagenomix exporter ===

* Getting softwares to export
Tools to export per role:
  [Role: binning]
        - metawrap_binning
        - metawrap_classify
        - metawrap_refine
        - metawrap_annotate
        - metawrap_reassemble
  [Role: MAG]
        - checkm2

* Getting files to export:
        141 files will be exported
* Getting output archive path
        > /Users/franck/programs/metagenomix/metagenomix/tests/example/agp/exports.tar.gz
* Getting copying commmands
* Getting archiving commmands

* To export, run the following:
sh <PIPELINE_FOLDER>/_exports/jobs/export_DD-MM-YYYY_HH-MM-SS.sh

* Then, copy(-edit)-paste to download from server to local:
scp <USERNAME>@<MACHINE>:/path/to/exports.tar.gz .



metagenomix export -i metagenomix/tests/example/agp/output -o metagenomix/tests/example/agp/exports.txt -s checkm* -s metawrap*
```

## Usage

```
Usage: metagenomix export [OPTIONS]

  Prepare an archive for specific pipeline outputs.

Options:
  -i, --folder TEXT         Path to pipeline output folder (`-o` for "create"
                            module)  [required]
  -o, --output TEXT         Output archive file  [required]
  -p, --pipeline TEXT       Path to the file containing the softwares to run
                            in order
  -s, --software TEXT       Software(s) to manage (or all in `-i/-p`)
  -e, --extensions TEXT     Files extensions to select
  -r, --regex TEXT          Regex for file names to select
  -l, --location TEXT       If not creating the tar locally, create it there
  --local / --no-local      Creates the tar locally, and not in the location
                            of `-l`
  -a, --account TEXT        User account for your HPC (in use for Slurm)
  --jobs / --no-jobs        Whether to make the archive in a job  [default:
                            jobs]
  --torque / --no-torque    Whether to prepare Torque jobs instead of Slurm
                            [default: no-torque]
  --verbose / --no-verbose  Whether to show input/outputs and other details
                            [default: verbose]
  --version                 Show the version and exit.
  --help                    Show this message and exit.

```