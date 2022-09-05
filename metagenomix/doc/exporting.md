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
    -i <PIPELINE_FOLDER> \
    -o /path/to/exports.txt \
    -s checkm2 \
    -s metawrap_*
```
where one export file will be created for all the outputs of `checkm2`, as 
well as `metawrap_annotate`, `metawrap_binning`, `metawrap_classify`,
`metawrap_reassemble`, and `metawrap_refine`.

### Extensions and regular expressions

It is possible to further subset the contents to export by targeting files 
ending with one or more specific extensions or matching regular expressions. 

#### Extensions

Targeting files with specific extension is done by using the `-e` (or
`--extension`) argument, possibly multiple times. 

For example:
```
metagenomix export \
    -i <PIPELINE_FOLDER> \
    -o /path/to/exports.txt \
    -s checkm2 \
    -e .tsv \
    -e .html
```
will export the `.tsv` and `.html` files present in the `checkm2` 
folder of the pipeline output.

#### Regular expressions

Targeting files matching regular expressions is done by using the `-r` (or
`--regex`) argument, possibly multiple times. 

For example:
```
metagenomix export \
    -i <PIPELINE_FOLDER> \
    -o /path/to/exports.txt \
    -s checkm2 \
    -r *table* \
    -r *ERR*
```
will export all the files that contains `table` or `ERR`.

## Output

Exporting using `metagenomix export` consists of copying and archiving target 
files while conserving the structure of the folders from the input locations.

### Archive

The path to the export file is **mandatory**, and must be provided using the 
`-o` (or `--output`) argument. If the filename does not have a `tar.gz` 
extension, it will be added. In the above example, `/path/to/exports.txt` 
will be replaced to `/path/to/exports.tar.gz`.  

By default (if no extension or regex is used to target files), then all the 
files, including the
[job output](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#job-outputs)
files will also be included in the export archive.

### Job script

Exporting many files can take a long time, and therefore, the default is to 
write the export commands into a job script, which is set to last for max. 2 
hours. It can be run by executing the bash script shown in the terminal, 
after `To export, run the following:`, which is composed of commands 
that will ensure that the actual job script outputs (`.o` and `.o`) will be 
written in the `_exports/jobs/output` sub-folders:
```
mkdir <PIPELINE_FOLDER>/_exports/jobs/output
cd <PIPELINE_FOLDER>/_exports/jobs/output
sbatch <PIPELINE_FOLDER>/_exports/jobs/export_DD-MM-YYYY_HH-MM-SS.slm
```

### Terminal display

For the example:
```
metagenomix export \
    -i <PIPELINE_FOLDER>  \
    -o /path/to/exports.txt \
    -s checkm* \
    -s metawrap*
```
the terminal will display explicit messages:
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
        > /path/to/exports.tar.gz
* Getting copying commmands
* Getting archiving commmands

* To export, run the following:
sh <PIPELINE_FOLDER>/_exports/jobs/export_DD-MM-YYYY_HH-MM-SS.sh

* Then, copy(-edit)-paste to download from server to local:
scp <USERNAME>@<MACHINE>:/path/to/exports.tar.gz .
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
  -e, --extension TEXT      Extension to select files
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