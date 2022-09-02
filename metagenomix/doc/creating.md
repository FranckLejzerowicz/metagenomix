# Creating

This explains how to create and run the analyses of a shotgun metagenomics 
pipeline `netagenomix create`.

## Input

### Fastq file folders

At least one folder containing .fastq (or .fastq.gz) files must be provided 
in command line, using one or several of the following arguments:

* `-i` (or `--fastq-dir-illumina`): Path to short Illumina reads fastq files folder(s).
* `-j` (or `--fastq-dir-pacbio`): Path to long PacBio reads fastq files folder(s).
* `-k` (or `--fastq-dir-nanopore`): Path to long MinION Nanopore reads fastq files folder-

It is perfectly possible to provide more than one of these arguments if the 
Illumina and PacBio and/or Nanopore data must be analysed altogether (the 
pipeline accounts for which software is able to handle which data type), and 
it is also possible to provide multiple folders per technology if the data is
located in multiple locations (whether for samples from different projects, 
or for the same re-sequenced samples), e.g.:

* inputs
    ```
    ├── illumina_fastqs_folder_1
    │   ├── sample1_1.fastq.gz
    │   ├── sample1_2.fastq.gz
    │   ├── sample2_1.fastq.gz
    │   └── sample2_2.fastq.gz
    ├── illumina_fastqs_folder_2
    │   ├── sample1_R1.fastq.gz
    │   ├── sample1_R2.fastq.gz
    │   ├── sample2_R1.fastq.gz
    │   ├── sample2_R2.fastq.gz
    │   ├── sample3_R1.fastq.gz
    │   └── sample3_R2.fastq.gz
    ├── pacbio_fastqs_folder_1
    │   └── sample1.fastq.gz
    └── nanopore_fastqs_folder_1
        ├── sample1.fastq.gz
        ├── sample2.fastq.gz
        └── sample3.fastq.gz
    ```
* command line
    ```
    metagenomix create \
        -i illumina_fastqs_folder_1 \
        -i illumina_fastqs_folder_2 \
        -j pacbio_fastqs_folder_1 \
        -k nanopore_fastqs_folder_1 \
    ```

### Sample names must match

Each fastq file name must start with the sample name **exactly** as it 
appear in the sample 
[metadata configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/metadata.md)
file. For Illumina, this is true whether the input files are paired-end or 
not (this is inferred as the sample names can be followed by `_[R]1.fastq[.gz]`
and `_[R]2.fastq[.gz]`).

Hence for the following metadata configuration file:
```
sample_name
sample1
sample2
```
would result in pipeline scripts being prepared to process this data:
```
illumina_fastqs_folder_1/sample1_1.fastq.gz
illumina_fastqs_folder_1/sample1_2.fastq.gz
illumina_fastqs_folder_1/sample2_1.fastq.gz
illumina_fastqs_folder_1/sample2_2.fastq.gz
illumina_fastqs_folder_2/sample1_1.fastq.gz
illumina_fastqs_folder_2/sample1_2.fastq.gz
illumina_fastqs_folder_2/sample2_1.fastq.gz
illumina_fastqs_folder_2/sample2_2.fastq.gz
pacbio_fastqs_folder_1/sample1.fastq.gz
nanopore_fastqs_folder_1/sample1.fastq.gz
nanopore_fastqs_folder_1/sample2.fastq.gz
```

### Sample merging

The sequences of sample for which there is more than on files located in 
multiple folder (e.g., re-sequencing or silent replicates) will be merged 
before running the first
[edit](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/preprocess.md#edit)
analysis so that all analysis outputs are performed on a per-sample basis. 
In the above example, this is the case for:

```
# sample1
illumina_fastqs_folder_1/sample1_1.fastq.gz
illumina_fastqs_folder_1/sample1_2.fastq.gz
illumina_fastqs_folder_2/sample1_1.fastq.gz
illumina_fastqs_folder_2/sample1_2.fastq.gz

# sample2
illumina_fastqs_folder_1/sample2_1.fastq.gz
illumina_fastqs_folder_1/sample2_2.fastq.gz
illumina_fastqs_folder_2/sample2_1.fastq.gz
illumina_fastqs_folder_2/sample2_2.fastq.gz
```
Note: in the case of Illumina it is possible that a sample appears with 
paired-end fastq files in one input folder and a single-end fastq file in 
another input folder. Then, the single-end files will be internally admitted 
as _another_ sample, associated with a copy of its metadata information. If a 
[paired-read merging](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#paired-read-merging)
or
[pooling](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#co-assembly-setup)
step is performed (whichever comes first), then the single-end reads will be 
merged into the file containing the non-extended reads resulting form
[merging analysis](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/alignment.md#paired-reads-merging).

## Output

The path to an output **folder** where all results will be written must be 
given using the `-o` (or `--output-dir`) argument.

### General structure

For every software, the output folder structure is always the same:
```
<current_sofware>
└── after_<previous_software(s)>
    ├── <technology/ies>
    #   # for "input unit" before pooling
    │   ├── <sample 1>
    │   ├── <sample 2>
    │   └── <sample n>
    #   # for "input unit" after pooling
    │   ├── <co-assembly pool 1>
    │   │   ├── <co-assembly pool group 1.1>
    │   │   ├── <co-assembly pool group 1.2>
    │   │   └── <co-assembly pool group 1.n>
    │   └── <co-assembly pool 2>
    │       ├── <co-assembly pool group 2.1>
    │       ├── <co-assembly pool group 2.2>
    │       └── <co-assembly pool group 2.n>
    ├── jobs
    │   ├── output
    │   ├── run_<input_unit_1 / chunk_1>.slm
    │   ├── run_<input_unit_2 / chunk_2>.slm
    │   └── run_<input_unit_n / chunk_n>.slm
    ├── provenance.txt
    └── run.sh
```

The created pipeline consists of one job script for each ([chunk](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#job-chunks)
of)
[input unit](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md),
which can be every sample or co-assembly pool group (depending on whether the 
analysis steps is before or after the 
[pooling](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#co-assembly-setup)
step), as explained
[here](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pooling.md).

### Software steps

A software step depends on the previous software step(s).

#### Limitation

Since a software can be set as output only **once** in the
[pipeline configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md#multiple-usage)
file, the user must edit this configuration to change the input to a given 
software. This will be reflected in the name of the first folder located 
inside a software output folder, which starts with `after_`. But, this only 
applies to the one previous software step. Indeed, the following two 
configurations will result in the same `alignment/after_trimming` output folder:

<table>
<tr>
<td>
Config without extra filtering
  <pre lang="courrier"> 
filtering
filtering          trimming
trimming           alignment
</pre>
</td>
<td>
Config with extra filtering
<pre lang="courrier">
filtering
filtering          extra_filtering
extra_filtering    trimming
trimming           alignment
</pre>
</td>
</tr>
</table>


#### Future solution

An internal mechanism to account that a software step is happening not 
only after one previous step but as a result of all previous steps will be 
established so that to separate outputs. This information is already 
available in the 
[provenance](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#provenance)
file (see below).


### Jobs

For each software step, the jobs that contain the actual command line 
scripts are located in the `jobs` folder. These scripts are executed by 
running the `run.sh` script located alongside this folder and displayed 
to the user as part of the
[terminal](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#terminal-output)
and
[file](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#run)
outputs of `metagenomix create`.


#### Job outputs

The `jobs` folder contains the `output` sub-folder where will be written the 
outputs of the Slurm/Torque jobs. These include: 
- stdout files (with extension `.o`): messages printed to the terminal 
  during the non-erroneous execution of the software.
- stderr files (with extension `.e`): error (or warnings) messages that 
  before software crash.

The `.o` and `.e` files can accumulate quickly and confuse users about the 
job status and notably to find useful stdout information. This is why the 
main `run.sh` scripts consist or three parts: 
```
# create to the jobs/output folder for the current software step
mkdir -p <output folder>/<software>/after_<previous software>/jobs/output
# move to this jobs/output folder before spawning the jobs
cd <output folder>/<software>/after_<previous software>/jobs/output
# spawn the jobs (here using Slurm)
sbatch <output folder>/<software>/after_<previous software>/jobs/run_<input_unit_1 / chunk_1>.slm
sbatch <output folder>/<software>/after_<previous software>/jobs/run_<input_unit_2 / chunk_2>.slm
sbatch <output folder>/<software>/after_<previous software>/jobs/run_<input_unit_n / chunk_n>.slm
```
As a result, all the jobs `.o` and `.e` output files will be written in this 
`jobs/output` location and nowhere else. This notably eases the
[monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
and
[management](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
for each step of the pipeline. The naming of these (slurm) job output files 
always has the following nomenclature:
`slurm-<software>.<project-name>.<input_unit_1 / chunk_1>_<JOBID>.e`


### Terminal output

The terminal will display information about the current pipeline, including 
whether the provided database paths point to existing and properly-formatted 
builds for the various aligners / formats that are requested (see docs on 
[databases]()
), as well as information about the jobs being written and their current 
analysis status. For more info on the progress of the software analyses, 
some command-line 
[arguments](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#monitoring-and-behaviour)
are available, whereas the
[monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
command will be much more useful.

After writing all the jobs for all softwares, the terminal display will 
first warn to `< PLEASE CONSIDER CHECKING THE COMMAND LINE SCRIPTS MANUALLY>` -
as it is worth checking the written commands before running them - and then 
provide the main `run.sh` bash script to run.

### Runs output

A folder named `_runs` will be created and filled with versionninf files 
containing the last terminal outputs as well as the dates at which the 
pipeline was created.  

### Running

In order to run the analyses for which all the necessary scripts have been 
generated, just run the `run.sh` script for the software you know needs to 
be run at this point. This can be done manually, or by copying and pasting 
the `sh /path/to/the/run.sh` commands displayed in the terminal.

For a software that is planned to run at the middle of the pipeline, it is 
obviously necessary that the previous softwares completed, and that their 
outputs are there to be used in input. Again, use
[monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
to make sure of that

#### Provenance

For each software step, a `provenance_HASH.txt` file is also located alongside 
the `jobs` folder. This file summarizes the pipeline analysis steps that led 
to the current results, including the names (and roles) of the previous 
softwares, as well as the full list of paramters used for each of these 
softwares.

The provenance filename contains a hash value in its names that is unique to 
all the softwares and their run parameters (except computing-resource params 
other than cpus/threads) that were used to result in the current outputs.

##### Example

For this short pipeline: 
```
cutadapt
cutadapt    bbmerge
```
The provenance file `bbmerge/after_cutadapt/provenance_HASH.txt` would contain: 

```
Pipeline steps to this output (and analysis type):
0	fastq   	: raw data
1	cutadapt	: preprocessing
2	bbmerge 	: paired read merging

Parameters for these steps' outputs:

 ===== 0: fastq (no param) =====

========== 1: cutadapt ==========
action: trim
...
zero_cap: false


========== 2: bbmerge ==========
branchlower: 3
...
ziplevel: 2
```

## Configuration files

Creating a pipeline relies on multiple
[configurations files](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configuration-files)
that all must/can be passed using the following commamd-line arguments: 

* Mandatory:
  * `-m` (or `--metadata`):
    [sample metadata](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/metadata.md)
    file with the names of the samples in the first column.
  * `-p` (or `--pipeline`): 
    [pipeline](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md)
    sequence of softwares to run one after another in order.
  * `-d` (or `--databases`): 
    [database](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/databases.md)
    names and their paths to be used by relevant pipeline softwares.
* Optional:
  * `-u` (or `--user-params`): 
    [user-defined parameters](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md)
    for the softwares and their computing resource needs.
  * `-M` (or `--modules`): 
    [modules](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/modules.md)
    to load (`module load <name>`) for the pipeline softwares. 
  * `-c` (or `--co-assembly`):
    [co-assembly](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/co-assembly.md)
    groups of samples based on metadata variable factors.
  * `-s` (or `--strains`):
    [strains](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/strains.md)
    foci for analyse of variation within species. 

## Computing and monitoring 

An established pipeline entails several softwares and each software is 
set to run with default (or user-edited)
[parameters](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md),
which can be either software-specific or related to computing resources 
and behaviour. 

### Computing

Computing-related parameters can be set for all softwares at once, directly at 
command line:
  * `-n` (or `--project-name`): `TEXT` that is **mandatory** as it will be 
    used in the name of the jobs and
    [job outputs](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#job-outputs).
  * `-x` (or `--chunks`): `INTEGER` for the maximum number of job
    [chunks](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#job-chunks)
    that are to be prepared per software, which depends on the number of
    [input unit](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md)
  

  * `--jobs/--no-jobs`: whether to prepare Slurm/Torque scripts.

If `--no-jobs` is used, then the main `run.sh` script for each software will 
contain the commands lines in `.sh` files instead of `.slm` (slurm) or `.pbs`
(torque) files. This is useful if softwares can run on your local machine. 
Thus, there will be not stdout/stderr outputs and the following arguments 
will have no effect:

  * `-a` (or `--account`): `TEXT` that is recognized by the 
    administration of some HPC environments (CPU hours accounting). 
  * `--torque/--no-torque`: this will change the format of the job scripts 
    so that they can be run on HPC using Torque as scheduler (instead of Slurm).
  * Arguments to define the [scratch location](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/scratch.md):
    * `-l,` `--localscratch`: `INTERER` for the amount of temporary storage 
      needed (in GB).
    * `--scratch/--no-scratch` 
    * `--userscratch/--no-userscratch` 
    * `--move-back/--no-move-back`: will prevent moving back the files from 
      the chosen scratch location (useful to keep files for manual work).

## Monitoring and behaviour

Although it is recommended to use the 
[monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
and
[management](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
commands of `metagenomix`, it is possible to gather information about all 
available parameters and about the status of every software in the 
configured pipeline, as well as to manage the content of the specific location 
used to cache
[Pfam-related](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/pfam_search.md)
analyses files:

* `--show-params / --no-show-params`: display all the parameters that 
  are available to use in the
  [parameters](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md)
  configuration file for the softwares of the pipeline.
* `--show-status / --no-show-status`: display a summary of the amnount of
  [input units](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md)
  that remain to be processed by `run.sh` as well as potential issues 
  related to missing outputs (incl., outputs from previous software steps).
* `--show-pfams / --no-show-pfams`: show which Pfam models where already 
  used and exist.
* `--purge-pfams / --no-purge-pfams`: clear out the Pfam models that where 
  already used.  
* `--verbose / --no-verbose`: display more information about data and jobs, 
  including the calls to [Xhpc](https://github.com/FranckLejzerowicz/Xhpc).
 
### Job chunks

The total number of jobs to run for software depends on the number of 
samples / co-assembly pool groups (read about
[input units](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md)), 
as well as - for relevant softwares - on the number of databases, models or
other parameters multiplying the amount of analyses that the pipeline is set 
to perform. 

By default and for each software, one job will be prepared for each
[input unit](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md)),
which for large projects can be hazardous, especially if each job requests 
large amounts of memory). Hence, it is possible for all softwares or for
specific ones to group commands lines into a given number of jobs ("chunks"):
* for all softwares: use `-x` (or `--chunks`) to provide an `INTEGER`
* for specific softwares: add `chunks: INTEGER` to the per-software sections 
  of the
  [parameters configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md) 
  file.

As a result, `metagenomix` will prepare a maximum of `INTEGER` jobs for the 
concerned softwares, which will be reflected in the names of the job files and 
stdour/stderr files. If there are less than `INTEGER` jobs to chunk, then no 
chunking will be operated and the names of the job files and stdour/stderr 
files will remain descriptive of the input units.

### Force re-writing

It is possible to ask the pipeline to re-write the command-lines for all 
softwares or specific ones, whether output was already been generated, but 
activating the **force** option:

* for all softwares: add `--force` to the `metagenomix create` command.
* for specific softwares: add `force: yes` to the per-software sections of the
  [parameters configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md) 
  file.

In this case, it might (currently) be preferable to rely on the pipeline 
[managing](https://github.com/FranckLejzerowicz/metagenomix/wiki/Running#Manage)
command and edit of remove outputs, so that there will be no already-existing
output and no need to **force** re-running/over-writing.

## Usage

```
Usage: metagenomix create [OPTIONS]

  Write jobs for your pipeline configuration.

Options:
  -i, --fastq-dir-illumina TEXT   Path to short Illumina reads fastq files
                                  folder(s)
  -j, --fastq-dir-pacbio TEXT     Path to long PacBio reads fastq files
                                  folder(s)
  -k, --fastq-dir-nanopore TEXT   Path to long MinION Nanopore reads fastq
                                  files folder
  -o, --output-dir TEXT           Path to pipeline output folder  [required]
  -m, --metadata TEXT             Path to the metadata file  [required]
  -p, --pipeline TEXT             Path to the file containing the softwares to
                                  run in order  [required]
  -d, --databases TEXT            Databases (yaml file)
  -u, --user-params TEXT          Parameters for the softwares of the pipeline
                                  (yaml file)
  -c, --co-assembly TEXT          Metadata column(s) defining the co-assembly
                                  groups (yaml file)
  -t, --strains TEXT              Species for strain-level analyses (yaml
                                  file)
  -n, --project-name TEXT         Name for your project  [required]
  -a, --account TEXT              User account for your HPC (in use for Slurm)
  -M, --modules TEXT              modules to use per software analyses (yaml
                                  file)
  -x, --chunks INTEGER            Number of jobs to split the commands into
                                  for each tool
  --force / --no-force            Force the re-writing of scripts for all
                                  commands(default is to not re-run if output
                                  file exists)  [default: no-force]
  --jobs / --no-jobs              Whether to prepare Torque jobs from scripts
                                  [default: jobs]
  --torque / --no-torque          Whether to prepare Torque jobs instead of
                                  Slurm  [default: no-torque]
  -l, --localscratch INTEGER      Use localscratch with the provided memory
                                  amount (in GB)
  --scratch / --no-scratch        Use the scratch folder to move files and
                                  compute  [default: no-scratch]
  --userscratch / --no-userscratch
                                  Use the userscratch folder to move files and
                                  compute  [default: no-userscratch]
  --move-back / --no-move-back    Do not move back from scrach (makes sense
                                  only for --userscratch)  [default: move-
                                  back]
  --show-params / --no-show-params
                                  Show all possible parameters for each
                                  software
  --show-status / --no-show-status
                                  Show status (needed inputs, done/to do
                                  outputs) for each software
  --show-pfams / --no-show-pfams  Show terms for which Pfam HMM models were
                                  already extracted before
  --purge-pfams / --no-purge-pfams
                                  Remove terms for Pfam HMM models that were
                                  already extracted before
  --verbose / --no-verbose        Whether to show input/outputs and other
                                  details  [default: no-verbose]
  --dev / --no-dev                For development...  [default: no-dev]
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```