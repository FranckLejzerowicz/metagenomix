# Managing

After 
[creating](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#creating)
a pipeline with
[configuration files](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configuration),
the user is provided with 
[job scripts](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#jobs)
to run the scheduled analyses, which write results in software outputs that 
are all organized according to a systematic 
[folder structure](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#general-structure).  

## Rationale

The main **advantage** is that a user can generate and run many jobs very
easily, which speeds-up the processing of shotgun metagenomics data.
But, one major **disadvantage** is that the user can rapidly be overwhelmed 
by the amount of outputs and their sizes. Indeed, very large files may 
accumulate in the file system for softwares such as 
[assembly](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#assembling), 
or
[alignment](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#alignment),
which can be a problem for users working with limited storage allocations 
(and notably on locations where compute nodes are mounted, which is the 
case for most HPC infrastructures). 

It is therefore necessary to **manage** the file system solicited for space, so
that it can host all input data, produce all output data and re-uses output 
data as input. HPC services often have _zero-storage policies_ on computing 
file systems, so that the final outputs or outputs that are not re-used as 
input for extended periods of time must be stored away on specific storage 
file system, but: 
* storage file systems are not necessarily accessible from jobs running on 
  compute nodes (trafficking files may not be possible as part of the jobs), 
* only a fraction of the project's input files may be relevant (i.e.,
  [input units](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md#inputsoutputs)
  for which no result was yet obtained).

Hence, it is **strongly recommended** that users very **regularly employ** 
the `metagenomix manage` command, which cannot be run as a job (must be run 
interactively).

## Inputs

Only one input argument is **mandatory**:
* `-i` / `--folder `: path to a 
[pipeline output folder](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#output)
(i.e., the `-o` argument to `metagenomix create`).

Yet, if no 
[management task]()
is specified, `metagenomix manage` will end without anything to manage. 

It is possible to restrict the management to the union of softwares that 
match what the user may provide as input to:
* `-p` / `--pipeline`: the softwares listed in the
  [pipeline configuration file](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md)
* `-s` / `--software`: regular expressions to match software name(s) (can be 
  used multiple times).

For example:
```
metagenomix manage \
  -i /path/on/cluster/for/some/output_folder \
  -k /path/on/storage/for/the_same/output_folder \
  -s checkm*
```

will manage any of `checkm`, `checkm_tetra`, `checkm_tree`, `checkm_treeqa`, 
`checkm_lineageset`, `checkm_analyze`, `checkm_qa`, `checkm_coverage`, 
`checkm_unbinned`, and `checkm2`.

### Management tasks

The `metagenomix manage` command allows performing three tasks, including 
[storage](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md#storage)
(`--store`),
[jobs](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md#jobs-cleansing)
output cleansing (`---jobs`), and file/folder 
[renaming](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md#renaming)
(`--rename`). It is also possible to
[remove](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md#removal)
jobs, files and folders during these various tasks (`--remove`). 

All management tasks will be performed interactively, i.e., by asking the 
user for input for each task and for each software, as follows:   

#### 1. Storage

If `--store` is activated (default is `--no-store`), this task will help the 
user collecting software outputs, in order to copy their files to a storage 
location, so that the original files from the computing location can be 
replaced by symlinks. These symlinks can later be followed by `metagenomix 
create` in order to create scripts necessary to fetch stored-away files 
before running new jobs (see section on
[stored inputs](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#stored-inputs)).

For this task to be possible, it is mandatory to provide the storage path:  
* `-o` / `--storage`: path to the storage location where the files and the 
  folder input folder structure will be copied (it is recommended to create a 
  folder name with the same name as the input folder), e.g. 

```
metagenomix manage \
  -i /path/on/cluster/for/some/output_folder \
  -o /path/on/storage/for/the_same/output_folder
```

##### Interaction

For each software and its previous software (see
[after](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md#limitation)
), the user will be interactively asked whether an entire folder should be 
stored away, e.g.:

```
  --------------------------
   Management task: Storage 
  --------------------------
  > 1 folder (./illumina=2.29 KB): Store? y/[n]: y
```

Some softwares may have several folders outputs (such as multiple 
co-assemblies for softwares running downstream of pooling in the pipeline, see
[here](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pooling.md)),
e.g.:

```
  --------------------------
   Management task: Storage 
  --------------------------
  > 3 folders (./pacbio=725 bytes; ./illumina=410 bytes; ./nanopore=600 bytes): Store all? y/n/[d]etails: 
```

In this case, the user can select which particular sub-folder should be 
stored or not:

```
  -------------------------
   Management task: Storage 
  --------------------------
  > 3 folders (./pacbio=725 bytes; ./illumina=410 bytes; ./nanopore=600 bytes): Store all? y/n/[d]etails: 
        * pacbio (725 bytes): Store? y/[n]: 
        * illumina (410 bytes): Store? y/[n]: 
        * nanopore (600 bytes): Store? y/[n]: 
```

##### Output

If the user enter `y` for at least one folder/sub-folder, a script (labelled 
with the data of creation) will be generated that must be executed:

```
==================================================
* Applying management decisions:
  - Storing
        -> please run the following script to spawn screen sessions
           sh <PATH>/output/_managed/DD-MM-YY_HHh/store.sh
```

### 2. Jobs cleansing


### 3. Renaming


### x. Removal


### Multipe screen sessions

If many/large files are to be stored away, the transfer can take a lot of 
time. Hence, it is possible to evenly split the total amount of transfers into 
multiple screen sessions, that will be detached to run in parallel:
* `-x` / `--chunks`: number of screen sessions to spawn

 


## Outputs

All management tasks consist of file system operations, that can be executed 
by running output
[screen]()
scripts. These scripts are written in a `_managed/DD-MM-YY_HHh` folder 
placed directly inside the main pipeline output folder (given using `-i`).

## Usage

```
Usage: metagenomix manage [OPTIONS]

  Deal with the contents of your pipeline output folder.

Options:
  -i, --folder TEXT         Path to pipeline output folder (`-o` for "create"
                            module)  [required]
  -p, --pipeline TEXT       Path to the file containing the softwares to run
                            in order
  -s, --software TEXT       Software(s) to manage (or all in `-i/-p`)
  -k, --disk TEXT           Path to a storage disk for `--store` task
  -x, --chunks INTEGER      Number of scripts for each `--store` task
  --jobs / --no-jobs        [Task] Enable job output management  [default: no-
                            jobs]
  --remove / --no-remove    [Task] Enable output removal  [default: no-remove]
  --rename / --no-rename    [Task] Enable output renaming  [default: no-
                            rename]
  --store / --no-store      [Task] Enable output storage  [default: no-store]
  --confirm / --no-confirm  Whether to ask for confirmation before applying
                            task  [default: confirm]
  --version                 Show the version and exit.
  --help                    Show this message and exit.
```