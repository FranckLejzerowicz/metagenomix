# metagenomix

`metagenomix` is a pipeline
[creator](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md),
[monitor](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
and 
[manager](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
that takes care of writing the command-lines for any shotgun metagenomics 
software, either as bash scripts or
[Slurm](https://slurm.schedmd.com/documentation.html)
/
[Torque](http://docs.adaptivecomputing.com/torque/4-0-2/help.htm)
jobs (incl. scratch space usage), based on user-defined
[configurations](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configurations)
for databases, co-assemblies, strain foci, as well as software-specific or 
computing resource parameters (incl. memory, scratch relocations, modules 
and conda environments). 

Outputs are scripts that the user needs to run sequentially, which typically 
can be handled by packages such as
[snakemake](https://snakemake.readthedocs.io/en/stable/): this is not (yet) 
used here as `metagenomix` is only meant to facilitate the creation, 
monitoring, management and access of shotgun metagenomic analyses results for 
personalized pipelines including any software.

Any software? Well, if not already in the
[softwares list](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares),
someone will need to add it to metagenomix, following the instructions to 
[contribute code](https://github.com/FranckLejzerowicz/metagenomix/wiki/Contributing).

Please read the documentation using the [Wiki pages](https://github.com/FranckLejzerowicz/metagenomix/wiki)

## Installation

```
pip install metagenomix
```

or

```
pip install --upgrade git+https://github.com/FranckLejzerowicz/metagenomix.git
```

### Depencencies

While a container solution with all the softwares and conda environments is 
being develop, it currently is the responsibility of the user to install all 
the databases and softwares that the pipeline will allow you to prepare 
command-lines for. Some softwares necessitate to be present either as a 
single binary file or with an entire folder (e.g., pre-trained models or 
scripts). Since some softwares require the user to edit configurations 
after install, some level of manual installation/tuning will be needed 
before usage.

## Usage

```
Usage: metagenomix [OPTIONS] COMMAND [ARGS]...

  Metagenomix command line manager

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  create   Write jobs for your pipeline configuration.
  export   Prepare an archive for specific pipeline outputs.
  manage   Edit the contents of your pipeline output folder.
  monitor  Check IO/job status of your pipeline configuration.
  merge    Combine the per-sample outputs into feature tables.
```

Detailed explanations for each command at
[Running](https://github.com/FranckLejzerowicz/metagenomix/wiki/Running)
Wiki page

