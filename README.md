:construction:

# metagenomix

## Description

metagenomix is a pipeline
[creator](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md),
[monitor](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
and 
[manager](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
that takes care of writing the command-lines for any shotgun metagenomics 
software, either as bash scripts or Slurm/Torque jobs (incl. scratch space 
usage), following  user-defined configurations for databases, co-assemblies, 
strain foci, as well as software-specific or computing resource parameters 
(incl. memory, scratch relocations, modules and conda environments).

Any software? Well, if not already in the
[softwares list](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares),
someone will need to add it to metagenomix, following the instructions to 
[contribute code](https://github.com/FranckLejzerowicz/metagenomix/wiki/Contributing).


This also eases importing issues and notably run this on a High-Performance 
Computer (HPC), either [Torque](http://docs.adaptivecomputing.com/torque/4-0-2/help.htm) 
or [Slurm](https://slurm.schedmd.com/documentation.html) scheduler, and it 
heavily relies on a you having installed and configured properly multiple conda
environments (explanations and environment.yml files provided in Wiki). 
This tool wrapping other tools could be done using 
[snakemake](https://snakemake.readthedocs.io/en/stable/) but this is meant to 
only create the scripts to run on your samples, not execute these scripts.

## Installation
```
pip install --upgrade git+https://github.com/FranckLejzerowicz/metagenomix.git
```

*_Note that python and pip should be python3_

### Depencencies

You must install all the databases and softwares (lists below) that the pipeline
will allow you to run. I recommend using conda environments for each tool and 
then specify the name of the environment for this tool in the run parameters
configuration file. I am working on an installer of conda environments to 
alleviate this time-consuming step, but since some tools may require you to edit
configurations for you system, this will not be fully integrated.  

#### Databases

* [PFAM](http://pfam.xfam.org/): allows searching HMM profiles per keywords from the Pfam-A catalogue,
(see [here](https://doi.org/10.1093/nar/gkp985) or 
[here](https://academic.oup.com/nar/article/38/suppl_1/D211/3112325?searchresult=1#64944278)) 
which you need to download: Pfam-A files "Pfam-A.hmm" and "Pfam-A.hmm.dat.gz"
[here](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases)
* [dbCAN](https://bcb.unl.edu/dbCAN/): allows searching HMM profiles per CAZy subfamily.
* [MAR](https://mmp2.sfb.uit.no/databases/): Allow querying the MarDB and MarRef databases for marine microbes. 
* [GTDB](https://gtdb.ecogenomic.org/)

## Input

IN CONSTRUCTION

## Outputs

IN CONSTRUCTION
