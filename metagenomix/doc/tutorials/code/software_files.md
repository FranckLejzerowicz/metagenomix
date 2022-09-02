# Code files

The code of the `metagenmix` project is organized into several folders.

### core

This
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/core)
contains the code that handle the generic pipeline mechanisms, 
which in logical order of execution include:

* `config.py`: Reading [command-line](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md) configuration ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/config.md)) 
* `databases.py`: Checking available databases and database builds for 
  specific sequence aligners ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/databases.md))
* `pipeline.py`: Reading pipeline configuration file, defining the 
  per-software and pipeline graph classes ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/pipeline.md))
* `parameters.py`: Checking the user-defined software parameters and 
  defining the defaults ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/parameters.md))
* `commands.py`: Collecting the command lines and data structures ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/commands.md))
* `jobs.py`: Writing pipeline command line scripts and other outputs ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/jobs.md))
* `output.py`: class collecting the current data and job outputs for a given 
  software ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/output.md))

### softwares

This
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/softwares)
contains the per-software code files:

* `alignment.py`: sequence alignment vs. databases and for merging ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/alignment.md))
* `annotation.py`: annotation of DNA and proteins sequence, contigs and 
  genomes ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/annotation.md))
* `args.py`: antibiotic resistance genes detection ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/args.md))
* `assembly.py`: DNA and protein-level sequence (hybrid) assembly ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/assembly.md))
* `binning.py`: binning of assembled contigs into per-sample MAGs ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/binning.md))
* `genomics.py`: delineation on MAGs across samples and genome quality check 
  and taxonomic annotation ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/genomics.md))
* `metamarker.py`: cross-sample detection of marker features (this tool is a 
  bit specific as it includes several commands) ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/metamarker.md))
* `phlans.py`: profilers from the BioBakery group - may be augmented with 
  more of the *PhlAn softwares ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/phlans.md))
* `plasmids.py`: detection of plasmids from reads and contigs ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/plasmids.md))
* `pooling.py`: preparation of the sample group for each of the 
  user-defined co-assembly pools ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pooling.md))
* `preprocess.py`: sequence filtering, adapter-trimming and quality checking 
  ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/preprocess.md))
* `profiling.py`: creation of taxonomic and functuional feature count tables 
  based on the alignment of reads to various databases ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/profiling.md))
* `simka.py`: across-samples ecological distance mesures based on reads k-mer
  ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/simka.md))
* `strains.py`: strain-level analyses including profiling (feature tables) 
  and more specific measures of inter/intra genome variation ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/strains.md))
* `viruses.py`: detection of virures from reads and contigs ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/viruses.md))

### resources

This
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/resources)
contains small files that ship with the `metagenomix` package 
and that could possibly be edited (not recommended):

* `scripts` (**DO NOT EDIT**): utilitarian python scripts called for 
  specific softwares (e.g., input formatting). 
* `pfam_terms`: diamond- and hmm- formatted files collected within 
  `metagenomix` for specific protein targets to be searched using the 
  softwares `search_diamond` and `search_hmmer`, respectively. This folder 
  can be managed either using `metagenomix create` and options 
  `--show-pfams/--no-show-pfams` and/or `--purge-pfams/--no-purge-pfams`, 
  or using `metagenomix manage` (in progress :construction:)
  ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/softwares/pfam_search.md))
* `wol`: key files from the [Web Of Life](https://biocore.github.io/wol/data)
  database and the [Woltka](https://github.com/qiyunzhu/woltka) classifier: 
  * `batch_down.sh`: [script](https://biocore.github.io/wol/data/genomes/) to setup the database locally
  * `make_down_list.py`: [script](https://biocore.github.io/wol/data/genomes/) to setup the database locally
  * `metadata.tsv`: Web of Life [metadata](https://biocore.github.io/wol/data/genomes/)
  * `tree.nwk`: Web of Life [tree](https://biocore.github.io/wol/data/trees/)
  * `genome_sizes.txt`: Web of Life [genome sizes](https://biocore.github.io/wol/data/)
  * `lineages.txt`: Web of Life [taxonomy](https://biocore.github.io/wol/data/)
* `run_params.yml`: default software parameters
* `softwares.txt`: list of available software names and their role  

### scripts

This
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/scripts)
contains the command line entry-point code, where all the command-line 
arguments are defined, including:

* `_standalone_metagenomix.py`: main standalone entry-point for `metagenomix`    
* `modules/`:
  * `create.py`: Write jobs for your pipeline configuration ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md))
  * `manage.py`: Edit the contents of your pipeline output folder ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md))
  * `monitor.py`: Check IO/job status of your pipeline configuration ([details](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md))

### tests

This
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/tests)
contains the unittests for all functions in the entire code (in progress...)

### doc

This 
[folder](https://github.com/FranckLejzerowicz/metagenomix/tree/main/metagenomix/doc)
documentation (all linked and accesible from the
[wiki](https://github.com/FranckLejzerowicz/metagenomix/wiki)).