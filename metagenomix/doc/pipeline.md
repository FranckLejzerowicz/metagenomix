:construction:

# Pipeline

The pipeline 
[configuration file](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configuration-files).
is mandatory as it establishes the sequence of softwares that are to be run 
one after the other.

## Format

The file must contain one or two [software names](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipline.md#software-names)
separated by a space (can be one or mulitple space(s), or a tab):
* if **only one software name** is given, then this software will use the raw 
  input 
  reads as input, those fastq files from the folders of the [input]() 
  options `-i` (or `--fastq-dir-illumina`), and/or `-j` (or 
  `--fastq-dir-pacbio`), and/or `-k` (or `--fastq-dir-nanopore`).


* if **two software names** are given, then the second space-separated 
  software will take as input the output of the first space-separated 
  software, for example:
  ```
  soft_1
  soft_1    soft_2
  soft_2    soft_3
  soft_2    soft_4
  ```
  will take the raw read and run this sequence of softwares: 
  ```
  soft_1 ──> soft_2 ──> soft_3
                  └───> soft_4
  ```
  Softwares can be used multiple times as inputs (as first space-separated 
  name in the pipeline configuration file), but currently, not as outputs. 
  So the following sequence of softwares will throw an error:
  ```
  soft_1    soft_3
  soft_2    soft_3
  ```
  The error will be explicit:
  ```
  Error in workflow "<path/to/pipeline/configuration/file.txt>":
          Can't run "soft_3" after "soft_2": "soft_3" already planned to use after "soft_1"
  Each tool must yield a single output re-used as input
  -> please re-run for each different workflow "graph"
  Exiting
  ```
  Read more on work-around to re-run `metagenomix` using different 
  configurations while
  [monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
  your progress and
  [managing](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
  your outputs.


* Lines starting with a hash comment character (`#`) will be ignored as well 
as blank lines. Hence, it can be useful to add comments to label your 
pipeline (for example, if sharing the configuration files with peers; see 
below).

### Content

Not all softwares can be run after any other:
  1. It is the user's responsibility to know what's going on and what to 
     analyse to answer a research question and obtain useful data and results.
  2. It may not make sense (e.g., running an assembler on the 
     output of a protein annotation software).
  3. A software can only be used once as output (solution in development)

To obtain a full list of the available softwares, please refer to the 
[softwares list](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares)
or look at the
[resources](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorial/code/software_files.md#resources)
file "`softwares.txt`".   

### Software names

Software names consist of one single word without underscore, except to call 
sub-commands for complex softwares, such as:
  * **metawrap**, which can either be called in the pipeline configuration 
    file using `metawrap` (which we call all the sub-commands), or using 
    `metawrap_binning`, `metawrap_refine`, etc to call specific sub-commands.
  * **checkm**, which can either be called using `checkm` (all sub-commands),
    or using `checkm_tree`, `checkm qa`, or other sub-commands. 
  * This is also the case for the custom **search_** softwares that 
    currently include `search_diamond` and `search_hmmer`. 
  
To learn more about the `metagenomix` mechanims of calling softwares, see 
the code [contributing](https://github.com/FranckLejzerowicz/metagenomix/wiki/Contributing) section.  

### Examples

Below are shown both a simple and an extensive example of valid pipeline 
configuration files. For more examples reflecting how previous research 
works could have been conducted using metagenomix, please explore examples 
of possible
[Research applications](https://github.com/FranckLejzerowicz/metagenomix/wiki/Research-applications).

#### Simple

```
# adapter trimmng (using raw reads as input) 
cutadapt
# kraken profiling (using cutadapt output as input)
cutadapt kraken2
# bracken profiling (using kraken2 output as input)
kraken2 bracken
```

#### Extensive

```
# ---------- start preprocessing ----------
# count raw reads
count

# adapter trimmng (using raw reads as input) 
hifiadapterfilt
cutadapt

# paired-reads merging (using cutadapt output as input) 
cutadapt bbmerge
# commented-out alternatives
# cutadapt pear
# cutadapt ngmerge
# ----------- end preprocessing -----------



# ---------- start profiling ----------
# align reads to bowtie2 databases (using bbmerge output as input)
bbmerge bowtie2
# profiling using woltka (using bowtie2 alignments as input)
bowtie2 woltka

# profiling using metaphlan (using bbmerge output as input)
bbmerge metaphlan

# kraken profiling (using bbmerge output as input)
bbmerge kraken2
# bracken profiling (using kraken2 output as input)
kraken2 bracken
# ------------ end profiling ------------



# ------------ start pooling ------------
# mandatory to start assembly workflow
bbmerge pooling (using bbmerge output as input)
# ------------- end pooling -------------



# ------------ start assembly workflow ------------
# assemble the merged fasta files that were pooled  
pooling spades (using co-assenbly pooled reads as input)
# other assemnblers can be used...
pooling canu
pooling flye

# check assembly quality  
spades quast (using spades assembly contigs as input)
# canu quast (quast can only be used once as output: change commenting & rerun)  
# flye quast (quast can only be used once as output: change commenting & rerun)

# bin contigs per sample
spades metawrap_binning
# refine binned contigs per sample
metawrap_binning metawrap_refine
# reassemble bins per sample
metawrap_refine metawrap_reassemble
# classify bins per sample
metawrap_refine metawrap_classify
# annotate genes of bins per sample
metawrap_refine metawrap_annotate
# ------------- end assembly workflow -------------



# ------------ start MAGs dereplication across samples ------------
metawrap_refine drep
# ------------- end MAGs dereplication across samples -------------



# ------------ start annotations ------------
spades prokka
spades barrnap
spades ccfind
spades viralverify
spades coconet
spades tiara
spades plasforest

# annotations can be done for each per-sample bins (e.g. plasmids) 
metawrap_refine plasmidfinder
metawrap_reassemble plasmidfinder

# or on the final MAGs
drep antismash

# predict proteins
spades prodigal
# drep prodigal (can be done on drep MAGs, but can only once as output)

# annotations of the predicted proteins
prodigal macsyfinder
prodigal integronfinder
prodigal deeparg_predict (<- modular tools, only calling a sub-command here)
prodigal search_diamond
prodigal search_hmmer

# ARGs annotations are possible on different inputs
spades abricate
drep metamarc
cutadapt karga
cutadapt kargva
cutadapt amrplusplus2

# ------------- end annotations -------------



# ------------ start contigs and genome quality checks ------------
# collecting sub-commands of a complex software (e.g., checkm)
spades checkm_tetra
drep checkm_coverage
metawrap_reassemble checkm_unbinned
metawrap_refine checkm_tree
# [note] sub-command may run in specific orders, here for "metawrap_refine"...
checkm_tree checkm_treeqa
checkm_tree checkm_lineageset
checkm_lineageset checkm_analyze
checkm_analyze checkm_qa
# [note] ...while running the software routine (all sub-command) avoids this
metawrap_refine checkm
```