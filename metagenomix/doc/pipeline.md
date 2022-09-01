# Pipeline

The pipeline 
[configuration file](https://github.com/FranckLejzerowicz/metagenomix/wiki/Configuration-files).
is mandatory as it establishes the sequence of softwares that are to be run 
one after the other.

### Format



### Content

Not all softwares can be run after any other:
  1. It is the user's responsibility to know what's going on and what to do 
     to answer research questions, and to get useful data and results.
  2. Many sequences of software in the pipeline do not make sense, such as 
     running an assembler software on the output of a protein annotation 
     software. Refer to the literature and to the list of software grouped 
     per category in the list of available [softwares](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares).

### Software names

Software names consist of one single word without underscore, except to call 
sub-commands for complex softwares, such as:
  * **metawrap**, which can either be called in the pipeline configuration 
    file using `metawrap` (which we call all the sub-commands), or using 
    `metawrap_binning`, `metawrap_refine`, etc to call specific sub-commands.
  * **checkm**, which can either be called using `checkm` (all sub-commands),
    or using `checkm_tree`, `checkm qa`, or other sub-commands. 

This is also the case for `search_diamond` and `search_hmmer`, which are 
both **search** functions (although)
