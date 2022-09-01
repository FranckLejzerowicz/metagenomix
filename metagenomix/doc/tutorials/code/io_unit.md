:construction:

# Inputs/Outputs

### Purpose

This page explains the types of input and output data structures that are 
being used for every software function, and how these are used and collected 
by various utility functions allowing the
[creation](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/creating.md),
[monitoring](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/monitoring.md)
and 
[management](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/managing.md)
of the pipeline jobs and outputs.

### Data structures

The code for every software function starts by iterating over the different 
samples or co-assembly pools (for non-
[holistic]()
softwares), using a `for` loop iterating over dictionary items consisting of 
the list of input files or folder (dict values) either per 
`(technology, sample_name)` or `(technology/ies, co-assembly group)` (dict 
keys).

The reason why the key of the input dictionary can be either per sample or 
per co-assembly pool group is because the essential input/output data 
structure in use by the pipeline will change depending on whether the 
software currently being processed is located upstream of downstream of the 
co-assembly [pooling]()
step.

### Scratch relocation

The
[io_update()](https://github.com/FranckLejzerowicz/metagenomix/blob/25a58495d2de7bdc282f22396a0d66db053f0f88/metagenomix/_io_utils.py#L393-L418)
function collects the input and output files and folders for each software 
class instance, so that when the jobs are written for this software, these 
files and folders can be moved to and from the
[scratch location](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/scratch.md)
if this is requested by the user. 

By default, jobs are set to compute from the file system location where they 
are queued (yet the main
[job scripting mechanism]()
first move to the job output folder) which is not a scratch location (slower 
computes). Hence, to leverage compute speed on scratch file systems, the 
user must request a specific [scratch location](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/scratch.md).
