# Pooling for co-assembly

## Overview

The pooling step is an important and **pivotal** step of the pipeline, since 
at this step occurs a major change in the way that `metagenomix` handles
[inputs/outputs](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/io_unit.md)
to collect and write the command-lines/jobs.

Essentially, the nature of the input to a software will depend on whether a 
software is set to run **upstream** or **downstream** of the
[pooling](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#pooling)
step according to the
[pipeline graph](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/graph.md)
resulting from a given
[pipeline configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md):
  * **upstream**: processing per sample (typical for
[profiling](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#profiling)
softwares)
  * **downstream**: processing per co-assembly pool group (typical on
[assembled](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#assembling)
data)

## Configuration

By default, there will be no co-assembly if the user does not provide a 
[co-assembly configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/co-assembly.md)
file. In fact (i.e., as per the
[core code](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/tutorials/code/core/commands.md))
there _will_ be a **per-sample** co-assembly, where each sample is co-assemble 
with itself only.

## Pipeline usage 

It is mandatory that the pipeline configuration file includes
[pooling](https://github.com/FranckLejzerowicz/metagenomix/wiki/Softwares#pooling)
as output software name and that it is used as input to an assembler. For 
example, these pipeline configurations are valid.

Is it perfectly fine to configure a pipeline that do not involve a pooling 
step, which means that this pipeline will not involve assembly-based 
analysis (as in the
[simple pipeline example](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md#simple)).