# Scratch locations

## Usage in `metagenomix`

Setting up scratch location for job to run faster there can be specified 
using the `metagenomix create` command line (applies to all softwares of the 
[pipeline configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/pipeline.md))
or for specific softwares in the 
[parameters configuration](https://github.com/FranckLejzerowicz/metagenomix/blob/main/metagenomix/doc/parameters.md)
file (takes precedence for this software).

If argument `-l` (or `--localscratch`) is used, it must be followed by an 
integer, indicating the amount of storage (in GB) that will be allocated to 
the job (e.g. `-l 100`). This will be translated in the HPC job script into 
the `--gres=localscratch:100G` directive (this is a slurm example).

When specifying the scratch location globally via the `metagenomix create` 
command, the user can provide either or all of these arguments:
* `-l <disk-size-in-GB>` (or `--localscratch <disk-size-in-GB>`)
* `--scratch/--no-scratch`
* `--userscratch/--no-userscratch`

If more that one of these scratch-related arguments is provided, then only 
one will be retained and used for all softwares, in the following order to 
priority:

High priority :arrow_up: localscratch > scratch > userscratch Low priority 
:arrow_down:

## Dependence on `Xphc`

### Overview

For `metagenomix` to create jobs, it relies on the companion software
[Xhpc](https://github.com/FranckLejzerowicz/Xhpc), which is automatically 
installed alongside and allows the preparation of HPC scripts based on basic 
bash scripts. For Xhpc to work, it is necessary that the user provide edit 
a `config.txt` file that is generated by this software (and for each 
`$USERNAME` independently) upon the first utilization.

### User input

#### Scratch locations

Note that if `--scratch` or `--userscratch` is activated for the first time, 
Xphc will query the ask to enter the absolute paths to these locations, or 
to stick with the default environment variables (if these exist). The user 
can choose to enter a path or an environment variable to replace this default. 
For example, on the NRIS's machine "SAGA", there are two scratch locations 
that already exist, which absolute paths are saved in the environment 
variables `SCRATCH` and `USERWORK` (see
[sigma2 documentation](https://documentation.sigma2.no/files_storage/clusters.html)).

#### Email

In addition, Xhpc will simply ask for the user's email address, to know 
where to redirect the HPC scheduler for send email at job completion/error (see
[here](https://github.com/FranckLejzerowicz/Xhpc#requisite)).

#### Changing

Upon subsequent usages, the user will not be asked again and `metagenomix` 
will just use the `Xhpc` user configs as defaults. Therefore, to change the 
email or scratch locations, one can run Xhpc with argument `--config-scratch`,
or directly edit the `config.txt` where Xhpc was installed.
