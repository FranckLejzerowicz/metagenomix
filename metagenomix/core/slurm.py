# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def set_tmpdir(self) -> list:
    """Define the path to a temporary folder based on the already existing
    TMPDIR folder.

    Parameters
    ----------
    self
        All arguments. Here only the following keys are of interest:
            torque: bool
                Adapt to Torque
    """
    # define the temporary folder
    d = self.config.directives
    tmpdir = '${TMPDIR}/%s_${%s}' % (self.job_name, d['job_id'])
    # set command to create the temporary folder
    tmp = [
        '# create and export the temporary directory',
        'mkdir -p %s' % tmpdir,
        'export TMPDIR="%s"' % tmpdir,
        'echo Temporary directory is ${TMPDIR}']
    return tmp


def set_scratching(self, params: dict) -> list:
    """Collect all the commands to set the scratch folder, for folder
    creation and deletion.

    Parameters
    ----------
    self
        All config attributes
            scratchs : str or int
                Path to scratch or userscratch, or localscratch int allocate
            job_id : str
                path to the stderr/stdout
    params : dict
        All arguments. Here only the following keys are of interest:
            scratch : str
                scratch space to use

    Returns
    -------
    scra : list
        Preamble commands
    """
    scra = []
    if params['scratch']:
        scratch = self.config.scratchs.get(params['scratch'], '/localscratch')
        # Define scratch directory
        scratch += '/${%s}' % self.config.directives['job_id']
        # Get commands to create and more to scratch directory
        scra.extend([
            '# Define and create a scratch directory',
            'SCRATCH_FOLDER="%s"' % scratch, 'export SCRATCH_FOLDER',
            'mkdir -p ${SCRATCH_FOLDER}', 'cd ${SCRATCH_FOLDER}',
            'echo Working directory is ${SCRATCH_FOLDER}'])
    return scra


def set_preamble(self, params: dict, job_script: str) -> list:
    """Collect the command for the job's preamble.

    Parameters
    ----------
    self
        All config attributes
            torque: bool
                Adapt to Torque
            std_path : str
                path to the stderr/stdout
    params : dict
        All arguments. Here only the following keys are of interest:
            env: str
                Conda environment
            machine: str
                Machine name, for specific conda usage
    job_script : str
        Path to the input job

    Returns
    -------
    prea : list
        Preamble commands
    """
    prea = ['# general environment info / behaviour', 'uname -a', 'set -e']
    if params['env']:
        env = ['\n# active conda environment',
               'echo "Conda environment is %s"' % params['env']]
        if params['machine'] == "saga":
            env.extend([
                'module load Anaconda3/2022.10',
                'export PS1=\\$',
                'source ${EBROOTANACONDA3}/etc/profile.d/conda.sh',
                'conda deactivate &>/dev/null'])
        env.append('conda activate %s' % params['env'])
        if params['machine'] == "saga":
            env.append('module purge')
        prea.extend(env)
    prea.extend(['\n# echo some info about the job',
                 'echo Running on host `hostname`',
                 'echo Time is `date`',
                 'echo Directory is `pwd`'])
    stdout = 'echo Job stdout is %s' % self.config.directives['std_path'],
    stderr = 'echo Job stderr is %s' % self.config.directives['std_path']
    prea.extend(['%s.o' % stdout, '%s.e' % stderr])
    prea.append('echo Job script: %s' % job_script)
    return prea


def get_nodes_cpus(self, params: dict) -> str:
    """Distribute the number of processors requested among the requested nodes.

    Parameters
    ----------
    self
        All config attributes
            torque: bool
                Adapt to Torque
    params : dict
        All arguments. Here only the following keys are of interest:
            nodes: int
                Number of nodes
            cpus: int
                Number of CPUs

    Returns
    -------
    directive : str
        Current directive
    """
    nnodes = int(params['nodes'])
    ncpus = int(params['cpus'])
    if self.config.torque:
        directive = '#PBS -l nodes=%s:ppn=%s' % (nnodes, ncpus)
    else:
        if nnodes > 1 and ncpus > 1:
            directive = '#SBATCH --nodes=%s' % nnodes
            directive += '\n#SBATCH --ntasks-per-node=1'
            directive += '\n#SBATCH --cpus-per-task=%s' % ncpus
        elif nnodes > 1:
            directive = '#SBATCH --nodes=%s' % nnodes
            directive += '\n#SBATCH --ntasks-per-node=1'
        elif ncpus > 1:
            directive = '#SBATCH --ntasks=1'
            directive += '\n#SBATCH --cpus-per-task=%s' % ncpus
        else:
            directive = '#SBATCH --ntasks=1'
            directive += '\n#SBATCH --cpus-per-task=1'
    return directive


def set_directives(self, params) -> list:
    d = self.config.directives
    dirs = []
    dirs.extend([d['shebang'], d['environment']])
    if d['account']:
        dirs.append(d['account'])
    if params['partition']:
        dirs.append('%s%s' % (d['partition'], params['partition']))
    else:
        dirs.append('%snormal' % d['partition'])
    if params['array_jobs']:
        dirs.append('#SBATCH --array=1-%s:' % d['array_jobs'])
    dirs.append('%s%s' % (d['job'], self.job_name))
    localscratch = '#SBATCH --gres=localscratch:'
    if isinstance(params['scratch'], int):
        dirs.append('%s%sG' % (localscratch, params['scratch']))
    elif d['localscratch']:
        dirs.append('%s%sG' % (localscratch, d['localscratch']))
    dirs.append(d['email'])
    dirs.append(d['oe'])
    dirs.append('%s%s:00:00' % (d['time'], params['time']))
    dirs.append(''.join([d['mem'], str(params['mem']), params['mem_dim'].upper()]))
    dirs.append(get_nodes_cpus(self, params))
    return dirs
