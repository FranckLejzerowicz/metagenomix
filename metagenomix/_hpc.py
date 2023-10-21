# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import isfile
from metagenomix._io_utils import get_first_line


def validate_email(email: str) -> bool:
    """A not-so-fancy that a string may be a valid email address.

    Parameters
    ----------
    email : str
        A string that maz look like an email address

    Returns
    -------
    bool
        True if the string looks like an email, else False
    """
    if email.count('@') != 1:
        return False
    user, domain = email.split('@')
    if not len(user):
        return False
    if '.' not in domain:
        return False
    domain_split = domain.split('.')
    if 0 in map(len, domain_split):
        return False
    return True


def get_email() -> str:
    """Collect the email address interactively from the user.

    Returns
    -------
    email : str
        email address of the user.
    """
    email = input('Please enter a valid email address: ')
    while not validate_email(email):
        email = input('Please enter a valid email address: ')
    return email


def create_config(log, email_fp: str) -> str:
    """Collect the email address interactively from the user
    and write it somewhere it can be reused.

    Parameters
    ----------
    email_fp : str
        Path to the email config file

    Returns
    -------
    email : str
        email address of the user.
    """
    email = get_email()
    write_email(log, email_fp, email)
    return email


def write_email(log, email_fp: str, email: str):
    """Write the email address somewhere it can be reused.

    Parameters
    ----------
    email_fp : str
        Path to the email config file
    email : str
        email address of the user.
    """
    with open(email_fp, 'w') as o:
        o.write('%s\n' % email)
    log.info('Written: %s' % email_fp)


def edit_config(self, log, email_fp: str):
    """
    Checks that the content of the existing config file
    is a valid email address.

    Parameters
    ----------
    self
        All config attributes
    log
        Logger
    email_fp : str
        Path to the file that may (or not) contains an email address

    Returns
    -------
    email : str
        email address of the user.
    """
    # parse the first line of the config file
    email = get_first_line(email_fp)
    if self.config_email:
        log.info("Registered email:", email)
        email = create_config(log, email_fp)
    return email


def get_scratch_area(scratch: str, default: str) -> str:
    """Collect the scratch/userscratch folder interactively from the user.

    Parameters
    ----------
    scratch : str
        "scratch" or "userscratch"
    default : str
        Default scratch folder on SAGA

    Returns
    -------
    scratch_folder : str
        Scratch folder
    """
    ret = input('(default [on SAGA]: "%s"): ' % default)
    if ret == 'y':
        scratch_folder = default
    else:
        if not isfile(ret):
            raise IOError('Enter a valid "%s" folder path...' % scratch)
        scratch_folder = ret
    return scratch_folder


def get_scratch(scratch: str) -> str:
    """Collect the scratch folders interactively from the user.

    Parameters
    ----------
    scratch : str
        "scratch" or "userscratch"

    Returns
    -------
    scratch_folder : str
        Scratch folder
    """
    print('Please enter "%s" location ("y" to confirm default)' % scratch)
    if scratch == 'userscratch':
        scratch_folder = get_scratch_area(scratch, '/cluster/work/users/$USER')
    else:
        scratch_folder = get_scratch_area(scratch, '/cluster/work/jobs')
    return scratch_folder


def write_scratches(scratch_fp: str, scratch_folder: str) -> None:
    """Write the scratch folder somewhere it can be reused.

    Parameters
    ----------
    scratch_fp : str
        Path to the config file
    scratch_folder : str
        Scratch folder
    """
    with open(scratch_fp, 'w') as o:
        o.write('%s\n' % scratch_folder)
    print('Written: %s' % scratch_fp)
    if scratch_folder == '${USERWORK}':
        print('Written automatically since your machine presets $USERNAME')


def create_scratch(self, scratch_fp: str, scratch: str) -> str:
    """Collect the scratch folders interactively from the user
    and write it somewhere it can be reused.

    Parameters
    ----------
    self
        config attributes
            config_scratch : bool
                Show current scratches folder and/or edit it
    scratch_fp : str
        Path to the file containing the scratch folder
    scratch : str
        "scratch" or "userscratch"

    Returns
    -------
    scratch_folder : str
        Scratch folder
    """
    if self.config_scratch:
        scratch_folder = get_scratch(scratch)
    elif scratch == 'userscratch' and 'USERWORK' in os.environ:
        scratch_folder = '${USERWORK}'
    else:
        scratch_folder = get_scratch(scratch)
    write_scratches(scratch_fp, scratch_folder)
    return scratch_folder


def edit_scratch(self, log, scratch_fp: str, scratch: str) -> str:
    """
    Checks that the content of the existing scratches file are valid paths.

    Parameters
    ----------
    self
        config attributes
            config_scratch : bool
                Show current scratches folder and/or edit it
    scratch_fp : str
        Path to the file containing the scratch folder
    scratch : str
        "scratch" or "userscratch"

    Returns
    -------
    scratch_folder : str
        Scratch folder
    """
    # parse the first line of the config file
    scratch_folder = get_first_line(scratch_fp)
    if self.config_scratch:
        log.info("Registered %s: %s" % (scratch, scratch_folder))
        scratch_folder = create_scratch(self, scratch_fp, scratch)
    return scratch_folder


def set_account(self) -> None:
    """Get the directives for the user account.

    Parameters
    ----------
    self
        All config attributes
            account : str
                Name of the account
            torque: bool
                Adapt to Torque
    """
    account = ''
    if self.account:
        if self.torque:
            account = '#PBS -A %s' % self.account
        else:
            account = '#SBATCH --account=%s' % self.account
    self.directives['account'] = account


def set_partition(self) -> None:
    """Get the directives for the allocation of a partition based on user input.

    Parameters
    ----------
    self
        All config attributes
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        self.directives['partition'] = '#PBS -q '
    else:
        self.directives['partition'] = '#SBATCH --partition='


def set_environment(self) -> None:
    """Get the directives for the user environment.

    Parameters
    ----------
    self
        All config attributes
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        self.directives['environment'] = '#PBS -V'
    else:
        self.directives['environment'] = '#SBATCH --export=ALL'


def set_job(self) -> None:
    """Get the directives for the user environment.

    Parameters
    ----------
    self
        All config attributes
            job: str
                Job name
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        self.directives['job'] = '#PBS -N '
    else:
        self.directives['job'] = '#SBATCH --job-name='


def set_email(self) -> None:
    """Get the directives for the emailing.

    Parameters
    ----------
    self
        All config attributes
            email: bool
                Send email (always if fail)
            torque: bool
                Adapt to Torque
            email_address : str
                email address of the user
    """
    if self.torque:
        if self.email:
            email = '#PBS -m ae'
        else:
            email = '#PBS -m a'
        email += '\n#PBS -M %s' % self.email
    else:
        if self.email:
            email = '#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80'
        else:
            email = '#SBATCH --mail-type=FAIL,TIME_LIMIT_80'
        email += '\n#SBATCH --mail-user=%s' % self.email
    self.directives['email'] = email


def set_stdout_stderr(self) -> None:
    """Get the directives for the working directory where the job stdout and
    stderr files will be written. The path (no extension) of the stdout and
    stderr files will also be collected.

    Parameters
    ----------
    self
        All config attributes
            workdir: str
                Working directory
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        work_dir = '${PBS_O_WORKDIR}/${PBS_JOBNAME}'
        job_id = 'PBS_JOBID'
        std_path = 'localhost:%s_${%s}' % (work_dir, job_id)
        oe = '#PBS -o %s.o\n#PBS -e %s.e' % (std_path, std_path)
    else:
        job_id = 'SLURM_JOB_ID'
        std_path_ = 'slurm-%sx_%sj' % ('%', '%')
        std_path = 'slurm-${SLURM_JOB_NAME}_${SLURM_JOB_ID}'
        oe = '#SBATCH -o %s.o\n#SBATCH -e %s.e' % (std_path_, std_path_)
    self.directives['std_path'] = std_path
    self.directives['job_id'] = job_id
    self.directives['oe'] = oe


def set_time(self) -> None:
    """Get the directives for the time.

    Parameters
    ----------
    self
        All config attributes
            time: str
                Wall time limit
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        self.directives['time'] = '#PBS -l walltime='
    else:
        self.directives['time'] = '#SBATCH --time='


def set_memory(self) -> None:
    """Get the directives for the memory.

    Parameters
    ----------
    self
        All config attributes
            torque: bool
                Adapt to Torque
    """
    if self.torque:
        self.directives['mem'] = '#PBS -l mem='
    else:
        self.directives['mem'] = '#SBATCH --mem='


def set_localscratch(self) -> None:
    """Get the directives for the localscratch --gres option.
    Only available for the SAGA cluster and on Slurm.

    Parameters
    ----------
    self
        All config attributes
            localscratch: tuple
                Use localscratch with the provided memory amount (in gb)
            torque: bool
                Adapt to Torque
    """
    self.directives['localscratch'] = None
    if not self.torque and self.localscratch:
        localscratch = '#SBATCH --gres=localscratch:'
        self.directives['localscratch'] = localscratch
        # args['localscratch'] = '/localscratch'
