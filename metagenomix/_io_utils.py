# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import sys
import yaml
import glob
import pandas as pd
from os.path import dirname, isdir, isfile


def read_yaml(
        file_path: str
) -> dict:
    """Simply reads a yaml and return its contents in a dictionary structure.

    Parameters
    ----------
    file_path: str
        Path to a yaml file.

    Returns
    -------
    yaml_dict : dict
        Dictionary returned by reading the yaml file (could be empty).
    """
    yaml_dict = {}
    if file_path:
        if not isfile(file_path):
            raise IOError('No yaml file "%s"' % file_path)
        with open(file_path) as yaml_handle:
            try:
                yaml_dict = yaml.load(yaml_handle, Loader=yaml.FullLoader)
            except AttributeError:
                yaml_dict = yaml.load(yaml_handle)
    return yaml_dict


def get_fastq_files(
        fastqs: list
) -> list:
    fastq_gz = [fastq for fastq in fastqs if '.gz' in fastq]
    if len(fastq_gz):
        return fastq_gz
    else:
        return fastqs


def get_fastq_paths(
        fastq_dirs: list
) -> list:
    fastqs = []
    for fastq_dir in fastq_dirs:
        fastqs.extend(glob.glob(fastq_dir + '/*.fastq*'))
    return fastqs


def mkdr(
        path: str,
        is_file: bool = False
) -> None:
    """Creates a folder is does not exist yet.

    Parameters
    ----------
    path : str
        Path of the file or folder to create.
    is_file : bool
        True is the path is a file path, False is it is a folder path.
    """
    if is_file:
        os.makedirs(dirname(path), exist_ok=True)
    else:
        os.makedirs(path, exist_ok=True)


def get_pfam_wget_cmd(
        pfam_dir: str
) -> str:
    """Get the validation (even through downloading) that the HMM file is there.

    Parameters
    ----------
    pfam_dir : str
        Path to the folder that should contain the "Pfam-A.*" files

    Returns
    -------
    cmd : str
        Command line to run to download the database
    """
    cmd = ''
    for ext in ['hmm', 'fasta', 'hmm.dat']:
        fp = '%s/Pfam-A.%s' % (pfam_dir, ext)
        if not isfile(fp):
            gz = 'pub/databases/Pfam/releases/Pfam35.0/Pfam-A.%s.gz' % ext
            cmd += 'wget -O %s http://ftp.ebi.ac.uk/%s\n' % (fp, gz)
            cmd += 'gunzip %s\n' % fp
    return cmd


def get_hmm_dat(
        pfam_dir: str
) -> pd.DataFrame:
    """Get as table the contents of
    http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.dat.gz

    Parameters
    ----------
    pfam_dir : str
        Path to the folder that should contain the "Pfam-A.*" files

    Returns
    -------
    pfam_pd : pd.DataFrame
        Pfam dat reformatted as table.
    """
    tsv = '%s/Pfam-A.hmm.dat.tsv' % pfam_dir
    if not isfile(tsv):
        records = []
        cur = {}
        with open('%s/Pfam-A.hmm.dat' % pfam_dir) as f:
            for line in f:
                line_decode = line.strip()
                if line_decode.startswith('//'):
                    records.append(cur)
                    cur = {}
                elif line_decode.startswith('#=GF'):
                    cur[line_decode[5:7]] = line_decode[7:].strip()
        pfam_pd = pd.DataFrame(records)
        pfam_pd.to_csv(tsv, index=False, sep='\t')
    else:
        pfam_pd = pd.read_table(tsv)
    return pfam_pd


def get_hmms_dias_cmd(
        hmm: str,
        term_pd: pd.DataFrame,
        term: str,
        odir: str
) -> tuple:
    """Collect the hmm profiles and diamond database for different term targets.

    Parameters
    ----------
    hmm : str
        Pfam-A.hmm file.
    term_pd : pd.DataFrame
        Data frame containing the Pfam-A.hmm.dat file, but reformatted above.
    odir : str
        Directory where all Pfam models will be kept.
    term : str
        Term to search among the Pfam hmms.

    Returns
    -------
    hmms_dias : dict
        [.hmm, .dmnd] files for each term.
    cmd : str
        Command to hmmfetch these term's hmms.
    """
    comp = re.compile("[\\ \-\",()%':&/.\[\]]")
    cmd = ''
    hmms_dias = {}
    pfam_out = '%s/%s' % (odir, comp.sub('_', term))
    mkdr(pfam_out)
    for r, row in term_pd.iterrows():
        acc, desc = row['AC'], comp.sub('_', row['DE'])
        name = '%s__%s' % (acc, desc)
        hmm_fp = '%s/%s.hmm' % (pfam_out, name)
        if not isfile(hmm_fp):
            cmd += 'hmmfetch %s %s >> %s\n' % (hmm, acc, hmm_fp)
        fa = '%s/fastas/%s.fa' % (odir, acc)
        fo = '%s/%s.fa' % (pfam_out, name)
        dia = '%s/%s.dmnd' % (pfam_out, name)
        if not isfile(dia):
            cmd += 'cp %s %s\n' % (fa, fo)
            cmd += 'diamond makedb --in %s -d %s\n' % (fo, dia)
            cmd += 'rm %s\n' % fo
        hmms_dias[name] = [hmm_fp, dia]
    return hmms_dias, cmd


def min_nlines(
        input_fp: str
) -> bool:
    """Check whether the number of lines in the file is above one.

    Parameters
    ----------
    input_fp : str
        Path to the input file.

    Returns
    -------
    ret : bool
        Whether the file contains at least one sequence.
    """
    ret = False
    with open(input_fp.replace('${SCRATCH_FOLDER}', '')) as f:
        for ldx, line in enumerate(f):
            if ldx > 1:
                ret = True
                break
    return ret


def not_paired(
        self,
        tech: str,
        fqs: list,
) -> bool:
    """Checks whether there are two input files, which is what is needed for
    merging. Otherwise, stop there and show a useful error message.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    fqs : list
        Paths to the input files

    Returns
    -------
    bool
        Whether the input files are not possibly pooled or not
    """
    nfiles = len(fqs)
    if nfiles != 2:
        self.outputs['outs'].setdefault((tech, self.sam_pool), []).extend(fqs)
        return True
    return False


def get_roundtrip(io) -> dict:
    roundtrip = {'to': inputs_to_scratch(io), 'from': outputs_back(io)}
    return roundtrip


def inputs_to_scratch(io) -> list:
    rsyncs, mkdirs = set(), set()
    # folders
    if ('I', 'd') in io:
        for folder_ in io[('I', 'd')]:
            folder = folder_.rstrip('/')
            src = folder_.rstrip('/').replace('${SCRATCH_FOLDER}', '')
            mkdirs.add('mkdir -p %s' % folder)
            rsyncs.add('rsync -aqruv %s/ %s' % (src, folder))
    # folders
    if ('O', 'd') in io:
        for folder in io[('O', 'd')]:
            mkdirs.add('mkdir -p %s' % folder.rstrip('/'))
    # files
    if ('I', 'f') in io:
        for file in io[('I', 'f')]:
            folder = dirname(file)
            src = file.replace('${SCRATCH_FOLDER}', '')
            mkdirs.add('mkdir -p %s' % folder)
            rsyncs.add('rsync -aqruv %s %s' % (src, file))
    return sorted(mkdirs) + sorted(rsyncs)


def outputs_back(io) -> list:
    outbound = set()
    if ('O', 'd') in io:
        # folders
        for folder_ in io[('O', 'd')]:
            folder = folder_.rstrip('/')
            src = folder_.rstrip('/').replace('${SCRATCH_FOLDER}', '')
            cmd = 'mkdir -p %s; rsync -aqruv %s/ %s' % (src, folder, src)
            cmd = 'if [ -d %s ]; then %s; fi' % (folder, cmd)
            outbound.add(cmd)
    if ('O', 'f') in io:
        # files
        for file in io[('O', 'f')]:
            src = file.replace('${SCRATCH_FOLDER}', '')
            folder = dirname(src)
            cmd = 'mkdir -p %s; rsync -aqruv %s %s' % (folder, file, src)
            cmd = 'if [ -f %s ]; then %s; fi' % (file, cmd)
            outbound.add(cmd)
    return sorted(outbound)


def get_scratch_cmds(
        self,
        key,
        soft,
        cur_cmds,
        cmds
) -> dict:
    if key in soft.io:
        roundtrip = get_roundtrip(soft.io[key])
        scratch_cmds = ['\n# Move to SCRATCH_FOLDER'] + roundtrip['to']
        scratch_cmds += ['\n# %s commands (%s)' % (soft.name, key)] + cur_cmds
        if self.config['move_back']:
            scratch_cmds += ['\n# Move from SCRATCH_FOLDER'] + roundtrip['from']
        cmds[key] = scratch_cmds
    else:
        cmds[key] = cur_cmds


def per_group_scratch(
        self,
        pool,
        soft,
        sam_cmds,
        cmds,
        commands
) -> None:
    if pool in sam_cmds:
        get_scratch_cmds(self, pool, soft, sam_cmds, cmds)
    else:
        for group in commands.pools[pool]:
            group_cmds = sam_cmds[group]
            get_scratch_cmds(self, (pool, group), soft, group_cmds, cmds)


def scratching(
        self,
        soft,
        commands
) -> dict:
    if soft.params['scratch'] and self.config.jobs:
        cmds = {}
        # Use commands.pools to unpack the commands well
        for sam, sam_cmds in soft.cmds.items():
            if isinstance(sam_cmds, list):
                get_scratch_cmds(self, sam, soft, sam_cmds, cmds)
            elif isinstance(sam_cmds, dict):
                per_group_scratch(self, sam, soft, sam_cmds, cmds, commands)
            else:
                sys.exit('The collected commands are neither list of dict!')
        return cmds
    else:
        return soft.cmds


def io_update(
        self,
        i_f=None,
        i_d=None,
        o_f=None,
        o_d=None,
        key=None
):
    for (IO_fd, val) in [
        (('I', 'f'), i_f),
        (('I', 'd'), i_d),
        (('O', 'f'), o_f),
        (('O', 'd'), o_d)
    ]:
        if not val:
            continue
        if isinstance(val, list):
            if key:
                self.outputs['io'][IO_fd].setdefault(key, set()).update(val)
            else:
                self.outputs['io'][IO_fd].update(val)
        elif isinstance(val, str):
            if key:
                self.outputs['io'][IO_fd].setdefault(key, set()).add(val)
            else:
                self.outputs['io'][IO_fd].add(val)


def to_do(
        file: str = None,
        folder: str = None
) -> bool:
    if file and isfile(file.replace('${SCRATCH_FOLDER}', '')):
        return False
    if folder and isdir(folder.replace('${SCRATCH_FOLDER}', '')):
        return False
    return True


def tech_specificity(
        self,
        data,
        tech: str,
        sam: str,
        specificity: list = []
) -> bool:
    """Returns a boolean that is True if the current technology can not be
    processed by the current tool, possibly because is has no input files.
    If there are files but the technology can not be processed, the files
    are passed as output so that they can be available to the next steps.

    Notes
    -----
    This behaviour is common to all software functions and therefore
    justifies this function, and some tools of a workflow ar not designed
    for some technologies so the pipeline should pass through these while
    retaining the files available for a next step.
    For example, the reads merging can only happen on the 'illumina' data.
    So the 'pacbio' and 'nanopore' data will just be passed as output to read
    merging (without merging) so that they can be available to the next step.

    Parameters
    ----------
    self : Commands class instance
        .outputs : dict
            All outputs
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    sam : str
        Sample name
    data : str
        Path to the input files or other input data structure
    specificity : list
        Technology that can be processed by the current software

    Returns
    -------
    bool
        Whether the technology is possibly processed using this tool
    """
    if not data or (specificity and tech not in specificity):
        self.outputs['outs'].setdefault((tech, sam), []).extend(data)
        return True
    return False


def caller(
        self,
        namespace
):
    """Calls as function the part of the software name that follows the first
    underscore.
    For example, software name "search_diamond" would call function`diamond()`.
    """
    module = sys.modules[namespace]
    for sub in ['', '_']:
        func = self.soft.name.split('_', 1)[-1] + sub
        if func == self.soft.name:
            continue
        if hasattr(module, func) and callable(getattr(module, func)):
            module_call = getattr(module, func)
            return module_call
    print('No function "%s" in module "%s"' % (func, namespace))
