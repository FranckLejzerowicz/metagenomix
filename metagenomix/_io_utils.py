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
import pandas as pd
from os.path import abspath, basename, dirname, isfile


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


def fill_fastq(fastq_paths: list, sams: set):
    """Populate a `fastq` dict with for each sample (keys) the list of
    fastq file paths (values), which would be either of length 1 if there
    is only one fastq file for the sample, or of length 2 if there are
    two, paired fastq file paths.

    Notes
    -----
    In fact, the values could have a length 4 since all *fastq* files are
    considered and there could be both `.fastq` and `.fastq.gz` file in the
    input folder. Only the `.fastq.gz` files would be used in a later stage.

    Parameters
    ----------
    fastq_paths : list
        All fastq files in the input folder
    sams : set
        All unique sample names that have a metadata

    Returns
    -------
    fastq : dict
        Fastq file path(s) per sample
    """
    fastq = {}
    for fastq_path in sorted(fastq_paths):
        for sam in sams:
            if basename(fastq_path).startswith(sam):
                if sam in fastq:
                    fastq[sam].append(abspath(fastq_path))
                else:
                    fastq[sam] = [abspath(fastq_path)]
                break
    return fastq


def mkdr(path: str, is_file: bool = False) -> None:
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


def get_pfam_wget_cmd(pfam_dir: str) -> str:
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


def get_hmm_dat(pfam_dir: str) -> pd.DataFrame:
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


def get_pfams_cmd(hmm: str, term_pd: pd.DataFrame,
                  term: str, odir: str) -> tuple:
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
    pfams : dict
        [.hmm, .dmnd] files for each term.
    cmd : str
        Command to hmmfetch these term's hmms.
    """
    comp = re.compile("[\\ \-\",()%':&/.\[\]]")
    cmd = ''
    pfams = {}
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
        pfams[name] = [hmm_fp, dia]
    return pfams, cmd


def reads_lines(file_fp: str) -> list:
    """Check if file exists or is empty.

    Parameters
    ----------
    file_fp : str
        Input file.

    Returns
    -------
    file_lines : list
        Line of the input file.
    """
    file_lines = []
    if isfile(file_fp):
        with open(file_fp) as f:
            file_lines = [x.strip() for x in f.readlines()]
    return file_lines


def check_min_lines_count(input_fp: str) -> bool:
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
    with open(input_fp) as f:
        for ldx, line in enumerate(f):
            if ldx > 1:
                ret = True
                break
    return ret


def get_out_dir_(out_dir: str, inputs: dict, sam_pool: str,
                group: str = None) -> tuple:
    """Get the output directory for either the

    Parameters
    ----------
    out_dir : str
    inputs : dict
    sam_pool : str
        Name of the current sample or pool.
    group : str

    Returns
    -------
    out_dir : str
    file_path : str
    """
    inputs = inputs[sam_pool]
    if group:
        inputs = inputs[group]
    if out_dir.endswith('after_plass'):
        file_path = inputs.replace('nuclassembly', 'assembly')
    elif out_dir.endswith('after_prodigal'):
        file_path = inputs[1]
        if 'after_spades' in file_path:
            out_dir = out_dir.replace('after_', 'after_spades_')
        elif 'after_plass' in file_path:
            out_dir = out_dir.replace('after_', 'after_plass_')
    else:
        file_path = inputs[1]
    out_dir = '%s/%s' % (out_dir, sam_pool)
    return out_dir, file_path


def get_out_dir(self, sam_pool, group: str = None) -> tuple:
    """Get the output directory for either the

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to pipeline output folder
        .inputs : dict
            Input files
    sam_pool : str
        Name of the current sample or pool.
    group : str

    Returns
    -------
    out_dir : str
    file_path : str
    """
    out_dir = self.dir
    inputs = self.inputs[sam_pool]
    if group:
        inputs = inputs[group]
    if self.dir.endswith('after_plass'):
        file_path = inputs.replace('nuclassembly', 'assembly')
    elif self.dir.endswith('after_prodigal'):
        file_path = inputs[1]
        if 'after_spades' in file_path:
            out_dir = out_dir.replace('after_', 'after_spades_')
        elif 'after_plass' in file_path:
            out_dir = out_dir.replace('after_', 'after_plass_')
    else:
        file_path = inputs[1]
    out_dir = '%s/%s' % (out_dir, sam_pool)
    return out_dir, file_path


def write_hmms(self) -> str:
    """Write a file to contain on each line, the paths to the .hmm
    files to search as part of integron_finder.

    Parameters
    ----------
    self : Commands class instance
        .dir : str
            Path to the output directory.
        .databases.hmms_dias : dict
            .hmm files per name of profile.

    Returns
    -------
    hmms_fp : str
        Path to the file containing the paths to the .hmm files to search.
    """
    hmms_fp = ''
    if self.databases.hmms_dias:
        mkdr(self.dir)
        hmms_fp = '%s/hmm.txt' % self.dir
        o = open(hmms_fp, 'w')
        for target, gene_hmms in self.databases.hmms_dias.items():
            for gene, (hmm, _) in gene_hmms.items():
                o.write('%s\n' % hmm)
        o.close()
    return hmms_fp


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
            cmd = 'if [ -d %s ]; then %s' % (folder, cmd)
            outbound.add(cmd)
    if ('O', 'f') in io:
        # files
        for file in io[('O', 'f')]:
            src = file.replace('${SCRATCH_FOLDER}', '')
            folder = dirname(src)
            cmd = 'mkdir -p %s; rsync -aqruv %s %s' % (folder, file, src)
            cmd = 'if [ -f %s ]; then %s' % (file, cmd)
            outbound.add(cmd)
    return sorted(outbound)


def get_scratch_cmds(key, soft, cur_cmds, cmds):
    if key in soft.io:
        roundtrip = get_roundtrip(soft.io[key])
        scratched_cmds = ['# Move to SCRATCH_FOLDER'] + roundtrip['to']
        scratched_cmds += ['\n# %s commands (%s)' % (soft.name, key)] + cur_cmds
        scratched_cmds += ['\n# Move from SCRATCH_FOLDER'] + roundtrip['from']
        cmds[key] = scratched_cmds
    else:
        cmds[key] = cur_cmds


def per_coassembly_scratch(pool, soft, sam_cmds, cmds, commands):
    if pool in sam_cmds:
        get_scratch_cmds(pool, soft, sam_cmds, cmds)
    else:
        for group in commands.pools[pool]:
            group_cmds = sam_cmds[group]
            get_scratch_cmds((pool, group), soft, group_cmds, cmds)


def scratching(soft, commands) -> dict:
    if soft.params['scratch']:
        cmds = {}
        # Use commands.pools to unpack the commands well
        for sam, sam_cmds in soft.cmds.items():
            if isinstance(sam_cmds, list):
                get_scratch_cmds(sam, soft, sam_cmds, cmds)
            elif isinstance(sam_cmds, dict):
                per_coassembly_scratch(sam, soft, sam_cmds, cmds, commands)
            else:
                sys.exit('The collected commands are neither list of dict!')
        return cmds
    else:
        return soft.cmds


def io_update(self, i_f=None, i_d=None, o_f=None, o_d=None, key=None):

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

