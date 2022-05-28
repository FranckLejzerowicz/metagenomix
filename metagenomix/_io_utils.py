# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import gzip
import yaml
import subprocess
import pandas as pd
from os.path import basename, dirname, isfile


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


def get_fastq_header(
        fastq_path: str
) -> str:
    """Get the first line of the fastq file.

    Parameters
    ----------
    fastq_path : str
        Path to a fastq.gz or fastq file.

    Returns
    -------
    fastq_line : str
        First line of the fastq file.
    """
    fastq_line = ''
    if fastq_path.endswith('fastq.gz'):
        with gzip.open(fastq_path) as f:
            for line in f:
                fastq_line = line.decode().strip()
                break
    elif fastq_path.endswith('fastq'):
        with open(fastq_path) as f:
            for line in f:
                fastq_line = line.strip()
                break
    return fastq_line


def get_cat_zcat(
        fastq_fp: str
) -> str:
    """Get the command `cat` or `gzcat` to run for fastq editing.

    Parameters
    ----------
    fastq_fp : str
        Fastq file path.

    Returns
    -------
    cat : str
        The command to use (`cat` or `gzcat`).
    """
    cat = "cat"
    if fastq_fp.endswith('.fastq.gz'):
        cat = "gzcat"
    return cat


def edit_fastq_cmd(
        fastq_fp: str,
        num: int,
        source: str
) -> str:
    """Get the unix command to run on each fastq file that
    needs editing, depending on the source of fastq file.

    Parameters
    ----------
    fastq_fp : str
        Fastq file path.
    num : int
        Can be 1 or 2 for the forward and reverse read, respectively.
    source : str
        Source of the fastq that defines format and editing type.

    Returns
    -------
    cmd : str
        Unix command to run to edit a fastq file.
    """
    cat = get_cat_zcat(fastq_fp)
    cmd = '%s %s | ' % (cat, fastq_fp)
    cmd += "awk '{ if (NR%s4==1) " % '%'
    if source == 'illumina':
        cmd += "{ print $1\"/%s\" } " % str(num)
    elif source == 'ebi':
        cmd += "{ gsub(/.%s .*/,\"/%s\",$0); print } " % (num, num)
    cmd += "else if (NR%s2 == 1) { print \"+\" } " % '%'
    cmd += "else { print } }' | gzip > %s_renamed\n" % fastq_fp
    cmd += "mv %s_renamed %s" % (fastq_fp, fastq_fp)
    if fastq_fp.endswith('.fastq'):
        cmd += ".gz\n"
    else:
        cmd += "\n"
    return cmd


def get_edit_fastq_cmd(
        reads: str,
        num: int
) -> str:
    """Get the unix command to run on each fastq file that needs editing.

    Parameters
    ----------
    reads : str
    num : int
        Can be 1 or 2 for the forward and reverse read, respectively.

    Returns
    -------
    cmd : str
        Unix command to run to edit a fastq file.
    """
    cmd = ''
    fastq_fp = reads[num - 1]
    line = get_fastq_header(fastq_fp)
    if not line.endswith('/%s' % num):
        line_split = line.split()
        if len(line_split) > 1:
            if line_split[1].startswith('%s:N:' % num):
                cmd = edit_fastq_cmd(fastq_fp, num, 'illumina')
            elif line_split[0].endswith('.%s' % num):
                cmd = edit_fastq_cmd(fastq_fp, num, 'ebi')
    return cmd


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


def get_pfam_file(
        pfam_hmm: str
) -> bool:
    """Get the validation (even through downloading) that the HMM file is there.

    Parameters
    ----------
    pfam_hmm : str
        "Pfam-A.hmm" or "Pfam-A.hmm.dat""

    Returns
    -------
    ret : bool
        Whether the target file is finally available or not.
    """
    ret = True
    if not isfile(pfam_hmm):
        user = input('"%s" not found...\nDownload? ([Y]/n)\n:' % pfam_hmm)
        if not user or user.lower()[0] == 'y':
            subprocess.call(
                ('wget -O %s http://ftp.ebi.ac.uk/pub/databases/Pfam/releases'
                 '/Pfam35.0/%s.gz' % (pfam_hmm, basename(pfam_hmm))).split())
        else:
            ret = False
    return ret


def get_hmm_dat(
        dat: str,
        tsv: str
) -> pd.DataFrame:
    """Get as table the contents of
    http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.dat.gz

    Parameters
    ----------
    dat : str
        Path to Pfam dat (Pfam-A.hmm.dat.gz)
    tsv : str
        Path to reformatted, output Pfam dat (Pfam-A.hmm.dat.tsv)

    Returns
    -------
    pfam_dat_pd : pd.DataFrame
        Pfam dat reformatted as table.
    """
    if not isfile(tsv):
        records = []
        cur = {}
        with open(dat) as f:
            for line in f:
                line_decode = line.strip()
                if line_decode.startswith('//'):
                    records.append(cur)
                    cur = {}
                elif line_decode.startswith('#=GF'):
                    cur[line_decode[5:7]] = line_decode[7:].strip()
        pfam_dat_pd = pd.DataFrame(records)
        pfam_dat_pd.to_csv(tsv, index=False, sep='\t')
    else:
        pfam_dat_pd = pd.read_table(tsv)
    return pfam_dat_pd


def get_pfams_cmd(
        hmm: str,
        term_pd: pd.DataFrame,
        pfam_dir: str,
        term: str
) -> tuple:
    """Collect the hmm profiles and diamond database for different term targets.

    Parameters
    ----------
    hmm : str
        Pfam-A.hmm file.
    term_pd : pd.DataFrame
        Data frame containing the Pfam-A.hmm.dat file, but reformatted above.
    pfam_dir : str
        Pfam database directory passed by the user.
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
    pfam_out = '%s/%s' % (pfam_dir, term)
    mkdr(pfam_out)
    for r, row in term_pd.iterrows():
        acc, desc = row['AC'], comp.sub('_', row['DE'])
        name = '%s__%s' % (acc, desc)
        hmm_fp = '%s/%s.hmm' % (pfam_out, name)
        if not isfile(hmm_fp):
            cmd += 'hmmfetch %s %s >> %s\n' % (hmm, acc, hmm_fp)
        fa = '%s/fastas/%s.fa' % (pfam_dir, acc)
        fo = '%s/fastas/%s.fa' % (pfam_dir, name)
        dia = '%s/%s.dmnd' % (pfam_out, name)
        if not isfile(dia):
            cmd += 'cp %s %s\n' % (fa, fo)
            cmd += 'diamond makedb --in %s -d %s\n' % (fo, dia)
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


def count_reads_cmd(idx: int, input_path: str, out: str, sam: str) -> str:
    """Get the command to count the reads in a fasta, fastq or fastq.gz file.

    Parameters
    ----------
    idx : int
        Unique, incemental numeric identifier.
    input_path : str
        Path to the input file name.
    out : str
        Path to the output file name.
    sam: str
        Sample name.

    Returns
    -------
    cmd : str
        Command to count reads and output the counts.
    """
    if input_path.endswith('fastq.gz'):
        cmd = "n%s=`gzcat %s | wc -l | " % (idx, input_path)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
    elif input_path.endswith('fasta'):
        cmd = "n%s=`wc -l %s | " % (idx, input_path)
        cmd += "sed 's/ //g' | awk '{x=$1/2; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/2; print x}'`"
    elif input_path.endswith('fastq'):
        cmd = "n%s=`wc -l %s | " % (idx, input_path)
        cmd += "sed 's/ //g' | awk '{x=$1/4; print x}'`\n"
        # cmd += "cut -d ' ' -f 1 | awk '{x=$1/4; print x}'`"
    else:
        raise IOError("Input sequence file invalid %s" % input_path)
    if idx:
        cmd += 'echo "%s,2,${n%s}" >> %s\n' % (sam, idx, out)
    else:
        cmd += 'echo "%s,1,${n%s}" > %s\n' % (sam, idx, out)
    return cmd


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


def get_out_dir(
        out_dir: str,
        inputs: dict,
        sam_pool: str,
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


def write_hmms(out_dir: str, hmms_dias: dict) -> str:
    """Write a file to contain on each line, the paths to the .hmm
    files to search as part of integron_finder.

    Parameters
    ----------
    out_dir : str
        Path to the output directory.
    hmms_dias : dict
        .hmm files per name of profile.

    Returns
    -------
    hmms_fp : str
        Path to the file containing the paths to the .hmm files to search.
    """
    hmms_fp = ''
    if hmms_dias:
        mkdr(out_dir)
        hmms_fp = '%s/hmm.txt' % out_dir
        o = open(hmms_fp, 'w')
        for target, gene_hmms in hmms_dias.items():
            for gene, (hmm, _) in gene_hmms.items():
                o.write('%s\n' % hmm)
        o.close()
    return hmms_fp