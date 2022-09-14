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
import hashlib
import itertools
import pandas as pd
from tabulate import tabulate
from os.path import basename, dirname, isdir, isfile



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
        sam: str,
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
    sam : str
        Name of the current sample
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
        self.soft.add_status(tech, sam, 1, message='unpaired reads')
        return True
    return False


def status_update(
        self,
        tech: str,
        inputs: list,
        pool: str = None,
        group: str = None,
        genome: str = None
) -> None:
    """Potentially add the fastq files to the status (files to generate).

    Parameters
    ----------
    self : Commands class instance
        .sam_pool : str
            Sample name
        .soft
            Software class instance
    tech : str
        Technology: 'illumina', 'pacbio', or 'nanopore'
    inputs : list
        Paths to input files
    pool : str
        Name of the current co-assembly
    group : str
        Name of the current co-assembly group
    genome : str
        MAGs/Genomes folder name or empty string (for assembly contigs)
    """
    to_dos = [x for x in inputs if to_do(x)]
    if to_dos:
        if pool:
            self.soft.add_status(tech, pool, to_dos,
                                 group=group, genome=genome)
        else:
            self.soft.add_status(tech, self.sam_pool, to_dos,
                                 group=group, genome=genome)


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
    if not isinstance(key, tuple):
        key = (key,)
    for (IO_fd, val) in [
        (('I', 'f'), i_f),
        (('I', 'd'), i_d),
        (('O', 'f'), o_f),
        (('O', 'd'), o_d)
    ]:
        if not val:
            continue
        if isinstance(val, list):
            self.outputs['io'][IO_fd].setdefault(key, set()).update(val)
        elif isinstance(val, str):
            self.outputs['io'][IO_fd].setdefault(key, set()).add(val)


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
    justifies this function, and some softwares of a workflow ar not designed
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
        self.soft.add_status(tech, sam, 0, message='technology incompatible')
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


def edit_sample_or_pool(not_done):
    sample_or_pool = 'pool'
    if not_done['group'].isnull().all():
        sample_or_pool = 'sample'
    not_done.rename(columns={'sample_or_pool': sample_or_pool}, inplace=True)


def show_not_done(not_done):
    """

    Parameters
    ----------
    not_done : pd.DataFrame

    Returns
    -------
    status_summary : str
        Summary of the analyses that remain to do
    """
    # replace "sample_or_pool" by "sample" or "pool" depending on the data type
    edit_sample_or_pool(not_done)
    # get the data that remains to analyse for the current software
    todo_pd = not_done.loc[not_done['status'] == 'To do'].copy()
    # make a one-liner summarizing what remains to analyse (per tech, pool etc)
    status_summary = summarize_status(todo_pd)
    print(status_summary)
    return status_summary


def print_vals(vals) -> str:
    vals = '; '.join(['%ss' % v if int(v.split()[0]) > 1 else v for v in
                        vals])
    return vals


def summarize_status(tab) -> str:
    """

    Parameters
    ----------
    tab

    Returns
    -------
    status_summary : str
        Summary of the analyses that remain to do
    """
    cols = tab.columns.tolist()[:2] + ['group', 'genome']
    dat = tab[cols].apply(lambda x: x.nunique()).astype(str).to_dict().items()
    vals = [' '.join(x[::-1]) for x in dat if int(x[1])]
    if vals:
        status_summary = 'to do:\t%s' % print_vals(vals)
    else:
        status_summary = 'nothing to do'
    return status_summary


def pretty_print(tab, header, typ) -> str:
    # stackoverflow.com/questions/41593793/
    # printtabulate-to-pretty-print-multiindex-pandas
    p_tab = tabulate(tab, headers=header, tablefmt='orgtbl')
    tab_width = [len(x) for x in p_tab.split('\n') if '|--' in x][0]
    dash = '-' * tab_width
    if not typ:
        typ = 'needed input'
    gap = (tab_width // 2) - (len(typ) // 2)
    header = '%s\n%s%s\n%s\n' % (dash, ' ' * gap, typ, dash)
    out_tab = header + p_tab + '\n' + dash
    p_tab = '\n\t\t' + re.sub('\\n', '\n%s' % ('\t' * 2), out_tab) + '\n'
    print(p_tab)
    return out_tab


def get_pivot_table_cols(
        tab: pd.DataFrame,
) -> list:
    """Get the columns of the status table that are to be used
    to report progress.

    Parameters
    ----------
    tab : pd.DataFrame

    Returns
    -------
    pivot_table_cols : list
    """
    tab_cols = tab.columns.tolist()
    cols = tab_cols[:2] + ['group']
    if 'genome' in tab_cols:
        cols += ['genome']
    pivot_table_cols = [col for col in cols if not tab[col].isnull().all()]
    return pivot_table_cols


def pivot_not_done_multi_header(t) -> tuple:
    # reorder column multiindex to get the tech above and group/genome below
    t.columns = t.columns.reorder_levels(range(len(t.columns.levels))[::-1])
    # sort table columns to suit the reordering
    t = t[[
        it for it in itertools.product(*[
            lev[::-1] if lx else lev for lx, lev in enumerate(t.columns.levels)
        ])
    ]]
    # edit the multi-index lie for pretty printing
    h = [t.index.names[0] + '/' + t.columns.names[0]] + list(
        map('\n'.join, t.columns.tolist()))
    return t, h


def print_tab_pv(
        tab_pv: pd.DataFrame,
        cols: list,
        typ: str = '',
) -> tuple:
    """

    Parameters
    ----------
    tab_pv : pd.DataFrame
    cols : list
    typ : str

    Returns
    -------
    out_tab : str
        Table to write out
    skip : bool
        Whether the print table is only for sample per tech
        (will break without printing content details since there is no content)
    """
    skip = False
    if len(cols[2:]) > 1:
        tab_pv_sorted, h = pivot_not_done_multi_header(tab_pv)
        out_tab = pretty_print(tab_pv_sorted, h, typ)
    elif len(cols[2:]):
        h = [tab_pv.index.names[0]] + ['\n'.join(x[::-1]) if x[0][:3] == 'Run'
                                       else x[-1] for x in tab_pv.columns]
        out_tab = pretty_print(tab_pv, h, typ)
    else:
        h = [tab_pv.index.names[0]] + tab_pv.columns.tolist()
        out_tab = pretty_print(tab_pv, h, typ)
        skip = True
    return out_tab, skip


def print_message(
        tab: pd.DataFrame,
) -> list:
    """Loop to print both the data remaining to analyze, both in terms of
    job counts and job contents (sample, co-assembly pool's group, genome).

    Parameters
    ----------
    tab : pd.DataFrame
        Status table for the

    Returns
    -------
    out_messages : list
        Pretty tables
    """
    out_messages = []
    cols = get_pivot_table_cols(tab)
    for (typ, aggfunc) in [
        ('counts', 'count'),
        ('contents', (lambda x: '\n'.join([y for y in x if y is not None]))),
    ]:
        tab_pv = tab[cols].pivot_table(
            values=cols[2:],
            columns=cols[0],
            index=cols[1],
            aggfunc=aggfunc)
        out_tab, break_here = print_tab_pv(tab_pv, cols, typ)
        out_messages.append(out_tab)
        if break_here:
            break
    return out_messages


def print_messages(not_done):
    out_messages = []
    for mess, mess_pd in not_done.groupby('message'):
        message = 'Attention -> "%s":' % mess
        print('\n\t\t%s ' % message)
        out_messages.append(message)
        out_messages.extend(print_message(mess_pd))
    return out_messages


def print_needed(path, not_done) -> str:
    tab = not_done.loc[not_done['status'] != 'To do'].copy()
    if len(tab):
        tab['status'] = ['\n'.join(map(basename, x)) for x in tab['status']]
        path = ['%s> %s' % (' ' * x, y) for x, y in enumerate(path[:-1])]
        status = 'Run:\n%s\n%s*' % ('\n'.join(path[:-1])[2:],
                                    path[-1].replace('> ', '=> *'))
        tab = tab.rename(columns={'status': status})
        cols = get_pivot_table_cols(tab) + [status]
        tab_pv = tab[cols].pivot_table(
            values=cols[2:],
            columns=cols[0],
            index=cols[1],
            aggfunc=(lambda x: '\n'.join(set(y for y in x if y is not None))))
        out_tab, _ = print_tab_pv(tab_pv, cols)
        return out_tab


def print_status_table(
        soft,
        show_status: bool = False
) -> None:
    if soft.status:
        status = pd.DataFrame(soft.status, columns=[
            'tech', 'sample_or_pool', 'status', 'group', 'message', 'genome'])
        status.drop_duplicates(inplace=True)
        not_done = status.loc[status['status'] != 'Done'].copy()
        if len(not_done):
            soft.tables.append(show_not_done(not_done))
            if show_status:
                # pretty table of data that remains to analyse when message
                soft.tables.extend(print_messages(not_done))
                # pretty table of data that remains to analyse upfront
                soft.tables.append(print_needed(soft.path, not_done))
        else:
            soft.tables.append('done')
            print('done')
    else:
        soft.tables.append('done')
        print('done')


def compute_hash(hash_string):
    h = hashlib.blake2b(digest_size=10)
    h.update(str(hash_string).encode('utf-8'))
    hashed = str(h.hexdigest())
    return hashed


def get_md5(fp):
    if isdir(fp):
        return ''
    md5 = hashlib.md5()
    with open(fp, "rb") as f:
        for chunk in f:
            md5.update(chunk)
    md5 = str(md5.hexdigest())
    return md5


def get_dates(run_fp) -> list:
    """Get the different dates at which the current pipeline configuration
    was run.

    Parameters
    ----------
    run_fp
        Path to the versioning file for this pipeline configuration
    Returns
    -------
    runs : list
        All the dates at which this pipeline configuration was run
    """
    runs = []
    if isfile(run_fp):
        with open(run_fp) as f:
            for line in f:
                if line.startswith('Date'):
                    runs.append(line.strip().split()[-1])
    return runs


def get_input_info(dat):
    m = '%s input unit' % dat.shape[0]
    if dat.shape[0] > 1:
        m += 's'
    info = '%s scripts (%s):' % (dat.name.nunique(), m)
    return info


def get_res_info(dat):
    d = [(len(v), [len(x) for x in v.values()]) for v in dat.values()]
    tech = '%s tech' % len(d)
    if len(d) > 1:
        tech += 's'
    fs = sum([x[0] for x in d])
    folders = '%s folder' % fs
    if fs > 1:
        folders += 's'
    av = round(sum([round(sum(x[1])/len(x[1]), 2) for x in d]) / len(d), 2)
    files = '%s file' % av
    if av > 1:
        files += 's'
    info = '%s; (%s and %s per folder)' % (tech, folders, files)
    return info


def get_size_info(dat):
    fs = '%s folder' % len(dat)
    if len(dat) > 1:
        fs += 's'
    sizes = {}
    for folder, size_ in dat.items():
        if size_ < 1024:
            size = "%s bytes" % size_
        elif size_ < 1024 * 1024:
            size = "%s KB" % round(size_ / 1024, 2)
        elif size_ < 1024 * 1024 * 1024:
            size = "%s MB" % round(size_ / (1024 * 1024), 2)
        elif size_ < 1024 * 1024 * 1024 * 1024:
            size = "%s GB" % round(size_ / (1024 * 1024 * 1024), 2)
        else:
            size = str(size_)
        sizes[folder] = size
    info = '%s (./%s)' % (fs, '; ./'.join(['='.join(x) for x in sizes.items()]))
    return info, sizes
