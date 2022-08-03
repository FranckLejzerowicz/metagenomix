#!/usr/bin/env python3
"""Fetch hierarchical mappings and descriptions of a given list of KEGG
Orthology (KO) entries from the KEGG server.

Usage:
    python me.py ko.list

Notes:
    This script utilizes the official KEGG API (https://www.kegg.jp/kegg/rest/
    keggapi.html) to query the KEGG database and retrieve relevant information
    in the following categories: orthology (ko), module, pathway, reaction,
    reaction class, compound, and disease.
    Tested and working with KEGG release 97.0+.

Restrictions:
    The official website states:
    "KEGG API is provided for academic use by academic users belonging to
    academic institutions." (https://www.kegg.jp/kegg/rest/)
    "The maximum number of identifiers that can be given is 10 (per query)."
    (https://www.kegg.jp/kegg/rest/keggapi.html)

References:
    The latest KEGG paper (Kanehisa et al., 2021):
    https://academic.oup.com/nar/article/49/D1/D545/5943834
"""

import sys
from time import sleep
from datetime import datetime
from functools import partial
from urllib.request import urlopen, HTTPError, URLError


__author__ = 'Qiyun Zhu'
__license__ = 'BSD-3-Clause'
__version__ = '0.0.1-dev'
__email__ = 'qiyunzhu@gmail.com'


# network connection parameters
server = 'http://rest.kegg.jp/'
step = 10      # no. of entries per query (max. 10 according to policy)
delay = 2      # time gap between two queries (sec)
retries = 5    # no. retries on failed query
timeout = 60   # waiting time before giving up (sec)


def fetch(api):
    """Fetch content from KEGG server.

    Parameters
    ----------
    api : str
        RESTful API command.

    Returns
    -------
    list of str
        Lines of text.
    """
    for i in range(retries):
        if i:
            print('Retrying...', end=' ', flush=True)
            sleep(delay)
        try:
            with urlopen(server + api, timeout=timeout) as response:
                return response.read().decode('utf-8').splitlines()
        except (HTTPError, URLError) as e:
            print(f'{e.code} {e.reason}.', end=' ', flush=True)


def kegg_info():
    """Return current KEGG release's version and statistics.

    Returns
    -------
    list of str
        KEGG release information.
    """
    return [x[17:] for x in fetch('info/kegg')]


def batch_query(ids, cmd, f, name='entries', step=10):
    """Perform batch query and retrieve results.

    Parameters
    ----------
    ids : list of str
        Entries to query.
    cmd : str
        API command ("get", "list", etc.)
    f : function
        Function to convert retrieved text into data.
    name : str, optional
        Task name to display.
    step : int, optional
        Number of queries to submit per time.

    Returns
    -------
    dict of dict
        Retrieved data.
    """
    data = {}
    print(f'Querying {len(ids)} {name}...', end=' ', flush=True)
    counter = 0
    ids = sorted(ids)
    for i in range(0, len(ids), step):
        batch = ids[i:i + step]
        text = fetch(cmd + '/' + '+'.join(batch))
        data = {**data, **f(text)}
        counter += len(batch)
        print(str(counter), end=' ', flush=True)
        sleep(delay)
    print('done.', flush=True)
    return data


def parse_list(text, code=None):
    """Parse KEGG list files.

    Parameters
    ----------
    text : list of str
        KEGG list text.
    code : str, optional
        Database code to strip from beginning.

    Returns
    -------
    dict
        Entry to description dictionary.

    Raises
    ------
    ValueError
        Entry has unexpected format (e.g., different code).
    """
    if code:
        ll = len(code) + 1
    res = {}
    for line in text:
        key, value = line.split('\t')
        if code:
            if not key.startswith(code + ':'):
                raise ValueError(f'Unexpected entry: {key}.')
            key = key[ll:]
        res[key] = value
    return res


def parse_flat(text, skeys=(), mkeys=()):
    """Parse KEGG flat files.

    Parameters
    ----------
    text : list of str
        KEGG flat text.
    skeys : tuple of str, optional
        Single keys to retrieve.
    mkeys : tuple of str, optional
        Multiple keys to retrieve.

    Returns
    -------
    dict of dict
        Processed data of each key under each entry.

    Examples
    --------
    ENTRY       K00699                      KO
    NAME        UGT
    DEFINITION  glucuronosyltransferase [EC:2.4.1.17]
    PATHWAY     ko00040  Pentose and glucuronate interconversions
                ko00053  Ascorbate and aldarate metabolism
                ko00140  Steroid hormone biosynthesis
    MODULE      M00014  Glucuronate pathway (uronate pathway)
                M00129  Ascorbate biosynthesis, animals, glucose-1P => ...
    DISEASE     H00208  Hyperbilirubinemia
                H01593  Osteoporosis
    """
    # data structure
    data = {}

    # current record
    entry = None  # current entry
    mkey = None   # current multi key

    # line heads
    sheads = tuple(x.upper().ljust(12) for x in skeys)
    mheads = tuple(x.upper().ljust(12) for x in mkeys)

    for line in text:

        # record starts
        if line.startswith('ENTRY       '):
            entry = line[12:].split()[0]
            data[entry] = {}
            continue

        # record ends
        if line == '///':
            entry, mkey = None, None

        # single keys
        if line.startswith(sheads):
            skey = skeys[sheads.index(line[:12])]
            data[entry][skey] = line[12:].rstrip(';')
            mkey = None
            continue

        # multi keys
        if line.startswith(mheads):
            mkey = mkeys[mheads.index(line[:12])]
            data[entry][mkey] = []

        # clear current key
        elif not line.startswith('            '):
            mkey = None

        # append targets to a multi key
        if mkey:
            data[entry][mkey].append(line[12:])

    return data


def extract_targets(data, keys):
    """Extract target entries and definitions from multiline terms.

    Parameters
    ----------
    data : dict of dict
        Main data structure.
    keys : list or tuple of str
        Keys under which targets will be extracted.

    Returns
    -------
    dict of dict
        Definitions of individual targets under each key.

    Examples
    --------
    K00699  glucuronosyltransferase [EC:2.4.1.17] [RN:R01383]
    K01195,K14756  beta-glucuronidase [EC:3.2.1.31] [RN:R01478]
    K00002  alcohol dehydrogenase (NADP+) [EC:1.1.1.2] [RN:R01481]
    """
    names = {x: {} for x in keys}
    for entry, datum in data.items():
        for key in keys:
            if key not in datum:
                continue
            res = []
            for line in datum[key]:

                # attempt to extract targets and names
                left, found, right = line.partition('  ')
                if not found:
                    continue
                targets = []
                for field in left.split(','):

                    # validate target entry format
                    if not field.isalnum():
                        targets = []
                        break
                    targets.append(field)

                # add one or multiple targets
                if targets:
                    res.extend(targets)
                    for target in targets:
                        names[key][target] = right

            datum[key] = sorted(set(res))
    return names


def extract_dblinks(data, dbs, key='dblinks'):
    """Extract database links from multiline terms.

    Parameters
    ----------
    data : dict of dict
        Main data structure.
    dbs : dict of str
        Map of database codes to names.
    key : str, optional
        Key of database link terms.

    Examples
    --------
    RN: R01478 R04979 R07818 R08127 R08260 R10830
    COG: COG3250
    GO: 0004566
    """
    for entry, datum in data.items():
        if key not in datum:
            continue
        for line in datum[key]:
            try:
                code, targets = line.split(': ', 1)
            except IndexError:
                continue
            if code in dbs:
                datum[dbs[code]] = targets.split()
        del(datum[key])


def write_smap(data, key, fname):
    """Write one-to-one mapping to file.

    Parameters
    ----------
    data : dict of dict
        Main data structure.
    key : str
        Key of data to write.
    fname : str
        Output file name.
    """
    with open(fname, 'w') as f:
        for entry, datum in data.items():
            if key in datum:
                print(entry, datum[key], sep='\t', file=f)


def write_mmap(data, key, fname):
    """Write one-to-many mapping to file.

    Parameters
    ----------
    data : dict of dict
        Main data structure.
    key : str
        Key of data to write.
    fname : str
        Output file name.
    """
    with open(fname, 'w') as f:
        for entry, datum in data.items():
            if key in datum:
                targets = []
                for value in datum[key]:
                    targets.extend(value.split(','))
                print(entry, '\t'.join(targets), sep='\t', file=f)


def write_names(names, fname):
    """Write names / descriptions of entries to file.

    Parameters
    ----------
    names : dict
        Name dictionary.
    fname : str
        Output file name.
    """
    with open(fname, 'w') as f:
        for key, name in sorted(names.items()):
            print(key, name, sep='\t', file=f)


def write_all(name, data, skeys=[], mkeys=[]):
    """Write all data to file.

    Parameters
    ----------
    name : str
        Name of current analysis.
    data : dict of dict
        Main data structure.
    skeys : iterable of str
        Single keys of data to write.
    mkeys : iterable of str
        Multiple keys of data to write.
    """
    for key in skeys:
        write_smap(data, key, f'{name}_{key}.txt')
    for key in mkeys:
        stem = 'ko' if key == 'orthology' else key
        write_mmap(data, key, f'{name}-to-{stem}.txt')


def rename_paths(code, data, names=None):
    """Convert pathway entries from "map", "rn" etc. to "ko".

    Parameters
    ----------
    code : str
        Expected pathway code.
    data : dict of dict
        Main data structure.
    names : dict of dict
        Also rename pathways in name dictionary.

    Raises
    ------
    ValueError
        Entry has unexpected format (e.g., different code).
    """
    ll = len(code) + 5
    for entry, datum in data.items():
        if 'pathway' in datum:
            newpaths = []
            for path in datum['pathway']:
                if len(path) != ll or not path.startswith(code):
                    raise ValueError(f'Unexpected pathway entry: {path}.')
                newpaths.append(f'ko{path[-5:]}')
            datum['pathway'] = newpaths
    if names and 'pathway' in names:
        for path in names['pathway']:
            if len(path) != ll or not path.startswith(code):
                raise ValueError(f'Unexpected pathway entry: {path}.')
        names['pathway'] = {f'ko{k[-5:]}': v for k, v in names[
            'pathway'].items()}


def get_ecs(definition):
    """Extract EC numbers from a KO definition.

    Parameters
    ----------
    definition : str
        KO definition.

    Returns
    -------
    list of str
        Extracted EC numbers.

    Examples
    --------
    K00930  acetylglutamate kinase [EC:2.7.2.8]
    K02618  oxepin-CoA hydrolase [EC:3.3.2.12 1.2.1.91]
    K09866  aquaporin-4
    """
    if definition.endswith(']'):
        idx = definition.find(' [EC:')
        if idx > 0:
            return definition[idx + 5:-1].split()


def get_compounds(equation):
    """Extract compound entries from an equation.

    Parameters
    ----------
    equation : str
        Equation string.

    Returns
    -------
    list of str
        Compounds extracted from the left side of equation.
    list of str
        Compounds extracted from the right side of equation.

    Examples
    --------
    C00068 + C00001 <=> C01081 + C00009
    C00029 + C00001 + 2 C00003 <=> C00167 + 2 C00004 + 2 C00080
    G10481(n+1) + G10620 <=> G10481(n) + G11108
    C17207(n) + (n-2) C00009 <=> C20861 + (n-2) C00636
    """
    res = []
    for side in equation.split(' <=> '):
        cpds = []
        for field in side.split(' + '):
            idx = field.find('C')
            if idx >= 0:
                cpd = field[idx:idx + 6]
                if len(cpd) == 6 and cpd[1:].isdigit():
                    cpds.append(cpd)
        res.append(sorted(set(cpds)))
    return res


def get_classes(data, key='class'):
    """Extract multiple classes from a single line.

    Parameters
    ----------
    data : dict of dict
        Main data structure.
    key : str, optional
        Key under which classes will be extracted.

    names
    -----
    A class line is delimited by "; ". Example:
        Pathway modules; Carbohydrate metabolism
    """
    for entry, datum in data.items():
        if key in datum:
            datum[key] = datum[key].split('; ')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    print(f'Task started at {datetime.now()}.')

    # get KEGG release info
    text = kegg_info()
    print('KEGG ' + text[1])
    with open('kegg_info.txt', 'w') as f:
        for line in text:
            print(line, file=f)

    # read query KOs
    with open(sys.argv[1], 'r') as f:
        kos = sorted(set(x.split('\t')[0] for x in f.read(
            ).splitlines() if not x.startswith('#')))
    print(f'KO entries to query: {len(kos)}.')

    # orthology (KO)
    skeys = ('name', 'definition')
    mkeys = ('module', 'pathway', 'disease', 'dblinks')
    f = partial(parse_flat, skeys=skeys, mkeys=mkeys)
    data = batch_query(kos, 'get', f, name='KOs')
    for ko, datum in data.items():
        if 'definition' in datum:
            ecs = get_ecs(datum['definition'])
            if ecs:
                datum['ec'] = ecs
    names = extract_targets(data, ('module', 'pathway', 'disease'))
    extract_dblinks(data, {'RN': 'reaction', 'COG': 'cog', 'GO': 'go'})
    mds = names['module'].keys()
    paths = names['pathway'].keys()
    dses = names['disease'].keys()
    rns = set().union(*[x['reaction'] for x in data.values()
                        if 'reaction' in x])
    mkeys = ('module', 'pathway', 'disease', 'ec', 'reaction', 'cog', 'go')
    write_all('ko', data, skeys, mkeys)

    # reaction
    skeys = ('name', 'definition', 'equation', 'enzyme')
    mkeys = ('orthology', 'module', 'pathway', 'rclass')
    f = partial(parse_flat, skeys=skeys, mkeys=mkeys)
    data = batch_query(rns, 'get', f, name='reactions')
    names = extract_targets(data, mkeys)
    for entry, datum in data.items():
        if 'enzyme' in datum:
            datum['enzyme'] = datum['enzyme'].split()
        if 'equation' in datum:
            left, right = get_compounds(datum['equation'])
            if left:
                datum['left_compound'] = left
            if right:
                datum['right_compound'] = right
            both = sorted(set(left + right))
            if both:
                datum['compound'] = both
    rename_paths('rn', data, names)
    skeys = ('name', 'definition', 'equation')
    mkeys = ('orthology', 'module', 'pathway', 'rclass', 'enzyme', 'compound',
             'left_compound', 'right_compound')
    write_all('reaction', data, skeys, mkeys)
    cpds = set().union(*[x['compound'] for x in data.values()
                         if 'compound' in x])
    rcs = names['rclass'].keys()
    mds = set(mds).union(names['module'].keys())
    paths = set(paths).union(names['pathway'].keys())

    # reaction class
    f = partial(parse_list, code='rc')
    names = batch_query(rcs, 'list', f, name='reaction classes')
    write_names(names, 'rclass_name.txt')

    # compound
    f = partial(parse_list, code='cpd')
    names = batch_query(cpds, 'list', f, name='compounds')
    write_names(names, 'compound_name.txt')

    # module
    skeys = ('name', 'definition', 'class')
    mkeys = ('orthology', 'pathway', 'reaction', 'compound')
    f = partial(parse_flat, skeys=skeys, mkeys=mkeys)
    data = batch_query(mds, 'get', f, name='modules')
    get_classes(data)
    names = extract_targets(data, mkeys)
    rename_paths('map', data, names)
    skeys = ('name', 'definition')
    mkeys = ('orthology', 'pathway', 'reaction', 'compound', 'class')
    write_all('module', data, skeys, mkeys)
    paths = set(paths).union(names['pathway'].keys())

    # pathway
    skeys = ('name', 'class')
    mkeys = ('orthology', 'module', 'disease', 'compound')
    f = partial(parse_flat, skeys=skeys, mkeys=mkeys)
    data = batch_query(paths, 'get', f, name='pathways')
    get_classes(data)
    names = extract_targets(data, mkeys)
    skeys = ('name',)
    mkeys = ('orthology', 'module', 'disease', 'compound', 'class')
    write_all('pathway', data, skeys, mkeys)
    dses = set(dses).union(names['disease'].keys())

    # disease
    f = partial(parse_list, code='ds')
    names = batch_query(dses, 'list', f, name='diseases')
    write_names(names, 'disease_name.txt')

    print(f'Task completed at {datetime.now()}.')


if __name__ == '__main__':
    main()