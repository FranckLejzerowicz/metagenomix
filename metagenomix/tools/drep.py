# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from metagenomix._io_utils import mkdr


def get_drep_bins(prev: str, pools: dict, inputs: dict) -> dict:
    """Get the genomes in an iterable format accommodating
    both the output of `metawrap_refine()` and of `metawrap_reassemble()`.

    Parameters
    ----------
    prev : str
        Name of the software preceding drep.
    pools : dict
        Samples per group as dispatched according to the pooling design.
    inputs : dict
        Outputs from the software preceding drep.

    Returns
    -------
    bins : dict
        Bins to dereplicate.
    """
    bins = {}
    for pool in pools:
        bins[pool] = {}
        for group, sam_paths in inputs[pool].items():
            if prev == 'metawrap_refine':
                bins[pool].setdefault('', []).append((group, sam_paths[1]))
            elif prev == 'metawrap_reassemble':
                for sam, paths in sam_paths.items():
                    for path in paths:
                        if path.endswith('strict'):
                            bins[pool].setdefault('strict', []).append(path)
                        elif path.endswith('permissive'):
                            bins[pool].setdefault('permissive', []).append(path)
            else:
                raise ValueError('No genome dereplication after %s' % prev)
    return bins


def get_drep_inputs(drep_dir: str, sam_paths: list):
    """Write the file containing the inputs bins to drep
    and the list of paths to these bins.

    Parameters
    ----------
    drep_dir :
        Path to the drep output folder.
    sam_paths : list
        List containing for each (group, bins_folder).

    Returns
    -------
    drep_in : str
        File containing the paths corresponding to each bin.
    paths : list
        List of paths corresponding to each bin.
    """
    paths = []
    mkdr(drep_dir)
    drep_in = '%s/input_genomes.txt' % drep_dir
    with open(drep_in, 'w') as o:
        for sam_path in sam_paths:
            for path in glob.glob('%s/*fa' % sam_path):
                paths.append(path)
                o.write('%s\n' % path)
    return drep_in, paths
