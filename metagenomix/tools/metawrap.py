# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from os.path import dirname, isfile


def metawrap_reassemble():
    pass


def get_sams_fastqs(mode: str, fastqs: dict):
    """Get the paths of the fastq files per sample or altogether.

    Parameters
    ----------
    mode : str
        Mode of processing for blobology:
        - "coassembly": all co-assembled samples.
        - "sample": each sample independently.
    fastqs : dict
        Path to the raw fastq files per sample.

    Returns
    -------
    sams_fastqs : dict
        Either same as `fastqs` of all fastqs under an empty key.
    """
    if mode == 'sample':
        sams_fastqs = fastqs
    elif mode == 'coassembly':
        sams_fastqs = {'': [y for x in fastqs.values() for y in x]}
    else:
        raise ValueError('No "blobology" param as "%s"' % mode)
    return sams_fastqs


def metawrap_quantify():
    pass

def metawrap_classify():
    pass

def metawrap_annotate():
    pass
