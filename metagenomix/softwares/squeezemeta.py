# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from os.path import dirname
from metagenomix._inputs import (sample_inputs, group_inputs,
                                 genome_key, genome_out_dir)
from metagenomix._io_utils import io_update, to_do, status_update


def squeezemeta(self) -> None:
    """SqueezeMeta is a full automatic pipeline for
    metagenomics/metatranscriptomics, covering all steps of the analysis.
    SqueezeMeta includes multi-metagenome support allowing the co-assembly
    of related metagenomes and the retrieval of individual genomes via
    binning procedures. Thus, SqueezeMeta features several unique
    characteristics.

    References
    ----------
    Tamames, J. and Puente-Sánchez, F., 2019. SqueezeMeta, a highly portable,
    fully automatic metagenomic analysis pipeline. Frontiers in microbiology,
    9, p.3349.

    Notes
    -----
    GitHub  : https://github.com/jtamames/SqueezeMeta
    Paper   : https://doi.org/10.3389/fmicb.2018.03349
    Docs    : https://github.com/jtamames/SqueezeMeta/wiki

    Parameters
    ----------
    self
    """
    pass
