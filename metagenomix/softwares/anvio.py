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


def anvio(self) -> None:
    """Anvi’o is a comprehensive platform that brings together many aspects
    of today’s cutting-edge computational strategies of data-enabled
    microbiology, including genomics, metagenomics, metatranscriptomics,
    pangenomics, metapangenomics, phylogenomics, and microbial population
    genetics in an integrated and easy-to-use fashion through extensive
    interactive visualization capabilities.

    References
    ----------
    Eren, A.M., Kiefl, E., Shaiber, A., Veseli, I., Miller, S.E., Schechter,
    M.S., Fink, I., Pan, J.N., Yousef, M., Fogarty, E.C. and Trigodet, F.,
    2021. Community-led, integrated, reproducible multi-omics with anvi’o.
    Nature microbiology, 6(1), pp.3-6.

    Notes
    -----
    GitHub  : https://github.com/merenlab/anvio
    Paper   : https://doi.org/10.1038/s41564-020-00834-3
    Docs    : https://anvio.org/

    Parameters
    ----------
    self
    """
    pass