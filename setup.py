# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import ast
from setuptools import find_packages, setup

classes = """
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

classifiers = [s.strip() for s in classes.split('\n') if s]

description = (
    "metagenomix is a pipeline of pipelines to conduct metagenomic analyses "
    "on Slurm/Torque"
)

with open("README.md") as f:
    long_description = f.read()

_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("metagenomix/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

standalone = ['metagenomix=metagenomix.scripts._standalone_metagenomix:standalone_metagenomix']

setup(
    name="metagenomix",
    version=version,
    license="BSD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Franck Lejzerowicz",
    author_email="franck.lejzerowicz@gmail.com",
    maintainer="Franck Lejzerowicz",
    maintainer_email="franck.lejzerowicz@gmail.com",
    url="https://github.com/FranckLejzerowicz/metagenomix",
    packages=find_packages(),
    install_requires=[
        "click",
        "biom-format",
        "scikit-bio",
        "scipy",
        "numpy",
        "pandas",
        "pyyaml",
        "seaborn",
        "Xhpc==2.21"
    ],
    classifiers=classifiers,
    entry_points={'console_scripts': standalone},
    package_data={
        'metagenomix': [
            'test/*/*/*',
            'resources/run_params.yml',
            'resources/wol_tree.nwk'
        ],
    },
    include_package_data=True,
    python_requires='>=3.6',
)
