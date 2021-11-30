# -*- coding: utf-8 -*-

# This file is part of polze, a package to locate poles and zeros of a
# meromorphic function with their multiplicities.
# Copyright (C) 2019  Benoit Nennig, benoit.nennig@isae-supmeca.fr

# polze is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# polze is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with polze. If not, see <https://www.gnu.org/licenses/>.

import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()


# Load version
def version():
    """ Get version from version.py."""
    v = None
    with open(os.path.join('./polze', 'version.py')) as f:
        for line in f:
            if line.lstrip().startswith('__version__'):
                v = line.split('=')[-1].strip().replace("'", "").replace('"', "")
                break
        return v


setuptools.setup(
    name="polze",
    version=version(),
    author="B. Nennig",
    author_email="benoit.nennig@isae-supmeca.fr",
    description="A package to locate poles and zeros of a meromorphic function with their multiplicities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nennigb/polze",
    packages=['polze'],
    install_requires=['numpy',
                      'scipy',
                      'matplotlib'],  # max version for python3.5 : matplotlib<=3.0.0
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics"
    ],
    python_requires='>=3.5',  # Tested from python > 3.5
)
