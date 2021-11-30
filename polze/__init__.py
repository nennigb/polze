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
# along with polze.  If not, see <https://www.gnu.org/licenses/>.

r"""
Polze  -- A package to locate poles and zeros of a meromorphic function with their multiplicities
==================================================================================================


.. include::../README.md
    :start-line:2
    :raw:
"""
from sys import modules as _modules
from ._polze import PZ, logger
from .version import __version__
# Append _polze.py docstring to top level doc.
from ._polze import __doc__ as _tutos

_modules[__name__].__doc__ += _tutos

# Add doc to logger object (pdoc hack)
logger = logger
"""Standard `logger` object from the standard logging library. It controls how
`polze` print on standard output. Verbosity can be limitted by setting
`polze.logger.setLevel(50)`. The default value is 20. For instance in your
scripts use,
```
import logging
logging.basicConfig(level=logging.ERROR)
logging.getLogger('polze').setLevel(logging.DEBUG)
```
"""


__all__ = ['PZ', '__version__', 'logger']
