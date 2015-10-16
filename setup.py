#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

major = 0
minor = 1
maintenance = 0

setup(
    name = "dolfinh5tools",
    version = "{:d}.{:d}.{:d}".format(major, minor, maintenance),
    description = "Tools that allow saving multiple timesteps of \
                    a field (i.e. dolfin function) to an hdf5 file, \
                    and to retrieve it again, using only the dolfin \
                    hdf5 commands.",
    author = "Marijan Beg",
    url = "https://github.com/fangohr/dolfinh5tools",
    packages = ["dolfinh5tools"]
)
