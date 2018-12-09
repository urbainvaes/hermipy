#!/usr/bin/env python
#
# Copyright (C) 2018 Urbain Vaes
#
# This file is part of hermipy, a python/C++ library for automating the
# Hermite Galerkin method.
#
# hermipy is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hermipy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup
from distutils.core import Extension
from distutils.command.build_ext import build_ext as build_ext_orig

import unittest
import os


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        cwd = os.getcwd()
        os.makedirs(self.build_temp, exist_ok=True)
        os.chdir(self.build_temp)
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY='
                      + cwd + '/' + self.build_lib]
        self.spawn(['cmake', cwd + '/cpp'] + cmake_args)
        if not self.dry_run:
            self.spawn(['make', '-j4'])
        os.chdir(str(cwd))


setup(name='Hermipy',
      version='v0.1',
      description='Library for the Hermite spectral method',
      author='Urbain Vaes',
      author_email='urbain@vaes.uk',
      url='https://github.com/urbainvaes/hermite',
      packages=['hermipy'],
      ext_modules=[Extension('hermite_cpp', sources=[])],
      cmdclass={'build_ext': build_ext},
      test_suite='tests',
      )
