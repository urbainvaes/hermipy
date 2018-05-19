#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext as build_ext_orig
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


setup(name='Hermite',
      version='v0.1',
      description='Library for the Hermite spectral method',
      author='Urbain Vaes',
      author_email='urbain@vaes.uk',
      url='https://github.com/urbainvaes/hermite',
      packages=['hermite'],
      ext_modules=[Extension('hermite_cpp', sources=[])],
      cmdclass={'build_ext': build_ext}
      )
