#!/usr/bin/env python

from distutils.core import setup, Extension

from distutils.command.build_ext import build_ext as build_ext_orig
import os
import pathlib


class build_ext(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = os.getcwd()
        cpp_dir = cwd + '/cpp'

        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        # example of cmake args
        config = 'Debug' if self.debug else 'Release'
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(extdir.parent.absolute()),
            '-DCMAKE_BUILD_TYPE=' + config
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', cpp_dir, ] + cmake_args)
        if not self.dry_run:
            self.spawn(['make', '-j4'])
        os.chdir(str(cwd))


setup(name='Distutils',
      version='v0.1',
      description='Library for the Hermite spectral method',
      author='Urbain Vaes',
      author_email='urbain@vaes.uk',
      url='https://github.com/urbainvaes/hermite',
      packages=['hermite'],
      ext_modules=[Extension('hermite_cpp', sources=[])],
      cmdclass={'build_ext': build_ext}
      )
