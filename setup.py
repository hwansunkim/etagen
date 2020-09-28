#!/usr/bin/env python

# Copyright (C) 2016 Whansun Kim and Edwin J. Son
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from setuptools import Extension
from numpy.distutils.core import setup
from numpy import get_include
from glob import glob
from sys import platform, version_info
from ctypes.util import find_library

def find_boost():
  if platform == 'darwin':
    if version_info[0] == 3:
      return "boost_python3-mt"
    else:
      return "boost_python-mt"
  else:
    boost_python = 'boost_python-py{}{}'.format(*version_info[:2])
    bpl = find_library(boost_python)
    if not bpl:
      boost_python = 'boost_python{}{}'.format(*version_info[:2])
      bpl = find_library(boost_python)
      if not bpl:
        print(boost_python)
        boost_python = 'boost_python'
        print(boost_python)
      else:
        print(bpl, 'found')
    else:
      print(bpl, 'found')
    return boost_python

def find_numpy():
  if version_info[0] == 3:
    boost_numpy = 'boost_numpy3'
  else:
    boost_numpy = 'boost_numpy'

  if platform == 'darwin':
    boost_numpy += '-mt'
  else:
    boost_numpy2 = boost_numpy + '-py{}{}'.format(*version_info[:2])
    boost_numpy3 = 'boost_numpy{}{}'.format(*version_info[:2])
    bnl2 = find_library(boost_numpy2)
    bnl3 = find_library(boost_numpy3)
    if bnl2:
      boost_numpy = boost_numpy2
      print(bnl2, 'found')
    elif bnl3:
      boost_numpy = boost_numpy3
      print(bnl3, 'found')

  return boost_numpy

def main():
  extra_options = dict(
    sources = [
      'emd.cpp',
      'cubic.cpp',
      'DHT.cpp',
      'trigger.cpp',
      'cluster.cpp',
#      'etagen.cpp',
    ],
    libraries = [
      find_boost(),
#      find_numpy()
    ],
    include_dirs = [get_include(), '/opt/local/include'],
  )

  boost_numpy = find_numpy()
  if not find_library(boost_numpy):
    print('Not found: Boost.Numpy')
    extra_options['sources'].append('etagen.numeric.cpp')
  else:
    extra_options['libraries'].append(boost_numpy)
    extra_options['sources'].append('etagen.numpy.cpp')

  extensions = [
    Extension(
      name = 'etagen._etagen',
#      sources = glob('src/*.c*'),
      extra_compile_args = ['-O3', '-fopenmp', '-std=c++11'],
      extra_link_args=["-fopenmp"],
#      libraries = [find_boost(), find_numpy()],
#      include_dirs = [get_include(), '/opt/local/include'],
      library_dirs = ['/opt/local/lib'],
      **extra_options
    ),
  ]

  setup(
    name = "etagen",
    version = "0.2",
    description = "python wapper of Etagen library",
    author = "Whansun Kim and Edwin J. Son",
    author_email = "hwansun.kim@gmail.com, eddy@nims.re.kr",
    packages = ['etagen'],
    package_dir = {'etagen': 'python'},
    ext_modules = extensions,
#    include_dirs = np.get_include()
    install_requires = ['numpy']
    )

if __name__ == '__main__':
  main()
