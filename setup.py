# Copyright (C) 2016 Whansun John Kim and Edwin J. Son
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

"""
Etagen - Trigger generator using Hilbert-Huang Transform implement by C++
This is Python wrapper of Etagen library.
"""

from distutils.core import setup, Extension
import numpy as np
import sys

if sys.platform == 'darwin':
	extra_options = dict(
		extra_compile_args = ['-O3', '-I/opt/local/include', '-I/opt/local/include/libomp', '-fopenmp', '-std=c++0x'],
		extra_link_args=['-L/opt/local/lib', '-L/opt/local/lib/libomp', "-fopenmp"],
		libraries=['boost_python-mt']
	)
elif sys.platform == 'linux2':
	extra_options = dict(
		extra_compile_args = ['-O3', '-fopenmp', '-I/usr/include/c++', '-std=c++0x'],
		extra_link_args=["-fopenmp"],
		libraries=['boost_python']
	)
else:
	extra_options = dict(
		extra_compile_args = ['-O3', '-fopenmp', '-std=c++11'],
		extra_link_args=["-fopenmp"],
		libraries=['boost_python']
	)

etagen_mod = Extension('etagen._etagen',
        sources = [
	'emd.cpp',
        'cubic.cpp',
        'DHT.cpp',
        'trigger.cpp',
        'cluster.cpp',
        'etagen.cpp',
	],
	**extra_options)
setup(name = "etagen",
        version = "0.1",
        description = "python wapper of Etagen library",
        author = "Whansun John Kim and Edwin J. Son",
        author_email = "hwansun.kim@gmail.com",
	packages=['etagen'],
	package_dir={'etagen': 'python'},
        ext_modules = [etagen_mod],
	include_dirs = np.get_include()
        )
