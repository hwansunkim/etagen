/* Copyright (C) 2016  Whansun Kim and Edwin J. Son
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef __TRIGGER_H__
#define __TRIGGER_H__

#include <boost/python.hpp>
using namespace boost::python;

object trigger_gen(int imfs, FLOAT* data, int n, FLOAT *amplitude, FLOAT *frequency, FLOAT m, FLOAT snr_th);
FLOAT* abs(int n, FLOAT* src);
FLOAT med_abs_dev(int n, FLOAT* data);
#endif
