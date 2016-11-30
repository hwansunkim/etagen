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

#ifndef __DHT_H__
#define __DHT_H__

#include "common.h"

void fir_filter(int n, FLOAT* w);
int dht(int n, FLOAT *src, FLOAT *des, FLOAT *w, int filter_length);
int hsa(int n, FLOAT *src, FLOAT *imag, FLOAT *amplitude, FLOAT *frequency, FLOAT fs, bool last);
FLOAT median(int n, FLOAT* src);

#endif
