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

#ifndef __EMD_H__
#define __EMD_H__

#include "common.h"

/* Weighted Side EMD Weight Function Type */
enum weightTpye
{
	exp_kernel,
	sin_kernel, 
	spline_kernel
};

int sift(const std::vector<FLOAT> *src, std::vector<FLOAT> *imf);
int stopcondition(const std::vector<FLOAT> *src);
void emd(int, const FLOAT*, FLOAT*, int&, const int, const int, const int);
int wsemd(int, int, int, FLOAT *const, FLOAT *const, FLOAT*, int&, const int, const int, const int);
int extrema(const FLOAT *src, const int n, std::vector<FLOAT> &x_pick, std::vector<FLOAT> &y_pick, std::vector<FLOAT> &x_bottom, std::vector<FLOAT> &y_bottom);
void weightfunction(FLOAT *w, int n, int type, FLOAT alpha);

#endif
