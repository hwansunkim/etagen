/* Copyright (C) 2016  Whansun Kim
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

#ifndef __CUBIC_H__
#define __CUBIC_H__

#include "common.h"

bool cubic_spline(std::vector<FLOAT>* x_series, std::vector<FLOAT>* y_series, std::vector<FLOAT> *destX, std::vector<FLOAT>* destY); 
bool monotonic_cubic_Hermite_spline(std::vector<FLOAT>*, std::vector<FLOAT>*, std::vector<FLOAT>*, std::vector<FLOAT>*);

#endif
