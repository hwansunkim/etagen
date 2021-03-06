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

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <algorithm>

#include "common.h"

const FLOAT PI = 3.141592653589793;

/*
+n must be even.
-64 <= n <= 64
https://en.wikipedia.org/wiki/Hilbert_transform
*/
    
void fir_filter(int n, FLOAT* w)
{
	FLOAT l = -n/2.0 + 1.0;
	for(int i=0; i < n; i++)
	{
		w[i] = -(2.0/(2.0*l-1.0)/PI);
		l += 1.0;
	}
}

FLOAT convol(int n, int offset, FLOAT *src, FLOAT *w)
{
	FLOAT ret = 0.0;
	for(int i=0; i < n; i++)
	{
		ret += src[2 * i + offset] * w[i]; 
	}	
	return ret;
}

int dht(int n, FLOAT *src, FLOAT *des, FLOAT *w, int filter_length)
{
#pragma omp parallel for
	for(int i=0; i < filter_length; i++)
	{
		int length = filter_length/2 + (i+1)/2;
		FLOAT tmp = convol(length, (i+1) % 2, src, w + filter_length/2 - (i+1)/2); 
		des[i] = tmp;
	}
#pragma omp parallel for
	for(int i = 0; i < n - filter_length * 2; i++)
	{
		FLOAT tmp = convol(filter_length, 1, src + i, w);
		des[i + filter_length] = tmp; 
	}
#pragma omp parallel for
	for(int i = 0 ; i < filter_length; i++)
	{
		int length = filter_length - (i+1)/2;
		FLOAT tmp = convol(length, 1, src + n - filter_length * 2 + i, w);
		des[i + n - filter_length] = tmp; 
	}

	return 0;
}

FLOAT median(int n, FLOAT* src)
{	
	std::vector<FLOAT> vec (src, src+n);
	if (vec.size() % 2 == 0)
    {
        std::nth_element (vec.begin(), vec.begin() + vec.size()/2 - 1, vec.end());
        double m = vec[vec.size()/2 - 1];
        m += *std::min_element (vec.begin() + vec.size()/2, vec.end());
        return m / 2;
    }
    else
    {
        std::nth_element (vec.begin(), vec.begin() + vec.size()/2, vec.end());
        return (vec[vec.size()/2]);
    }
}

int hsa(int n, FLOAT *src, FLOAT * imag, FLOAT *amplitude, FLOAT *frequency, FLOAT fs, bool last)
{
	const FLOAT pi2 = PI *2.0, pi4 = PI *4.0, pi24 = PI *24.0;
/*
	int Nf = last ? n-1 : n;
	FLOAT *fi = (FLOAT*)malloc(sizeof(FLOAT)*(Nf+1));
	if (!last) fi[Nf] = std::atan2(imag[Nf], src[Nf]);
*/
	FLOAT *fi = (FLOAT*)malloc(sizeof(FLOAT)*(n));
	FLOAT dfi;
/* 
http://en.cppreference.com/w/cpp/numeric/math/atan2 
-PI <= atan2 <= PI
*/

#pragma omp parallel for
	for(int i=0; i < n; i++)
	{
		amplitude[i] = sqrt(src[i] * src[i] + imag[i] * imag[i]);
		fi[i] = std::atan2(imag[i], src[i]);
	}

	// unwrap routine (should be done sequentially) - added by Eddy
	for(int i=1; i < n; i++)
	{
		dfi = std::fmod(fi[i] - fi[i-1] + PI, pi2);
		if (dfi < 0.0) dfi += pi2;
		fi[i] = fi[i-1] + dfi - PI;
	}

#pragma omp parallel for
	//for(int i=0; i< n-1; i++)
	//for(int i=0; i< Nf; i++)
	for(int i=0; i< n; i++)
	{	
        if(i==0) frequency[i] = (fi[i+1] - fi[i])*fs/pi2;
        else if(i==n-1) frequency[i] = (fi[i] - fi[i-1])*fs/pi2;
        else if(i==1 || i==n-2) frequency[i] = (fi[i+1] - fi[i-1])*fs/pi4;
        else frequency[i] = (-fi[i+2] + 8*fi[i+1] - 8*fi[i-1] + fi[i-2])*fs/pi24;
/*
		frequency[i] = (fi[i+1] - fi[i])*fs/pi2;
		if(frequency[i] > fs/2) frequency[i] -= fs/2;
		if(frequency[i] < -fs/2) frequency[i] += fs/2;
		if(frequency[i] < 0) frequency[i] += fs/2;
*/
	}
	free(fi);
/*
	FLOAT m = median(Nf, frequency);
	if( m <= fs/4.0)
	{
		for(int i=0; i< Nf; i++)
		{	
			if( frequency[i] > std::max(m*3, fs/4.0)) frequency[i] -= fs/2;
		}
	}
*/
	return 0;
}
