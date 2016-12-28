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

#include <boost/python.hpp>
#include <iostream>
#include <numpy/ndarrayobject.h>
#include <vector>
#include <cmath>
#include "DHT.h"

using namespace boost::python;

object DoubleToNumpyArray( int nd, int s, FLOAT *data );


int peak(const FLOAT *src, const int n, FLOAT m, std::vector<FLOAT> &v)
{
	v.push_back(0.0);

	int p_sign = src[1] - src[0] > 0 ? 1 : 0 ;
	int sign;
	for(int i=1; i < n-1; i++)
	{
		sign = (src[i+1] - src[i] > 0)?1:0;
		if(sign < p_sign)
		{	
			if(src[i] > m)
			{
				v.push_back(src[i]);
			}
			else
			{
				v.push_back(0.0);
			}
		}
		else
		{
			v.push_back(0.0);
		}
		p_sign = sign;
	}
	v.push_back(0.0);

	return 0;	
}

int peak(const FLOAT *src, const int n, std::vector<FLOAT> &x, std::vector<FLOAT> &y)
{
	x.push_back(0.0);
	y.push_back(src[0]);

	int p_sign = src[1] - src[0] > 0 ? 1 : 0 ;
	int sign;
	for(int i=1; i < n-1; i++)
	{
		sign = (src[i+1] - src[i] > 0)?1:0;
		if(sign < p_sign)
		{	
			x.push_back((FLOAT)i);
			y.push_back(src[i]);
		}
		p_sign = sign;
	}
	x.push_back((FLOAT)n-1.0);
	y.push_back(src[n-1]);

	return 0;	
}

int valley(std::vector<FLOAT> index, std::vector<FLOAT> src, std::vector<FLOAT> &x, std::vector<FLOAT> &y)
{
	int n = index.size();
	x.push_back(index[0]);
	y.push_back(src[0]);

	int p_sign = src[1] - src[0] > 0 ? 1 : 0 ;
	int sign;
	for(int i=1; i < n-1; i++)
	{
		sign = (src[i+1] - src[i] > 0)?1:0;
		if(sign > p_sign)
		{					
			x.push_back(index[i]);
			y.push_back(src[i]);
		}
		p_sign = sign;
	}
	x.push_back(index[n-1]);
	y.push_back(src[n-1]);

	return 0;	
}

FLOAT* abs(int n, FLOAT* src)
{
	FLOAT* ret = (FLOAT*)malloc(sizeof(FLOAT)*n);

	for(int i=0; i< n; i++)
	{
		ret[i] = std::abs(src[i]);
	}
	return ret;
}

FLOAT snr(FLOAT* src, int n, FLOAT m)
{
	FLOAT ret = 0.0;

	for(int i=0; i< n; i++)
	{
		ret += src[i]*src[i];
	}
	return sqrt(ret)/m;
}

void quartile(FLOAT* pF, int n, FLOAT *fmax, FLOAT *fmin)
{
	std::vector<FLOAT> vec(pF, pF + n - 1);

	if (vec.size() % 4 < 2)
	{
		std::nth_element (vec.begin(), vec.begin() + vec.size()/4  , vec.end());
		*fmin = vec[vec.size()/4];
		*fmin += *std::max_element (vec.begin(), vec.begin() + vec.size()/4);

		std::nth_element (vec.begin(), vec.begin() + vec.size()/4*3 + vec.size()%2 - 1, vec.end());
		*fmax = vec[vec.size()/4*3 + vec.size()%2 - 1 ];
		*fmax += *std::min_element (vec.begin() + vec.size()/4*3 + vec.size()%2 , vec.end());

		*fmin /= 2.0;
		*fmax /= 2.0;	
	}
	else
	{
		std::nth_element (vec.begin(), vec.begin() + vec.size()/4, vec.end());
		*fmin = vec[vec.size()/4];	
		std::nth_element (vec.begin(), vec.begin() + vec.size()/4*3 + vec.size()%2 + 1, vec.end());
		*fmax = vec[vec.size()/4*3 + vec.size()%2 + 1 ];				
	}
}

FLOAT med_abs_dev(int n, FLOAT* data)
{
	FLOAT* abs_dev = (FLOAT*)malloc(sizeof(FLOAT)*n);
	FLOAT data_med = median(n, data);

	for(int i=0; i< n; i++)
	{
		abs_dev[i] = std::abs(data[i] - data_med);
	}

	return median(n, abs_dev);
}

object trigger_gen(int imfs, FLOAT* data, int n, FLOAT *amplitude, FLOAT *frequency, FLOAT m, FLOAT snr_th)
{
	FLOAT fac_med = 1.4;
	FLOAT* data_abs = abs(n, data);
	FLOAT data_med = med_abs_dev(n, data);
	//FLOAT data_med = median(n, data_abs);
	//fac_med *= (1.0 + data_med / m);
	data_med *= fac_med;
	m *= fac_med;

	std::vector<FLOAT> am_peak;
	am_peak.reserve(n);
	std::vector<FLOAT> tmp_index;
	std::vector<FLOAT> tmp_peak;

	std::vector<FLOAT> d_index;
	std::vector<FLOAT> d_valley;
	
	peak(data_abs, n, tmp_index, tmp_peak);
	valley(tmp_index,tmp_peak, d_index, d_valley);
	//peak(amplitude, n, std::min(m/2.0, data_med), am_peak);
	peak(amplitude, n, data_med, am_peak);

	int info_size = 9;
	std::vector<int> trg_index;
	std::vector<int> index;
	std::vector<FLOAT> trg_snr;
	std::vector<FLOAT>::iterator start = am_peak.begin();
	for(int i=0; i < (int)d_index.size() - 1; i++)
	{
		std::vector<FLOAT>::iterator max_am = std::max_element(start + d_index[i], start + d_index[i+1]);
		FLOAT _snr = snr(&data[(int)d_index[i]], d_index[i+1]-d_index[i], m);
		
		if((max_am[0] > 0) && (_snr >= snr_th))
		{	
			trg_index.push_back(i);
			index.push_back(std::distance(start, max_am));
			trg_snr.push_back(_snr);
		}
	}
	if (index.size() == 0) return DoubleToNumpyArray(0, info_size, (FLOAT*)0);
	FLOAT *trigger_info = (FLOAT*)malloc(sizeof(FLOAT)*index.size()*info_size); 	
#pragma omp parallel for
	for(int j=0; j < (int)index.size(); j++)
	{
		int i = trg_index[j];
		FLOAT *pF = frequency + (int)d_index[i];
		int nF = d_index[i+1] - d_index[i];		
		FLOAT fmin, fmax;
		quartile(pF, nF, &fmax, &fmin);
		trigger_info[info_size*j + 0] = imfs;
		trigger_info[info_size*j + 1] = d_index[i];
		trigger_info[info_size*j + 2] = d_index[i+1];
		trigger_info[info_size*j + 3] = index[j];
		trigger_info[info_size*j + 4] = amplitude[index[j]];
		trigger_info[info_size*j + 5] = frequency[std::max(0,index[j]-1)];
		trigger_info[info_size*j + 6] = fmin;		
		trigger_info[info_size*j + 7] = fmax;
		trigger_info[info_size*j + 8] = trg_snr[j];
	}
	free(data_abs);
	return DoubleToNumpyArray(index.size(), info_size, trigger_info);
}
