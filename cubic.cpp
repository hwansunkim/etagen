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

#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <math.h>

#include "cubic.hpp"
#define eps 0

FLOAT h00(FLOAT t)
{
    return 2.0*t*t*t - 3.0*t*t +1;
}
FLOAT h10(FLOAT t)
{
    return t*(1.0-t)*(1.0-t);
}
FLOAT h01(FLOAT t)
{
    return t*t*(3.0-2.0*t);
}
FLOAT h11(FLOAT t)
{
    return t*t*(t-1.0);
}

bool Linear_Interpolation(const std::vector<FLOAT> *x_src, const std::vector<FLOAT> *y_src, std::vector<FLOAT>* destX, std::vector<FLOAT>* destY)
{
    FLOAT x1 = (*x_src)[0];   
    FLOAT x2 = (*x_src)[1];   
    FLOAT y1 = (*y_src)[0];   
    FLOAT y2 = (*y_src)[1];  

    destX->push_back((*x_src)[0]);
    destY->push_back((*y_src)[0]);

    for(int i = (*x_src)[0] + 1; i < (*x_src)[1]; i++)
    {
        destX->push_back((FLOAT)i);
        destY->push_back((y1-y2)/(x1-x2)*(i -x1) + y1);
    }

    destX->push_back((*x_src)[1]);
    destY->push_back((*y_src)[1]);

    return true;
}

bool monotonic_cubic_Hermite_spline(std::vector<FLOAT>* x_src, std::vector<FLOAT>* y_src, std::vector<FLOAT>* destX, std::vector<FLOAT>* destY)
{

	if((int)x_src->size() == 2)
	{
		Linear_Interpolation(x_src, y_src, destX, destY);
		return true;
	}
	// 0-based index 사용.
	int n = (int)x_src->size();
	int k = 0;
	FLOAT *m = new FLOAT[n];
	FLOAT ak, bk, akbk;
	FLOAT *delta_k = new FLOAT[n-1];
	m[0] = ((*y_src)[1] - (*y_src)[0])/((*x_src)[1] - (*x_src)[0]);
	m[n-1] = ((*y_src)[n-1] - (*y_src)[n-2])/((*x_src)[n-1]-(*x_src)[n-2]);

	for(k=0; k < n-1; k++)
	{
		delta_k[k] = ((*y_src)[k+1] - (*y_src)[k])/((*x_src)[k+1] - (*x_src)[k]);
        // if(isnan(delta_k[k]) || delta_k[k] == 0) printf("%lf %lf %lf %lf\n", (*x_src)[k+1], (*x_src)[k], (*y_src)[k], (*y_src)[k+1] );
	}
	for(k = 1; k<n-1; k++)
	{
		m[k] = (delta_k[k-1] + delta_k[k])/2;
	}
	m[0] = delta_k[0];
	m[n-1] = delta_k[n-2];

	for(k = 0; k<n-1; k++)
	{
		if(fabs(delta_k[k]) == 0.0)
		{
			m[k] = m[k+1] = 0.0;
		}
		else
		{
			ak = m[k]/delta_k[k];
			bk = m[k+1]/delta_k[k];
			akbk = ak*ak + bk*bk;
			if(akbk > 9)
			{
				m[k] = 3/(sqrt(akbk)) * ak * delta_k[k];
				m[k+1] = 3/(sqrt(akbk)) * bk * delta_k[k];
			}
            if(isnan(m[k]) || isnan(m[k+1])) 
            {
                m[k] = m[k+1] = 0.0;
            }
		}

	}

	FLOAT cur_x = 0.0;
	FLOAT next_x = 0.0;
	FLOAT cur_y = 0.0;
	FLOAT next_y = 0.0;
	FLOAT h = 0.0;
	FLOAT x = 0.0;
	FLOAT t = 0.0, y = 0.0;

	for(k = 0; k<n-1; k++)
	{
		cur_x = (FLOAT)((int)(0.5 + (*x_src)[k]));
		next_x = (FLOAT)((int)((*x_src)[k+1]));
		cur_y = (*y_src)[k];
		next_y = (*y_src)[k+1];
		h = next_x - cur_x;

		for(x = cur_x; x<next_x; x+=1.0)
		{
			t = (x-cur_x)/h;
			destX->push_back(x);

			y = cur_y*h00(t) + h*m[k]*h10(t) + next_y*h01(t) + h*m[k+1]*h11(t);
            // if(isnan(y)) printf("t %lf %lf %lf %lf \n", t, h, m[k], m[k+1]);
			destY->push_back(y);
		}
	}
	t = (x-cur_x)/h;
	destX->push_back(x);
	y = cur_y*h00(t) + h*m[k]*h10(t) + next_y*h01(t) + h*m[k+1]*h11(t);
	destY->push_back(y);

	delete[] m;
	delete[] delta_k;

	return true;
}

bool cubic_spline(std::vector<FLOAT>* x_series, std::vector<FLOAT>* y_series, std::vector<FLOAT> *destX, std::vector<FLOAT>* destY)
{   
    if((int)x_series->size() == 2)
    {
        Linear_Interpolation(x_series, y_series, destX, destY);
        return true;
    }
    int n = (int)x_series->size()-1;
    // Step 1.
    FLOAT *h = new FLOAT[n];
    FLOAT *alpha = new FLOAT[n];
    
    /* h size is n - 1 */
    for(int i = 0; i < n; i++){
        h[i] = (*x_series)[i+1] - (*x_series)[i];
    }
 
    // Step 2.
    for(int i = 1; i < n;i++){
        alpha[i]= 3*((*y_series)[i+1]-(*y_series)[i])/h[i]-3*((*y_series)[i]-(*y_series)[i-1])/h[i-1];
    }
 
    // Step 3.
    FLOAT *l = new FLOAT[n+1];
    FLOAT *u = new FLOAT[n];
    FLOAT *z = new FLOAT[n+1];
    FLOAT *c = new FLOAT[n+1];
    FLOAT *b = new FLOAT[n];
    FLOAT *d = new FLOAT[n];
 
    l[0] = 1; u[0] = 0; z[0] = 0;
 
    // Step 4.
    for(int i = 1; i < n; i++){
        l[i] = 2*((*x_series)[i+1] - (*x_series)[i-1]) - h[i-1]*u[i-1];
        u[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i];
    }
 
    // Step 5.
    l[n] = 1;     z[n] = 0;     c[n] = 0;
 
    // Step 6.
    for(int i = n-1; i>=0; i--){
        c[i] = z[i] - u[i]*c[i+1];
        b[i] = ((*y_series)[i+1] - (*y_series)[i])/h[i] - h[i]*(c[i+1] + 2*c[i])/3;
        d[i] = (c[i+1] - c[i]) / (3*h[i]);
    }

	FLOAT x, x_offset, Sx;
    int t;
    for(t = 0; t < n; t++)
    {
        x = (*x_series)[t];
        for(; x < (*x_series)[t+1]; x+=1.0)
        {
            x_offset = x - (*x_series)[t];
            Sx = (*y_series)[t] + b[t]*x_offset + c[t]*x_offset*x_offset + d[t]*x_offset*x_offset*x_offset;
            
            destX->push_back(x);
            destY->push_back(Sx);

        }
    }
    /*
    x_offset = x - (*x_series)[t];
    Sx = (*y_series)[t] + b[t]*x_offset + c[t]*x_offset*x_offset + d[t]*x_offset*x_offset*x_offset;

    destX->push_back(x);
    destY->push_back(Sx);
    */
    if ((FLOAT)destX->size() < (*x_series)[n] - (*x_series)[0])
    {
        return false;
    }
    else
    {
        destX->push_back((*x_series)[n]);
        destY->push_back((*y_series)[n]);
    }
		
    delete [] h;
    delete [] alpha;
    delete [] l;
    delete [] u;
    delete [] z;
    delete [] c;
    delete [] b;
    delete [] d;
 
    return true;
}
