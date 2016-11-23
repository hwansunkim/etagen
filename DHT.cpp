#include <vector>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <algorithm>
#define FLOAT double

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
	const FLOAT pi2 = PI *2.0;
	int Nf = last ? n-1 : n;
	FLOAT *fi = (FLOAT*)malloc(sizeof(FLOAT)*(Nf+1));
	if (!last) fi[Nf] = std::atan2(imag[Nf], src[Nf]);
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

#pragma omp parallel for
	//for(int i=0; i< n-1; i++)
	for(int i=0; i< Nf; i++)
	{	
		frequency[i] = (fi[i+1] - fi[i])*fs/pi2;
		if(frequency[i] > fs/2) frequency[i] -= fs/2;
		if(frequency[i] < -fs/2) frequency[i] += fs/2;
		if(frequency[i] < 0) frequency[i] += fs/2;
	}
	free(fi);
	FLOAT m = median(Nf, frequency);
	if( m <= fs/4.0)
	{
		for(int i=0; i< Nf; i++)
		{	
			if( frequency[i] > std::max(m*3, fs/4.0)) frequency[i] -= fs/2;
		}
	}
	return 0;
}
