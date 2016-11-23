#include <iostream>
#include <vector>
#include <cstring>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include "cubic.hpp"
#include "emd.hpp"

#define mono1
#define DEBUG1

#define sgn(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

int extrema(const FLOAT *src, const int n, std::vector<FLOAT> &x_peak, std::vector<FLOAT> &y_peak, std::vector<FLOAT> &x_bottom, std::vector<FLOAT> &y_bottom)
{
	x_peak.push_back(0.0);
	y_peak.push_back(src[0]);
	x_bottom.push_back(0.0);
	y_bottom.push_back(src[0]);

	//int p_sign = (src[1] - src[0] > 0) ? 1 : 0 ;
	FLOAT diff = src[1] - src[0];
	int p_sign = sgn(diff);
	int p_x = 0;
	int sign;
	for(int i=1; i < n-1; i++)
	{
		//sign = (src[i+1] - src[i] > 0)?1:0;
		diff = src[i+1] - src[i];
		sign = sgn(diff);
		if(sign == 0)
			continue;
		if(sign < p_sign)
		{
			//x_peak.push_back((FLOAT)i);
			x_peak.push_back(round((FLOAT)(i+1+p_x)/2.0));
			y_peak.push_back(src[i]);
		}
		else if(sign > p_sign)
		{
			//x_bottom.push_back((FLOAT)i);
			x_bottom.push_back(round((FLOAT)(i+1+p_x)/2.0));
			y_bottom.push_back(src[i]);
		}
		p_sign = sign;
		p_x = i;
	}
	x_peak.push_back((FLOAT)n-1.0);
	y_peak.push_back(src[n-1]);
	x_bottom.push_back((FLOAT)n-1.0);
	y_bottom.push_back(src[n-1]);


	int ret = fmin(x_peak.size(), x_bottom.size());
	/*
	if (ret < 3)
		return ret;
	*/

	if (x_peak.size() > 4)
	{
		FLOAT x1 = x_peak[1];
		FLOAT y1 = y_peak[1];
		FLOAT x2 = x_peak[2];
		FLOAT y2 = y_peak[2];
		FLOAT t = x_peak[0];
		y_peak[0] = fmax(y_peak[0], (y1-y2)/(x1-x2)*(t -x1) + y1);
	
		int index = x_peak.size()-1;

		x1 = x_peak[index - 2];
		y1 = y_peak[index - 2];
		x2 = x_peak[index - 1];
		y2 = y_peak[index - 1];
		t = x_peak[index];
		y_peak[index] = fmax(y_peak[index], (y1-y2)/(x1-x2)*(t -x1) + y1);
	}

	if (x_bottom.size() > 4)
	{
		FLOAT x1 = x_bottom[1];
		FLOAT y1 = y_bottom[1];
		FLOAT x2 = x_bottom[2];
		FLOAT y2 = y_bottom[2];
		FLOAT t = x_bottom[0];

		y_bottom[0] = fmin(y_bottom[0], (y1-y2)/(x1-x2)*(t -x1) + y1);

	
		int index = x_bottom.size()-1;

		x1 = x_bottom[index - 2];
		y1 = y_bottom[index - 2];
		x2 = x_bottom[index - 1];
		y2 = y_bottom[index - 1];
		t = x_bottom[index];
		y_bottom[index] = fmin(y_bottom[index], (y1-y2)/(x1-x2)*(t -x1) + y1);
	}

	return ret;	
}

int sift(const int N, const FLOAT *src, std::vector<FLOAT> *imf, const int mono)
{
	std::vector<FLOAT> x_series;
	std::vector<FLOAT> y_series;
	std::vector<FLOAT> destX_peak;
	std::vector<FLOAT> destY_peak;
	std::vector<FLOAT> destX_bottom;
	std::vector<FLOAT> destY_bottom;
	
	std::vector<FLOAT> x_peak;
	std::vector<FLOAT> y_peak;

	std::vector<FLOAT> x_bottom;
	std::vector<FLOAT> y_bottom;

	std::vector<FLOAT> destY_middle;

	//peak & bottom search
	// if (extrema(src, N, x_peak, y_peak, x_bottom, y_bottom) < 3) return -1;
	extrema(src, N, x_peak, y_peak, x_bottom, y_bottom);

	if(mono == 1)
	{
		monotonic_cubic_Hermite_spline(&x_bottom, &y_bottom, &destX_bottom, &destY_bottom);
		monotonic_cubic_Hermite_spline(&x_peak, &y_peak, &destX_peak, &destY_peak);
	}
	else
	{
		cubic_spline(&x_bottom, &y_bottom, &destX_bottom, &destY_bottom);
		cubic_spline(&x_peak, &y_peak, &destX_peak, &destY_peak);
	}
	if (((int)destY_bottom.size() != N) || ((int)destY_peak.size() != N))
	{
		std::cerr << "spline error" << std::endl;
		return -1;
	}
	for(int i=0; i < N; i++)
	{
		imf->push_back(src[i] - (destY_peak[i] + destY_bottom[i])/2.0);
	}
	return 0;
}

int sift(const std::vector<FLOAT> *src, std::vector<FLOAT> *imf, const int mono)
{
	std::vector<FLOAT> x_series;
	std::vector<FLOAT> y_series;
	std::vector<FLOAT> destX_peak;
	std::vector<FLOAT> destY_peak;
	std::vector<FLOAT> destX_bottom;
	std::vector<FLOAT> destY_bottom;
	
	std::vector<FLOAT> x_peak;
	std::vector<FLOAT> y_peak;

	std::vector<FLOAT> x_bottom;
	std::vector<FLOAT> y_bottom;

	std::vector<FLOAT> destY_middle;

	//peak & bottom search
	// if (extrema(src->data(), src->size(), x_peak, y_peak, x_bottom, y_bottom) < 3) return -1;
	extrema(src->data(), src->size(), x_peak, y_peak, x_bottom, y_bottom);
	
	if( mono == 1)
	{
		monotonic_cubic_Hermite_spline(&x_bottom, &y_bottom, &destX_bottom, &destY_bottom);
		monotonic_cubic_Hermite_spline(&x_peak, &y_peak, &destX_peak, &destY_peak);
	}
	else
	{
		cubic_spline(&x_bottom, &y_bottom, &destX_bottom, &destY_bottom); 	
		cubic_spline(&x_peak, &y_peak, &destX_peak, &destY_peak);
	}
	int N = src->size();
	if (((int)destY_bottom.size() != N) || ((int)destY_peak.size() != N))
	{
		std::cerr << "spline error" << std::endl;
		return -1;
	}
	for(int i=0; i < N; i++)
	{
		imf->push_back((*src)[i] - (destY_peak[i] + destY_bottom[i])/2.0);
	}
	return 0;
}

int stopcondition(const std::vector<FLOAT> *src)
{
	int i;
	//int p_sign = (*src)[1] - (*src)[0] > 0 ? 1 : 0;
	FLOAT diff = (*src)[1] - (*src)[0];
	int p_sign = sgn(diff);
	int sign;
	//int p_zero = (*src)[0] > 0 ? 1 : 0;
	int p_zero = sgn((*src)[0]);
	int zero;
	int Next = 0;
	int Nzero = 0;
	for(i=1; i < (int)src->size()-1; i++)
	{
		//sign = ((*src)[i+1] - (*src)[i] > 0)?1:0;
		diff = (*src)[i+1] - (*src)[i];
		sign = sgn(diff);
		if(sign != 0)
		{
			if(sign != p_sign)
			{
				Next++;
			}
			p_sign = sign;
		}

		zero = sgn((*src)[i]);
		if(zero != 0)
		{
			if(zero != p_zero)
			{
				Nzero++;
			}
			p_zero = zero;
		}
	}   

	//  printf("%d - %d = %d\n", Nzero, Next, Nzero-Next);
	return Nzero - Next > 0? Nzero - Next : Next - Nzero;
}

void emd(int N, const FLOAT *origin, FLOAT *imf_series, int &MAX_IMF, const int MAX_ITERATION, const int mono, const int s_number)
{
	std::vector<FLOAT> vec1, vec2;
	std::vector<FLOAT> *src = &vec1;
	std::vector<FLOAT> *des = &vec2;
	std::vector<FLOAT> *vswap;

#ifdef DEBUG	
	std::cout << "EMD Processing..." << N << std::endl;
	std::cout << "MAX_IMF : " << MAX_IMF << std::endl;
	std::cout << "MAX_ITERATION : " << MAX_ITERATION << ", s_number : " << s_number << std::endl;
#endif

	FLOAT time_series[N];
	memcpy(time_series, origin, sizeof(FLOAT)*N);
	for(int index =0; index < MAX_IMF; index++)
	{
		(*src).clear();
		sift(N, time_series, src, mono);
		for(int k = 1; k < MAX_ITERATION; k++)
		{
			(*des).clear();
			if(sift(src, des, mono) < 0)
			{
				break;
			}
			
			if (stopcondition(des) <= s_number)
			{
				break;
			}
			vswap = src;
			src = des;
			des = vswap;
		}
		for(int i=0; i < N; i++)
		{
			if(std::isnan((*des)[i])) 
			{
				std::cerr << "imf [" << index <<"], " << i << std::endl;
				MAX_IMF = index;
				break;
			}
			time_series[i] -= (*des)[i];
			imf_series[index*N + i] = (*des)[i];
		}	
	}
}

void weightfunction_exp(FLOAT *w, int n, FLOAT alpha)
{
	//#pragma omp parallel for
	for(int i = 0; i < n; i++)
	{
		//w[i] = std::exp(-1.0/2.0*std::pow(alpha*(2.0*i/n -1.0),2.0));
		FLOAT x = alpha*(i/(n-1) - 0.5);
		w[i] = std::exp(-2.0*x*x);
	}
}

void weightfunction_sin(FLOAT *w, int n)
{
	const FLOAT pi = 3.141592653;
	//#pragma omp parallel for
	for(int i = 0; i < n; i++)
	{
		//w[i] = std::pow(std::sin(i*pi/n),2.0);
		w[i] = std::sin(i*pi/(n-1));
		w[i] *= w[i];
	}
}


void weightfunction(FLOAT *w, int n, int type, FLOAT alpha)
{
	switch(type)
	{
		case exp_kernel:
			weightfunction_exp(w, n, alpha);
			break;
		case sin_kernel:
			weightfunction_sin(w, n);
			break;
		case spline_kernel:
			break;
		default:
			break;
			//fprintf(stderr, "Not suppoted Weighted Function Type\n");
	}

}

int wsemd(int N, int segment, int bufnum, FLOAT *const weight, FLOAT *const origin, FLOAT *imf_series, 
	int &numberofimf, const int iteration, const int mono, const int s_number)
{
	const int fragment = segment/bufnum;
	segment = fragment * bufnum;
	const int imfsize = (N/fragment - 2*(bufnum-1));

	int total_buf_size = N / fragment - bufnum + 1;
	FLOAT **emd_buf;

	if(N % fragment > 0)
	{
		std::cout<<"data_size can't divide by segment size"<<std::endl;
		return -1;  
	}


	emd_buf = (FLOAT**)malloc(sizeof(FLOAT*)*total_buf_size);

	for(int i =0 ; i < total_buf_size; i++)
	{
		emd_buf[i] = (FLOAT*)malloc(sizeof(FLOAT)*segment*numberofimf);
	}

	/* initialize buffer */
	#pragma omp parallel for
	for(int i=0; i < total_buf_size; i++)			
	{
		FLOAT* pdata = origin + fragment*i;
		emd(segment, pdata, emd_buf[i], numberofimf, iteration, mono, s_number);	
	}
	
	/* Front Part of imf series */
	for(int frag = 0; frag < bufnum - 1; frag++)
	{
		for (int i = 0; i < numberofimf; i++)							//calcuate each imfs
		{
			FLOAT *imftmp = imf_series + N * i + fragment * frag;		//target imf start point
			for(int j = 0; j < fragment; j++)							//calcuating fragment
			{
				FLOAT wtmp, wsum = 0.0;
				FLOAT sum = 0.0;
				for(int k =0; k < frag +1 ; k++)		//sum of buffer for buffer length
				{
					wtmp = (frag == 0) ? 1 : weight[fragment * (frag - k) + j];
					sum += wtmp * emd_buf[k][ fragment*(frag - k) + segment * i + j ];
					wsum += wtmp;
					//wtmp += weight[fragment * k + j];
				}
				*(imftmp + j) = sum / wsum;
			}
		}

	}

	/* Middle of imf series */
	for(int frag = 0; frag < imfsize; frag++)					//calculate each fragments
	{
		for (int i = 0; i < numberofimf; i++)							//calcuate each imfs
		{
			FLOAT *imftmp = imf_series + N * i + fragment * (bufnum - 1 + frag);		//target imf start point
			for(int j = 0; j < fragment; j++)							//calcuating fragment
			{
				FLOAT wtmp, wsum = 0.0;
				FLOAT sum = 0.0;
				for(int k =0; k < bufnum; k++)		//sum of buffer for buffer length
				{
					wtmp = weight[fragment * (bufnum - 1 - k) + j];
					sum += wtmp * emd_buf[frag + k][ fragment*(bufnum - 1 - k) + segment * i + j ];
					wsum += wtmp;
				}
				*(imftmp + j) = sum / wsum;
			}
		}
	}

	/*Back-End part of imf series */
	for(int frag = 0; frag < bufnum - 1; frag++)
	{
		for (int i = 0; i < numberofimf; i++)							//calcuate each imfs
		{
			FLOAT *imftmp = imf_series + N * i + fragment * (imfsize + bufnum - 1 + frag);		//target imf start point
			for(int j = 0; j < fragment; j++)							//calcuating fragment
			{
				FLOAT wtmp, wsum = 0.0;
				FLOAT sum = 0.0;
				for(int k = bufnum - 1; k > frag; k--)		//sum of buffer for buffer length
				{
					wtmp = (frag == bufnum - 2) ? 1 : weight[fragment * k + j];
					sum += wtmp * emd_buf[imfsize + frag + bufnum - 1 - k][ fragment * k + segment * i + j ];
					wsum += wtmp;
				}
				*(imftmp + j) = sum / wsum;
			}
		}

	}

	for(int i = 0; i < total_buf_size; i++)
	{
		free(emd_buf[i]);
	}
	free(emd_buf);
	return 0;

}
