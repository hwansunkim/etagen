#ifndef __DHT_H__
#define __DHT_H__

#include "common.h"

void fir_filter(int n, FLOAT* w);
int dht(int n, FLOAT *src, FLOAT *des, FLOAT *w, int filter_length);
int hsa(int n, FLOAT *src, FLOAT *imag, FLOAT *amplitude, FLOAT *frequency, FLOAT fs, bool last);
FLOAT median(int n, FLOAT* src);

#endif