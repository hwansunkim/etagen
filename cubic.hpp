#ifndef __CUBIC_H__
#define __CUBIC_H__

#include "common.h"

bool cubic_spline(std::vector<FLOAT>* x_series, std::vector<FLOAT>* y_series, std::vector<FLOAT> *destX, std::vector<FLOAT>* destY); 
bool monotonic_cubic_Hermite_spline(std::vector<FLOAT>*, std::vector<FLOAT>*, std::vector<FLOAT>*, std::vector<FLOAT>*);

#endif