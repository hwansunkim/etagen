#ifndef __TRIGGER_H__
#define __TRIGGER_H__

#include <boost/python.hpp>
using namespace boost::python;

object trigger_gen(int imfs, FLOAT* data, int n, FLOAT *amplitude, FLOAT *frequency, FLOAT m, FLOAT snr_th);
FLOAT* abs(int n, FLOAT* src);
#endif
