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
#include <numpy/ndarrayobject.h>
#include <vector>
#include <cmath>
#include "emd.hpp"
#include "cubic.hpp"
#include "DHT.h"
#include "omp.h"
#include "trigger.h"
#include "cluster.hpp"
#include <iostream>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define CFLOAT std::complex<FLOAT>
#define NPFLOAT NPY_DOUBLE
#define NPCFLOAT NPY_CDOUBLE

#define MAX_SIFT 15
#define MIN_DHT_FILTER_LENGTH 128
#define MIN_STRIDE 256
#define DEFAULT_SNR_THRESHOLD 5.5
#define MIN_GEN_TRG_LENGTH 2048

#define getArray(xarr) (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) \
			                        PyArray_EnsureArray(xarr))

using namespace boost::python;

handle<> handle_none(Py_None);
static object none(handle_none);

inline numeric::array pyArray_FromData(int ndim, npy_intp* shape, FLOAT* data)
{
	PyObject* pyObj = PyArray_SimpleNewFromData(ndim, shape, NPFLOAT, data);
	handle<> handle(pyObj);
	numeric::array arr(handle);
	return arr;
}

inline numeric::array pyArray_FromData(int ndim, npy_intp* shape, CFLOAT* data)
{
	PyObject* pyObj = PyArray_SimpleNewFromData(ndim, shape, NPCFLOAT, data);
	handle<> handle(pyObj);
	numeric::array arr(handle);
	return arr;
}

numeric::array stdVecToNumpyArray( std::vector<FLOAT> const& vec )
{
	npy_intp size = vec.size();
	FLOAT *data = size ? const_cast<FLOAT *>(&vec[0]) 
	: static_cast<FLOAT *>(NULL);
	return pyArray_FromData(1, &size, data);
}

numeric::array DoubleToNumpyArray( int s, FLOAT *data )
{
	npy_intp size = s;
	return pyArray_FromData(1, &size, data);
}

numeric::array DoubleToNumpyArray( int nd, int s, FLOAT *data )
{
	npy_intp size[2];
	size[0]= nd;
	size[1]= s;
	return pyArray_FromData(2, size, data);
}

numeric::array DoubleToNumpyArray( int nd, int s, CFLOAT *data )
{
	npy_intp size[2];
	size[0]= nd;
	size[1]= s;
	return pyArray_FromData(2, size, data);
}

class etagen
{
	private:
		int MAX_ITERATION;
		int monotonic;
		int weight_type;
		FLOAT alpha;
		int s_number;
		int numberofbuf;
		int bufoffset;
		int imfbuf_size;
		FLOAT *weight;
		int filter_length;
		FLOAT *w;
		triggerCluster tc;
	public:
		object data;
		object imf;
		object res;
		object hht;
		object insa;
		object insf;
		object utrgs;
		object trgs;

		int numberofimf;
		int data_size;

		FLOAT start_time;
		FLOAT fsr;

		etagen(object & _data = none, FLOAT _fsr = 1.0,
		       FLOAT _start_time = 0.0)
		{
			//public variables initialization
			imf       = none;
			res       = none;
			hht       = none;
			insa      = none;
			insf      = none;
			utrgs     = none;
			trgs      = none;

			//private variables initialization
			imfbuf_size = 0;
			bufoffset = 0;
			weight = NULL;
			filter_length = 0;
			w = NULL;

			//using member functions
			set_start_time(_start_time);
			set_fsr(_fsr);
			set_emd_param();
			set_data(_data);
		}
		~etagen()
		{
			if (weight != NULL) free(weight);
			if (w != NULL) free(w);
		}
		void set_data(object&);
		void set_start_time(FLOAT _st)
		{
			start_time = _st;
		}
		void set_fsr(FLOAT _fs)
		{
			fsr = _fs;
		}
		//void add_data(object&);
		//object& get_data();
		void set_emd_param(int num_imfs=0, int num_sifts=MAX_SIFT,
				   int S_number=0, int monotonic_spline=0,
				   int emd_size=2048, int num_seg=8,
				   int w_type=sin_kernel);
		void show_emd_param();
		void emd_();
		void wsemd_();
		void set_imf(object&);
		//object& get_imf;
		void set_res(object&);
		void hilbert_(int _filter_len=MIN_DHT_FILTER_LENGTH, int _stride=MIN_STRIDE);
		void set_hht(object&);
		//object& get_hht;
		void set_insa(object&);
		//object& get_insa;
		void set_insf(object&);
		//object& get_insf;
		list gen_utrgs(FLOAT, int, int);
		numeric::array get_waveform(int);
		//object& get_utrgs;
		numeric::array gen_trgs(numeric::array, FLOAT, FLOAT, FLOAT);
		//object& get_trgs;
};

void etagen::set_emd_param(int num_imfs, int num_sifts, int S_number,
		   int monotonic_spline, int emd_size, int num_seg, int w_type)
{
	int max_imfs = (int) log2(emd_size) - 1;
	numberofimf = num_imfs ? std::min(max_imfs, num_imfs) : max_imfs;
	MAX_ITERATION = num_sifts > 0 ? num_sifts : MAX_SIFT;
	monotonic = monotonic_spline;
	weight_type = w_type;
	alpha = (FLOAT) numberofbuf / 2.0;
	s_number = S_number;
	numberofbuf = num_seg > 0 ? num_seg : 8;
	bufoffset = (int) (emd_size / numberofbuf);
	if (weight != NULL) free(weight);
	weight = (FLOAT*)malloc(sizeof(FLOAT)*numberofbuf * bufoffset);
	weightfunction(weight, numberofbuf * bufoffset, weight_type, alpha);
}

void etagen::show_emd_param()
{
	std::stringstream msg;
	msg << "Parameters are set as follows." << std::endl;
	msg << "[General]" << std::endl;
	msg << "Number of IMFs:\t" << numberofimf << std::endl;
	msg << "Number of siftings:\t" << MAX_ITERATION << std::endl;
	msg << "S-Number:\t" << s_number << std::endl;
	msg << "[wSEMD settings]" << std::endl;
	msg << "EMD size:\t" << numberofbuf * bufoffset << std::endl;
	msg << "Number of segments:\t" << numberofbuf << std::endl;
	PySys_WriteStdout("%s", msg.str().c_str());
}

void etagen::set_data(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	data = a;
	data_size = PyArray_Size(getArray(a.ptr()));
}

void etagen::set_imf(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"imf should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a.ptr())[1] != data_size)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of imf and data does not match.");
		throw_error_already_set();
	}
	imf = a;
	numberofimf = PyArray_DIMS(a.ptr())[0];
	imfbuf_size = data_size * numberofimf;
}

void etagen::set_res(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"res should be a 1D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_Size(a.ptr()) != data_size)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of res and data does not match.");
		throw_error_already_set();
	}
	res = a;
}

void etagen::set_hht(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"hht should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a.ptr())[1] != PyArray_DIMS(data.ptr())[0])
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of hht and data does not match.");
		throw_error_already_set();
	}
	hht = a;
}

void etagen::set_insa(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"insa should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a.ptr())[1] != PyArray_DIMS(data.ptr())[0])
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of insa and data does not match.");
		throw_error_already_set();
	}
	insa = a;
}

void etagen::set_insf(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"insf should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a.ptr())[1] != PyArray_DIMS(data.ptr())[0]-1)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of insf and data does not match.");
		throw_error_already_set();
	}
	insf = a;
}

void etagen::emd_()
{
	FLOAT* _time_series = static_cast<FLOAT *>(PyArray_DATA(data.ptr()));
	FLOAT* _imf_series = (FLOAT*)malloc(sizeof(FLOAT) * data_size * \
					    numberofimf);
	emd(data_size, _time_series, _imf_series, numberofimf,
	    MAX_ITERATION, monotonic, s_number); 
	imf = DoubleToNumpyArray(numberofimf, data_size, _imf_series);
	free(_imf_series);
	imfbuf_size = data_size * numberofimf;
}

void etagen::wsemd_()
{
	int emd_size = numberofbuf * bufoffset;	//emd slide size

	FLOAT* _time_series = static_cast<FLOAT *>(PyArray_DATA(data.ptr()));
	FLOAT* _imf_series = (FLOAT*)malloc(sizeof(FLOAT) * data_size * \
					    numberofimf);
	wsemd(data_size, emd_size, numberofbuf, weight,
	      _time_series, _imf_series, numberofimf,
	      MAX_ITERATION, monotonic, s_number); 
	imf = DoubleToNumpyArray(numberofimf, data_size, _imf_series);
	free(_imf_series);
	imfbuf_size = data_size * numberofimf;
}

void etagen::hilbert_(int _filter_len, int _stride)
{
	npy_intp* shape = PyArray_DIMS(imf.ptr());
	npy_intp fshape[2] = {shape[0], shape[1] - 1};

	int _fl = (_filter_len > MIN_DHT_FILTER_LENGTH) ? _filter_len
		  : MIN_DHT_FILTER_LENGTH;
	int _strd = (_stride >= MIN_STRIDE) ? _stride : data_size;
	if (filter_length != _fl)
	{
		filter_length = _fl;
		if (w != NULL) free(w);
		w = (FLOAT*)malloc(sizeof(FLOAT) * filter_length);
		fir_filter(filter_length, w);
	}

	FLOAT* _imf_series = static_cast<FLOAT *>(PyArray_DATA(imf.ptr()));
	FLOAT* _dht_result = (FLOAT*)malloc(sizeof(FLOAT) * imfbuf_size);
	CFLOAT* _hht_series = (CFLOAT*)malloc(sizeof(CFLOAT) * imfbuf_size);
	FLOAT* _dht_amplitude = (FLOAT*)malloc(sizeof(FLOAT) * imfbuf_size);
	FLOAT* _dht_frequency = (FLOAT*)malloc(sizeof(FLOAT) * fshape[0] * fshape[1]);

#pragma omp parallel for			
	for(int i=0; i < shape[0]; i++)
	{
		dht(shape[1], _imf_series + shape[1]*i, _dht_result + shape[1]*i, w, filter_length);
#pragma omp parallel for			
		for(int j=0; j < shape[1]; j += _strd)
		{
			bool last = (j + _strd < shape[1]) ? false : true;
			hsa(std::min(_strd, (int)shape[1] - j), _imf_series + shape[1]*i + j, _dht_result + shape[1]*i + j, _dht_amplitude + shape[1]*i + j, _dht_frequency + fshape[1]*i + j, fsr, last);
		}
	}
#pragma omp parallel for			
	for(int i=0; i < shape[0] * shape[1]; i++)
	{
		*(_hht_series + i) = CFLOAT (*(_imf_series + i), *(_dht_result + i));
	}

	hht = (object) DoubleToNumpyArray(shape[0], shape[1], _hht_series);
	insa = (object) DoubleToNumpyArray(shape[0], shape[1], _dht_amplitude);
	insf = (object) DoubleToNumpyArray(fshape[0], fshape[1], _dht_frequency);
	free(_dht_result);
	free(_hht_series);
	free(_dht_amplitude);
	free(_dht_frequency);
}

list etagen::gen_utrgs(FLOAT snr_th, int sidx, int len)
{
	list ret;
	if ((len < MIN_GEN_TRG_LENGTH) || (len > data_size))
	{
		len = data_size;
		std::stringstream msg;
		msg << "*length* was reset: " << len << std::endl;
		PySys_WriteStdout("%s", msg.str().c_str());
	}
	if (sidx < 0)
	{
		sidx = 0;
		std::stringstream msg;
		msg << "*start* was reset: " << sidx << std::endl;
		PySys_WriteStdout("%s", msg.str().c_str());
	}
	if (sidx > data_size - len)
	{
		sidx = data_size - len;
		std::stringstream msg;
		msg << "*start* was reset: " << sidx << std::endl;
		PySys_WriteStdout("%s", msg.str().c_str());
	}
	FLOAT* _time_series = static_cast<FLOAT *>(PyArray_DATA(data.ptr()));
	FLOAT* _imf_series = static_cast<FLOAT *>(PyArray_DATA(imf.ptr()));
	FLOAT* _dht_amplitude = static_cast<FLOAT *>(PyArray_DATA(insa.ptr()));
	FLOAT* _dht_frequency = static_cast<FLOAT *>(PyArray_DATA(insf.ptr()));
	for(int i=0; i < numberofimf; i++)
	{
		ret.append(trigger_gen(i, _imf_series + data_size*i + sidx, len, _dht_amplitude + data_size*i + sidx, _dht_frequency + (data_size-1)*i + sidx, med_abs_dev(data_size, _time_series), snr_th));
		//ret.append(trigger_gen(i, _imf_series + data_size*i + sidx, len, _dht_amplitude + data_size*i + sidx, _dht_frequency + (data_size-1)*i + sidx, median(data_size, abs(data_size, _time_series)), snr_th));
	}
	return ret;
}

numeric::array etagen::gen_trgs(numeric::array _utrgs, FLOAT snr_th, FLOAT ttol, FLOAT ftol)
{
	npy_intp _ntrgs = PyArray_SIZE(_utrgs.ptr());
	if (_ntrgs == 0) return DoubleToNumpyArray(0, 14, (FLOAT*)0);
	PyObject *pyarr = PyArray_FromAny(PyArray_ToList((PyArrayObject*) _utrgs.ptr()), NULL, 0, 0, 0, NULL);
	trginfo *_trg = (trginfo*) PyArray_DATA(pyarr);
	cltinfo *clt;
	npy_intp* shape = PyArray_DIMS(imf.ptr());
	FLOAT* _imf_series = static_cast<FLOAT *>(PyArray_DATA(imf.ptr()));
	FLOAT* _time_series = static_cast<FLOAT *>(PyArray_DATA(data.ptr()));
	
	tc.set_param(ttol, ftol, med_abs_dev(data_size, _time_series), _imf_series, shape[0], shape[1], start_time, fsr);
	for(int i=0; i < _ntrgs; i++)
	{
		tc.feed(&_trg[i]);
	}

	clt = tc.getClusteredTrigger(snr_th);

	return DoubleToNumpyArray(tc.numofCluster, 14, (FLOAT*)clt);
}

numeric::array etagen::get_waveform(int index)
{
	int len;
	FLOAT *wave;
	len = tc.getWaveform(index, &wave);
	return DoubleToNumpyArray(len, wave);
}
object extrema_(object & a)
{
	if (PyArray_NDIM(a.ptr()) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	PyObject* pyarr = getArray(a.ptr());
	FLOAT* data = (FLOAT*) PyArray_DATA(pyarr);
	int n = PyArray_Size(pyarr);
	std::vector<FLOAT> x_peak, y_peak, x_bot, y_bot;
	extrema(data, n, x_peak, y_peak, x_bot, y_bot);
	return make_tuple(stdVecToNumpyArray(x_peak),
			  stdVecToNumpyArray(y_peak),
			  stdVecToNumpyArray(x_bot),
			  stdVecToNumpyArray(y_bot));
}

object spline_(object & _x, object & _y, bool mono=false)
{
	if ((PyArray_NDIM(_x.ptr()) != 1) || (PyArray_NDIM(_y.ptr()) != 1))
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_Size(_x.ptr()) != PyArray_Size(_y.ptr()))
	{
		PyErr_SetString(PyExc_ValueError,
				"x and y should have same size.");
		throw_error_already_set();
	}
	PyObject* xarr = getArray(_x.ptr());
	PyObject* yarr = getArray(_y.ptr());
	FLOAT* xdata = (FLOAT*) PyArray_DATA(xarr);
	FLOAT* ydata = (FLOAT*) PyArray_DATA(yarr);
	int n = PyArray_Size(xarr);
	std::vector<FLOAT> x_pts (xdata, xdata+n);
	std::vector<FLOAT> y_pts (ydata, ydata+n);
	std::vector<FLOAT> x_data, y_data;
	mono ? monotonic_cubic_Hermite_spline(&x_pts, &y_pts, &x_data, &y_data)
	     : cubic_spline(&x_pts, &y_pts, &x_data, &y_data);
	return (object) stdVecToNumpyArray(y_data);
}

BOOST_PYTHON_MODULE(_etagen)
{
	numeric::array::set_module_and_type("numpy", "ndarray");
	import_array();
	enum_<weightTpye>("weightType")
	.value("EXP_KERNEL", exp_kernel)
	.value("SIN_KERNEL", sin_kernel)
	//.value("SPLINE_KERNEL", spline_kernel)
	;
	def("find_extrema", &extrema_,
	    arg("data"),
	    ""
	    );
	def("evaluate_spline", &spline_,
	    (arg("x"), arg("y"), arg("monotonic")=false),
	    ""
	    );
	/*
	def("emd", &emd_,
	    ""
	    );
	def("wsemd", &wsemd_,
	    ""
	    );
	*/
	class_<etagen>("etagen",
	               "EtaGen class",
		       init<object &, FLOAT, FLOAT>((arg("data")=none,
			       arg("fsr")=1.0, arg("start_time")=0.0),
			   "construct an EtaGen instance"
			   )
	    )
	.def_readonly("start_time", &etagen::start_time,
	    "start time of input data"
	    )
	.def_readonly("fsr",        &etagen::fsr,
	    "sampling rate"
	    )
	.def_readonly("data_len",       &etagen::data_size,
	    "Length of the input data"
	    )
	.def_readonly("nimfs",       &etagen::numberofimf,
	    "Number of IMFs"
	    )
	.def_readonly("data",       &etagen::data,
	    "input data"
	    )
	.def_readonly("imfs",       &etagen::imf,
	    "IMF time series (EMD)"
	    )
	.def_readonly("res",        &etagen::res,
	    "residual time series"
	    )
	.def_readonly("hht",        &etagen::hht,
	    "HHT time series (EMD + HSA)"
	    )
	.def_readonly("insa",       &etagen::insa,
	    "instantaneous amplitude of each IMFs"
	    )
	.def_readonly("insf",       &etagen::insf,
	    "instantaneous frequency of each IMFs"
	    )
	.def_readwrite("utrgs",      &etagen::utrgs,
	    "preliminary (unclustered) triggers"
	    )
	.def_readwrite("trgs",       &etagen::trgs,
	    "(clustered) triggers"
	    )
	.def("set_start_time",      &etagen::set_start_time,
	    ""
	    )
	.def("set_fsr",             &etagen::set_fsr,
	    ""
	    )
	.def("set_emd_param",       &etagen::set_emd_param,
	    (arg("num_imfs")=0, arg("num_sifts")=MAX_SIFT,
	       arg("S_number")=0, arg("monotonic_spline")=0,
	       arg("emd_size")=2048, arg("num_seg")=8,
	       arg("w_type")=sin_kernel),
	    ""
	    )
	.def("show_emd_param",      &etagen::show_emd_param,
	    ""
	    )
	.def("set_data",            &etagen::set_data,
	    ""
	    )
	//.def("add_data",            &etagen::add_data,
	//    ""
	//    )
	//.def("get_data",            &etagen::get_data,
	//    ""
	//    )
	.def("emd",                 &etagen::emd_,
	    ""
	    )
	.def("wsemd",               &etagen::wsemd_,
	    ""
	    )
	.def("set_imf",             &etagen::set_imf,
	    ""
	    )
	//.def("get_imf",             &etagen::get_imf,
	//    ""
	//    )
	.def("set_res",             &etagen::set_res,
	    ""
	    )
	.def("hilbert",             &etagen::hilbert_,
	    (arg("filter_length")=MIN_DHT_FILTER_LENGTH, arg("stride")=MIN_STRIDE),
	    ""
	    )
	.def("set_hht",             &etagen::set_hht,
	    ""
	    )
	//.def("get_hht",             &etagen::get_hht,
	//    ""
	//    )
	.def("set_insa",            &etagen::set_insa,
	    ""
	    )
	//.def("get_insa",            &etagen::get_insa,
	//    ""
	//    )
	.def("set_insf",            &etagen::set_insf,
	    ""
	    )
	//.def("get_insf",            &etagen::get_insf,
	//    ""
	//    )
	.def("gen_utrgs",           &etagen::gen_utrgs,
	    (arg("snr_threshold")=DEFAULT_SNR_THRESHOLD, arg("start")=0, arg("length")=0),
	    ""
	    )
	//.def("get_utrgs",           &etagen::get_utrgs,
	//    ""
	//    )
	.def("gen_trgs",            &etagen::gen_trgs,
	    (arg("triggers"), arg("snr_threshold")=DEFAULT_SNR_THRESHOLD, arg("t_tolerance")=0, arg("f_tolerance")=0),
	    ""
	    )
	.def("get_waveform",            &etagen::get_waveform,
	   ""
	   )
	;
}
