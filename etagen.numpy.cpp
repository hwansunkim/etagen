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

/**
 * @package EtaGen
 * # EtaGen
 * An event trigger generator based on Hilbert-Huang Transform<br>
 * It is a Python module but the core library of EtaGen is built in C/C++
 *
 * ## Sample USAGE:
 * * This shows how to generate event triggers from *data* by EtaGen<br>
 *     >>> from etagen import etagen, kernel
 *
 *     >>> h = etagen(data, fsr=1024)
 *
 *     >>> h.set_emd_param(num_imfs=8, num_sifts=15, S_number=1, emd_size=1024, num_seg=4, w_type=kernel.SIN_KERNEL)
 *
 *     >>> h.show_emd_param()<br>
 *     Parameters are set as follows.<br>
 *     [General]<br>
 *     Number of IMFs: 8<br>
 *     Number of siftings:     15<br>
 *     S-Number:       1<br>
 *     [wSEMD settings]<br>
 *     EMD size:       1024<br>
 *     Number of segments:     4
 *
 *     >>> h.wsemd()
 *
 *     >>> h.hilbert(filter_length=128, stride=1024)
 *
 *     >>> h.get_utriggers(snr_th=5, stride=4*fsr, overlap=2*fsr)<br>
 *     Generating triggers with 5-snr threshold in segments of length 4096, overlapping 2048 samples and skipping 0 samples from boundaries...<br>
 *     ... generated 25 trigger event(s)
 *
 *     >>> h.get_triggers(t_tolerance=0.001, f_tolerance=0.5, snr_th=5.5)<br>
 *     u_snr_th should be larger than the one used to generate triggers: assuming u_snr_th = 5<br>
 *     Clustering triggers of u_snr > 5 with time tolerance=0.001, frequency tolerance=0.5 and dropping clusters of snr < 5.5...<br>
 *     total Clusters : 25 > 5.500000<br>
 *     total Clusters : 12 > 5.500000<br>
 *     ... generated 12 trigger cluster(s)
 *
 *     >>> h.trgs[['c_time','c_freq','p_time','p_freq','npts','snr']]<br>
 *     array([ (0.041015625, 297.67912076702106, 0.041015625, 297.67912076702106, 1, 5.546819634000152),<br>
 *     (0.6464843749999999, 344.6262005898434, 0.646484375, 344.6262005898434, 1, 5.797973757254617),<br>
 *     (1.2744140625, 260.1084208607948, 1.2744140625, 260.1084208607948, 1, 5.546288486345601),<br>
 *     (1.5664062499999998, 142.61211763115602, 1.56640625, 142.61211763115602, 1, 6.513913449375204),<br>
 *     (1.8193359374999998, 333.3958253278954, 1.8193359375, 333.39582532789547, 1, 6.632835601773412),<br>
 *     (2.7138671875, 341.5032411355364, 2.7138671875, 341.5032411355364, 1, 6.609969624271628),<br>
 *     (4.1552734375, 220.52834605967874, 4.1552734375, 220.5283460596787, 1, 5.8852032771220015),<br>
 *     (4.493164062499999, 201.227897198888, 4.4931640625, 201.22789719888803, 1, 5.890950833813454),<br>
 *     (5.766601562500001, 229.12911997930595, 5.7666015625, 229.12911997930595, 1, 5.557695477749056),<br>
 *     (6.659179687499999, 245.8274380780258, 6.6591796875, 245.8274380780258, 1, 5.8990177018284875),<br>
 *     (6.328125, 34.86349503817338, 6.328125, 34.86349503817338, 1, 5.834112766561207),<br>
 *     (6.0, 67.39566743927662, 6.0, 67.39566743927662, 1, 5.9916237288573955)], <br>
 *     type=[('c_time', '<f8'), ('c_freq', '<f8'), ('p_time', '<f8'), ('p_freq', '<f8'), ('npts', '<i8'), ('snr', '<f8')])
 *
 */

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <vector>
#include <cmath>
#include "emd.hpp"
#include "cubic.hpp"
#include "DHT.h"
#include "omp.h"
#include "trigger.h"
#include "cluster.hpp"
#include <iostream>
#include <stdio.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define CFLOAT std::complex<FLOAT>
#define NPFLOAT NPY_DOUBLE
#define NPCFLOAT NPY_CDOUBLE

#define MAX_SIFT 15
#define MIN_DHT_FILTER_LENGTH 128
#define MIN_STRIDE 256
#define DEFAULT_SNR_THRESHOLD 5.5
#define MIN_GEN_TRG_LENGTH 2048

namespace p = boost::python;
namespace np = boost::python::numpy;

#if PY_VERSION_HEX >= 0x03000000
#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#else
#define NUMPY_IMPORT_ARRAY_RETVAL
#endif

#define import_array() {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import"); return NUMPY_IMPORT_ARRAY_RETVAL; } }

handle<> handle_none(Py_None);
static p::object none(handle_none);

inline int const PyArray_NDIM(p::object obj)
{
	return np::from_object(obj).get_nd();
}

inline Py_intptr_t const* PyArray_DIMS(p::object obj)
{
	return np::from_object(obj).get_shape();
}

inline int PyArray_Size(np::ndarray arr)
{
	return arr.shape(0);
}
inline int PyArray_Size(p::object obj)
{
	return np::from_object(obj).shape(0);
}
inline FLOAT* PyArray_DATA(p::object obj)
{
	return (FLOAT*)(np::from_object(obj).get_data());
}
inline np::ndarray pyArray_FromData(int ndim, int* shape, FLOAT* data)
{
	p::tuple shape_;
	p::tuple stride_;
	if(ndim == 1)
	{
		shape_ = make_tuple(shape[0]);
		stride_ = make_tuple(sizeof(FLOAT));
	}
	else
	{
		shape_ = make_tuple(shape[0], shape[1]);
		stride_ = make_tuple((int)(sizeof(FLOAT) * shape[1]), (int)sizeof(FLOAT));
	}
	np::ndarray arr = np::from_data(data, np::dtype::get_builtin<FLOAT>(),
			shape_,
			stride_,
			p::object());
	return arr;
}

inline np::ndarray pyArray_FromData(int ndim, int* shape, CFLOAT* data)
{
	p::tuple shape_;
	p::tuple stride_;
	if(ndim == 1)
	{
		shape_ = make_tuple(shape[0]);
		stride_ = make_tuple((int)sizeof(CFLOAT));
	}
	else
	{
		shape_ = make_tuple(shape[0], shape[1]);
		stride_ = make_tuple((int)(sizeof(CFLOAT) * shape[1]), (int)sizeof(CFLOAT));
	}
	np::ndarray arr = np::from_data(data, np::dtype::get_builtin<CFLOAT>(),
			shape_,
			stride_,
			p::object());
	return arr;
}

np::ndarray stdVecToNumpyArray( std::vector<FLOAT> const& vec )
{
	int size = vec.size();
	FLOAT *data = size ? const_cast<FLOAT *>(&vec[0]) 
	: static_cast<FLOAT *>(NULL);
	return pyArray_FromData(1, &size, data);
}

np::ndarray DoubleToNumpyArray( int s, FLOAT *data )
{
	int size = s;
	return pyArray_FromData(1, &size, data);
}

np::ndarray DoubleToNumpyArray( int nd, int s, FLOAT *data )
{
	int size[2];
	size[0]= nd;
	size[1]= s;
	return pyArray_FromData(2, size, data);
}

np::ndarray DoubleToNumpyArray( int nd, int s, CFLOAT *data )
{
	int size[2];
	size[0]= nd;
	size[1]= s;
	return pyArray_FromData(2, size, data);
}

/**@class etagen
 * EtaGen class
 *
 * USAGE (an example):<br>
 *     h = etagen.etagen(data, fsr, start_time)
 *
 * arguments:<br>
 *     data          input data which is assumed to be whitened<br>
 *     fsr           sampling frequency<br>
 *     start_time    start time of data in sec.
 *
 * return:<br>
 *     h             an instance of etagen class<br>
 *
 * attribute stored:<br>
 *     h.data<br>
 *     h.fsr<br>
 *     h.start_time
 *
 * example:<br>
 *     >>> import numpy as np
 *
 *     >>> from etagen import etagen
 *
 *     >>> h = etagen(np.random.randn(1024),32)
 *
 *     >>> print h.fsr<br>
 *     32.0
 *
 *     >>> print h.start_time<br>
 *     0.0
 *
 */
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
		p::object data;
		p::object imf;
		p::object res;
		p::object hht;
		p::object insa;
		p::object insf;
		p::object utrgs;
		p::object trgs;

		int numberofimf;
		int data_size;

		FLOAT start_time;
		FLOAT fsr;

		etagen(p::object & _data = none, FLOAT _fsr = 1.0,
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
		void set_data(p::object&);
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
		void set_imf(p::object&);
		//object& get_imf;
		void set_res(p::object&);
		void hilbert_(int _filter_len=MIN_DHT_FILTER_LENGTH, int _stride=MIN_STRIDE);
		void set_hht(p::object&);
		//object& get_hht;
		void set_insa(p::object&);
		//object& get_insa;
		void set_insf(p::object&);
		//object& get_insf;
		list gen_utrgs(FLOAT, int, int);
		np::ndarray get_waveform(int);
		np::ndarray get_waveform(long);
		//object& get_utrgs;
		np::ndarray gen_trgs(np::ndarray, FLOAT, FLOAT, FLOAT);
		//object& get_trgs;
		void show();
};

/** (re)set the parameters for (wS)EMD
 *
 * arguments:<br>
 *     num_imfs     number of IMFs to be obtained<br>
 *     num_sifts    maximum number of siftings to obtain each IMF<br>
 *     S_number     criterion to stop sifting,<br>
 *                  (number of extrema) - (number of zero-crossings) <= S_number<br>
 *     emd_size     (for wSEMD) number of data samples in each EMD segments<br>
 *     num_seg      (for wSEMD) number of segments to average out the IMFs<br>
 *     w_type       (for wSEMD) type of weighting kernel<br>
 *                  (see docstring of etagen.kernel)<br>
 *
 */
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

/** show the parameters for (wS)EMD
 *
 */
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

/** (re)set the input data
 *
 * example:<br>
 *     >>> import numpy as np
 *
 *     >>> from etagen import etagen
 *
 *     >>> h = etagen(np.random.randn(1024))
 *
 *     >>> print h.data.shape<br>
 *     (1024,)
 *
 *     >>> h.set_data(np.random.randn(32))
 *
 *     >>> print h.data.shape<br>
 *     (32,)
 *
 */
void etagen::set_data(p::object & a)
{
	if (PyArray_NDIM(a) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	data = a;
	data_size = PyArray_Size(a);
}

/** set or replace IMFs
 *
 */
void etagen::set_imf(object & a)
{
	if (PyArray_NDIM(a) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"imf should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a)[1] != data_size)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of imf and data does not match.");
		throw_error_already_set();
	}
	imf = a;
	numberofimf = PyArray_DIMS(a)[0];
	imfbuf_size = data_size * numberofimf;
}

/** (re)set the residual
 *
 */
void etagen::set_res(object & a)
{
	if (PyArray_NDIM(a) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"res should be a 1D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_Size(a) != data_size)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of res and data does not match.");
		throw_error_already_set();
	}
	res = a;
}

/** (re)set the Hilbert-transformed IMFs
 *
 */
void etagen::set_hht(object & a)
{
	if (PyArray_NDIM(a) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"hht should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a)[1] != PyArray_DIMS(data)[0])
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of hht and data does not match.");
		throw_error_already_set();
	}
	hht = a;
}

/** (re)set the instantaneous amplitudes of IMFs
 *
 */
void etagen::set_insa(object & a)
{
	if (PyArray_NDIM(a) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"insa should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a)[1] != PyArray_DIMS(data)[0])
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of insa and data does not match.");
		throw_error_already_set();
	}
	insa = a;
}

/** (re)set the instantaneous frequencies of IMFs
 *
 */
void etagen::set_insf(object & a)
{
	if (PyArray_NDIM(a) != 2)
	{
		PyErr_SetString(PyExc_ValueError,
				"insf should be a 2D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_DIMS(a)[1] != PyArray_DIMS(data)[0]-1)
	{
		PyErr_SetString(PyExc_ValueError,
				"Lengths of insf and data does not match.");
		throw_error_already_set();
	}
	insf = a;
}

/** decompose self.data into IMFs by EMD<br>
 * (parameters are set by self.set_emd_parameter)
 *
 */
void etagen::emd_()
{
	FLOAT* _time_series = PyArray_DATA(data);
	FLOAT* _imf_series = (FLOAT*)malloc(sizeof(FLOAT) * data_size * \
					    numberofimf);
	emd(data_size, _time_series, _imf_series, numberofimf,
	    MAX_ITERATION, monotonic, s_number); 
	imf = (object)DoubleToNumpyArray(numberofimf, data_size, _imf_series);
	imfbuf_size = data_size * numberofimf;
}

/** decompose self.data into IMFs by wSEMD<br>
 * (parameters are set by self.set_emd_parameter)
 *
 */
void etagen::wsemd_()
{
	int emd_size = numberofbuf * bufoffset;	//emd slide size

	FLOAT* _time_series = PyArray_DATA(data);
	FLOAT* _imf_series = (FLOAT*)malloc(sizeof(FLOAT) * data_size * \
					    numberofimf);

	wsemd(data_size, emd_size, numberofbuf, weight,
	      _time_series, _imf_series, numberofimf,
	      MAX_ITERATION, monotonic, s_number); 
	imf = (object)DoubleToNumpyArray(numberofimf, data_size, _imf_series);
	imfbuf_size = data_size * numberofimf;
}

/** Hilbert spectral analysis
 *
 * arguments:<br>
 *     filter_length    length of FIR filter for DHT<br>
 *     stride           number of samples to estimate instantaneous frequency
 *
 * attribute stored:<br>
 *     h.hht            Hilbert-transformed IMFs<br>
 *     h.insa           instantaneous amplitudes of IMFs<br>
 *     h.insf           instantaneous frequnecies of IMFs
 *
 */
void etagen::hilbert_(int _filter_len, int _stride)
{
	Py_intptr_t const* shape = PyArray_DIMS(imf);
	Py_intptr_t fshape[2] = {shape[0], shape[1] - 1};
	
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

	FLOAT* _imf_series = PyArray_DATA(imf);
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
	//free(_dht_result);
	//free(_hht_series);
	//free(_dht_amplitude);
	//free(_dht_frequency);
}

/** METHOD FOR THE INTERNAL USE
 *
 */
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
	FLOAT* _time_series = PyArray_DATA(data);
	FLOAT* _imf_series = PyArray_DATA(imf);
	FLOAT* _dht_amplitude = PyArray_DATA(insa);
	FLOAT* _dht_frequency = PyArray_DATA(insf);

	for(int i=0; i < numberofimf; i++)
	{
		ret.append(trigger_gen(i, _imf_series + data_size*i + sidx, len, _dht_amplitude + data_size*i + sidx, _dht_frequency + (data_size-1)*i + sidx, med_abs_dev(data_size, _time_series), snr_th));
		//ret.append(trigger_gen(i, _imf_series + data_size*i + sidx, len, _dht_amplitude + data_size*i + sidx, _dht_frequency + (data_size-1)*i + sidx, median(data_size, abs(data_size, _time_series)), snr_th));
	}
	return ret;
}

/** METHOD FOR THE INTERNAL USE
 *
 */
np::ndarray etagen::gen_trgs(np::ndarray _utrgs, FLOAT snr_th, FLOAT ttol, FLOAT ftol)
{
	int _ntrgs = PyArray_Size(_utrgs);
	if (_ntrgs == 0) return DoubleToNumpyArray(0, cltinfo_num, (FLOAT*)0);
	np::ndarray pyarr =  _utrgs.copy();
	trginfo *_trg = (trginfo*) pyarr.get_data();
	cltinfo *clt;
	Py_intptr_t const* shape = PyArray_DIMS(imf);
	FLOAT* _imf_series = PyArray_DATA(imf);
	FLOAT* _time_series = PyArray_DATA(data);
	
	tc.set_param(ttol, ftol, med_abs_dev(data_size, _time_series), _imf_series, shape[0], shape[1], start_time, fsr);
	for(int i=0; i < _ntrgs; i++)
	{
		trginfo *tmp = new trginfo[1];
		tmp->id = _trg[i].id;
		tmp->start_index = _trg[i].start_index;
		tmp->end_index = _trg[i].end_index;
		tmp->peak_index = _trg[i].peak_index;
		tmp->amplitude = _trg[i].amplitude;
		tmp->frequency = _trg[i].frequency;
		tmp->fmin = _trg[i].fmin;
		tmp->fmax = _trg[i].fmax;
		tmp->snr = _trg[i].snr;
		tc.feed(tmp);
	}
	clt = tc.getClusteredTrigger(snr_th);
	return DoubleToNumpyArray(tc.numofCluster, cltinfo_num, (FLOAT*)clt);
}

/** reconstruct the waveform of the *indx*-th trigger
 *
 * arguments:<br>
 *     indx    the index of the trigger to reconstruct the waveform
 *
 * return type:<br>
 *     Numpy ndarray type
 *
 */
np::ndarray etagen::get_waveform(int index)
{
	long indx = (long) index;

	return get_waveform(indx);
}

np::ndarray etagen::get_waveform(long index)
{
	int len;
	FLOAT *wave = NULL;
	try
	{
		len = tc.getWaveform(index, &wave);
	} catch (std::bad_alloc& e) {
		len = tc.getWaveform(index, &wave);
	}
	return DoubleToNumpyArray(len, wave);
}
void etagen::show()
{
	tc.show();
}

object extrema_(object & a)
{
	if (PyArray_NDIM(a) != 1)
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	FLOAT* data = PyArray_DATA(a);
	int n = PyArray_Size(a);
	std::vector<FLOAT> x_peak, y_peak, x_bot, y_bot;
	extrema(data, n, x_peak, y_peak, x_bot, y_bot);
	return make_tuple(stdVecToNumpyArray(x_peak),
			  stdVecToNumpyArray(y_peak),
			  stdVecToNumpyArray(x_bot),
			  stdVecToNumpyArray(y_bot));
}

object spline_(object & _x, object & _y, bool mono=false)
{
	if ((PyArray_NDIM(_x) != 1) || (PyArray_NDIM(_y) != 1))
	{
		PyErr_SetString(PyExc_ValueError,
				"data should be a 1D arraylike object.");
		throw_error_already_set();
	}
	if (PyArray_Size(_x) != PyArray_Size(_y))
	{
		PyErr_SetString(PyExc_ValueError,
				"x and y should have same size.");
		throw_error_already_set();
	}
	FLOAT* xdata = PyArray_DATA(_x);
	FLOAT* ydata = PyArray_DATA(_y);
	int n = PyArray_Size(_x);
	std::vector<FLOAT> x_pts (xdata, xdata+n);
	std::vector<FLOAT> y_pts (ydata, ydata+n);
	std::vector<FLOAT> x_data, y_data;
	mono ? monotonic_cubic_Hermite_spline(&x_pts, &y_pts, &x_data, &y_data)
	     : cubic_spline(&x_pts, &y_pts, &x_data, &y_data);
	return (object) stdVecToNumpyArray(y_data);
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_waveform_overloads, get_waveform, 1, 1)

np::ndarray (etagen::*get_waveform_int)(int)   = &etagen::get_waveform;
np::ndarray (etagen::*get_waveform_long)(long) = &etagen::get_waveform;

BOOST_PYTHON_MODULE(_etagen)
{
	np::initialize();
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
			   "construct an EtaGen instance\n"
			   "USAGE (an example):\n"
			   "    h = etagen.etagen(data, fsr, start_time)\n"
			   "arguments:\n"
			   "    data          input data which is assumed to be whitened\n"
			   "    fsr           sampling frequency\n"
			   "    start_time    start time of data in sec.\n"
			   "return:\n"
			   "    h             an instance of etagen class\n"
			   "attribute stored:\n"
			   "    h.data\n"
			   "    h.fsr\n"
			   "    h.start_time\n"
			   "example:\n"
			   ">>> import numpy as np\n"
			   ">>> from etagen import etagen\n"
			   ">>> h = etagen(np.random.randn(1024),32)\n"
			   ">>> print h.fsr\n"
			   "32.0\n"
			   ">>> print h.start_time\n"
			   "0.0\n"
			   )
	    )
	.def_readonly("start_time", &etagen::start_time,
	    "start time of input data\n"
	    )
	.def_readonly("fsr",        &etagen::fsr,
	    "sampling rate\n"
	    )
	.def_readonly("data_len",       &etagen::data_size,
	    "Length of the input data\n"
	    )
	.def_readonly("nimfs",       &etagen::numberofimf,
	    "Number of IMFs\n"
	    )
	.def_readonly("data",       &etagen::data,
	    "input data\n"
	    )
	.def_readonly("imfs",       &etagen::imf,
	    "IMF time series (EMD)\n"
	    )
	.def_readonly("res",        &etagen::res,
	    "residual time series\n"
	    )
	.def_readonly("hht",        &etagen::hht,
	    "HHT time series (EMD + HSA)\n"
	    )
	.def_readonly("insa",       &etagen::insa,
	    "instantaneous amplitude of each IMFs\n"
	    )
	.def_readonly("insf",       &etagen::insf,
	    "instantaneous frequency of each IMFs\n"
	    )
	.def_readwrite("utrgs",      &etagen::utrgs,
	    "preliminary (unclustered) triggers\n"
	    )
	.def_readwrite("trgs",       &etagen::trgs,
	    "(clustered) triggers\n"
	    )
	.def("set_start_time",      &etagen::set_start_time,
	    "(re)set the start time of self.data\n"
	    "example:\n"
	    ">>> import numpy as np\n"
	    ">>> from etagen import etagen\n"
	    ">>> h = etagen(np.random.randn(1024))\n"
	    ">>> print h.start_time\n"
	    "0.0\n"
	    ">>> h.set_start_time(10)\n"
	    ">>> print h.start_time\n"
	    "10.0\n"
	    )
	.def("set_fsr",             &etagen::set_fsr,
	    "(re)set the sampling frequency of self.data\n"
	    "example:\n"
	    ">>> import numpy as np\n"
	    ">>> from etagen import etagen\n"
	    ">>> h = etagen(np.random.randn(1024))\n"
	    ">>> print h.fsr\n"
	    "1.0\n"
	    ">>> h.set_fsr(32)\n"
	    ">>> print h.fsr\n"
	    "32.0\n"
	    )
	.def("set_emd_param",       &etagen::set_emd_param,
	    (arg("num_imfs")=0, arg("num_sifts")=MAX_SIFT,
	       arg("S_number")=0, arg("monotonic_spline")=0,
	       arg("emd_size")=2048, arg("num_seg")=8,
	       arg("w_type")=sin_kernel),
	    "(re)set the parameters for (wS)EMD\n"
	    "arguments:\n"
	    "    num_imfs     number of IMFs to be obtained\n"
	    "    num_sifts    maximum number of siftings to obtain each IMF\n"
	    "    S_number     criterion to stop sifting,\n"
	    "                 (number of extrema) - (number of zero-crossings) <= S_number\n"
	    "    emd_size     (for wSEMD) number of data samples in each EMD segments\n"
	    "    num_seg      (for wSEMD) number of segments to average out the IMFs\n"
	    "    w_type       (for wSEMD) type of weighting kernel\n"
	    "                 (see docstring of etagen.kernel)\n"
	    )
	.def("show_emd_param",      &etagen::show_emd_param,
	    "show the parameters for (wS)EMD\n"
	    )
	.def("set_data",            &etagen::set_data,
	    "(re)set the input data"
	    "example:\n"
	    ">>> import numpy as np\n"
	    ">>> from etagen import etagen\n"
	    ">>> h = etagen(np.random.randn(1024))\n"
	    ">>> print h.data.shape\n"
	    "(1024,)\n"
	    ">>> h.set_data(np.random.randn(32))\n"
	    ">>> print h.data.shape\n"
	    "(32,)\n"
	    )
	//.def("add_data",            &etagen::add_data,
	//    ""
	//    )
	//.def("get_data",            &etagen::get_data,
	//    ""
	//    )
	.def("emd",                 &etagen::emd_,
	    "decompose self.data into IMFs by EMD\n"
	    "(parameters are set by self.set_emd_parameter)\n"
	    )
	.def("wsemd",               &etagen::wsemd_,
	    "decompose self.data into IMFs by wSEMD\n"
	    "(parameters are set by self.set_emd_parameter)\n"
	    )
	.def("set_imf",             &etagen::set_imf,
	    "set or replace IMFs\n"
	    )
	//.def("get_imf",             &etagen::get_imf,
	//    ""
	//    )
	.def("set_res",             &etagen::set_res,
	    "(re)set the residual\n"
	    )
	.def("hilbert",             &etagen::hilbert_,
	    (arg("filter_length")=MIN_DHT_FILTER_LENGTH, arg("stride")=MIN_STRIDE),
	    "Hilbert spectral analysis\n"
	    "arguments:\n"
	    "    filter_length    length of FIR filter for DHT\n"
	    "    stride           number of samples to estimate instantaneous frequency\n"
	    "attribute stored:\n"
	    "    h.hht            Hilbert-transformed IMFs\n"
	    "    h.insa           instantaneous amplitudes of IMFs\n"
	    "    h.insf           instantaneous frequnecies of IMFs\n"
	    )
	.def("set_hht",             &etagen::set_hht,
	    "(re)set the Hilbert-transformed IMFs\n"
	    )
	//.def("get_hht",             &etagen::get_hht,
	//    ""
	//    )
	.def("set_insa",            &etagen::set_insa,
	    "(re)set the instantaneous amplitudes of IMFs\n"
	    )
	//.def("get_insa",            &etagen::get_insa,
	//    ""
	//    )
	.def("set_insf",            &etagen::set_insf,
	    "(re)set the instantaneous frequencies of IMFs\n"
	    )
	//.def("get_insf",            &etagen::get_insf,
	//    ""
	//    )
	.def("gen_utrgs",           &etagen::gen_utrgs,
	    (arg("snr_threshold")=DEFAULT_SNR_THRESHOLD, arg("start")=0, arg("length")=0),
	    "METHOD FOR THE INTERNAL USE\n"
	    )
	//.def("get_utrgs",           &etagen::get_utrgs,
	//    ""
	//    )
	.def("gen_trgs",            &etagen::gen_trgs,
	    (arg("triggers"), arg("snr_threshold")=DEFAULT_SNR_THRESHOLD, arg("t_tolerance")=0, arg("f_tolerance")=0),
	    "METHOD FOR THE INTERNAL USE\n"
	    )
	.def("get_waveform",            get_waveform_int,
	   get_waveform_overloads(
	   "reconstruct the waveform of the *indx*-th trigger\n"
	   "arguments:\n"
	   "    indx    the index of the trigger to reconstruct the waveform\n"
	   "return type:\n"
	   "    Numpy ndarray type"
	   ))
	.def("get_waveform",            get_waveform_long,
	   "reconstruct the waveform of the *indx*-th trigger\n"
	   "arguments:\n"
	   "    indx    the index of the trigger to reconstruct the waveform\n"
	   "return type:\n"
	   "    Numpy ndarray type"
	   )
	.def("show",             &etagen::show,
			""
		)
	;
}
