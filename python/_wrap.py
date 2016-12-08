#!/usr/bin/env python2

# Copyright (C) 2016  Edwin J. Son
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import _etagen
from numpy import array, arange, hstack, vstack, loadtxt, median
from matplotlib.pyplot import scatter

class etagen(_etagen.etagen):
	# dtype for trigger informations
	utrg_dtype = [('imf_idx',int), ('s_time',float), ('e_time',float), ('p_time',float), ('p_amp',float), ('p_freq',float), ('f_min',float), ('f_max',float), ('snr',float)]

	trg_dtype = [('s_time',float), ('e_time',float), ('c_time',float), ('c_freq',float), ('c_energy',float), ('p_time',float), ('p_freq',float), ('p_amp',float), ('p_imf_idx',int), ('p_snr',float), ('f_min',float), ('f_max',float), ('npts',int), ('snr',float)]

	def get_utriggers(self, snr_th=3, stride=0, overlap=0, skip=0):
		"""
		self.get_utriggers(snr_th=3, stride=0, overlap=0, skip=0):

		Generate triggers, provided that Intrinsic Mode Functions(IMFs) are
		obtained before self.get_utriggers() is called.

		Required item:

			self.imfs	IMFs obtained by self.get_imfs()

		Input parameters:

			s		A threshold parameter. Default is 3.
			stride		The length of segment in which triggers are generated.
			overlap		The length of overlap of two neighboring segments.
			skip		The length of skipped samples at boundaries of IMFs.
	
		Output:

			self.utrgs	The trigger information generated.
					The information is [IMF index, start time, end time,
					peak time, peak amplitude, peak frequency, snr].
		"""
		if self.data is None:
			raise RuntimeError("data should be loaded before get_utriggers()")
		if self.imfs is None:
			self.wsemd()
		if (self.insa is None) or (self.insf is None):
			self.hilbert()
		if overlap < 0:
			overlap = 0
		if skip < 0:
			skip = 0
		if stride <= 0:
			stride = self.data_len - 2*skip
		self.min_snr = snr_th
		#self.std = self.data.std()
		self.std = median(abs(self.data))
		self.trg_skip = skip
		self.trg_stride = stride
		self.trg_overlap = overlap
		print "Generating triggers with {0}-snr threshold in segments of length {1}, overlapping {2} samples and skipping {3} samples from boundaries...".format(snr_th, stride, overlap, skip)
		# imf_idx, s_time, e_time, p_time, p_amp, p_freq, f_min, f_max, snr
		self.utrgs = array([], dtype=etagen.utrg_dtype)#.reshape(0,1)
		for j in range(skip, self.data_len-skip-overlap, stride-overlap):
			_trgs = self.gen_utrgs(self.min_snr, j, stride)
			if sum([len(i) for i in _trgs]) > 0:
				_trgs = vstack(_trgs)
			else:
				continue
			_trgs = _trgs[_trgs[:,-4] > 0]	  # p_freq > 0
			if _trgs.shape[0] == 0:
				continue
			if j > skip:
				_trgs = _trgs[_trgs[:,3] >= overlap/2]
			if _trgs.shape[0] == 0:
				continue
			if j+stride < self.data_len-skip:
				_trgs = _trgs[_trgs[:,3] < stride - overlap/2]
			_trgs[:,1:4] /= self.fsr
			_trgs[:,1:4] += self.start_time + float(j)/self.fsr
			self.utrgs = hstack((self.utrgs, array(zip(_trgs[:,0].astype(int), _trgs[:,1], _trgs[:,2], _trgs[:,3], _trgs[:,4], _trgs[:,5], _trgs[:,6], _trgs[:,7], _trgs[:,8]), dtype=etagen.utrg_dtype)))
		self.utrgs.sort(order=['imf_idx', 'p_time'])
		print "... generated %i trigger event(s)" % (self.utrgs.shape[0])

	def get_utriggers_from_file(self, trg_file, **kwargs):
		"""
		self.get_utriggers_from_file(trg_file, **kwargs):

		Read triggers from file

		Input parameters:

			trg_file	Filename that contains trigger informations in ASCII format
			**kwargs	kwargs for numpy.loadtxt or numpy.fromfile
	
		Output:

			self.utrgs	The trigger information generated.
					The information is [IMF index, start time, end time,
					peak time, peak amplitude, peak frequency, snr].
		"""
		# imf_idx, s_time, e_time, p_time, p_amp, p_freq, f_min, f_max, snr
		self.utrgs = loadtxt(trg_file, dtype=etagen.utrg_dtype, **kwargs)
		self.min_snr = int(self.utrgs['snr'].min())
		print "... read %i trigger event(s)" % (self.utrgs.shape[0])

	def get_triggers(self, snr_th=10, t_tolerance=0., f_tolerance=0., u_snr_th=3):
		"""
		self.get_triggers(maxDist=2):

		Cluster trigger events generated by self.get_utriggers().

		Required item:

			self.utrgs	triggers generated by self.get_utriggers()

		Input parameters:

			maxDist		maximum distance to cluster tiggers. Default is 2.
			unit_time	unit time for distance. Default is 0.01(=10ms).
			u_snr_th	A threshold parameter for each trigger event. Default is 3.
			snr_th	A threshold parameter for a cluster of triggers.
					Default is 10.
	
		Output:

			self.utrgs	The trigger information generated if not calculated before.
					The information is [IMF index, start time, end time,
					peak time, peak amplitude, peak frequency, snr].
			self.trgs	The trigger cluster information generated.
		"""
		if self.data is None:
			raise RuntimeError("data should be loaded before get_triggers()")
		if self.utrgs is None:
			self.get_utriggers(snr_th=u_snr_th)
		if u_snr_th < self.min_snr:
			print "u_snr_th should be larger than the one used to generate triggers: assuming u_snr_th =", self.min_snr
			self.u_snr_threshold = self.min_snr
		else:
			self.u_snr_threshold = u_snr_th
		self.snr_threshold = snr_th
		#self.waveforms = []
		print "Clustering triggers of u_snr > {0} with time tolerance={1}, frequency tolerance={2} and dropping clusters of snr < {3}...".format(self.u_snr_threshold, t_tolerance, f_tolerance, snr_th)
		# s_time, e_time, c_time, c_freq, c_energy, p_time, p_freq, p_amp, p_imf_idx, p_snr, npts, snr

		carg = self.utrgs['snr']>=self.u_snr_threshold
		utrgs = array(self.utrgs[carg])

		_trgs = self.gen_trgs(utrgs, snr_threshold=self.snr_threshold, t_tolerance=t_tolerance, f_tolerance=f_tolerance)
		self.trgs = array(zip(_trgs[:,0], _trgs[:,1], _trgs[:,4], _trgs[:,6], _trgs[:,5], _trgs[:,7], _trgs[:,9], _trgs[:,8], _trgs[:,10].astype(int), _trgs[:,11], _trgs[:,2], _trgs[:,3], _trgs[:,12].astype(int), _trgs[:,13]), dtype=etagen.trg_dtype)
		print "... generated %i trigger cluster(s)" % (self.trgs.shape[0])

	def get_triggers_from_file(self, trg_file, **kwargs):
		"""
		self.get_triggers_from_file(trg_file, **kwargs):

		Read clustered trigger events from file

		Input parameters:

			trg_file	Filename that contains clustered trigger informations
					in ASCII format
			**kwargs	kwargs for numpy.loadtxt or numpy.fromfile
	
		Output:

			self.trgs	The trigger cluster information generated.
		"""
		# s_time, e_time, c_time, c_freq, c_energy, p_time, p_freq, p_amp, p_imf_idx, p_snr, npts, snr
		self.trgs = loadtxt(trg_file, dtype=etagen.trg_dtype, **kwargs)
		self.snr_threshold = int(self.trgs['snr'].min())
		print "... read %i trigger cluster(s)" % (self.trgs.shape[0])

	def copy(self, start=None, end=None):
		"""
		hhtout = self.copy(start=None, end=None):

		Copy HHT instance [from *start* second to *end* second, if specified]
		*start*/*end* should be greater than self.start_time

		Input parameters:

			start		start time to copy the data
			end		end time to copy the data

		Output parameters:

			hhtout		copied HHT instance
		"""
		if start is None:
			start_idx = 0
		elif start < self.start_time:
			raise ValueError("start should be greater than self.start_time: %f" % (self.start_time))
		else:
			start_idx = int((start-self.start_time) * self.fsr)
		if end is None:
			end_idx = self.data_len
		elif end < self.start_time:
			raise ValueError("end should be greater than self.start_time: %f" % (self.start_time))
		else:
			end_idx = int((end-self.start_time) * self.fsr)
		hhtout = etagen(data=self.data[start_idx:end_idx], start_time=self.start_time+start_idx/float(self.fsr), fsr=self.fsr)
		if self.imfs is not None:
			hhtout.set_imfs(self.imfs[:,start_idx:end_idx])
			hhtout.set_res(self.res[start_idx:end_idx])
			hhtout.nimfs = self.nimfs
		if self.hht is not None:
			hhtout.set_hht(self.hht[:,start_idx:end_idx])
		if self.insa is not None:
			hhtout.set_insa(self.insa[:,start_idx:end_idx])
		if self.insf is not None:
			hhtout.set_insf(self.insf[:,start_idx:end_idx-1])
		end_time = hhtout.start_time + float(hhtout.data_len)/hhtout.fsr
		if self.utrgs is not None:
			hhtout.utrgs  = self.utrgs[(self.utrgs['p_time'] > hhtout.start_time) & (self.utrgs['p_time'] < end_time)]
			hhtout.utrgs[hhtout.utrgs['s_time'] < hhtout.start_time]['s_time'] = hhtout.start_time
			hhtout.utrgs[hhtout.utrgs['e_time'] > end_time]['e_time'] = end_time
			hhtout.min_snr = self.min_snr
			try:
				if start is None and end is None:
					hhtout.trg_skip = self.trg_skip
				hhtout.trg_stride = self.trg_stride
				hhtout.trg_overlap = self.trg_overlap
			except:
				pass
		if self.trgs is not None:
			idx  = where((self.trgs['c_time'] > hhtout.start_time) & (self.trgs['c_time'] < end_time))[0]
			hhtout.trgs  = self.trgs[idx]
			hhtout.trgs[hhtout.trgs['s_time'] < hhtout.start_time]['s_time'] = hhtout.start_time
			hhtout.trgs[hhtout.trgs['e_time'] > end_time]['e_time'] = end_time
			#hhtout.waveforms = []
			#for i in idx:
				#sidx = max(0, int((hhtout.start_time - self.trgs[i]['s_time']) * self.fsr))
				#eidx = len(self.waveforms[i]) - max(0, int((end_time - self.trgs[i]['e_time']) * self.fsr))
				#hhtout.waveforms.append(self.waveforms[i][sidx:-eidx])
			try:
				hhtout.u_snr_threshold = self.u_snr_threshold
			except:
				pass
			hhtout.snr_threshold = self.snr_threshold
		return hhtout

	def plot_tfd(self, indices=None, **kwargs):
		"""
		self.plot_tfd(indices=None, **kwargs):

		Plot time-frequency-distribution (TFD)

		Input parameters:

			indices		List of indices of IMFs used for plotting TFD
					All IMFs will be used by default.
			**kwargs	kwargs for matplotlib.pyplot.scatter
		"""
		if indices is None:
			indices = range(self.nimfs)
			print "plotting TFD for all IMFs"
		else:
			print "plotting TFD for IMF", indices
		timef = self.start_time + arange(1,self.insa.shape[-1]).reshape((1,-1)).repeat(len(indices),axis=0)/float(self.fsr)
		insf = self.insf[indices]
		args = insf > 0
		insa = self.insa[indices,1:][args]
		idx = insa.argsort()
		scatter(timef[args][idx], insf[args][idx], c=insa[idx], **kwargs)
		#for i in id:
			#scatter(timef[self.insf[i]>0], self.insf[i,self.insf[i]>0], c=self.insa[i,self.insf[i]>0], **kwargs)
		print "... done"

	def plot_trgs(self, indices=None, clustered=True, scale=True, **kwargs):
		"""
		self.plot_trgs(indices=None, **kwargs):

		Plot trigger-gram

		Input parameters:

			indices		List of indices of triggers used for plotting trigger-gram
					All triggers will be used by default.
			**kwargs	kwargs for matplotlib.pyplot.scatter
		"""
		if clustered:
			trgs = self.trgs
		else:
			trgs = self.utrgs
		if indices is None:
			print "plotting all triggers"
		else:
			trgs = trgs[indices]
			print "plotting given triggers"
		idx = trgs['snr'].argsort()
		if 's' not in kwargs.keys():
			kwargs['s'] = 20
		if scale:
			kwargs['s'] = trgs[idx]['snr'] * kwargs['s'] / median(trgs[idx]['snr'])
		scatter(trgs[idx]['p_time'], trgs[idx]['p_freq'], c=trgs[idx]['snr'], **kwargs)
		print "... done"
