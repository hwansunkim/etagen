#!/usr/bin/env python

from __future__ import division, print_function
from etagen import etagen, kernel
import numpy as np
from time import time
# to check run time
t0 = time()

dur   = 32                      # Duration of data in second
fsr   = 16 * 1024               # Sampling frequency in Hz

nsift = 15                      # Maximum number of sifting
snum  = 0                       # S-number (a stop criterion)
mono  = 0                       # Set 1 to use monotonic cubic spline
wtype = kernel.EXP_KERNEL       # Kernel type for wSEMD
nblks = 8                       # Number of segments for wSEMD
bsize = 128                     # segment size = EMD size / Number of segments
esize = nblks * bsize           # EMD size for wSEMD
nimfs = min(8, int(np.log2(esize)-1))   # Number of IMFs

fleng = 128                     # Filter length for DHT
hsast = 1024                    # Length of segment used to calculate frequency

min_snr = 5                     # Minimum SNR to generate trigger
trgst   = 16 * fsr              # Length of segment used to generate triggers
trgov   = 8 * fsr               # Overlap between neighbor segments

snr_th  = 5.5                   # SNR theshold for trigger clusters
t_tol   = 0.001                 # Time tolerance [s] to cluster triggers
f_tol   = 0.5                   # Frequency tolerance (ratio)

# set random seed
np.random.seed(0)

# make etagen instance with Gaussian random data
t1 = time()
h = etagen(np.random.randn(dur * fsr), fsr)
print("{0}s elapsed to generate data".format(time() - t1))

# set wSEMD parameters
t1 = time()
h.set_emd_param(nimfs, nsift, snum, mono, esize, nblks, wtype)

# get IMFs by wSEMD
h.wsemd()
print("{0}s elapsed to get {1} IMFs".format(time() - t1, h.nimfs))

# run HSA
t1 = time()
h.hilbert(filter_length=fleng, stride=hsast)
print("{0}s elapsed to get instantaneous amplitude and frequency".format(time() - t1))

# generate unclustered triggers
t1 = time()
h.get_utriggers(snr_th=min_snr, stride=trgst, overlap=trgov)
print("{0}s elapsed to generate {1} triggers".format(time() - t1, h.utrgs.shape[0]))
ntrg = min(10, h.utrgs.shape[0])
print("Top {0} SNR of triggers: {1}".format(ntrg, h.utrgs[h.utrgs['snr'].argsort()[::-1]][:ntrg]['snr']))

# cluster triggers
t1 = time()
h.get_triggers(snr_th=snr_th, t_tolerance=t_tol, f_tolerance=f_tol)
print("{0}s elapsed to get {1} trigger clusters".format(time() - t1, h.trgs.shape[0]))
ntrg = min(10, h.trgs.shape[0])
print("Top {0} SNR of trigger clusters: {1}".format(ntrg, h.trgs[h.trgs['snr'].argsort()[::-1]][:ntrg]['snr']))

print("{0}s elapsed in total".format(time() - t0))

for i in range(100):
	wave = h.get_waveform(i)
	x = np.arange(wave.shape[0])

	from pylab import *
	plt.plot(x,wave)
	plt.show()


