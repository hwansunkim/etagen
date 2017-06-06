# EtaGen
An event trigger generator based on Hilbert-Huang Transform<br>
It is a Python module but the core library of EtaGen is built in C/C++

## Sample USAGE:
* This shows how to generate event triggers from *data* by EtaGen<br>
```
    >>> from etagen import etagen, kernel

    >>> h = etagen(data, fsr=1024)

    >>> h.set_emd_param(num_imfs=8, num_sifts=15, S_number=1, emd_size=1024, num_seg=4, w_type=kernel.SIN_KERNEL)

    >>> h.show_emd_param()<br>
    Parameters are set as follows.<br>
    [General]<br>
    Number of IMFs: 8<br>
    Number of siftings:     15<br>
    S-Number:       1<br>
    [wSEMD settings]<br>
    EMD size:       1024<br>
    Number of segments:     4

    >>> h.wsemd()

    >>> h.hilbert(filter_length=128, stride=1024)

    >>> h.get_utriggers(snr_th=5, stride=4*fsr, overlap=2*fsr)<br>
    Generating triggers with 5-snr threshold in segments of length 4096, overlapping 2048 samples and skipping 0 samples from boundaries...<br>
    ... generated 25 trigger event(s)

    >>> h.get_triggers(t_tolerance=0.001, f_tolerance=0.5, snr_th=5.5)<br>
    u_snr_th should be larger than the one used to generate triggers: assuming u_snr_th = 5<br>
    Clustering triggers of u_snr > 5 with time tolerance=0.001, frequency tolerance=0.5 and dropping clusters of snr < 5.5...<br>
    total Clusters : 25 > 5.500000<br>
    total Clusters : 12 > 5.500000<br>
    ... generated 12 trigger cluster(s)

    >>> h.trgs[['c_time','c_freq','p_time','p_freq','npts','snr']]<br>
    array([ (0.041015625, 297.67912076702106, 0.041015625, 297.67912076702106, 1, 5.546819634000152),<br>
       (0.6464843749999999, 344.6262005898434, 0.646484375, 344.6262005898434, 1, 5.797973757254617),<br>
       (1.2744140625, 260.1084208607948, 1.2744140625, 260.1084208607948, 1, 5.546288486345601),<br>
       (1.5664062499999998, 142.61211763115602, 1.56640625, 142.61211763115602, 1, 6.513913449375204),<br>
       (1.8193359374999998, 333.3958253278954, 1.8193359375, 333.39582532789547, 1, 6.632835601773412),<br>
       (2.7138671875, 341.5032411355364, 2.7138671875, 341.5032411355364, 1, 6.609969624271628),<br>
       (4.1552734375, 220.52834605967874, 4.1552734375, 220.5283460596787, 1, 5.8852032771220015),<br>
       (4.493164062499999, 201.227897198888, 4.4931640625, 201.22789719888803, 1, 5.890950833813454),<br>
       (5.766601562500001, 229.12911997930595, 5.7666015625, 229.12911997930595, 1, 5.557695477749056),<br>
       (6.659179687499999, 245.8274380780258, 6.6591796875, 245.8274380780258, 1, 5.8990177018284875),<br>
       (6.328125, 34.86349503817338, 6.328125, 34.86349503817338, 1, 5.834112766561207),<br>
       (6.0, 67.39566743927662, 6.0, 67.39566743927662, 1, 5.9916237288573955)], <br>
      dtype=[('c_time', '<f8'), ('c_freq', '<f8'), ('p_time', '<f8'), ('p_freq', '<f8'), ('npts', '<i8'), ('snr', '<f8')])
```
