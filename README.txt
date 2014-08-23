Version 1.0:

The initial version of the script takes the .FITLC files as input and
outputs a list of CID, Ia_prob, current_r, peak_r, photo z, phase,
followed by an output of priorities for different limiting magnitudes
(from mag 20 to 25 in increments of 0.5).

Example:

>python2.7 spec_priority.py DESY1_FITOPT000_01313859.FITLC_v1 

1313859 1 23.1965040612 22.9791513473 0.5 3
[  1.34753949e-04   1.34753949e-04   1.34753949e-04   1.34753949e-04
   1.34753949e-04   1.34753949e-04   1.34753949e-04   1.22484767e+01
   1.02271674e+01   8.20585817e+00   6.18454893e+00]
