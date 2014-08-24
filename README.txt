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

Version 1.1:

After Rollin's fixes to give reasonable output, Chris updated the data
files.  Now we have everything that we need other than t_max and
t_max_err.  I also have the current_mjd to be set to 5 days past peak
since these data are from last year.



Remaining bugs:

Currently grabbing t_max from the models, which will be okay, but
should be improved to get max in rest-frame B or similar

Need some sort of error on the phase, which is an error on the time of
max.

There's a hack where we add 0.5 days to the current time to determine
the phase.  That's to avoid the lag between when the program is run
and when we get data.  I'm assuming that spectroscopy is about half a
day after the last time things were run.

We're assuming the Bayesian probability, but we haven't tested that.

Need host mass from some file



Wanted features:

Need to address when we don't have dm15 (or velocity?) for some
spectroscopic objects, but we do have a redshift (and
classification).  Probably just need to do something like the
host_flag.

Should allow an arbitrary MJD input as current_mjd.  Then we can
easily run things in the past.

Make some of the hard-coded variables more flexible (total SNe, dz,
phase limit, etc).

Make the variable splits true variables.

