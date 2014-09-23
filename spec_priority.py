def gauss_function(x, a, x0, sigma):

    return a*np.exp(-(x-x0)**2/(2*sigma**2))




def find_mag_weight(lim_mag, current_r):

   weight = np.minimum(np.maximum(1 - (lim_mag - current_r) * 0.3, 1e-3), 1)
   if (current_r > lim_mag):
      weight = 1e-5

   return weight




def find_phase_weight(phase_lim, phase, phase_err):

   if (phase - phase_err < phase_lim):
       #To account for the difference between when the program is run and when the spectrum is taken:
       phase = phase + 0.5
       phase_err = (phase_err**2 + 0.5**2)**0.5
   
       phase_arr = np.arange(-20, 20, 0.1)
   
       weight_func = np.zeros(len(phase_arr))
   
       index = np.where(phase_arr < -2)
       weight_func[index] = np.maximum(1 - 0.007 * (phase_arr[index] - 2)**2, 1e-3)
   
       index = np.where(phase_arr > 0)
       weight_func[index] = np.maximum(1 - 0.01 * phase_arr[index]**2, 1e-3)
   
       index = np.where((phase_arr >= -2) & (phase_arr <= 0))
       weight_func[index] = 1
   
       index = np.where(phase_arr >= phase_lim)
       weight_func[index] = 1e-3

       phase_dist = gauss_function(phase_arr, 1, phase, phase_err)
       phase_dist = phase_dist / phase_dist.sum()

       weight = (phase_dist * weight_func).sum()

   else:
       weight = 1e-3
   
   return weight
   



def find_z_weight(photoz, photoz_err, dz):

   n_params = 3
   
   n_bin = 1 / dz
   
   weight = np.zeros((n_bin, 2, 2, 2))

   zmin = np.arange(0, 1, dz)
   zmax = np.arange(0+dz, 1+dz, dz)
   
   n = 501

   z_arr = np.arange(0, 1.25, 0.0025)
   
   z_dist = gauss_function(z_arr, 1, photoz, photoz_err)
   z_dist = z_dist / z_dist.sum()
   
   for i in np.arange(0, n_bin):
      index = np.where((z_arr >= zmin[i]) & (z_arr < zmax[i]))
      weight[i,:,:,:] = (z_dist[index]).sum()
   
   return weight
   



def find_host_flag(host_split, host_mass):

   n_params = 2
   
   absurd = 20
   
   if (host_mass > absurd):   
      host_flag = 0
   if (host_mass < host_split):   
      host_flag = 1
   if ((host_mass >= host_split) & (host_mass < absurd)):   
      host_flag = 2
   
   return host_flag
   



def find_priority(peak_r, current_r, ia_prob, photoz, photoz_err, phase, phase_err, phase_first, host_mass, lim_mag_arr, spec_dist):

#peak_r is the peak magnitude in the r band, either measured or
#expected.

#current_r is the current magnitude in the r band, hopefully from a
#model fit.

#ia_prob is the Ia probability from exterior sources.

#photoz is the host-galaxy photo-z

#photoz_err is the photo-z error.

#phase is the model rest-frame phase relative to maximum

#phase_err is the uncertainty in the phase


#spec_dist is a 6x2x2x2 array to indicate the number of SNe in each
#redshift bin with each of three properties, dm15, velocity, and host
#mass.


#Now the more complicated algorithm to get the probability
   dz = 0.2
   
   n_bin = 1 / dz
   tot_sn = 160
   
   ideal_dist = np.zeros((n_bin, 2, 2, 2)) + tot_sn / 8 / float(n_bin)

   
   
   #Define different weights.  These will need to be updated for
   #different telescopes, conditions, etc
   phase_lim = 10
   phase_weight = find_phase_weight(phase_lim, phase, phase_err)
   
   dz = 0.2
   z_weight = find_z_weight(photoz, photoz_err, dz)
   
   host_split = 10
   host_flag = find_host_flag(host_split, host_mass)
   
   diff_dist = np.maximum((ideal_dist - spec_dist) * z_weight, 1e-3)
   if (host_flag == 0):   
      diff_dist = diff_dist / 2.
   if (host_flag == 1):   
      diff_dist[0,:,:,:] = 0
   if (host_flag == 2):   
      diff_dist[1,:,:,:] = 0
   
   dist_weight = diff_dist.sum()

   if (phase_first <= 0):
       first_weight = 1
   else:
       first_weight = 1e-5

   priorities = np.zeros([len(lim_mag_arr)])
   for i in np.arange(0, len(lim_mag_arr)):
      mag_weight = find_mag_weight(lim_mag_arr[i], current_r)

      priority = mag_weight * phase_weight * first_weight * dist_weight * ia_prob
      
      #Do the magnitude complete sample
      mag_complete = 20.5
      if ((peak_r < mag_complete) & (current_r <= lim_mag_arr[i])) :   
         priority = priority + 10
      
      
      #Do the volume complete sample
      z_complete = 0.20
      if ((photoz <= z_complete) or ((photoz - photoz_err <= z_complete) & (photoz_err < 0.2)) & (current_r <= lim_mag_arr[i])):   
         priority = priority + 10
      
      
      #Add another bump for really bright things
      mag_bright = 19
      if (current_r < mag_bright):   
         priority = priority + 10
      
      if ((current_r > lim_mag_arr[i]) & (priority > 1e-15)):
          priority = 1e-15
          
      priorities[i] = priority
   
   return priorities

   


def create_spec_dist(file):

   dm15_split = 1.1
   vel_split = 11.8
   host_split = 10
   
   data = open(file, "r")
   lines = data.readlines()
   data.close

   z = []
   dm15 = []
   vel = []
   host_mass = []

   for line in lines:
      p = line.split()
      z.append(float(p[1]))
      dm15.append(float(p[2]))
      vel.append(float(p[3]))
      host_mass.append(float(p[4]))

   dz = 0.2
   n_bin = 1 / dz
   zmin = np.arange(0, 1, dz)
   zmax = np.arange(0+dz, 1+dz, dz)

   spec_dist = np.zeros((n_bin, 2, 2, 2))

   for i in np.arange(0, n_bin):
      spec_dist[i,0,0,0] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 < dm15_split) & (vel < vel_split) & (host_mass < host_split))[0])
      spec_dist[i,0,0,1] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 < dm15_split) & (vel < vel_split) & (host_mass >= host_split))[0])
      spec_dist[i,0,1,0] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 < dm15_split) & (vel >= vel_split) & (host_mass < host_split))[0])
      spec_dist[i,0,1,1] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 < dm15_split) & (vel >= vel_split) & (host_mass >= host_split))[0])
      spec_dist[i,1,0,0] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 >= dm15_split) & (vel < vel_split) & (host_mass < host_split))[0])
      spec_dist[i,1,0,1] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 >= dm15_split) & (vel < vel_split) & (host_mass >= host_split))[0])
      spec_dist[i,1,1,0] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 >= dm15_split) & (vel >= vel_split) & (host_mass < host_split))[0])
      spec_dist[i,1,1,1] = len(np.where((z >= zmin[i]) & (z < zmax[i]) & (dm15 >= dm15_split) & (vel >= vel_split) & (host_mass >= host_split))[0])

   return spec_dist




def find_first(mjd, flux, flux_err, dataflag):

    good1 = np.where((dataflag == 1) & (flux_err > 0))
    mjd1 = mjd[good1]
    flux1 = flux[good1]
    flux_err1 = flux_err[good1]

    good2 = np.where(flux1/flux_err1 >= 3)

    tfirst = min(mjd1[good2])

    return tfirst




def parse_file(file):

   current_time = Time.now()
   current_mjd = current_time.mjd
      
   #Read in data
   data = open(file, "r")
   lines = data.readlines()
   data.close

   name = lines[0].split()[2]
   cid = int( lines[1].split()[2] )
   ia_prob = float( lines[10].split()[2] )
   tmax = float( lines[6].split()[2] )
   tmax_err = float( lines[6].split()[4] )
   if (abs(tmax_err) > 10):
       tmax_err = 3

   #Get spectroscopic redshift if available; if not, get the photo z
   z = float( lines[4].split()[2] )
   if (z >= 0):
       z_err = float( lines[4].split()[4] )
   else:
       z = float( lines[5].split()[2] )
       z_err = float( lines[5].split()[4] )

   if (z_err <= 0):
       z_err = 0.2
       
   mjd = []
   flux = []
   flux_err = []
   dataflag = []
   band = []
   ifit = []

   for line in lines[15:]:
       p = line.split()
       mjd.append(float(p[1]))
       flux.append(float(p[3]))
       flux_err.append(float(p[4]))
       dataflag.append(int(p[5]))
       band.append(p[6])
       ifit.append(int(p[8]))

   mjd = np.array(mjd)
   flux = np.array(flux)
   flux_err = np.array(flux_err)
   dataflag = np.array(dataflag)
   band = np.array(band)
   ifit = np.array(ifit)
   

   #Hack to run on old data.  Remove if running on current date.
   #current_mjd = tmax + 5
   ##############

   #Determine the current phase
   current_phase = (current_mjd - tmax) / (1. + z)

   #Determine the phase for each point
   phase = (mjd - tmax) / (1. + z)
   
   #Determine MJD and phase of first S/N > 3 data point
   t_first = find_first(mjd, flux, flux_err, dataflag)
   phase_first = (t_first - tmax) / (1. + z)

   #Convert from flux to magnitudes
   mag = flux[:]
   #Need to fix negative flux; then convert
   for element in range(len(flux)):
      if flux[element] < 1e-5:
         flux[element] = 1e-5
      mag[element] =  -2.5 * math.log10(flux[element]) + 27.5

   mag = np.array(mag)
   
   #Determine the peak and current r mag
   r_model_index = np.where((dataflag == 0) & (band == 'r') & (ifit == 1))
   r_model_mag = mag[r_model_index]

   current_phase_r_index = np.where(abs(phase[r_model_index] - current_phase) == min(abs(phase[r_model_index] - current_phase)))
   peak_r_index = np.where(abs(phase[r_model_index]) == min(abs(phase[r_model_index])))

   current_r = np.mean(r_model_mag[current_phase_r_index])
   peak_r = np.mean(r_model_mag[peak_r_index])

   return name, cid, peak_r, current_r, ia_prob, z, z_err, current_phase, tmax_err, current_mjd, phase_first




import sys
import math
import json
import numpy as np
from numpy import random
from astropy.time import Time
import re

fitlc_output = sys.argv[1]

name, cid, peak_r, current_r, ia_prob, z, z_err, phase, phase_err, current_mjd, phase_first = parse_file(fitlc_output)

minlim = 20
maxlim = 25
dlim = 0.5

lim_mag_arr = np.arange(minlim, maxlim+dlim, dlim)

host_mass = 50

spec_dist = create_spec_dist('des_spec')

priorities = find_priority(peak_r, current_r, ia_prob, z, z_err, phase, phase_err, phase_first, host_mass, lim_mag_arr, spec_dist)

doc = {"version"        :   "20140920",
       "mag_limits"     :   list( lim_mag_arr ),
       "candidate_name" :   name        ,
       "candidate_id"   :   int( cid )  ,
       "ia_prob"        :   ia_prob     ,
       "current_r"      :   current_r   ,
       "peak_r"         :   peak_r      ,
       "redshift"       :   z           ,
       "redshift_err"   :   z_err       ,
       "phase"          :   phase       ,
       "mjd"            :   current_mjd ,
       "priorities"     :   list( priorities )  } 

print json.dumps( doc )

