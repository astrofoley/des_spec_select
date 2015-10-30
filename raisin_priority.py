#!/usr/bin/env python

import  math

from    astropy.time import Time
import  numpy as np


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
        
        index = np.where(phase_arr < -5)
        weight_func[index] = np.maximum(1 - 0.007 * (phase_arr[index] - 2)**2, 1e-3)
        
        index = np.where(phase_arr > 0)
        weight_func[index] = np.maximum(1 - 0.01 * phase_arr[index]**2, 1e-3)
        
        index = np.where((phase_arr >= -5) & (phase_arr <= 0))
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
    
    if (photoz >= 0):
    
        weight = np.zeros(n_bin)
    
        zmin = np.arange(0, 1, dz)
        zmax = np.arange(0+dz, 1+dz, dz)
    
        n = 501
    
        z_arr = np.arange(0, 1.25, 0.0025)
       
        z_dist = gauss_function(z_arr, 1, photoz, photoz_err)
        z_dist = z_dist / z_dist.sum()
    
        for i in np.arange(0, n_bin):
            index = np.where((z_arr >= zmin[i]) & (z_arr < zmax[i]))
            weight[i] = (z_dist[index]).sum()
    
    else:
        weight = np.ones(n_bin)
    
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
   



def find_priority(name, peak_r, current_r, ia_prob, photoz, photoz_err, phase, phase_err, phase_first, host_mass, lim_mag_arr, spec_dist):

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
    dz = 0.01
    
    n_bin = 1 / dz
    tot_sn = 160
    
    ideal_dist = np.zeros((n_bin, 2, 2, 2)) + tot_sn / 8 / float(n_bin)
    
    
    
    #Define different weights.  These will need to be updated for
    #different telescopes, conditions, etc
    phase_lim = 0
    phase_weight = find_phase_weight(phase_lim, phase, phase_err)
    
    dz = 0.01
    z_weight = find_z_weight(photoz, photoz_err, dz)
    
    z_arr = np.arange(0, 1, dz)
    ideal_z = 0.49
    ideal_dz = 0.05
    
    ideal_dist = gauss_function(z_arr, 1, ideal_z, ideal_dz)
    ideal_dist = ideal_dist / ideal_dist.sum()
    redshift_weight = (ideal_dist * z_weight).sum()
    
    if (phase_first <= 0):
        first_weight = 1
    else:
        first_weight = 1e-5
    
    priorities = np.zeros([len(lim_mag_arr)])
    for i in np.arange(0, len(lim_mag_arr)):
        mag_weight = find_mag_weight(lim_mag_arr[i], current_r)
    
        priority = mag_weight * phase_weight * first_weight * redshift_weight * ia_prob
    
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
    if (ia_prob >= 0.5):
        r_model_index = np.where((dataflag == 0) & (band == 'r') & (ifit == 1))
    else:
        r_model_index = np.where((dataflag == 0) & (band == 'r') & (ifit == 3))
 
    r_model_mag = mag[r_model_index]
 
    current_phase_r_index = np.where(abs(phase[r_model_index] - current_phase) == min(abs(phase[r_model_index] - current_phase)))
    peak_r_index = np.where(abs(phase[r_model_index]) == min(abs(phase[r_model_index])))
 
    current_r = np.mean(r_model_mag[current_phase_r_index])
    peak_r = np.mean(r_model_mag[peak_r_index])
 
 
    return name, cid, peak_r, current_r, ia_prob, z, z_err, current_phase, tmax_err, current_mjd, phase_first


def priorities_from_fitlc( fitlc_output ) :

    name, cid, peak_r, current_r, ia_prob, z, z_err, phase, phase_err, current_mjd, phase_first = parse_file(fitlc_output)

    minlim = 20
    maxlim = 25
    dlim = 0.5
    
    lim_mag_arr = np.arange(minlim, maxlim+dlim, dlim)
    
    host_mass = 50
    
    spec_dist = create_spec_dist('des_spec')
    
    priorities = find_priority(name, peak_r, current_r, ia_prob, z, z_err, phase, phase_err, phase_first, host_mass, lim_mag_arr, spec_dist)

    return {
        "target_id" : name,
        "rf_raisin_priority" : {
            "version"       :   "20141029"         ,
            "mag_limits"    :   list( lim_mag_arr ),
            "ia_prob"       :   ia_prob            ,
            "current_r"     :   current_r          ,
            "peak_r"        :   peak_r             ,
            "redshift"      :   z                  ,
            "redshift_err"  :   z_err              ,
            "phase"         :   phase              ,
            "priorities"    :   list( priorities ) ,
        }
    }

if __name__ == "__main__" :

    import  argparse
    import  json
    import  pprint
    import  sys
    
    import  atc_tools
    
    # Parse command line.
    
    parser = argparse.ArgumentParser()
    parser.add_argument( "fitlc_outputs", 
            help  = "FITLC file paths.", 
            nargs = "+" )
    parser.add_argument( "--ignore-errors", "-i", 
            help   = "Ignore FITLC files with errors. Still disalllow inserts though.",
            action = "store_true" )
    parser.add_argument( "--force", "-f",
            help   = "Validate/insert priorities even if there were FITLC file errors.",
            action = "store_true" )
    parser.add_argument( "--json", "-j", 
            help   = "Just compute and output JSON.", 
            action = "store_true" )
    parser.add_argument( "--chunk-size", "-c",
            help    = "Validate/insert this many documents at a time [%(default)s].",
            default = 100,
            type    = int )
    parser.add_argument( "--quiet", "-q", 
            help   = "Output warnings and errors only.",
            action = "store_true" )
    parser.add_argument( "--test", "-t", 
            help   = "Test mode: Server only validates, no inserts.", 
            action = "store_true" )
    args = parser.parse_args()

    # Generate priority documents.

    docs = list()
    for path in args.fitlc_outputs :
        try :
            doc = priorities_from_fitlc( path )
        except :
            if args.ignore_errors :
                print "Error. Skipping:", path
                continue
            raise
        docs.append( doc )

    # Halt if errors were encountered, unless force.

    if len( docs ) != len( args.fitlc_outputs ) :
        message = "Expected {} results, have only {}.".format( len( args.fitlc_outputs ), len( docs ) )
        if args.force :
            print message
        else :
            raise ValueError( message )

    # In JSON mode, output JSON and quit.

    if args.json :
        print json.dumps( docs, sort_keys = True, indent = 4 )
        sys.exit( 0 )

    # Otherwise, validate or insert in chunks.

    if not args.quiet :
        pprint.pprint( docs )

    posts = atc_tools.client().posts
    for start in range( 0, len( docs ), args.chunk_size ) :
        end = start + args.chunk_size
        result = posts.post( docs[ start : end ], test_mode = args.test )
        if not args.quiet :
            print start, "to", end, "of", len( docs )
            pprint.pprint( result )
