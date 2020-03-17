import matplotlib
matplotlib.use('Agg')
import csv
import pandas as pd
import numpy as np
#import filters  # THIS LIBRARY IS NOT INCLUDED IN THE REPOSITORY
import matplotlib.pyplot as pyp
import scipy.signal as signal
import numpy.fft as fft
import os         # ADD THIS SO WE CAN USE ENVIRONMENT VARIABLES
#import matplotlib
from scipy.optimize import curve_fit
import scipy.optimize
from scipy.interpolate import interp1d
import sys
import argparse
import skimage
from skimage.filters import gaussian
#%matplotlib inline


#print os.environ
print os.environ.get('TAU_ACC_ZHAIRES_DIR')
#print os.environ['TAU_ACC_ZHAIRES_DIR']

############################################
# Usage: generate electric fields from ZHAireS Pulses. There are several steps
#	1. Define two environment variables $TAU_ACC_ZHAIRES_DIR = the location of the ZHAireS pulse files in npz format
#					    $TAU_ACC_DETECTOR_PLOTS_DIR = the location to store plots from this analysis
#	2. Make csv file directories in each of Harm's sim folders (Zenith angle / Energy folder (where the out.root files are). ) 
#		using the ~/Dropbox/MountainTop/harms_sims/makeDirs.sh script			
#	3. Convert ZHAireS root files (from Harm) to npz files. This is a two step process.
#		a. In the top level folder where the ZHAireS ROOT files are, start the root interactive shell with the dumpData.C 
#			script, which converts the ROOT files to csv files. #
#			root
#			.L dumpData.C
#			dumpData(35)
#					.... where 35 is the altitude of the detector in km.
#
#		b. In the same top level folder, run python convertCsvToNpz.py, checking first that you are running the correct
#			altitudes and zenith angles
#	3. run this script to generate the lookup tables. They will be stored in TAU_ACC_ZHAIRES_DIR 
#			as interpolator_efields_%dkm.npz'%h, 'w')
	

pastel_colors = ['#FF6666','#FFCC66','#CCFF66','#66FF66','#66FFCC','#66FFFF','#66CCFF','#6666FF','#CC66FF','#FF66FF','#FF6FCF']
bright_colors = ['#FF0000','#FF8000','#FFFF00','#80FF00','#00FF00','#00FF80','#00FFFF','#0080FF','#0000FF','#8000FF','#FF00FF']
bright_colors2 = bright_colors = ['#FF0000','#FF8000','#80FF00','#00FF80','#0080FF','#0000FF','#8000FF','#FF00FF','#FF0080']
pyp.rcParams['font.size']=12
pyp.rcParams['legend.labelspacing'] = 0.1
cmap = matplotlib.cm.get_cmap('inferno')

###########################################################
def read_npz_files(antenna_height, decay, zenith):
    #npzfile = np.load(os.environ['TAU_ACC_ZHAIRES_DIR'] + "/" + "altitude_%dkm_decay_%dkm_zenith_%d.npz"%(antenna_height, decay, zenith))
    npzfile = np.load(os.environ['TAU_ACC_ZHAIRES_DIR'] + "/" + "altitude_%2.1fkm_decay_%2.1fkm_zenith_%d.npz"%(antenna_height, decay, zenith))

    efield_td = npzfile['efield_td'][()] ## for some reason just reading the arrays gives a 0-dim array, so skip this.
    efield_fd = npzfile['efield_fd'][()]
    return efield_td, efield_fd

###########################################################
# WANT TO INTERPOLATE THE EFIELDS AT
# ALL ZENITH ANGLES, PSI ANGLES, and STARTING FREQUENCYS
# FOR 10-MHz SUBBANDS

def band_peak_efield(efield, freq, f_LO, bandwidth):
    df = freq[1]-freq[0]
    cut = np.logical_and(freq>=f_LO, freq<=f_LO + bandwidth)
    # factor of two comes from spectrum
    return 2*np.sum(np.abs(efield[cut]))*df

def efield_vs_angle(h, decay, z, f_Lo_list, bandwidth):

    off_angle_array = np.arange(0.0,80*0.04, 0.04)
    i_off_angle_array = range(0, len(off_angle_array))
    epeak_list = np.zeros((len(f_Lo_list), len(i_off_angle_array)))
    cc = 0
    for i_offangle in i_off_angle_array:
        #if(i_offangle%10==0): print i_offangle
        offangle = off_angle_array[i_offangle]
        efield_td, efield_fd = read_npz_files(h,decay, z)
        efield = np.sqrt(pow(efield_td[i_offangle]['x_v_per_m'],2) + pow(efield_td[i_offangle]['y_v_per_m'],2))
        time = efield_td[i_offangle]['time_s']
        time -= time[0]
        dt = np.abs(time[1]-time[0])
        E_fft = np.fft.rfft(efield)*dt*1.e6 # use dt in \mirco-s so that fft is in units of V/m/MHz
        fr = np.fft.rfftfreq(len(efield), dt*1.e6) # frequencies in MHz
        for i_f_Lo in range(0,len(f_Lo_list)):
            f_Lo = f_Lo_list[i_f_Lo]
            epeak_list[i_f_Lo, cc] = band_peak_efield(E_fft, fr, f_Lo, bandwidth)
        cc+=1
    return off_angle_array, epeak_list

###########################################################
# Construct an array of the peak efield for a range of starting frequencies
# and bandwidths and psi angles
def construct_epeak_array(h):
	zenith_list = np.array([55, 60, 65, 70, 75, 80, 85, 87, 89])
	f_Lo_list = np.arange(10., 1610., 10.)
	if( h < 37): 
		decay_list = np.arange(0., h, 0.5)
		if( h == 0.5):
			decay_list = np.array([0, 0.5])
	else:
		decay_list = np.arange(0,10,1)
	bandwidth = 10.

	###
	# Construct an array of the peak efield for a range of starting frequencies
	# and bandwidths and psi angles
	## Shape of the array is 
	##  [zenith angle [55,65,60,70,75,80,85,87,89 starting_frequency (10-1600 MHz in 10 MHz steps), bandwidth (10 MHz), off_angles (0.04-3.2 in 0.04 deg steps)]

	#epeak_array = []
	#for z in zenith_list:
	#    psi_list, epeak_list = efield_vs_angle(h, decay, z, f_Lo_list, bandwidth)
	#    epeak_array.append(epeak_list)
	#    
	#epeak_array = np.array(epeak_array)

	epeak_array = []
	for d in decay_list:
    		z_epeak_array = []
    		for z in zenith_list:
        		if( d == h):
				psi_list, epeak_list = efield_vs_angle(h, d-0.5, z, f_Lo_list, bandwidth)
			else:
				psi_list, epeak_list = efield_vs_angle(h, d, z, f_Lo_list, bandwidth)
        		z_epeak_array.append(epeak_list)
   		epeak_array.append(z_epeak_array)
	epeak_array = np.array(epeak_array)

	return decay_list, zenith_list, psi_list, f_Lo_list, epeak_array

###########################################################
def epeak_zenith_angle_slice(i_d, i_ze, epeak_array):
    return epeak_array[i_d, i_ze,:, :] 
def epeak_start_freq_slice(i_d,i_f_Lo, epeak_array):
    return epeak_array[i_d, :,i_f_Lo,:] 
def epeak_psi_angle_slice(i_d, i_psi, epeak_array):
    return epeak_array[i_d, :, :, i_psi] 

###########################################################
def  clean_numerical_noise_zenith(h, epeak_array, decay_list, zenith_list, psi_list, f_Lo_list, i_d, i_ze, gauss_blur_sigma = 2.5, ncontours = 10, plot_suffix=""):
	#print "clean_numerical_noise_zenith", h, decay_list[i_d], zenith_list[i_ze], gauss_blur_sigma, ncontours,plot_suffix
	# check if gauss_blur_sigma and ncountors are numbers or arrays
	g = gauss_blur_sigma
	if np.isscalar(gauss_blur_sigma) == False:
	   	print i_ze, gauss_blur_sigma
		g = gauss_blur_sigma[i_ze]
	n = ncontours
	if np.isscalar(ncontours) == False:
	   n = ncontours[i_ze]
	### get the decay
	decay = decay_list[i_d]

	### Use gaussian blur to set mask
	ze = zenith_list[i_ze]
	titles = ["Emergence Angle %d deg."%(90.-int(z)) for z in zenith_list]
	fig = pyp.figure(figsize=(12,8))
	pyp.subplots_adjust(hspace=0.3, wspace=0.3)
	pyp.suptitle(titles[i_ze], fontsize=18)

	# set up the 2-d histogram
	p, start_freq = np.meshgrid(psi_list, f_Lo_list)
	#epeak_array[i_ze,:, :] 
	H = epeak_zenith_angle_slice(i_d, i_ze, epeak_array)
	Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

	# plot the 2-d histogram
	pyp.subplot(2, 2, 1)
	pyp.pcolormesh( p,start_freq, np.log10(Hmasked), cmap='inferno', )
	pyp.xlim(psi_list.min(), psi_list.max())
	pyp.ylim(f_Lo_list.min(), f_Lo_list.max())
	pyp.colorbar(label="log$_{10}$(Peak E-field (V/m))")
	pyp.xlabel(r'$\theta_{view}$ (deg.)')
	pyp.ylabel('Start Frequency (MHz)')
	pyp.grid(True, which='both')
	pyp.title("Original")

	# construct a Gaussian blur
	pyp.subplot(2, 2, 2)
	Hgauss_blur = gaussian(Hmasked, sigma=g)
	pyp.pcolormesh( p,start_freq, np.log10(Hgauss_blur), cmap='inferno', )
	pyp.xlim(psi_list.min(), psi_list.max())
	pyp.ylim(f_Lo_list.min(), f_Lo_list.max())
	pyp.colorbar(label="log$_{10}$(Peak E-field (V/m))")
	pyp.xlabel(r'$\theta_{view}$ (deg.)')
	pyp.ylabel('Start Frequency (MHz)')
	pyp.grid(True, which='both')
	pyp.title("Gaussian Blurred Mask")

	# find the contours:
	CS = pyp.contour(p,start_freq, np.log10(Hgauss_blur),n)

	# Decide which contour to keep
	num_contour_segs = np.array([len(CS.allsegs[i]) for i in range(len(CS.allsegs))])
	len_contour_segs = np.zeros(len(num_contour_segs))
	for i in range(len(CS.allsegs)):
	    for j in range(len(CS.allsegs[i])):
		len_contour_segs[i] += len(CS.allsegs[i][j])

	# if the mean is within 2*rms on both sides then keep up to 2 segments
	# if not keep only one segment
	max_num_contour_segments = 2
	upper_right = np.log10(Hmasked[-20:-1, -20:-1] )
	upper_left = np.log10(Hmasked[-20:-1, 0:20])
	#print "upper right", np.abs(np.mean(upper_right)), np.std(upper_right)
	#print "upper left", np.abs(np.mean(upper_left)), np.std(upper_left)
	#print np.abs(np.mean(upper_left)), " <?", np.abs(np.mean(upper_right)) - 5*np.std(upper_right)

	# if the upper left quadrant is higher than the upper right quadrant
	if np.abs(np.mean(upper_left)) + 2*np.std(upper_left)< np.abs(np.mean(upper_right)) :
	    max_num_contour_segments = 1
	    
	# reject contours that are 50% lower than or 200% higher than the mean 
	# to get rid of small islands and really jagged contours
	i=0
	c=0
	while( (len_contour_segs[c] < 0.5*np.mean(len_contour_segs) or
		len_contour_segs[c] > 2.*np.mean(len_contour_segs)) or
		num_contour_segs[i] > max_num_contour_segments):
	    #print "in while", c, np.where(num_contour_segs < max_num_contour_segments)[0][c]
	    i = np.where(num_contour_segs <= max_num_contour_segments)[0][c]
	    c+=1
	i = np.where(num_contour_segs <= max_num_contour_segments)[0][c]
	level = CS.levels[i]

	# plot the cleaned 2-d histogram
	pyp.subplot(2, 2, 3)

	# Mask pixels less than the average of the noise
	Hmasked2 = np.ma.masked_where(np.log10(Hgauss_blur)<level,Hmasked)
	Hmasked2[Hmasked2.mask] = 0.
	pyp.pcolormesh( p,start_freq, np.log10(Hmasked2), cmap='inferno', )
	pyp.xlim(psi_list.min(), psi_list.max())
	pyp.ylim(f_Lo_list.min(), f_Lo_list.max())
	pyp.colorbar(label="log$_{10}$(Peak E-field (V/m))")
	pyp.xlabel(r'$\theta_{view}$ (deg.)')
	pyp.ylabel('Start Frequency (MHz)')
	pyp.grid(True, which='both')
	pyp.title("Cleaned log scale")

	# Mask pixels less than the average of the noise
	pyp.subplot(2, 2, 4)
	pyp.pcolormesh( p,start_freq, Hmasked2*1e6, cmap='inferno', )
	pyp.xlim(psi_list.min(), psi_list.max())
	pyp.ylim(f_Lo_list.min(), f_Lo_list.max())
	pyp.colorbar(label="Peak E-field ($\mu$V/m)")
	pyp.xlabel(r'$\theta_{view}$ (deg.)')
	pyp.ylabel('Start Frequency (MHz)')
	pyp.grid(True, which='both')
	pyp.title("Cleaned linear scale")
	
	# Save the cleaned image
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_maps_altitude%2.1fkm_decay%2.1fkm_ze%d_cleaning%s.png"
			%(h, decay, int(ze), plot_suffix ) )

	epeak_array[i_d, i_ze,:, :] = Hmasked2
	#print "Shape at end of clean_numerical_noise_zenith: ", np.shape(epeak_array)
	return epeak_array

def clean_numerical_noise(h, epeak_array, decay_list, zenith_list, psi_list, f_Lo_list, gauss_blur_sigma = 2.5, ncontours = 10, plot_suffix=""):	
	#print "clean_numerical_noise",h, gauss_blur_sigma, ncontours
	#epeak_array = []
	#for d in decay_list:
    	#	z_epeak_array = []
    	#	for z in zenith_list:
        #		psi_list, epeak_list = efield_vs_angle(h, d, z, f_Lo_list, bandwidth)
        #		z_epeak_array.append(epeak_list)
   	#	epeak_array.append(z_epeak_array)
	#epeak_array = np.array(epeak_array)
	clean_epeak_array = []
	for i_d in range(len(decay_list)):
		z_epeak_array = []
		# check if gauss_blur_sigma and ncountors are numbers or arrays
		g = gauss_blur_sigma
		if np.isscalar(gauss_blur_sigma) == False:
		   print np.shape(gauss_blur_sigma), i_d, gauss_blur_sigma
		   g = gauss_blur_sigma[i_d]
		n = ncontours
		if np.isscalar(ncontours) == False:
		   n = ncontours[i_d]

		for i_ze in range(len(zenith_list)):

			c_epeak_array = clean_numerical_noise_zenith(h, epeak_array, decay_list, zenith_list, psi_list, f_Lo_list, i_d, i_ze, 
					gauss_blur_sigma=g, ncontours=n,plot_suffix=plot_suffix)
			z_epeak_array.append(c_epeak_array[i_d, i_ze, :, :])
		clean_epeak_array.append(z_epeak_array)
	#print "Shape at end of clean_numerical_noise: ", np.shape(clean_epeak_array)
	return np.array(clean_epeak_array)

###########################################################
def plot_epeak_zenith_psi(h,decay, epeak_array, zenith_list, psi_list, choose_start_freq, i_d, plot_suffix="", log=False):
	### Plot epeak vs zenith angle vs psi angle for fixed starting frequ(ency
	#choose_start_freq = [10., 30., 50., 200., 300., 1000.]
	#ind = [1, 3, 5, 8, 20, 30]
	titles = ["%d MHz"%(freq) for freq in choose_start_freq]

	choose_start_freq = [10., 30., 50., 200., 300., 1000.]
	titles = ["%d MHz"%(freq) for freq in choose_start_freq]
	fig = pyp.figure( figsize=(12,8))
	pyp.subplots_adjust(hspace=0.3, wspace=0.4)

	for i,freq in enumerate(choose_start_freq):
	    i_f_Lo = int(freq/10)
	    ax = pyp.subplot(2,3,i+1)
	    #ax.set_xscale('log')
	    #ax.set_yscale('log')
	    p, zenith = np.meshgrid( psi_list, zenith_list)
	    #H=epeak_array[:, ind[i_f_Lo], 0, :]*1e6
	    #print np.shape(epeak_array), i_d, i_f_Lo
	    H = epeak_start_freq_slice(i_d, i_f_Lo, epeak_array)*1e6
	    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
	    if log:
	    	Hmasked = np.log10(Hmasked)
	    pyp.pcolormesh(p, zenith, Hmasked, cmap='inferno', )
	    pyp.ylim(zenith_list.min(), zenith_list.max()+1)
	    pyp.xlim(psi_list.min(), psi_list.max())
	    if( log):
            	pyp.colorbar(label="log(Peak E-field ($\mu V/m$))")
            else:
            	pyp.colorbar(label="Peak E-field ($\mu V/m$)")
	    pyp.ylabel('Zenith Angle (deg.)')
	    pyp.xlabel('$\psi$ (deg.)')
	    pyp.grid(True, which='both')
	    pyp.title(titles[i])
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_maps_altitude%2.1fkm_decay_%2.1fkm_ze_vs_psi%s.png"%(h,decay, plot_suffix) )

###########################################################
def plot_epeak_zenith_start_freq(h,decay, epeak_array, zenith_list, f_Lo_list, choose_psi, i_d, plot_suffix="", log=False):
	### Second plot will be epeak vs zenith angle vs starting frequency for fixed psi angle
	ind = [0, 11, 27, 35, 47, 63]
	#choose_psi = psi_list[ind]
	titles = ["$\psi$=%2.2f deg."%(p) for p in choose_psi]
	fig = pyp.figure( figsize=(12,8))
	pyp.subplots_adjust(hspace=0.3, wspace=0.4)

	for i, p in enumerate(choose_psi):
	    i_psi = ind[i]
	    ax = pyp.subplot(2,3,i+1)
	    ax.set_xscale('log')
	    #ax.set_yscale('log')
	    start_freq, zenith = np.meshgrid(  f_Lo_list, zenith_list)
	    #H=epeak_array[:, :, 0, i_psi]*1e6
	    H = epeak_psi_angle_slice(i_d, i_psi, epeak_array)*1e6
	    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
	    if log:
	    	Hmasked = np.log10(Hmasked)
	    pyp.pcolormesh(start_freq, zenith, Hmasked, cmap='inferno', )
	    pyp.ylim(zenith_list.min(), zenith_list.max())
	    pyp.xlim(f_Lo_list.min(), f_Lo_list.max())
	    if( log):
            	pyp.colorbar(label="log(Peak E-field ($\mu V/m$))")
            else:
            	pyp.colorbar(label="Peak E-field ($\mu V/m$)")
	    pyp.ylabel('Zenith Angle (deg.)')
	    pyp.xlabel('Starting Frequency (MHz)')
	    pyp.grid(True, which='both')
	    pyp.title(titles[i])
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_maps_altitude%2.1fkm_decay%2.1fkm_startfreq_vs_ze%s.png"%(h, decay, plot_suffix) )

###########################################################
def plot_epeak_start_freq_psi(h,decay,epeak_array, f_Lo_list, psi_list, choose_zenith,i_d, plot_suffix="", log=False):
	### Third plot will be epeak vs starting frequency vs psi_list for fixed psi angles
	#ind = [0,1,2,3,4,5,6,7,8]
	#choose_zenith = zenith_list[ind]
	titles = ["Zenith %d deg."%(int(z)) for z in choose_zenith]
	fig = pyp.figure( figsize=(12,8))
	pyp.subplots_adjust(hspace=0.5, wspace=0.5)

	for i_ze, ze in enumerate(choose_zenith):
	    ax = pyp.subplot(3,3,i_ze+1)
	    #ax.set_xscale('log')
	    ax.set_yscale('log')
	    p, start_freq = np.meshgrid(psi_list, f_Lo_list)
	    #H=epeak_array[i_ze, :, 0, :]*1e6
	    H = epeak_zenith_angle_slice(i_d, i_ze, epeak_array)*1e6
	    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
	    if log:
	    	Hmasked = np.log10(Hmasked)
	    pyp.pcolormesh( p,start_freq, Hmasked, cmap='inferno', )
	    pyp.xlim(psi_list.min(), psi_list.max())
	    pyp.ylim(f_Lo_list.min(), f_Lo_list.max())
	    if( log):
            	pyp.colorbar(label="log(Peak E-field ($\mu V/m$))")
            else:
            	pyp.colorbar(label="Peak E-field ($\mu V/m$)")
	    pyp.xlabel('$\psi$ (deg.)')
	    pyp.ylabel('Start Frequency (MHz)')
	    pyp.grid(True, which='both')
	    pyp.title(titles[i_ze])
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_maps_altitude%2.1fkm_decay%2.1fkm_startfreq_psi%s.png"%(h,decay, plot_suffix) )

###########################################################
def write_interpolator_file(h, epeak_array, f_Lo_list, psi_list, zenith_list, decay_list ):
	#############
	# Save interpolators for ever frequency step. The interpolation variables are zenith angle (55-89 degrees) and psi/off-angle (0-3.2 deg)
	interp_file = open(os.environ['TAU_ACC_ZHAIRES_DIR']\
			       +'/interpolator_efields_%2.1fkm.npz'%(h), 'w')

	#interpolator_list = []
	#for i_f_Lo, f_Lo in enumerate(f_Lo_list):
	#    zz,  pp = np.meshgrid(zenith_list, psi_list)
	#    interpolator = scipy.interpolate.LinearNDInterpolator(
	#	np.vstack((zz.flatten(), pp.flatten())).T, 
	#		  epeak_array[:,i_f_Lo,:].T.flatten()
	#		     )
	#    interpolator_list.append(interpolator)
	interpolator_list = []
	for i_f_Lo, f_Lo in enumerate(f_Lo_list):
		dd, zz, pp = np.meshgrid( zenith_list, decay_list, psi_list)
		interpolator = scipy.interpolate.LinearNDInterpolator(
              			zip(dd.flatten(), zz.flatten(), pp.flatten()),
              			epeak_array[:,:,i_f_Lo,:].flatten()) 
		# call the interpolator
		# interpolator(zenith, decay, psi) like in meshgrid
		interpolator_list.append(interpolator)

	np.savez(interp_file, efield_interpolator_list=interpolator_list)
	interp_file.close()
	print "Wrote interpolator file ", interp_file

###########################################################
def clean_full_set(args, gauss_array=None, ncontours_array=None):
	# read CSV pulse files and calculate the frequency-domain peak electric field
	decay_list, zenith_list, psi_list, f_Lo_list, epeak_array = construct_epeak_array(args.altitude)

	# plot slices
	choose_start_freq = [10., 30., 50., 200., 300., 1000.]
	ind = [0, 11, 27, 35, 47, 63]
	choose_psi = psi_list[ind]
	ind = [0,1,2,3,4,5,6,7,8]
	choose_zenith = zenith_list[ind]
	print("plotting the uncleaned versions")
	for i_d, decay in enumerate(decay_list):
		plot_epeak_zenith_psi(args.altitude,decay,epeak_array, zenith_list, psi_list, choose_start_freq, i_d)
		plot_epeak_zenith_start_freq(args.altitude,decay,epeak_array, zenith_list, f_Lo_list, choose_psi, i_d)
		plot_epeak_start_freq_psi(args.altitude,decay,epeak_array, f_Lo_list, psi_list, choose_zenith,i_d)

		plot_epeak_zenith_psi(args.altitude,decay,epeak_array, zenith_list, psi_list, choose_start_freq, i_d, log=True, plot_suffix="_log")
		plot_epeak_zenith_start_freq(args.altitude,decay,epeak_array, zenith_list, f_Lo_list, choose_psi, i_d, log=True, plot_suffix="_log")
		plot_epeak_start_freq_psi(args.altitude,decay,epeak_array, f_Lo_list, psi_list, choose_zenith,i_d, log=True, plot_suffix="_log")

	# clean up the numerical noise from the simulations		
	print "clean_numerical_noise", args.altitude, args.gauss_blur_sigma, args.ncontours
	epeak_array = clean_numerical_noise(args.altitude, epeak_array, decay_list, zenith_list, psi_list, f_Lo_list, 
		gauss_blur_sigma = args.gauss_blur_sigma, ncontours = args.ncontours )

	# plot cleaned slices
	print "plotting cleaned ones"
	for i_d, decay in enumerate(decay_list):
 		plot_epeak_zenith_psi(args.altitude,decay,epeak_array, zenith_list, psi_list, choose_start_freq, i_d, plot_suffix="_clean")
		plot_epeak_zenith_start_freq(args.altitude,decay,epeak_array, zenith_list, f_Lo_list, choose_psi, i_d, plot_suffix="_clean")
		plot_epeak_start_freq_psi(args.altitude,decay,epeak_array, f_Lo_list, psi_list, choose_zenith, i_d, plot_suffix="_clean")

 		plot_epeak_zenith_psi(args.altitude,decay,epeak_array, zenith_list, psi_list, choose_start_freq, i_d, log=True,plot_suffix="_clean_log")
		plot_epeak_zenith_start_freq(args.altitude,decay,epeak_array, zenith_list, f_Lo_list, choose_psi, i_d, log=True,plot_suffix="_clean_log")
		plot_epeak_start_freq_psi(args.altitude,decay,epeak_array, f_Lo_list, psi_list, choose_zenith, i_d, log=True,plot_suffix="_clean_log")

	# save the interpolators
	write_interpolator_file(args.altitude, epeak_array, f_Lo_list, psi_list, zenith_list, decay_list )

###########################################################
def test_zenith(args):
	# read CSV pulse files and calculate the frequency-domain peak electric field
	zenith_list, psi_list, f_Lo_list, epeak_array = construct_epeak_array(args.altitude, args.decay)
	i_ze = np.where(zenith_list == args.zenith)[0][0]
	
	# clean up the numerical noise from the simulations
	epeak_array = clean_numerical_noise_zenith(args.altitude, args.decay, epeak_array, zenith_list, psi_list, f_Lo_list, i_ze,
	gauss_blur_sigma = args.gauss_blur_sigma, ncontours = args.ncontours )

####################################################################################
def load_efield_interpolator(EFIELD_LUT_file_name):
    # the interpolator is for 10-MHz subbands and 
    # is called as interpolator(zenith_angle, starting_frequency,psi_angle)
    # zenith angle is the shower zenith angle in deg.
    # staring_frequency is the lowest frequency in the band in MHz
    # psi_angle is the angel off the shower axis, 
    # equivalent to view angle in deg.
    interp_file = np.load(EFIELD_LUT_file_name)
    return interp_file['efield_interpolator_list'][()]
####################################################################################

def plot_interp_zenith_psi(h, decay, interpolator,choose_start_freq, plot_suffix="", log=False):
    #choose_start_freq = [15., 35., 55., 205., 305., 1005.]
    #ind = [np.where(f_Lo_list == choose_start_freq[i])[0][0] for i in range(len(choose_start_freq))]
    #ind = [1, 3, 5, 8, 20, 30]
    titles = ["%d MHz"%(freq) for freq in choose_start_freq]
    fig = pyp.figure( figsize=(12,8))
    pyp.subplots_adjust(hspace=0.3, wspace=0.4)

    interp_psi_list = np.arange(0.05, 3.2, 0.05)
    interp_zenith_list = np.arange(60, 90, 2.5)
    df = 10.
    for i_choose_freq, freq in enumerate(choose_start_freq):
        i_f_Lo = int(round(freq / df - 1))
        ax = pyp.subplot(2,3,i_choose_freq+1)
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        p, zenith = np.meshgrid( interp_psi_list, interp_zenith_list)
        H = interpolator[i_f_Lo](zenith, decay, p)*1e6
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
	if log:
	    	Hmasked = np.log10(Hmasked)
        pyp.pcolormesh(p, zenith, Hmasked, cmap='inferno', )
        pyp.ylim(interp_zenith_list.min(), interp_zenith_list.max())
        pyp.xlim(interp_psi_list.min(), interp_psi_list.max())
        if( log):
            pyp.colorbar(label="log(Peak E-field ($\mu V/m$))")
        else:
            pyp.colorbar(label="Peak E-field ($\mu V/m$)")
        pyp.ylabel('Zenith Angle (deg.)')
        pyp.xlabel('$\psi$ (deg.)')
        pyp.grid(True, which='both')
        pyp.title(titles[i_choose_freq])
    pyp.suptitle("Interpolated Peak E-fields")
    pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/interp_efield_maps_altitude%2.1fkm_decay%2.1fkm_ze_vs_psi%s.png"%(h, decay, plot_suffix) )

def compare_1d_plots(h, epeak_array,efield_interpolator_list,decay_list, 
                     zenith_list, psi_list, f_Lo_list, 
                     choose_decay,  choose_ze, choose_f_Lo, log=False, plot_suffix="" ):

    pyp.figure(figsize=(24,4))
    pyp.subplots_adjust(hspace=0.3, wspace=0.5)
    pyp.subplot(1,3,1)

    #zenith_list = np.array([55,60,65,70,75,80,85,87,89])
    i_choose_ze = np.where(zenith_list==choose_ze)[0][0]
    i_choose_d  = np.where(decay_list==choose_decay)[0][0]
    epeak_zen = epeak_zenith_angle_slice(i_choose_d,i_choose_ze, epeak_array)*1e6
    norm = matplotlib.colors.Normalize(vmin=-400, vmax=1600)
    for i_f_Lo in [1,3,7,9,26,34,67,100, 150]:
    	if( log):
		pyp.plot(psi_list, np.log10(epeak_zen[ i_f_Lo, :]), color = cmap(norm(f_Lo_list[i_f_Lo])), linewidth=2,
                	label="%d MHz"%(int(f_Lo_list[i_f_Lo])))
	else:	
        	pyp.plot(psi_list,epeak_zen[ i_f_Lo, :], color = cmap(norm(f_Lo_list[i_f_Lo])), linewidth=2,
                	label="%d MHz"%(int(f_Lo_list[i_f_Lo])))
        #pyp.plot(psi_list, interpolator(87, f_Lo_list[i_f_Lo], psi_list)*1e6, color = cmap(norm(i_f_Lo)), linewidth=2)

    for i_f_Lo in np.array([1,3,7,9,26,34,67,100, 150]):
        if( log):
		pyp.plot(psi_list, np.log10(efield_interpolator_list[i_f_Lo](zenith_list[i_choose_ze], decay_list[i_choose_d], psi_list)*1e6), 'o',
                 color=cmap(norm(f_Lo_list[i_f_Lo])), label="%d MHz"%(int(f_Lo_list[i_f_Lo])))
	else:
		pyp.plot(psi_list, efield_interpolator_list[i_f_Lo](zenith_list[i_choose_ze],  decay_list[i_choose_d], psi_list)*1e6, 'o',
                 color=cmap(norm(f_Lo_list[i_f_Lo])), label="%d MHz"%(int(f_Lo_list[i_f_Lo])))

    pyp.legend(loc=[1.01,0])
    pyp.xlabel(" $\psi$ (deg.)")
    if( log ):
    	pyp.ylabel("log(Peak E-field ($\mu V/m$))")
    else:
    	pyp.ylabel("Peak E-field ($\mu V/m$)")
        pyp.ylim(0,)


    pyp.xlim(psi_list.min(), psi_list.max())
    pyp.title("Zenith %d$^{\circ}$, %2.1f km"%(int(zenith_list[i_choose_ze]), decay_list[i_choose_d]))
    #pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/3dinterp_efield_maps_altitude%dkm_psi_ze%d.png"%(altitude, int(zenith_list[i_ze], )) )

    pyp.subplot(1,3,2)
    df = 10.
    i_choose_f_Lo = int(round(choose_f_Lo / df - 1))
    #print choose_f_Lo, i_choose_f_Lo, f_Lo_list[i_choose_f_Lo]
    epeak_zen = epeak_start_freq_slice(i_choose_d, i_choose_f_Lo, epeak_array)*1e6
    norm = matplotlib.colors.Normalize(vmin=55, vmax=100)
    for i_ze, ze in enumerate(zenith_list):
    	if( log ):
		pyp.plot(psi_list, np.log10(epeak_zen[ i_ze, :]), color = cmap(norm(ze)), linewidth=3,
                 label="%d$^{\circ}$"%(int(ze)))
        else:
		pyp.plot(psi_list, epeak_zen[ i_ze, :], color = cmap(norm(ze)), linewidth=3,
                 label="%d$^{\circ}$"%(int(ze)))

    for ze in zenith_list:#np.arange(55, 90, 5):
        if( log):
		pyp.plot(psi_list,  np.log10(efield_interpolator_list[i_choose_f_Lo]( ze, decay, psi_list,)*1e6), marker='o',
                 color=cmap(norm(ze)), linestyle='--', label="%d$^{\circ}$"%(int(ze)))
	else:
		pyp.plot(psi_list,  efield_interpolator_list[i_choose_f_Lo]( ze, decay, psi_list,)*1e6, marker='o',
                 color=cmap(norm(ze)), linestyle='--', label="%d$^{\circ}$"%(int(ze)))
    pyp.xlabel(" $\psi$ (deg.)")
    if( log ):
    	pyp.ylabel("log(Peak E-field ($\mu V/m$))")
    else:
    	pyp.ylabel("Peak E-field ($\mu V/m$)")
        pyp.ylim(0,)

    pyp.xlim(psi_list.min(), psi_list.max())
    pyp.legend(loc=[1.01, 0.1])
    pyp.title("%d MHz, %2.1f km"%((int(f_Lo_list[i_choose_f_Lo]), decay_list[i_choose_d])))
    
    pyp.subplot(1,3,3)
    df = 10.
    epeak_decay = epeak_array[:, i_choose_ze, i_choose_f_Lo, :]*1e6
    ze = zenith_list[i_choose_ze]
    norm = matplotlib.colors.Normalize(vmin=0., vmax=altitude)
    for i_d, d in enumerate(decay_list):
    	if( log ):
            pyp.plot(psi_list, np.log10(epeak_decay[ i_d, :]), color = cmap(norm(d)), linewidth=3,
                 label="%2.1f km"%(d))
        else:
            pyp.plot(psi_list, epeak_decay[ i_d, :], color = cmap(norm(d)), linewidth=3,
                 label="%2.1f km"%(d))

    for d in decay_list:
        if( log):
            pyp.plot(psi_list,  np.log10(efield_interpolator_list[i_choose_f_Lo]( ze, d, psi_list,)*1e6), marker='o',
                 color=cmap(norm(d)), linestyle='--', label="%2.1f km"%(d))
        else:
            pyp.plot(psi_list,  efield_interpolator_list[i_choose_f_Lo]( ze, d, psi_list,)*1e6, marker='o',
                 color=cmap(norm(d)), linestyle='--', label="%2.1f km"%(d))
    pyp.xlabel(" $\psi$ (deg.)")
    if( log ):
    	pyp.ylabel("log(Peak E-field ($\mu V/m$))")
    else:
    	pyp.ylabel("Peak E-field ($\mu V/m$)")
        pyp.ylim(0,)

    pyp.xlim(psi_list.min(), psi_list.max())
    pyp.legend(loc=[1.01, 0.1])
    pyp.title("%d MHz, %d$^{\circ}$"%((int(f_Lo_list[i_choose_f_Lo]), int(zenith_list[i_choose_ze]))))
    pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + \
                "/3dinterp_efield_maps_altitude%2.1fkm_decay%2.1fkm_ze%d_psi_fLo%dMHz.png"%(h, choose_decay, choose_ze, choose_f_Lo)) 


###########################################################
def test_interpolator(args):
	# read CSV pulse files and calculate the frequency-domain peak electric field
	zenith_list, psi_list, f_Lo_list, epeak_array = construct_epeak_array(args.altitude, args.decay)
	#i_ze = np.where(zenith_list == args.zenith)[0][0]

	# read the interpolator file
	EFIELD_LUT_filename = os.environ['TAU_ACC_ZHAIRES_DIR']\
			       +'/interpolator_efields_%dkm.npz'%args.altitude
	load_efield_interpolator(EFIELD_LUT_file_name)
	efield_interpolator_list = load_efield_interpolator(EFIELD_LUT_file_name)

	# compare 2d histograms of the original pulses to the interpolated pulses
	# for epeak vs. zenith vs. psi and for several starting frequencies
	# note that the interpolated parameters are zenith angle and view angle / psi angle
	choose_start_freq = [10., 30., 50., 200., 300., 1000.] # the interpolator was made with 10-MHz subband spacing
	plot_epeak_zenith_psi(altitude, decay, epeak_array, 
                                zenith_list, psi_list, choose_start_freq)
	choose_start_freq = [15., 35., 55., 205., 305., 1005.] # the interpolator was made with 10-MHz subband spacing. These are in between
	plot_interp_zenith_psi(altitude, decay, efield_interpolator_list, choose_start_freq)
	

	# compare 1-d distributions from the interpolator for several different frequencies
	choose_ze=80
	choose_f_Lo=300
	choose_decay=0.5
	compare_1d_plots(args.altitude, epeak_array,efield_interpolator_list, zenith_list, psi_list, f_Lo_list, choose_decay, choose_ze, choose_f_Lo )

if __name__ == "__main__":
	
	parser=argparse.ArgumentParser(description='Build peak electric field lookup tables')
  	parser.add_argument("-a", "--altitude",  default=3, type=float)
	parser.add_argument("-d", "--decay", default=0, type=float)
	parser.add_argument("-g", "--gauss_blur_sigma", default = 3.0, type=float)
	parser.add_argument("-n", "--ncontours", default=10, type=int)
	parser.add_argument("-z", "--zenith", default=-999, type=int)
	parser.add_argument("-t", "--test", default=False, type=bool)
	args=parser.parse_args()

	if args.test:
		test_interpolator(args)
	else:
		# use a spreadsheet to set the parameters
		clean_parms = pd.read_csv("clean_params_electricfields_mountaintop.csv")	      
		#zenith_list = [55, 60, 65, 70, 75, 80, 85, 87, 89]
		zenith_list = [55, 60, 65, 70, 75, 80, 85, 87, 89 ]
		for a in [0.5]:#[1.0, 2.0, 3.0, 4.0]:
			glist = []
			nclist = []
			for d in np.arange(0, 4.0, 0.5):
				if (d<a)and(a>0.5) or (d<=a)and(a==0.5):	
					print "MAKING INTERPOLATORS FOR ALTITUDE : ", a, " DECAY HEIGHT ", d
					args.altitude = a
					#args.decay = d
					zglist = []
					znclist = []
					for ize, ze in enumerate(zenith_list):
						cut = (clean_parms.altitude == a) & (clean_parms.decay == d) & (clean_parms.zenith == ze)
						zglist.append( float(clean_parms[cut].gauss_blur.values[0] ) )
						znclist.append( int(clean_parms[cut].ncontours.values[0] ) )
					glist.append(zglist)
					nclist.append(znclist)
			args.gauss_blur_sigma = glist
			args.ncontours = nclist

			clean_full_set(args)

	pyp.show()
'''	if args.zenith < 0.:
		if( args.test ):
			test_interpolator(args)
		else:
			clean_full_set(args)
	else:
		test_zenith(args)

'''

