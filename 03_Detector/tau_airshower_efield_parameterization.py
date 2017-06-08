import csv
import pandas as pd
import numpy as np
#import filters  # THIS LIBRARY IS NOT INCLUDED IN THE REPOSITORY
import matplotlib.pyplot as pyp
import scipy.signal as signal
import numpy.fft as fft
import os         # ADD THIS SO WE CAN USE ENVIRONMENT VARIABLES
import matplotlib
from scipy.optimize import curve_fit
import scipy.optimize
from scipy.interpolate import interp1d
pastel_colors = ['#FF6666','#FFCC66','#CCFF66','#66FF66','#66FFCC','#66FFFF','#66CCFF','#6666FF','#CC66FF','#FF66FF','#FF6FCF']
bright_colors = ['#FF0000','#FF8000','#FFFF00','#80FF00','#00FF00','#00FF80','#00FFFF','#0080FF','#0000FF','#8000FF','#FF00FF']
bright_colors2 = bright_colors = ['#FF0000','#FF8000','#80FF00','#00FF80','#0080FF','#0000FF','#8000FF','#FF00FF','#FF0080']
pyp.rcParams['font.size']=10
pyp.rcParams['legend.labelspacing'] = 0.1
print os.environ.get('TAU_ACC_ZHAIRES_DIR')


def read_npz_files(antenna_height, zenith):
    npzfile = np.load(os.environ['TAU_ACC_ZHAIRES_DIR'] + "/" + "altitude_%dkm_zenith_%d.npz"%(antenna_height, zenith))
    efield_td = npzfile['efield_td'][()] ## for some reason just reading the arrays gives a 0-dim array, so skip this.
    efield_fd = npzfile['efield_fd'][()]
    return efield_td, efield_fd

# WANT TO PARAMETERIZE THE EFIELDS AT ALL ZENITH ANGLES 
# FOLLOWING TECHNIQUE USED FOR ANITA
# WE WILL HAVE TO CHOOSE A FEW BANDS OF INTEREST BASED ON
# RESULTS FROM EFIELD_PULSES_ANGULAR_ACCEPTANCE
# this is for E=1.e17

def band_peak_efield(efield, freq, f_LO, bandwidth):
    df = freq[1]-freq[0]
    cut = np.logical_and(freq>=f_LO, freq<=f_LO + bandwidth)
    return np.sum(np.abs(efield[cut]))*df

def efield_vs_angle(h, z, f_Lo_list, bandwidth_list):
    # calculate the coherent sum in the frequency domain for simulation at different
    # psi angles (angle relative to the shower axis or decay axis)
    # perform this calculation for each band defined by the starting frequency (f_Lo_list)
    # and bandwidth (bandwidth list)
    off_angle_array = np.arange(0.04,40*0.04, 0.04)
    i_off_angle_array = range(1, len(off_angle_array)+1)
    epeak_list = np.zeros((len(f_Lo_list), len(bandwidth_list), len(i_off_angle_array)))
    cc = 0
    for i_offangle in i_off_angle_array:
        offangle = off_angle_array[i_offangle-1]
        efield_td, efield_fd = read_npz_files(h,z)
        efield = efield_td[i_offangle]['y_v_per_m']
        time = efield_td[i_offangle]['time_s']
        time -= time[0]
        dt = time[1]-time[0]
        E_fft = np.fft.rfft(efield)*dt*1.e6 # use dt in \micros so that fft is in units of V/m/MHz
        fr = np.fft.rfftfreq(len(efield), dt*1.e6) # frequencies in MHz
        for i_f_Lo in range(0,len(f_Lo_list)):
            f_Lo = f_Lo_list[i_f_Lo]
            for i_bandwidth in range(0,len(bandwidth_list)):
                bandwidth = bandwidth_list[i_bandwidth]
                epeak_list[i_f_Lo, i_bandwidth, cc] = band_peak_efield(E_fft, fr, f_Lo, bandwidth)
        cc+=1
    return off_angle_array, epeak_list


def lorentzian_gaussian_background_func(psi, E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2):
    # function that has a floating lorentzian peak added to a background gaussian centered on zero
    v = (psi - gauss_peak) / gauss_width
    E_field = E_0*(frac_gauss*np.exp(-v**2/2.)  + (1-frac_gauss)/(1+v**2)) + + E_1*np.exp(-psi**2/2./width2**2)
    return E_field


###############
# Full width half-max calculation of an Y vs X plot
# use for finding the inital parameters for the fit
def FWHM(X,Y):
    half_max = max(Y) / 2.
    max_ind = np.where(Y==max(Y))[0][0]
    under_half_max_ind = np.where(Y<half_max)[0]

    if( len(under_half_max_ind) == 0):
        return abs((X[-1] - X[0]))
    elif( under_half_max_ind.all() < max_ind):
        return 2*abs((X[max_ind] - X[under_half_max_ind[-1]]))
    
   
    
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    pyp.figure(1)
    pyp.plot(X,Y)
    pyp.figure(2)
    pyp.plot(X[:-1],d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    print left_idx, right_idx
    if( len(left_idx) == 0):
        left_idx = 0
    else:
        left_idx = left_idx[0]
    if( len(right_idx)== 0):
        right_idx = len(d) - 1
    else:
        right_idx = right_idx[0]
    print left_idx, right_idx
    return X[right_idx] - X[left_idx] #return the difference (full width)


###############
def epeak_zenith_fLo_bandwidth(h, zenith_list, f_Lo_list,bandwidth_list):
	# load the epeak pulses from the files and calculate the coherent sum in band
	epeak_array = []
	for z in zenith_list:
    		psi, epeak_list = efield_vs_angle(h, z, f_Lo_list, bandwidth_list)
    		epeak_array.append(epeak_list)
    
	epeak_array = np.array(epeak_array)
	return psi, bandwidth_list, zenith_list, f_Lo_list, epeak_array

################
def save_curve_fit_parameters(curve_fit_parameters, h):
	fit_file = open(os.environ['TAU_ACC_ZHAIRES_DIR']+'/curve_fits_%dkm.npz'%h, 'w')
	np.savez(fit_file, curve_fit_parms=curve_fit_parameters)
	fit_file.close()

if __name__ == '__main__':
	#####
	# Fit curves using a starting frequency of 25 MHz and varying the bandwidth.
	#### Each altitude requires different inital parameters
	#### Should come up with a more robust universal algorithm for setting the initial parameters


	#### This set of initial parameters works for the h=37km files.	
	# load the efield pulses from the files
	h=37
	psi, bandwidth_list, zenith_list, f_Lo_list, epeak_array = epeak_zenith_fLo_bandwidth(h=h,
								zenith_list=np.array([60, 70, 75, 80, 85, 87]), f_Lo_list=np.arange(10.,1000., 10.), bandwidth_list=np.array([10])   )
	i_bandwidth = 0
	bandwidth = bandwidth_list[i_bandwidth]

	curve_fit_parameters = np.zeros((len(zenith_list), len(f_Lo_list), 6))
	for i_ze in range(0, len(zenith_list)):
	    for i_f_Lo in range(0,len(f_Lo_list)):
		f_Lo = f_Lo_list[i_f_Lo]
		E_0  =  max(epeak_array[i_ze][i_f_Lo][i_bandwidth])# V/m at ground level for 10^17 eV tau lepton.
		E_1 = epeak_array[i_ze][i_f_Lo][i_bandwidth][0]   # V/m at ground level for 10^17 eV tau lepton.
		frac_gauss = E_1/E_0 # fraction of Gaussian peak ( 1-frac_gauss is lorentzian)
		gauss_peak = psi[np.where(max(epeak_array[i_ze][i_f_Lo][i_bandwidth]) == epeak_array[i_ze][i_f_Lo][i_bandwidth][:])[0][0]]      # degree. 
		if( i_f_Lo == 0):
		    gauss_width = FWHM(psi, epeak_array[i_ze][i_f_Lo][i_bandwidth][:])
		    width2 = 0.05
		    p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] # initial parameters for the curve fit.
		else:
		    # Use the last gaussian and lorentzian widths to seed the fit.
		    # The fit is sensitive to the individual frac_gauss, E_0, and E_1, but the widths are changing slowly.
		    gauss_width = FWHM(psi, epeak_array[i_ze][i_f_Lo][i_bandwidth][:])
		    p0 = popt
		    if(  frac_gauss > 0.5 ):
			# If the psi = 0 Efield strength is roughly equivalent to the peak E-field strength, 
			# supress the second Gaussian
			# Consider doing a single gaussian fit.
			p0[0] = E_0
			p0[2] = gauss_peak
			p0[4] = E_1*0.01
		    else: 
			p0[0] = E_0
			p0[2] = gauss_peak
			p0[4] = E_1
	   
	        popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
	        curve_fit_parameters[i_ze, i_f_Lo, :] = popt
		if( i_f_Lo % 10 == 0):
		    pyp.figure(1)
	   	    pyp.subplot(2, 3, i_ze+1)
		    pyp.plot(psi, epeak_array[i_ze,i_f_Lo, i_bandwidth,:]*1e6, 'o', color = bright_colors2[i_f_Lo%9], label="%d-%d MHz"%(int(f_Lo), int(f_Lo + bandwidth)))
		    pyp.plot(psi, lorentzian_gaussian_background_func(psi, *popt)*1e6,  'k')
	    pyp.figure(1)
	    pyp.subplot(2, 3, i_ze+1)
	    pyp.title("Zenith %d$^{\circ}$"%(int(zenith_list[i_ze])), fontsize=12)
	    if( i_ze == 0 or i_ze == 3):
	    	pyp.ylabel("Peak E-Field ($\mu$V/m))")
	    if( i_ze > 2):
		pyp.xlabel("$\psi$ (deg.)")
	    pyp.xlim(0)
	    pyp.ylim(0, 1.2*epeak_array[i_ze,:, i_bandwidth,:].max()*1e6)
	    if( i_ze == 5):
	    	pyp.legend(loc=[1.01,0.01])
	pyp.suptitle("Fit Parameterization $10^{17}$ eV %d km"%h)
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_eparms_altitude%dkm.png"%(h)) 
	save_curve_fit_parameters(curve_fit_parameters, h)

	#####
	#### This set works for the 1 km altitude
	h=1
	psi, bandwidth_list, zenith_list, f_Lo_list, epeak_array = epeak_zenith_fLo_bandwidth(h, 
								zenith_list=np.array([60, 70, 75, 80, 85, 87]), f_Lo_list=np.arange(10.,1000., 10.), bandwidth_list=np.array([10])   )
	i_bandwidth = 0
	bandwidth = bandwidth_list[i_bandwidth]

	curve_fit_parameters = np.zeros((len(zenith_list), len(f_Lo_list), 6))
	for i_ze in range(0, len(zenith_list)):
	    pyp.figure(i_ze)
	    for i_f_Lo in range(0,len(f_Lo_list)):
		f_Lo = f_Lo_list[i_f_Lo]
		E_0  =  max(epeak_array[i_ze][i_f_Lo][i_bandwidth])# V/m at ground level for 10^17 eV tau lepton.
		E_1 = epeak_array[i_ze][i_f_Lo][i_bandwidth][-1]   # V/m at ground level for 10^17 eV tau lepton.
		frac_gauss = E_1/E_0 # fraction of Gaussian peak ( 1-frac_gauss is lorentzian)
		gauss_peak = psi[np.where(max(epeak_array[i_ze][i_f_Lo][i_bandwidth]) == epeak_array[i_ze][i_f_Lo][i_bandwidth][:])[0][0]]      # degree. 
		if( i_f_Lo == 0):
		    gauss_width = FWHM(psi, epeak_array[i_ze][i_f_Lo][i_bandwidth][:])
		    width2 = 0.05
		    p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] # initial parameters for the curve fit.
		else:
		    # Use the last gaussian and lorentzian widths to seed the fit.
		    # The fit is sensitive to the individual frac_gauss, E_0, and E_1, but the widths are changing slowly.
		    p0 = popt
		if( frac_gauss < 0.05):
		    p0[4] = E_1
	   
		try:
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		except:
		    print "Fit failed, trying again"
		    p0 = curve_fit_parameters[i_ze-1, i_f_Lo, :]
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		
		curve_fit_parameters[i_ze, i_f_Lo, :] = popt
		if( i_f_Lo % 10 == 0):
		    pyp.figure( 100)
		    pyp.subplot(2, 3, i_ze+1)
		    pyp.plot(psi, epeak_array[i_ze,i_f_Lo, i_bandwidth,:]*1e6, 'o', color = bright_colors2[i_f_Lo%9], label="%d-%d MHz"%(int(f_Lo), int(f_Lo + bandwidth)))
		    #popt = abs(popt)
		    pyp.plot(psi, lorentzian_gaussian_background_func(psi, *popt)*1e6,  'k')
	    pyp.figure( 100)
 	    pyp.subplot(2, 3, i_ze+1)
	    pyp.title("Zenith %d$^{\circ}$"%(int(zenith_list[i_ze])), fontsize=12)
	    if( i_ze == 0 or i_ze == 3):
	    	pyp.ylabel("Peak E-Field ($\mu$V/m))")
	    if( i_ze > 2):
		pyp.xlabel("$\psi$ (deg.)")
	    pyp.xlim(0)
	    pyp.ylim(0, 1.2*epeak_array[i_ze,:, i_bandwidth,:].max()*1e6)
	    if( i_ze == 5):
	    	pyp.legend(loc=[1.01,0.01])
	pyp.suptitle("Fit Parameterization $10^{17}$ eV %d km"%h)
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_eparms_altitude%dkm.png"%(h)) 
	save_curve_fit_parameters(curve_fit_parameters, h)
	#####
	#### This set works for the 2 km altitude
	h=2
	psi, bandwidth_list, zenith_list, f_Lo_list, epeak_array = epeak_zenith_fLo_bandwidth(h, 
								zenith_list=np.array([60, 70, 75, 80, 85, 87]), f_Lo_list=np.arange(10.,1000., 10.), bandwidth_list=np.array([10])   )
	i_bandwidth = 0
	bandwidth = bandwidth_list[i_bandwidth]

	curve_fit_parameters = np.zeros((len(zenith_list), len(f_Lo_list), 6))
	for i_ze in range(0, len(zenith_list)):

	    for i_f_Lo in range(0,len(f_Lo_list)):
		f_Lo = f_Lo_list[i_f_Lo]
		E_0  =  max(epeak_array[i_ze][i_f_Lo][i_bandwidth])# V/m at ground level for 10^17 eV tau lepton.
		E_1 = epeak_array[i_ze][i_f_Lo][i_bandwidth][-1]   # V/m at ground level for 10^17 eV tau lepton.
		frac_gauss = E_1/E_0 # fraction of Gaussian peak ( 1-frac_gauss is lorentzian)
		gauss_peak = psi[np.where(max(epeak_array[i_ze][i_f_Lo][i_bandwidth]) == epeak_array[i_ze][i_f_Lo][i_bandwidth][:])[0][0]]      # degree. 
		if( i_f_Lo == 0):
		    gauss_width = FWHM(psi, epeak_array[i_ze][i_f_Lo][i_bandwidth][:])
		    width2 = 0.05
		    p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] # initial parameters for the curve fit.
		else:
		    # Use the last gaussian and lorentzian widths to seed the fit.
		    # The fit is sensitive to the individual frac_gauss, E_0, and E_1, but the widths are changing slowly.
		    p0 = popt
		    p0[2] = gauss_peak
		try:
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		except:
		    print "Fit failed, trying again", i_ze, i_f_Lo
		    if( i_ze > 0):
			p0 = curve_fit_parameters[i_ze-1, i_f_Lo, :]
		    else:
			p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] 
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		curve_fit_parameters[i_ze, i_f_Lo, :] = popt
		if( i_f_Lo % 10 == 0):
		    pyp.figure(200)
		    pyp.subplot(2, 3, i_ze+1) 
		    pyp.plot(psi, epeak_array[i_ze,i_f_Lo, i_bandwidth,:]*1e6, 'o', color = bright_colors2[i_f_Lo%9], label="%d-%d MHz"%(int(f_Lo), int(f_Lo + bandwidth)))
		    pyp.plot(psi, lorentzian_gaussian_background_func(psi, *popt)*1e6,  'k')
	    pyp.figure( 200)
	    pyp.subplot(2, 3, i_ze+1)
	    pyp.title("Zenith %d$^{\circ}$"%(int(zenith_list[i_ze])), fontsize=12)
	    if( i_ze == 0 or i_ze == 3):
	    	pyp.ylabel("Peak E-Field ($\mu$V/m))")
	    if( i_ze > 2):
		pyp.xlabel("$\psi$ (deg.)")
	    pyp.xlim(0)
	    pyp.ylim(0, 1.2*epeak_array[i_ze,:, i_bandwidth,:].max()*1e6)
	    if( i_ze == 5):
	    	pyp.legend(loc=[1.01,0.01])
	pyp.suptitle("Fit Parameterization $10^{17}$ eV %d km"%h)
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_eparms_altitude%dkm.png"%(h)) 
	save_curve_fit_parameters(curve_fit_parameters, h)

	#####
	#### This set works for the 3 km altitude
	h=3
	psi, bandwidth_list, zenith_list, f_Lo_list, epeak_array = epeak_zenith_fLo_bandwidth(h, 
								zenith_list=np.array([60, 70, 75, 80, 85, 87]), f_Lo_list=np.arange(10.,1000., 10.), bandwidth_list=np.array([10])   )
	i_bandwidth = 0
	bandwidth = bandwidth_list[i_bandwidth]

	pyp.figure(300)
	pyp.subplots(2,3)
	curve_fit_parameters = np.zeros((len(zenith_list), len(f_Lo_list), 6))
	for i_ze in range(0, len(zenith_list)):

	    for i_f_Lo in range(0,len(f_Lo_list)):
		f_Lo = f_Lo_list[i_f_Lo]
		E_0  =  max(epeak_array[i_ze][i_f_Lo][i_bandwidth])# V/m at ground level for 10^17 eV tau lepton.
		E_1 = epeak_array[i_ze][i_f_Lo][i_bandwidth][-1]   # V/m at ground level for 10^17 eV tau lepton.
		frac_gauss = E_1/E_0 # fraction of Gaussian peak ( 1-frac_gauss is lorentzian)
		gauss_peak = psi[np.where(max(epeak_array[i_ze][i_f_Lo][i_bandwidth]) == epeak_array[i_ze][i_f_Lo][i_bandwidth][:])[0][0]]      # degree. 
		if( i_f_Lo == 0):
		    gauss_width = FWHM(psi, epeak_array[i_ze][i_f_Lo][i_bandwidth][:])
		    width2 = 0.05
		    p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] # initial parameters for the curve fit.
		else:
		    # Use the last gaussian and lorentzian widths to seed the fit.
		    # The fit is sensitive to the individual frac_gauss, E_0, and E_1, but the widths are changing slowly.
		    p0 = popt
		    p0[2] = gauss_peak
		try:
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		except:
		    print "Fit failed, trying again", i_ze, i_f_Lo
		    if( i_ze > 0):
			# if there are older curve fits for the same frequency band, try using those initial parameters
			p0 = curve_fit_parameters[i_ze-1, i_f_Lo, :]
		    else:
			# if not, go with the best guess
			p0 = [E_0, frac_gauss, gauss_peak, gauss_width, E_1, width2] 
		    popt, pcov = curve_fit(lorentzian_gaussian_background_func, psi, epeak_array[i_ze, i_f_Lo,i_bandwidth,:], p0=p0)
		curve_fit_parameters[i_ze, i_f_Lo, :] = popt
		if( i_f_Lo % 10 == 0):
		    pyp.figure(300)
		    pyp.subplot(2, 3, i_ze+1)
		    pyp.plot(psi, epeak_array[i_ze,i_f_Lo, i_bandwidth,:]*1e6, 'o', color = bright_colors2[i_f_Lo%9], label="%d-%d MHz"%(int(f_Lo), int(f_Lo + bandwidth)))
		    pyp.plot(psi, lorentzian_gaussian_background_func(psi, *popt)*1e6,  'k')
	    pyp.subplot(2, 3, i_ze+1)
	    pyp.title("Zenith %d$^{\circ}$"%(int(zenith_list[i_ze])), fontsize=12)
	    if( i_ze == 0 or i_ze == 3):
	    	pyp.ylabel("Peak E-Field ($\mu$V/m))")
	    if( i_ze > 2):
		pyp.xlabel("$\psi$ (deg.)")
	    pyp.xlim(0)
	    pyp.ylim(0, 1.2*epeak_array[i_ze,:, i_bandwidth,:].max()*1e6)
	    if( i_ze == 5):
	    	pyp.legend(loc=[1.01,0.01])
	pyp.suptitle("Fit Parameterization $10^{17}$ eV %d km"%h)	
	pyp.savefig(os.environ['TAU_ACC_DETECTOR_PLOTS_DIR'] + "/efield_eparms_altitude%dkm.png"%(h) )
	save_curve_fit_parameters(curve_fit_parameters, h)
	
	pyp.show()
