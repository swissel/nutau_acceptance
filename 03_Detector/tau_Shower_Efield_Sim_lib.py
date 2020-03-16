import numpy as np
import os
import sys 
import re
from Tau_Decay_Pythia_Simulator import Tau_Decay_Simulator
# tau_exit_num_events = 1.e6

# Define constants
Earth_radius   = 6371.      # km
tau_life_time  = 2.906e-13  # seconds
tau_mass       = 1776.82e6  # eV/c2
speed_of_light = 299792.458 # km / second
kB_W_Hz_K = 1.38064852e-23 # Watts / Hz / K

# Threshold voltage based on galactic noise and system temperature:
T_ice  = 290     # Kelvin (water temperature)
T_sys  = 140.     # System temperature at the lower part of the ANITA band from elog 119
kB     = 1.38e-23 # in Volts^2/Ohm/Hz/Kelvin
Z_0 = 377. # Ohms impedance of free space
Z_L = 50. # Ohms 50-Ohm characteristic impedance
Z_A = 50. # Ohms antenna impedance
frac_sky = 0.5 # fraction of sky visible to the antenna
#Gain_dB = 10. # antenna gain
#BW     = 400.e6   # bandwidth for first two ANITA bands.
#threshold_voltage = 8.9e-06 # V, this is 8.9 microVolts for sky temp 290 Kelvin T_sys
#threshold_voltage =np.sqrt(kB*(T_ant+T_sys)*Z_load*BW) # V, this is 8.22 microVolts for the parameters given above

#threshold_voltage = 72.e-6 # V, this is based on taking half the peak field for the weakest ANITA-1 EAS candidate. (0.466 mV/m convolved with antenna effective height at 300 MHz). 
#threshold_voltage_snr = 5.0
#print >> sys.stderr,  'threshold_voltage_snr', threshold_voltage_snr

####################################################################################

def galactic_temperature(f_MHz):
    # Dulk 2001
    nu = f_MHz # Hz
    tau = 5.0 * pow(nu, -2.1)
    Ieg = 1.0610e-20
    Ig  = 2.48e-20
    # Iv in  W/m^2/Hz/sr
    Iv = Ig * pow(nu, -0.52) * (1-np.exp(-tau))/tau + Ieg * pow(nu, -0.80) * np.exp(-tau)
    
    kB = 1.38064852e-23 # Watts / Hz / K
    c = speed_of_light * 1e3 # m/s
    temp = Iv * c**2 / (2*(nu*1e6)**2)/kB
    ## IV is the intensity
    ## temp is the galactic noise temperature in Kelvin
    return Iv, temp # W/m^2/Hz/sr, K

####################################################################################

#def noise_v_sq(freq_MHz, Z_L, Z_A, Gain_dB, Nphased=1):
#    # TODO: Include reflection losses at the antennas as (1-Gamma^2)
#    val = kB_W_Hz_K*(frac_sky*galactic_temperature(freq_MHz)[1] + (1.0-frac_sky)*T_ice)*Z_L*Nphased
#    # This is in Volts^2/Hz
#    return val

####################################################################################
def reflection_coefficient(Z_A, Z_L):
    # Z_A is the antenna impedance
    # Z_L is the load impedance characteristic to the system
    return (Z_A - np.conj(Z_L))/(Z_A + Z_L)


def noise_voltage(freq_min_MHz, freq_max_MHz, df,  Z_L, Z_A, Gain_dB, Nphased=1.):
    # Note that Gain_dB isn't used here because the noise is integrated over the full beam of the antenna
    # TODO: Update with more detailed antenna model

    # calculate the reflection coefficient
    Gamma = reflection_coefficient(Z_A, Z_L)
    eff_load = 1. - np.abs(Gamma)**2

    # Radiation resistance of the antennas
    R_A = np.real(Z_A)

    # this is in V^2/Hz
    T_gal = galactic_temperature(np.arange(freq_min_MHz, freq_max_MHz + df, df))[1] 
    gal_noise = np.sqrt(Nphased * np.trapz(T_gal) * frac_sky * df * 1e6 * kB_W_Hz_K * R_A * eff_load  )
    sys_noise = np.sqrt(Nphased * (T_ice*(1.-frac_sky) * eff_load + T_sys) * (freq_max_MHz - freq_min_MHz) * kB_W_Hz_K * Z_L  ) 
    
    # assuming we're phasing after the amplifier rather than before.
    # if before, then the system temperature would also decrease as 1/N 
    # note that temperature ~ power, so this is decreasing the noise voltage by 1/sqrt(N)
    # as you would expect from incoherent noise
    combined_temp = Nphased*(T_gal*frac_sky * eff_load + T_ice*(1.-frac_sky) * eff_load + T_sys)

    #print >> sys.stderr,  np.sum(T_gal)
    #print >> sys.stderr,  np.sum(combined_temp)
    #print >> sys.stderr,  np.sqrt(np.sum(T_gal) * df * 1e6 * kB_W_Hz_K * Z_L )
    
    # this is in V (rms)
    combined_noise = np.sqrt( np.trapz(combined_temp) * df * 1e6 * kB_W_Hz_K * Z_L  )

    return combined_noise, gal_noise, sys_noise

####################################################################################

def tau_lepton_decay_range(tau_lepton_energy):
    return speed_of_light * 10**(tau_lepton_energy)  / tau_mass * tau_life_time

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

'''def E_field_interp(efield_interpolator_list, view_angle_deg, zenith_angle_deg, f_Lo, f_High, log10_tau_energy, distance_exit_km, distance_decay_km):
    # Lorentzian beam pattern based on 10-MHz filtered subbands of Harm's results
    # Returns electric field peak in V/m
  
    # Since the efields are stored in 10-MHz subbands
    # integrate over the range from f_Lo to f_High in 10-MHz bands
    E_field = np.zeros(len(distance_decay_km))
    #E_field = 0.
    # TODO: Right now forcing the parameters outside the interpolation range to the edges
    # shoudl replace with extrapolation
    z = zenith_angle_deg.copy()
    v = view_angle_deg.copy()
    z[z>89.] = 89.
    z[z<55.] = 55.
    v[v<0.04] = 0.04
    v[v>3.16] = 3.16

    df = 10.
    for freq in np.arange(f_Lo, f_High, df):
        i_f_Lo = int(round(freq / df - 1))
        E_field += efield_interpolator_list[i_f_Lo](z,v) 
	
    # account for ZHAIReS sims only extending to 3.16 deg 
    #TODO: right now sampling a gaussian centered at zero to extrapolate to the wider angles. Should verify with ZHAireS sims out to wider angles
    E_field[view_angle_deg>3.16] = 0.

    E_field *= distance_exit_km/distance_decay_km   # distance to tau decay point correction
    E_field *= 10**(log10_tau_energy - 17.) # Energy scaling
    return E_field'''

def E_field_interp(efield_interpolator_list, view_angle_deg, zenith_angle_deg, 
                   altitude, decay_altitude_km, f_Lo, f_High, 
                   log10_tau_energy, distance_exit_km, distance_decay_km):
    # Lorentzian beam pattern based on 10-MHz filtered subbands of Harm's results
    # Returns electric field peak in V/m
  
    # Since the efields are stored in 10-MHz subbands
    # integrate over the range from f_Lo to f_High in 10-MHz bands
    E_field = np.zeros(len(distance_decay_km))

    # Forcing the parameters outside the interpolation range to the edges
    z = zenith_angle_deg.copy()
    v = view_angle_deg.copy()
    d = decay_altitude_km.copy()
    z[z>89.] = 89.
    z[z<55.] = 55.
    v[v<0.04] = 0.04
    v[v>3.16] = 3.16
    d[d>altitude-0.5] = altitude-0.5

    df = 10.
    for freq in np.arange(f_Lo, f_High, df):
        i_f_Lo = int(round(freq / df - 1))
        E_field += efield_interpolator_list[i_f_Lo](z,d,v) 
	
    # account for ZHAIReS sims only extending to 3.16 deg 
    #TODO: right now sampling a gaussian centered at zero to extrapolate to the wider angles. Should verify with ZHAireS sims out to wider angles
    E_field[view_angle_deg>3.16] = 0.

    E_field *= distance_exit_km/distance_decay_km   # distance to tau decay point correction
    E_field *= 10**(log10_tau_energy - 17.) # Energy scaling
    return E_field
####################################################################################

def Voltage_interp(efield_interpolator_list, view_angle_deg, zenith_angle_deg, 
                   altitude, decay_altitude_km, f_Lo, f_High, 
                   log10_tau_energy, distance_decay_km,
		   Gain_dB,Z_A, Z_L, Nphased=1):
    # Lorentzian beam pattern based on 10-MHz filtered subbands of Harm's results
    # Returns electric field peak in V/m
  
    # Since the efields are stored in 10-MHz subbands
    # integrate over the range from f_Lo to f_High in 10-MHz bands
    Voltage = np.zeros(len(distance_decay_km))
    #E_field = 0.
    # TODO: Right now forcing the parameters outside the interpolation range to the edges
    # shoudl replace with extrapolation
    z = zenith_angle_deg.copy()
    v = view_angle_deg.copy()
    d = decay_altitude_km.copy()

    z[z>89.] = 89.
    z[z<55.] = 55.
    v[v<0.04] = 0.04
    v[v>3.16] = 3.16
    d[d>altitude-0.5] = altitude-0.5
    d[d<0] = 0.

    # these are the simulation arrays used to generate the lookup tables
    zenith_list = np.array([50, 55, 60, 65, 70, 75, 80, 85, 87, 89])
    decay_altitude_list = np.arange(0., altitude,0.5)

    # find the nearest neighbor for both the zenith angle at the exit point and the decay alttidue
    zhaires_sim_icethick = 0.0
    zhaires_sim_detector_altitude = altitude
    e_zhaires_tau_shower = 1e17
    r_zhaires_tau_shower = np.zeros(len(v))
    for i in range(len(v)): # loop over all the sims
    	i_ze = find_nearest(zenith_list, zenith_angle_deg[i])[0]
    	i_d  = find_nearest(decay_altitude_list, decay_altitude_km[i], lower_bound = 0)[0]
    	nearest_zenith_angle = zenith_list[i_ze]
    	nearest_decay_altitude = decay_altitude_list[i_d]
    	r_zhaires_tau_shower[i] = get_distance_decay_to_detector_zenith_exit(zhaires_sim_icethick , nearest_decay_altitude,
									   zhaires_sim_detector_altitude, nearest_zenith_angle)

    df = 10.
    for freq in np.arange(f_Lo, f_High, df):
        i_f_Lo = int(round(freq / df - 1))
        # using the average frequency in the bin to calculate the voltage
        Voltage += E_to_V_signal(efield_interpolator_list[i_f_Lo](z,d,v), Gain_dB, (freq+df)/2., Z_A, Z_L, Nphased)

    # account for ZHAIReS sims only extending to 3.16 deg 
    Voltage[view_angle_deg>3.16] = Voltage[view_angle_deg>3.16]*np.exp( -(view_angle_deg[view_angle_deg>3.16]-0.)**2 / (2*3.16)**2)
    
    Voltage *= (r_zhaires_tau_shower / distance_decay_km )  # distance to tau decay point correction
    Voltage *= (pow(10, log10_tau_energy)/ e_zhaires_tau_shower) # Energy scaling
    return Voltage
####################################################################################
# Parameterization of efield simulations
#
#########
def lorentzian_gaussian_background_func(psi, E_0, frac_gauss, gauss_peak,
                                        gauss_width, E_1, width2):
    v = (psi - gauss_peak) / gauss_width
    # background gaussian centered at psi=0
    # lorentizian + Gaussian to model the beam (takes care of some asymmetries inside vs. outside cone)
    E_field = E_0*(frac_gauss*np.exp(-v**2/2.)  + (1-frac_gauss)/(1+v**2)) + E_1*np.exp(-psi**2/2./width2**2)
    
    return E_field
    
def efield_anita_generic_parameterization(energy, distance_shower_to_detector, theta_view, parm_decay_altitude=0): 
    if( parm_decay_altitude!= 5):
    	if parm_decay_altitude != 0:
    	     print >> sys.stderr,  "Warning: Using parameterization for 0 km tau decay altitude, but requested ", parm_decay_altitude
    	# Parameterization for a 10^17 eV tau shower at zenith angle of 65deg / emergence angle 25deg
	# decay altitude of 0 km, ground elevation of 3 km, detector altitude of 37 km
    	parms = [  1.41662197e-04,   8.10616766e-01,   9.35347123e-01,
             1.45537805e-01,1.13492535e-06,   1.00000000e+00]
    	# Evaluates to 86.417407396098497 km
    	r_zhaires_tau_shower = get_distance_decay_to_detector_zenith_exit(2, 0, 37, 65,)
    	e_zhaires_tau_shower = 1e17
    else:
   	# Parameterization for a 10^17 eV tau shower at zenith angle of 65 deg / emergence angle 25deg
   	# decay altitude of 5 km, ground elevation of 3 km, detector altitude of 37 km
    	parms = [  2.85992061e-04,   7.76446976e-01,   4.78891355e-01,   1.46997140e-01,
   	3.19229984e-06,  2.90189325e-01]
    	# Evaluates to 67.682759132837418 km
    	r_zhaires_tau_shower = get_distance_decay_to_detector_zenith_exit(2, 5, 37, 65,)
    	e_zhaires_tau_shower = 1e17

    epeak = lorentzian_gaussian_background_func(theta_view, *parms)
    escaled = epeak * (energy / e_zhaires_tau_shower) * (r_zhaires_tau_shower / distance_shower_to_detector )
    return escaled
   
############
# Now test that this parameterization will work for several decay altitudes
def efield_anita_generic_parameterization_decay(energy, decay_altitude, distance_shower_to_detector, theta_view):
    
    # Parameterization for a 10^17 eV tau shower at zenith angle of 60deg / emergence angle 30deg
    # ground elevation of 3 km, detector altitude of 37 km and different decay altitudes 

    # The first row in the array sis for decay altitude of 0 km; the last is 9 km.
    parm_array = np.array([[  1.74467904e-04,   7.63866727e-01,   8.77653517e-01,  1.55970631e-01, 1.66957500e-06,   1.00000000e+00],
                           [  1.99696451e-04,   6.21403627e-01,   7.68656948e-01,  1.62155531e-01, 3.59185823e-07,   3.44511913e-07],
                           [  2.31546719e-04,   8.19862178e-01,   6.72984202e-01,  1.63150020e-01, 2.60513430e-06,   9.64407268e-02],
                           [  2.63833633e-04,   6.00079563e-01,   5.70486369e-01,  1.63871081e-01, 1.39574219e-06,   1.00000000e+00],
                           [  2.80864266e-04,   6.78254444e-01,   4.72776446e-01,  1.71672019e-01, 5.77921580e-06,   4.34553854e-01],
                           [  3.17318501e-04,   8.11459190e-01,   4.02198021e-01,  1.70212480e-01, 9.96559946e-06,   2.21691882e-01],
                           [  2.86854486e-04,   8.16986311e-01,   3.08890072e-01,  1.77961387e-01, 9.99999999e-06,   1.06078872e-01],
                           [  2.89818471e-04,   8.15486187e-01,   2.45378525e-01,  1.80607523e-01, 9.99999536e-06,   4.35679248e-02],
                           [  2.58729422e-04,   9.12020646e-01,   1.82649564e-01,  1.88871931e-01, 9.99999785e-06,   2.59300117e-02],
                           [  1.93973823e-04,   9.99999455e-01,   1.03194543e-01,  1.86353972e-01, 4.35152861e-06,   9.99999999e-01]])
    
    decay_altitude_list = np.array([0,1,2,3,4,5,6,7,8,9])
    escaled = np.zeros(len(energy))
    for i in range(len(theta_view)):
    	i_d  = find_nearest(decay_altitude_list, decay_altitude[i], lower_bound = 0)[0]
    	if( i_d > -1): # leave the electric field zero for the ones that are not above ground level
    		nearest_decay_altitude = decay_altitude_list[i_d]
    		parms = parm_array[nearest_decay_altitude]
    		epeak = lorentzian_gaussian_background_func(theta_view[i], *parms)
     
    		# r_parm_array = np.array([ 86.4174073961 ,  65.4677078729 ,  63.4691177949 ,  61.4714655076 ,  
    		#                59.4747492519 ,  57.4789672742 ,  55.4841178256 ,  53.4901991629 ,
    		#               51.4972095476 ,  49.5051472466  ])
    		#r_zhaires_tau_shower = r_parm_array[nearest_decay_altitude]
    		r_zhaires_tau_shower = get_distance_decay_to_detector_zenith_exit(2, nearest_decay_altitude, 37, 60)
    		e_zhaires_tau_shower = 1e17
    		escaled[i] = epeak * (energy[i] / e_zhaires_tau_shower) * (r_zhaires_tau_shower / distance_shower_to_detector[i] )
    return escaled
############
# useful to find the nearest zenith angle and/or decay altitude for the lookup table of parameterization
##########
def find_nearest(array, values, lower_bound=None, upper_bound=None):
    # finds the nearest values in the arrays
    # if the values are outside the desire range, sets the index to -1
    values = np.atleast_1d(values)
    indices = np.abs(np.round(np.subtract.outer(array, values))).argmin(0)
    out = array[indices]
    if( lower_bound != None):
        bound_ind = np.where(values < lower_bound)
        indices[bound_ind] = -1
    if( upper_bound != None):
        bound_ind = np.where(values > upper_bound)
        indices[bound_ind] = -1
    return indices

############
# Now test that this parameterization will work for several decay altitudes and zenith angles
def efield_anita_generic_parameterization_decay_zenith(energy, decay_altitude, zenith_exit_deg, distance_shower_to_detector, theta_view, parm_2d):
    
    # Parameterization for a 10^17 eV tau shower at zenith angle of 60deg / emergence angle 30deg
    # ground elevation of 3 km, detector altitude of 37 km
    #
    # The first row in the array sis for decay altitude of 0 km; the last is 9 km.
    
    # these are the simulation arrays used to generate the parameterizations
    zenith_list = np.array([50, 55, 60, 65, 70, 75, 80, 85, 87, 89])
    decay_altitude_list = np.array([0,1,2,3,4,5,6,7,8,9])

    escaled = np.zeros(len(energy))
    for i in range(len(theta_view)):
         if decay_altitude[i] > 0: # if the decay altitude < 0, leave the electric field at 0

		# find the nearest neighbor for both the zenith angle at the exit point and the decay alttidue
                i_ze = find_nearest(zenith_list, zenith_exit_deg[i])[0]
                i_d  = find_nearest(decay_altitude_list, decay_altitude[i], lower_bound = 0)[0]
	 
		# if the decay altitude is < 0, then throw this event out
                if i_d >= 0:	
                	nearest_zenith_angle = zenith_list[i_ze]
                	nearest_decay_altitude = decay_altitude_list[i_d]
                	parms = parm_2d[i_ze, i_d]

                	epeak = lorentzian_gaussian_background_func(theta_view[i], *parms)
                	# Distance from the shower to the detector for the parameterized LDFs 
                	# at different decay altitudes and decay zenith angles
                	zhaires_sim_icethick = 2.0
                	zhaires_sim_detector_altitude = 37.
                	e_zhaires_tau_shower = 1e17
                	r_zhaires_tau_shower = get_distance_decay_to_detector_zenith_exit(zhaires_sim_icethick , nearest_decay_altitude,
									   zhaires_sim_detector_altitude, nearest_zenith_angle)
                	escaled[i] = epeak * (energy[i] / e_zhaires_tau_shower) * (r_zhaires_tau_shower / distance_shower_to_detector[i] )
    return escaled

####################################################################################

def E_to_V_signal(E_pk, Gain_dB, freq_MHz, Z_A, Z_L, Nphased=1):
    # Derived assuming:
    #		P_r   : the power received at the antenna P_r = A P_inc = |V_A|^2/(8 R_A) = V_A  Z_L / (Z_L + Z_A) 
    #			Note: V_A = 2 V_L for a perfectly matched antenna -- which we assume here
    #		P_inc : incident power at the antenna P_inc = |E|^2 / (2 Z_0) 
    #		A : antenna aperture = \lambda^2/(4pi) G eff_load eff_pol 
    #		eff_load: load mismatch factor eff_load = (1 - Gamma^2), where Gamma is the reflection coefficient = (Z_L - Z_A*)/(Z_L + A_Z);	assuming that we have a perfect match
    #		eff_pol : polarization mismatch factor; assuming that this is built into ZHAireS pulses
    # calculate the reflection coefficient
    Gamma = reflection_coefficient(Z_A, Z_L)
    eff_load = 1. - np.abs(Gamma)**2

    # Radiation resistance of the antennas
    R_A = np.real(Z_A)

    #print >> sys.stderr,  Gamma, R_A, Z_L / (Z_A + Z_L), eff_load, Z_A, Z_L 
    V_A = 2. * E_pk * (speed_of_light*1.e3)/(freq_MHz * 1.e6) * np.sqrt(R_A/Z_0 * pow(10., Gain_dB/10.)/4./np.pi * eff_load) * Nphased
    V_L = V_A * Z_L / (Z_A + Z_L) # V_L = 1/2 * V_A for a perfectly matched antenna
    return V_L

####################################################################################

def get_random_interaction_point(cos_th_min, num_ev):
  cos_th_E = np.random.uniform(cos_th_min,1., num_ev)
  sin_th_E = np.sqrt(1.-cos_th_E**2)
  x = Earth_radius*sin_th_E
  y = 0.*sin_th_E  # keep it in the x-z plane
  z = Earth_radius*cos_th_E
  return x,y,z, cos_th_E, sin_th_E
  
####################################################################################

def get_random_direction_vector(cos_th_E, sin_th_E, num_events):
  # for an area, the direction is randomly sampled in cos^2. 
  # this is due to the dot product in the integral
  cos_th_particle = np.sqrt(np.random.uniform(0.,1., num_events) )
  ph_particle     = np.random.uniform(0.,2.*np.pi, num_events)
  k_x = np.sqrt(1.-cos_th_particle**2)*np.cos(ph_particle)
  k_y = np.sqrt(1.-cos_th_particle**2)*np.sin(ph_particle)
  k_z = cos_th_particle
  #k_x_rot = k_x
  #k_y_rot = k_y
  #k_z_rot = k_z
  k_x_rot = +k_x*cos_th_E + k_z*sin_th_E
  k_y_rot =  k_y
  k_z_rot = -k_x*sin_th_E + k_z*cos_th_E
  return k_x_rot, k_y_rot, k_z_rot
  
####################################################################################

def get_tau_lepton_energy(inouts, tau_neutrino_energy, cos_theta_exit):
    # dummy function until you can sort out the look-up table.
    # this may need to be converted into a class that gets initialized once and runs lookups quickly.
    P_exit = np.sqrt(1-cos_theta_exit**2)**40
    log10_tau_lepton_energy = np.log10(tau_neutrino_energy*np.sqrt(1-cos_theta_exit**2))
    log10_tau_lepton_energy = np.min([np.random.normal(log10_tau_lepton_energy, log10_tau_lepton_energy/10.), log10_tau_lepton_energy])
    return 10**log10_tau_lepton_energy, P_exit

####################################################################################

def get_decay_range(tau_lepton_energy):
    c = 2.998e5 # km/s
    tau_rest_frame_life_time = 2.9e-13 # seconds
    tau_rest_frame_mass = 1.7768e9 # eV/c^2
    gamma = tau_lepton_energy / tau_rest_frame_mass
    L = gamma * tau_rest_frame_life_time * c # in km
    #print >> sys.stderr,  '**L',L
    return np.random.exponential(L)

####################################################################################

def get_position(lam, x_exit, y_exit, z_exit, k_x, k_y, k_z):
    x = x_exit + lam * k_x
    y = y_exit + lam * k_y
    z = z_exit + lam * k_z
    return x,y,z

####################################################################################

def get_view_angle(k_x, k_y, k_z, x_pos, y_pos, z_pos, x_det, y_det, z_det):
    #x_pd = x_pos - x_det
    #y_pd = y_pos - y_det
    #z_pd = z_pos - z_det
    #cos_view_angle = -( k_x * x_pd + k_y * y_pd + k_z * z_pd ) / np.sqrt(x_pd**2 + y_pd**2 + z_pd**2)
    #return np.arccos(cos_view_angle) # radians
    x_pd = x_det - x_pos
    y_pd = y_det - y_pos
    z_pd = z_det - z_pos
    cos_view_angle = ( k_x * x_pd + k_y * y_pd + k_z * z_pd ) / np.sqrt(x_pd**2 + y_pd**2 + z_pd**2)
    return np.arccos(cos_view_angle) # radians

####################################################################################

def get_distance_to_detector(x_pos, y_pos, z_pos, x_det, y_det, z_det):
    x_pd = x_pos - x_det
    y_pd = y_pos - y_det
    z_pd = z_pos - z_det

    return np.sqrt( x_pd**2 + y_pd**2 + z_pd**2 )

####################################################################################
def update_exit_point(k_x, k_y, k_z, x_exit, y_exit, z_exit, icethick_geom, icethick_now):
    # Assuming that the geometric sims are run with one ice thickness
    # and these sims are run with a different ice thickness" 
    r_exit_geom = np.sqrt(x_exit**2 + y_exit**2 + z_exit**2)
    r_exit_new  = r_exit_geom + (icethick_now - icethick_geom)
    
    # define a new vector along the r^hat direction
    # the dot product of the unit vecto rhat and k-hat gives the angle 
    # of the cord through the ice
    cosTh = 1./r_exit_geom *( x_exit * k_x + y_exit * k_y + z_exit * k_z )
    
    a = 1.0
    b = 2 * r_exit_geom * cosTh 
    c = (r_exit_geom**2 - r_exit_new**2)
    d = (-b + np.sqrt(b**2 - 4.*a*c) )/(2.*a)
    
    x_exit = x_exit +  d * k_x
    y_exit = y_exit +  d * k_y
    z_exit = z_exit +  d * k_z
    return x_exit, y_exit, z_exit 
    

def decay_point_geom(k_x, k_y, k_z, x_exit, y_exit, z_exit, X0_decay, x_det, y_det, z_det):
    #def decay_point_geom(k_x, k_y, k_z, x_exit, y_exit, z_exit, X0_decay, x_det, y_det, z_det, emergence_angle):
    # all (input and output) angles in radians
    
    # This is a little complicated, but it works.
    #emergence_phi_angle = np.arctan(k_y/k_x);
    #z_decay = X0_decay * np.sin(emergence_angle) + z_exit
    #y_decay = -(X0_decay * np.cos(emergence_angle))*np.sin(emergence_phi_angle) + y_exit
    #x_decay = -(X0_decay * np.cos(emergence_angle))*np.cos(emergence_phi_angle) + x_exit
    
    # Using euler angles is much cleaner.
    x_decay = X0_decay * k_x + x_exit
    y_decay = X0_decay * k_y + y_exit
    z_decay = X0_decay * k_z + z_exit
    
    decay_view_angle = get_view_angle(k_x, k_y, k_z, x_decay, y_decay, z_decay, x_det, y_det, z_det)
    decay_dist_to_detector = get_distance_to_detector(x_decay, y_decay, z_decay, x_det, y_det, z_det )
    
    return x_decay, y_decay, z_decay, decay_view_angle, decay_dist_to_detector

def get_altitude(x, y, z, ground_altitude=Earth_radius):
   # finds the altitude of the decay point, given a ground altitude relative to the center of the
   # Earth and the decay position in Cartesian coordinates
   r = np.sqrt(x**2 + y**2 + z**2)
   return r - ground_altitude

def get_zenith_angle(k_x, k_y, k_z, x, y, z):
	# finds the zenith angle relative to the Earth normal for a vector
	# defined by a Cartesian point (x,y,z) and it's propagation direction (k-hat)
	r = np.sqrt(x*x + y*y + z*z)
	# the zenith angle is the angle between n-hat (the normal to the surface)
	# and the propagation direction of the shower k-hat
	cos_zenith = k_x * x / r + k_y * y / r + k_z * z / r
	return np.arccos(cos_zenith) # radians

####################################################################################

def distance_detector_to_ground(zenith_angle_rad, ground_elevation, detector_altitude):
	''' Calculates the distance from the exit point to the detector assuming that the Earth has an average radius up to sea level
	    plus the ground_elevation (i.e. the ice or water thickness) .The distance, d, is the given by the spherical earth distance
	    between the two points, which can be calculated from the law of cosines with three legs of a triangle:
	    	1. d sin theta
		2. R_Earth + z_g  + d cos theta 
		3. R_Earth + h_det
		Note that legs 1 & 2 are orthogonal

	Inputs: zenith_angle_rad: Zenith angle in radians of exit (theta above)
		ground_elevation: thickness of the ground above sea level (z_g above, in km because the Earth radius is in km)
		detector_altitude: location of the detector in km above sea level.
	'''
	a = 1.
	b = 2.* np.cos(zenith_angle_rad) * ( Earth_radius + ground_elevation)
	c = (Earth_radius + ground_elevation)**2 - (Earth_radius + detector_altitude)**2

	d = (-b + np.sqrt(b**2 - 4.*a*c) )/(2.*a)
	return d
	
####################################################################################
def get_X0( ground_elevation, decay_altitude, zenith_exit_deg, R_e = Earth_radius):
    a = 1
    b = 2*np.cos(zenith_exit_deg*np.pi/180.)*(R_e + ground_elevation)
    c = (R_e + ground_elevation)**2 - (R_e + ground_elevation + decay_altitude)**2
    X0 =( -b + np.sqrt(b**2 - 4*a*c) )/(2*a)
    return X0

def get_decay_zenith_angle(ground_elevation, decay_altitude, X0, R_e=Earth_radius):
    A = X0
    B = R_e + ground_elevation + decay_altitude
    C = R_e + ground_elevation
    
    cosZenithDecay = (A**2 + B**2 - C**2) / (2*A*B)
    return 180./np.pi * np.arccos(cosZenithDecay)

def get_distance_decay_to_detector(ground_elevation, 
                                   decay_altitude, detector_altitude,
                                   zenith_decay_deg, R_e=Earth_radius):
    a = 1
    b = 2*np.cos(zenith_decay_deg*np.pi/180.)*(R_e + ground_elevation + decay_altitude)
    c = (R_e + ground_elevation + decay_altitude)**2 - (R_e + detector_altitude)**2
    d =( -b + np.sqrt(b**2 - 4*a*c) )/(2*a)
    return d

def get_distance_decay_to_detector_zenith_exit(ground_elevation, 
                                   decay_altitude, detector_altitude,
                                   zenith_exit_deg, R_e=Earth_radius):
    # finds the distance between the two points defined by the
    # decay altitude, detector altitude, and the zenith angle at the point
    if( decay_altitude ==0 ):
        return get_X0( ground_elevation, detector_altitude, zenith_exit_deg)
    
    X0 = get_X0(ground_elevation, decay_altitude, zenith_exit_deg)
    zenith_decay_deg = get_decay_zenith_angle(ground_elevation, decay_altitude, X0)
    return get_distance_decay_to_detector(ground_elevation, decay_altitude, 
                                          detector_altitude, zenith_decay_deg)

def get_geometric_zenith_angle(x_det, y_det, z_det, x_exit, y_exit, z_exit):
    # geometric zenith angle at the exit point 
    # comes from the dot product between the r-vector to the exit point and
    # the vector pointing from the exit point ot the detector
    r_det_exit = np.sqrt(pow(x_exit-x_det, 2) + pow(y_exit-y_det,2) + pow(z_exit-z_det, 2))
    r_exit = np.sqrt(x_exit*x_exit + y_exit*y_exit + z_exit*z_exit)
    r_det = np.sqrt(x_det*x_det + y_det*y_det + z_det*z_det)

    n_x = x_exit/r_exit
    n_y = y_exit/r_exit
    n_z = z_exit/r_exit
    # this is the geometry k-vector!
    kg_x = (x_det - x_exit) / r_det_exit
    kg_y = (y_det - y_exit) / r_det_exit
    kg_z = (z_det - z_exit) / r_det_exit

    cosZe_geom = n_x * kg_x + n_y * kg_y + n_z * kg_z
    zenith = np.arccos(cosZe_geom)
    return zenith

####################################################################################

def GEOM(altitude, num_events, view_angle_cut_deg = 5.):
    ##############################################
    # STEPS
    ##############################################

    # 0. Calculate the Earth angle of the horizon.
    cos_th_min = 1./(1 + altitude/Earth_radius)

    # 1. Pick random location in the field of view
    x_exit, y_exit, z_exit, cos_th_E, sin_th_E = get_random_interaction_point(cos_th_min, num_events)

    # 2. Pick a random direction with respect to the normal of the exit point
    k_x, k_y, k_z = get_random_direction_vector(cos_th_E, sin_th_E, num_events)

    # 3. Calculate tau exit angle
    cos_theta_exit = (k_x*x_exit + k_y*y_exit + k_z*z_exit)/Earth_radius

    # 4. Get view angle from the ground
    view_angle = get_view_angle(k_x, k_y, k_z, x_exit, y_exit, z_exit, 0., 0., Earth_radius + altitude)

    # 5. Reject events below that view angle
    cut = view_angle < view_angle_cut_deg * np.pi/180.

    return x_exit[cut], z_exit[cut], k_x[cut], k_y[cut], k_z[cut], cos_theta_exit[cut], view_angle[cut]


####################################################################################

def load_tau_LUTs(LUT_file_name):
    f = np.load(LUT_file_name)
    data_array    = f['data_array']
    th_exit_array = f['th_exit_array']
    num_sim       = f['num_sim']
    f.close()
    P_exit = []
    for k in range(0, len(th_exit_array)):
        P_exit.append(float(len(data_array[k]))/float(num_sim[k]))
    return np.array(90. + th_exit_array), np.array(P_exit), np.array(data_array)

####################################################################################

def parse_input_args(input_arg_string):
    vals = str(input_arg_string)[10:-1]
    vals = vals.split(',')
    keys = []
    values = []
    for val in vals:
        keys.append(val.split('=')[0].split()[-1])
        values.append(val.split('=')[1])
        #print >> sys.stderr,  keys[-1], values[-1]
    input_dict = { keys[i] : values[i] for i in range(0,len(keys)) }
    return input_dict

####################################################################################
def A_OMEGA_tau_exit(geom_file_name, LUT_file_name, EFIELD_LUT_file_name, cut_ang, f_Lo, f_High, outTag='test', 
			N=-1, noise='default', Gain_dB=10.0, Nphased=1, LUT=True, icethick_geom = 0.0, threshold_voltage_snr=5.0,
			start_event=0): 
    
    print >> sys.stderr,  "Inputs to A_OMEGA_tau_exit:\n=============================="
    print >> sys.stderr,  "geom_file_name", geom_file_name
    print >> sys.stderr,  "LUT_file_name", LUT_file_name
    print >> sys.stderr,  "EFIELD_LUT_file_name", EFIELD_LUT_file_name
    print >> sys.stderr,  "cut_ang", cut_ang
    print >> sys.stderr,  "Bandwidth f_Lo, f_High", f_Lo, f_High
    print >> sys.stderr,  "Output Tag: ", outTag
    print >> sys.stderr,  "N", N
    print >> sys.stderr,  "Gain_dB", Gain_dB
    print >> sys.stderr,  "Nphased", Nphased
    print >> sys.stderr,  "LUT", LUT
    print >> sys.stderr,  "icethick_geom ", icethick_geom

    # 1. Load geometry file
    GEOM = np.load(geom_file_name)
    GEOM_inputs = parse_input_args(GEOM['input_args'])
    print >> sys.stderr,  '\nGEOMETRY FILE INPUT ARGUMENTS\n', 
    for key in GEOM_inputs.keys():
        print >> sys.stderr,  '\t',key.ljust(20), GEOM_inputs[key]

    if( N == -1):
    	N_tot = len(GEOM['cos_theta_exit'])
    else:
    	N_tot = int(N)
    print >> sys.stderr,  'Number of Geometry Events to scan', N_tot
    print >> sys.stderr,  ''

    cos_theta_exit  = GEOM['cos_theta_exit'][start_event:start_event + N_tot]
    exit_view_angle = GEOM['exit_view_angle'][start_event:start_event + N_tot]
    x_exit          = GEOM['x_exit'][start_event:start_event + N_tot]
    y_exit          = np.zeros(len(x_exit))
    #y_exit          = GEOM['y_exit']
    z_exit          = GEOM['z_exit'][start_event:start_event + N_tot]
    k_x             = GEOM['k_x'][start_event:start_event + N_tot]
    k_y             = GEOM['k_y'][start_event:start_event + N_tot]
    k_z             = GEOM['k_z'][start_event:start_event + N_tot]
    A_geom          = GEOM['A_geom']
    altitude = float(GEOM_inputs['altitude'])
    GEOM.close()
   
    # 2. Set up the detector position
    x_det = 0.
    y_det = 0.
    z_det = altitude + Earth_radius
    
    # 3. Impose a viewing angle cut
    view_cut = exit_view_angle*180./np.pi<cut_ang
    N_cut = len(cos_theta_exit[view_cut])
    print >> sys.stderr, "Number of events: ", N_cut, "Start event:", start_event, "Maximum number of events:", N_tot, "Max - start:", N_tot-start_event
    A_Omega = A_geom* float(N_cut) / float(N_tot)
    # TODO: fix to calculate dist_to_detector from decay point rather than from
    # the exit point which is what it is now.
    dist_exit_to_detector = get_distance_to_detector(x_exit, y_exit, z_exit, x_det, y_det, z_det)[view_cut] 
    GEOM_theta_exit = np.arccos(cos_theta_exit[view_cut])*180./np.pi
    exit_view_angle = exit_view_angle[view_cut] #radians
    zenith_angle_exit = np.arccos(cos_theta_exit) # radians
    theta_emergence = np.pi/2. - zenith_angle_exit # radians
    
    # 7. set up some arrays to be calculated on the fly
    dist_decay_to_detector = np.zeros(len(dist_exit_to_detector))
    x_decay = np.zeros(len(x_exit[view_cut]))
    y_decay = np.zeros(len(y_exit[view_cut]))
    z_decay = np.zeros(len(z_exit[view_cut]))
    decay_view_angle = np.zeros(len(exit_view_angle))
    X0_dist = np.zeros(len(dist_exit_to_detector))
    log10_tau_energy = np.zeros(len(exit_view_angle))
    log10_shower_energy = np.zeros(len(exit_view_angle))
    Peak_Voltage_SNR = np.zeros(len(exit_view_angle))
    zenith_angle_geom = np.zeros(len(exit_view_angle))

    # 8. Set up histogram for the differential acceptance
    binwidth=0.3 # 0.3deg wide bins
    emerge_angle_bins = np.arange(0., 90.+ binwidth, binwidth)
    nsim_emerge_angle = np.zeros(len(emerge_angle_bins))
    
    # 9. Compute the noise in this band
    #       integrating in 10-MHz steps, to match the frequency bins of the peak voltage
    df = 10.# MHz
    noises = noise_voltage(f_Lo, f_High, df,  Z_L, Z_A, Gain_dB, Nphased) # factor of 1e6 because noise temperature is in V/Hz.
    print >> sys.stderr,  "Galactic Noise voltage ", noises[1]*1e6, " micro-Volts"
    print >> sys.stderr,  "Galactic noise temperature K ", galactic_temperature(np.arange(f_Lo, f_High+df, df))[1], sum(galactic_temperature(np.arange(f_Lo, f_High+df, df))[1])
    print >> sys.stderr,  "System noise ", noises[2]
    if( noise == 'sys'):
    	Noise_Voltage = noises[2]
    	print >> sys.stderr,  "Noise voltage ", Noise_Voltage*1e6, " micro-Volts, due to T_sys + T_ant"
    elif( noise == 'gal'):
    	Noise_Voltage = noises[1]
    	print >> sys.stderr,  "Noise voltage ", Noise_Voltage*1e6, " micro-Volts, due to G*T_gal"
    else: # default is the combination
    	Noise_Voltage = noises[0]
    	print >> sys.stderr,  "Noise voltage ", Noise_Voltage*1e6, " micro-Volts, due to G*T_gal + T_sys + T_ant"
    
    # 5. Load Energy Look-up Table
    print >> sys.stderr,  "Loading energy look-up table: ", LUT_file_name
    LUT_th_exit, LUT_P_exit, LUT_log10_E_tau = load_tau_LUTs(LUT_file_name)
    P_LUT = np.zeros(len(x_exit[view_cut]))   # zero until it is
    P_range = np.zeros(len(x_exit[view_cut])) # zero until it is proven to decay before it passes the detector
    P_det = np.zeros(len(x_exit[view_cut]))   # zero until it is proven to be detectable
    ice_thick_match = re.match( r'.*/(\d+\.\d+)km\_ice.*', LUT_file_name)
    ice_thick = 2.0
    if ice_thick_match:
    	ice_thick = float(ice_thick_match.group(1))
    print >> sys.stderr, "Ice thickness ", ice_thick, " km"

    # 6. Load the Efield interpolator for this altitude 
    #    N. B.: EField_LUT_file_name should have the altitude in the file name
    if LUT:
    	global efield_interpolator_list
    	efield_interpolator_list = load_efield_interpolator(EFIELD_LUT_file_name)
    else:
    	global parm_2d 
    	parm_2d = load_efield_parameterization()
    
    # 7. Load the Tau Decay Simulator
    ###
    TDS = Tau_Decay_Simulator()

    # 11. Loop Through Geometry Events
    sum_P_exit = 0.
    sum_P_exit_P_range = 0.	     # zero until it is proven to decay before it passes the detector
    sum_P_exit_P_range_P_det = 0.    # zero until it is proven to be detectable
    triggered_events = []
    #ranged_events = []
    energy_frac =TDS.sample_energy_fraction(num_events=len(GEOM_theta_exit))
    for k in range(0,len(GEOM_theta_exit)):

	# 7.1 Get LUT exit angle closest to geometry exit angle. Note these are both zenith angles. The GEOM_theta_exit is the zenith angle of the particle
        idx = np.argmin(np.abs(LUT_th_exit - GEOM_theta_exit[k]))
        # 7.2 Get tau lepton energy and decay position
        P_LUT[k] = LUT_P_exit[idx]
        if( P_LUT[k] > 1.e-15):  # make sure the probability of this event is non-zero  
	    # sample the pexit lookup tables for the exiting tau energy
            log10_tau_energy[k] = LUT_log10_E_tau[idx][np.random.randint(0,len(LUT_log10_E_tau[idx]))] # random tau energy
            # estimated decay range based on tau energy
            decay_range = tau_lepton_decay_range(log10_tau_energy[k])                     
	    # sample the energy that goes into the shower
            log10_shower_energy[k] = np.log10( energy_frac[k] )  + log10_tau_energy[k]
	    
            x_exit[k], y_exit[k], z_exit[k] = update_exit_point(k_x[k], k_y[k], k_z[k], x_exit[k], y_exit[k], z_exit[k], icethick_geom, ice_thick)
	
	    # the exit point and the detector position define the exit point zenith angle and emergence angle
            zenith_angle_geom[k] = get_geometric_zenith_angle(x_det, y_det, z_det, x_exit[k], y_exit[k], z_exit[k]) # radians, zenith angle in the geometry tables at the exit point

	    # sample exponentially distributed decay positions	
            X0_dist[k] = np.random.exponential(scale=decay_range)                          
            x_decay[k], y_decay[k], z_decay[k], decay_view_angle[k], dist_decay_to_detector[k] = decay_point_geom(k_x[k], k_y[k], k_z[k], x_exit[k], y_exit[k], z_exit[k], X0_dist[k], x_det, y_det, z_det)
            # If the event is contained within the range, then the probability is 1.
            # If the shower decays beyond the detector, then it has a negative x-position. 
	    # The range probability in that case is zero.
            if((X0_dist[k] < dist_exit_to_detector[k]) and (x_decay[k] > 0.)):
                P_range[k] = 1.
	
        if(k%100000 == 0 and k>0):
            print >> sys.stderr, 'Progress: %d events '%k, log10_tau_energy[k]
            
    # 12. Calculate the electric field and voltage at the detector 
    # 
    # Electric field is calculated by summing the interpolated electric fields in 10-MHz subbands
    # perform array-wise interpolations (outside of the event loop) to minimize number of computations. 
    # Voltage is calculated in three steps:
    #	1. Interpolate electric field from ZHAireS simulations in 10 MHz bins
    #   2. Convert electric field to voltage based on detector model (gain, nphased, 
    #   3. Sum over all 10-MHz bins in the desired frequency band (f_Lo to f_High)
    
    # 13.1 calculate the decay altitude and zenith angles
    decay_altitude = get_altitude(x_decay, y_decay, z_decay, ground_altitude=Earth_radius+ice_thick)
    zenith_angle_decay =  get_zenith_angle(k_x, k_y, k_z, x_decay, y_decay, z_decay) # zenith angle of the shower at the decay point

    if( LUT ):
        Peak_Voltage = Voltage_interp( efield_interpolator_list, decay_view_angle*180./np.pi, zenith_angle_decay*180./np.pi,
				       altitude, decay_altitude, 
				       f_Lo, f_High, log10_shower_energy, dist_decay_to_detector, Gain_dB, Z_A, Z_L, Nphased)
    else:
	# 0-km decay parameterization
	#Peak_Efield = efield_anita_generic_parameterization(pow(10, log10_tau_energy), dist_decay_to_detector, decay_view_angle*180./np.pi, parm_decay_altitude=0)
	# 5-km decay parameterization
	#Peak_Efield   = efield_anita_generic_parameterization(pow(10, log10_tau_energy),dist_decay_to_detector, decay_view_angle*180./np.pi, parm_decay_altitude=5)
	# multiple decay altitudes
	#Peak_Efield   = efield_anita_generic_parameterization_decay(pow(10, log10_tau_energy), zhs_decay_altitude, dist_decay_to_detector, decay_view_angle*180./np.pi)
	# multiple decay altitudes and zenith angles
	
	# the ZHAireS simulations ere all run at ice thicknesses of  2 km, 
    	# so you assume that the altitude is calculated from the Earth radius + 2 km
    	zhs_decay_altitude = get_altitude(x_decay, y_decay, z_decay, ground_altitude=Earth_radius+2.0)
    	Peak_Efield   = efield_anita_generic_parameterization_decay_zenith(pow(10, log10_shower_energy), zhs_decay_altitude, zenith_angle_decay*180./np.pi,  
								dist_decay_to_detector,  decay_view_angle*180./np.pi, parm_2d)
 	# 13.2 generalizing to 300 MHz
    	Peak_Voltage = E_to_V_signal(Peak_Efield, Gain_dB, 300., Z_A, Z_L, Nphased) 
    	Peak_Voltage_Threshold = E_to_V_signal(Epk_to_pk_threshold, Gain_dB, 300., Z_A, Z_L, Nphased) / Vpk_to_Vpkpk_conversion
    #print >> sys.stderr,  "Thresholds : Peak Voltage (not pk-to-pk) :", Peak_Voltage_Threshold, "V/m, Epeak-to-peak :", Epk_to_pk_threshold, " V/m, pk2pk to pk conversion ", Vpk_to_Vpkpk_conversion

    # 14. Check for trigger at the detector
    for k in range(0,len(GEOM_theta_exit)):
        if( P_LUT[k] > 1.e-15):
            if( P_range[k] == 1.):
                Peak_Voltage_SNR[k] = Peak_Voltage[k] / Noise_Voltage
		#ranged_events.append(np.array( [ log10_tau_energy[k], dist_exit_to_detector[k], X0_dist[k], dist_decay_to_detector[k], Peak_Voltage[k], exit_view_angle[k]*180./np.pi, decay_view_angle[k]*180./np.pi,  zenith_angle[k]*180./np.pi ]))

                if(Peak_Voltage_SNR[k] > threshold_voltage_snr):
                    P_det[k] = 1.
                    triggered_events.append(np.array( [ log10_tau_energy[k], dist_exit_to_detector[k], X0_dist[k], dist_decay_to_detector[k], Peak_Voltage[k], exit_view_angle[k]*180./np.pi, decay_view_angle[k]*180./np.pi,  zenith_angle_exit[k]*180./np.pi, zenith_angle_decay[k]*180./np.pi, zenith_angle_geom[k]*180./np.pi, decay_altitude[k], P_LUT[k], P_range[k], P_det[k], log10_shower_energy[k] ]))
        
        sum_P_exit                += P_LUT[k]
        sum_P_exit_P_range        += P_LUT[k] * P_range[k]
        sum_P_exit_P_range_P_det  += P_LUT[k] * P_range[k] * P_det[k]
	
	# Each sample should report
	#1. the emergence angle corresponding to the exit point.
	#2. the emergence angle of the particle.
	#3. P_exit of the emergence angle of the sampled particle.
	#4. P_range of the sampled particle.
	#5. P_trig of the sampled particle.

        if(k%100000 == 0 and k>0):
            #if(k%1== 0 and k>0):
            print >> sys.stderr,  'After %d events: %d events triggered: '%(len(GEOM_theta_exit[0:k]), len(triggered_events))
            print >> sys.stderr,  '\t %1.3e km^2 sr'%A_Omega
            print >> sys.stderr,  '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit / float(k))
            print >> sys.stderr,  '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range / float(k) )
            print >> sys.stderr,  '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range_P_det / float(k) )
            print >> sys.stderr,  '\t ',np.array(triggered_events).shape
            print >> sys.stderr,  ''
        #all_events.append(np.array( [ log10_tau_energy[k], dist_exit_to_detector[k], X0_dist[k], dist_decay_to_detector[k], Peak_Voltage[k], exit_view_angle[k]*180./np.pi, decay_view_angle[k]*180./np.pi,  zenith_angle_exit[k]*180./np.pi, zenith_angle_decay[k]*180./np.pi, zenith_angle_geom[k]*180./np.pi, decay_altitude[k], P_LUT[k], P_range[k], P_det[k]]))

    print >> sys.stderr,  'After all %d events, %d events triggered: '%(len(GEOM_theta_exit), len(triggered_events))
    print >> sys.stderr,  '\t %1.3e km^2 sr'%(A_Omega)
    print >> sys.stderr,  '\t %1.3e km^2 sr, %d events exited Earth'%(A_Omega*sum_P_exit/float(N_cut),N_tot*sum_P_exit/float(N_cut) )
    print >> sys.stderr,  '\t %1.3e km^2 sr, %d events decayed before detector'%(A_Omega*sum_P_exit_P_range/float(N_cut), N_tot*sum_P_exit_P_range/float(N_cut))
    print >> sys.stderr,  '\t %1.3e km^2 sr, %d events triggered'%(A_Omega*sum_P_exit_P_range_P_det/float(N_cut), N_tot*sum_P_exit_P_range_P_det/float(N_cut))
    np.savez("energyfraction.npz", shower_energy_fraction=energy_frac)
    np.savez(outTag+'.npz', A_Omega_start    = A_Omega,
                            A_Omega_exit     = A_Omega*sum_P_exit/float(N_cut),
                            A_Omega_range    = A_Omega*sum_P_exit_P_range/float(N_cut),
                            A_Omega_trig     = A_Omega*sum_P_exit_P_range_P_det/float(N_cut),
                            N_events_start   = N_tot,
                            N_triggered      = len(triggered_events))


    np.savez(outTag+'_events.npz', A_Omega_start    = A_Omega,
                            A_Omega_exit     = A_Omega*sum_P_exit/float(N_cut),
                            A_Omega_range    = A_Omega*sum_P_exit_P_range/float(N_cut),
                            A_Omega_trig     = A_Omega*sum_P_exit_P_range_P_det/float(N_cut),
			    noise_voltage    = Noise_Voltage,
                            triggered_events = np.array(triggered_events),
			    emerge_angle_bins = emerge_angle_bins,
			    nsim_emerge_angle = nsim_emerge_angle,
                            #ranged_events = np.array(ranged_events),
			    #all_events = np.array(all_events))
			    )
    print >> sys.stderr,  "Wrote ", outTag+'.npz and ', outTag+'_events.npz'

    exit()
