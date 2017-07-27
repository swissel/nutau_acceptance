import numpy as np
import os
tau_exit_num_events = 1.e6

# Define constants
Earth_radius   = 6371.      # km
tau_life_time  = 2.906e-13  # seconds
tau_mass       = 1776.82e6  # eV/c2
speed_of_light = 299792.458 # km / second
kB_W_Hz_K = 1.38064852e-23 # Watts / Hz / K

# Threshold voltage based on galactic noise and system temperature:
T_ant  = 105.     # Kelvin (based on elog 655 entry by P. Gorham)
T_sys  = 140.     # System temperature at the lower part of the ANITA band from elog 119
kB     = 1.38e-23 # in Volts^2/Ohm/Hz/Kelvin
Z_0 = 377. # Ohms impedance of free space
Z_L = 50. # Ohms 50-Ohm characteristic impedance
R_L = 50. # Ohms real part of antenna impedance
Gain_dB = 10. # antenna gain
#BW     = 400.e6   # bandwidth for first two ANITA bands.
#threshold_voltage = 8.9e-06 # V, this is 8.9 microVolts for sky temp 290 Kelvin T_sys
#threshold_voltage =np.sqrt(kB*(T_ant+T_sys)*Z_load*BW) # V, this is 8.22 microVolts for the parameters given above

#threshold_voltage = 72.e-6 # V, this is based on taking half the peak field for the weakest ANITA-1 EAS candidate. (0.466 mV/m convolved with antenna effective height at 300 MHz). 
threshold_voltage_snr = 5.0
print 'threshold_voltage_snr', threshold_voltage_snr
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

def noise_v_sq(freq_MHz, Z_L, Z_A, Gain_dB):
    # TODO: Include reflection losses at the antennas as (1-Gamma^2)
    val = kB_W_Hz_K*galactic_temperature(freq_MHz)[1]*Z_L*pow(10., Gain_dB/10.)
    # This is in Volts^2/Hz
    return val

####################################################################################

def noise_voltage(freq_min_MHz, freq_max_MHz, df,  Z_L, Z_A, Gain_dB):
    # TODO: Include reflection losses at the antennas as (1-Gamma^2)
    gal_noise = np.sqrt(np.sum(noise_v_sq(np.arange(freq_min_MHz, freq_max_MHz + df, df), Z_L, R_L, Gain_dB)*1e6 * df ))
    sys_noise = np.sqrt( kB_W_Hz_K * (T_sys + T_ant) * (freq_max_MHz - freq_min_MHz)*1e6 * Z_L)
    
    # this is in V^2/Hz
    T_gal = galactic_temperature(np.arange(freq_min_MHz, freq_max_MHz + df, df))[1]
    combined_temp = pow(10., Gain_dB/10.)*T_gal + T_sys + T_ant

    #print np.sum(T_gal)
    #print np.sum(combined_temp)
    #print np.sqrt(np.sum(T_gal) * df * 1e6 * kB_W_Hz_K * Z_L )
    
    # this is in V
    combined_noise = np.sqrt( np.sum(combined_temp) * df * 1e6 * kB_W_Hz_K * Z_L)

    return combined_noise, gal_noise, sys_noise

####################################################################################

def tau_lepton_decay_range(tau_lepton_energy):
    return speed_of_light * 10**(tau_lepton_energy)  / tau_mass * tau_life_time

####################################################################################
def load_efield_interpolator(h, EFIELD_LUT_file_name=os.environ['TAU_ACC_ZHAIRES_DIR'] + '/interpolator_efields_%dkm.npz'):
    # the interpolator is for 10-MHz subbands and 
    # is called as interpolator(zenith_angle, starting_frequency,psi_angle)
    # zenith angle is the shower zenith angle in deg.
    # staring_frequency is the lowest frequency in the band in MHz
    # psi_angle is the angel off the shower axis, 
    # equivalent to view angle in deg.
    interp_file = np.load(EFIELD_LUT_file_name)
    return interp_file['efield_interpolator_list'][()]
####################################################################################

def E_field_interp(view_angle_deg, zenith_angle_deg, f_Lo, f_High, log10_tau_energy, distance_km):
    # Lorentzian beam pattern based on 10-MHz filtered subbands of Harm's results
    # Returns electric field peak in V/m
  
    # Since the efields are stored in 10-MHz subbands
    # integrate over the range from f_Lo to f_High in 10-MHz bands
    E_field = np.zeros(len(distance_km))
    #E_field = 0.
    # TODO: Right now forcing the parameters outside the interpolation range to the edges
    # shoudl replace with extrapolation
    z = zenith_angle_deg.copy()
    v = view_angle_deg.copy()
    z[z>87.] = 87.
    z[z<60.] = 60.
    v[v<0.04] = 0.04
    v[v>3.16] = 3.16

    #zenith_angle_deg[zenith_angle_deg > 87.] = 87.
    #zenith_angle_deg[zenith_angle_deg < 60.] = 60.
    #view_angle_deg[view_angle_deg < 0.04] = 0.04
    #view_angle_deg[view_angle_deg > 3.16 ] = 3.16
    #print E_field, zenith_angle_deg, view_angle_deg
    df = 10.
    for freq in np.arange(f_Lo, f_High + df, df):
        i_f_Lo = int(round(freq / df - 1))
        E_field += efield_interpolator_list[i_f_Lo](z,v) * df
	
    # account for ZHAIReS sims only extending to 3.16 deg 
    #TODO: right now sampling a gaussian centered at zero to extrapolate to the wider angles, should verify this parameterization
    E_field[view_angle_deg>3.16] = E_field[view_angle_deg>3.16]*np.exp( -(view_angle_deg[view_angle_deg>3.16]-0.)**2 / (2*3.16)**2)
    
    E_field *= 86.4/distance_km   # distance to tau decay point correction
    E_field *= 10**(log10_tau_energy - 17.) # Energy scaling
    return E_field

####################################################################################

def E_to_V_signal(E_pk, Gain_dB):
    return E_pk * (speed_of_light*1.e3)/300.e6 * np.sqrt(50./377. * pow(10., Gain_dB/10.)/4/np.pi)

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
    #print '**L',L
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
    data_array = f['data_array']
    th_exit_array = f['th_exit_array']
    f.close()
    P_exit = []
    for k in range(0, len(th_exit_array)):
        P_exit.append(float(len(data_array[k]))/tau_exit_num_events)
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
        #print keys[-1], values[-1]
    input_dict = { keys[i] : values[i] for i in range(0,len(keys)) }
    return input_dict

####################################################################################

def A_OMEGA_tau_exit(geom_file_name, LUT_file_name, EFIELD_LUT_file_name, cut_ang, f_Lo, f_High, outTag='test', N=-1):

    # 1. Load geometry file
    GEOM = np.load(geom_file_name)
    GEOM_inputs = parse_input_args(GEOM['input_args'])
    print '\nGEOMETRY FILE INPUT ARGUMENTS\n', 
    for key in GEOM_inputs.keys():
        print '\t',key.ljust(20), GEOM_inputs[key]

    if( N == -1):
    	N_tot = len(GEOM['cos_theta_exit'])
    else:
    	N_tot = int(N)
    print 'Number of Geometry Events to scan', N_tot
    print ''

    cos_theta_exit  = GEOM['cos_theta_exit'][0:N_tot]
    exit_view_angle = GEOM['exit_view_angle'][0:N_tot]
    x_exit          = GEOM['x_exit'][0:N_tot]
    y_exit          = np.zeros(len(x_exit))
    #y_exit          = GEOM['y_exit']
    z_exit          = GEOM['z_exit'][0:N_tot]
    k_x             = GEOM['k_x'][0:N_tot]
    k_y             = GEOM['k_y'][0:N_tot]
    k_z             = GEOM['k_z'][0:N_tot]
    A_geom          = GEOM['A_geom']
    altitude = float(GEOM_inputs['altitude'])
    GEOM.close()
   
    # Set up the detector position
    x_det = 0.
    y_det = 0.
    z_det = altitude + Earth_radius

    # 2. Impose a viewing angle cut
    view_cut = exit_view_angle*180./np.pi<cut_ang
    N_cut = len(cos_theta_exit[view_cut])
    A_Omega = A_geom* float(N_cut) / float(N_tot)
    # TODO: fix to calculate dist_to_detector from decay point rather than from
    # the exit point which is what it is now.
    dist_exit_to_detector = get_distance_to_detector(x_exit, y_exit, z_exit, x_det, y_det, z_det)[view_cut]  #np.sqrt( x_exit**2 + (z_exit - altitude - Earth_radius)**2 ) [view_cut]
    GEOM_theta_exit = np.arccos(cos_theta_exit[view_cut])*180./np.pi
    exit_view_angle = exit_view_angle[view_cut] #radians
    zenith_angle = np.arccos(cos_theta_exit) # radians
    theta_emergence = np.pi/2. - zenith_angle # radians
    
    # 3. set up some arrays to be calculated on the fly
    dist_decay_to_detector = np.zeros(len(dist_exit_to_detector))
    x_decay = np.zeros(len(x_exit[view_cut]))
    y_decay = np.zeros(len(y_exit[view_cut]))
    z_decay = np.zeros(len(z_exit[view_cut]))
    decay_view_angle = np.zeros(len(exit_view_angle))
    X0_dist = np.zeros(len(dist_exit_to_detector))
    log10_tau_energy = np.zeros(len(exit_view_angle))
    
    # 4. Compute the noise in this band
    #       integrating in 10-MHz steps, to match the frequency bins of the peak voltage
    df = 10.# MHz
    noises = noise_voltage(f_Lo, f_High, df,  Z_L, R_L, Gain_dB) # factor of 1e6 because noise temperature is in V/Hz.
    print "Galactic Noise voltage ", noises[1]*1e6, " micro-Volts"
    print "Galactic noise temperature K ", galactic_temperature(np.arange(f_Lo, f_High+df, df))[1], sum(galactic_temperature(np.arange(f_Lo, f_High+df, df))[1])
    print "System noise ", noises[2]
    Noise_Voltage = noises[0]
    print "Noise voltage ", Noise_Voltage*1e6, " micro-Volts, due to G*T_gal + T_sys + T_ant"
    
    # 5. Load Energy Look-up Table
    print "Loading energy look-up table: ", LUT_file_name
    LUT_th_exit, LUT_P_exit, LUT_log10_E_tau = load_tau_LUTs(LUT_file_name)
    P_LUT = np.zeros(len(x_exit[view_cut]))   # zero until it is
    P_range = np.zeros(len(x_exit[view_cut])) # zero until it is proven to decay before it passes the detector
    P_det = np.zeros(len(x_exit[view_cut]))   # zero until it is proven to be detectable


    # 6. Load the Efield interpolator for this altitude
    global efield_interpolator_list
    efield_interpolator_list = load_efield_interpolator(altitude, EFIELD_LUT_file_name)
    
    # 7. Loop Through Geometry Events
    sum_P_exit = 0.
    sum_P_exit_P_range = 0.
    sum_P_exit_P_range_P_det = 0.    
    triggered_events = []
    all_events = []
    for k in range(0,len(GEOM_theta_exit)):
        # 5.1 Get LUT exit angle closest to geometry exit angle
        idx = np.argmin(np.abs(LUT_th_exit - GEOM_theta_exit[k]))
        #print idx, LUT_P_exit[idx], GEOM_theta_exit[k]
        # 5.2 Get tau lepton energy and decay position
        #P_range = 0.                  # zero until it is proven to decay before it passes the detector
        #P_det = 0.                    # zero until it is proven to be detectable
        P_LUT[k] = LUT_P_exit[idx]
        if( P_LUT[k] > 1.e-15):  # make sure the probability of this event is non-zero  
            log10_tau_energy[k] = LUT_log10_E_tau[idx][np.random.randint(0,len(LUT_log10_E_tau[idx]))] # random tau energy
            decay_range = tau_lepton_decay_range(log10_tau_energy[k])                      # estimated decay range
            X0_dist[k] = np.random.exponential(scale=decay_range)                          # sample exponentially distributed decay positions
            x_decay[k], y_decay[k], z_decay[k], decay_view_angle[k], dist_decay_to_detector[k] = decay_point_geom(k_x[k], k_y[k], k_z[k], x_exit[k], y_exit[k], z_exit[k], X0_dist[k], x_det, y_det, z_det)
            # If the event is contained within the range, then the probability is 1.
            # If the shower decays beyond the detector, then it has a negative x-position. The range probability in that case is zero.
            #if( (x_decay[k] < 0.) and (X0_dist[k] < dist_exit_to_detector[k])):
            #    print "Decay beyond detector", x_decay[k], x_exit[k], X0_dist[k], dist_decay_to_detector[k], dist_exit_to_detector[k]
            if((X0_dist[k] < dist_exit_to_detector[k]) and (x_decay[k] > 0.)):
                P_range[k] = 1.
              	# E_field_interp(view_angle_deg, zenith_angle_deg, f_Lo, f_High, log10_tau_energy, distance_km)
	    #elif( (x_decay[k] <= 0.) ):
	    #	print k, " decayed past the detector ", P_range[k], " ", decay_view_angle[k], " ", x_decay[k]
	if(k%100000 == 0 and k>0):
	    print 'Progress: %d events '%k, log10_tau_energy[k]

    # Calculate the electric field by summing the interpolated electric fields in 10-MHz subbands
    # perform array-wise interpolations (outside of the event loop) to minimize number of computations. 
    Peak_E_field = E_field_interp(decay_view_angle*180./np.pi, zenith_angle*180./np.pi, f_Lo, f_High, log10_tau_energy, dist_decay_to_detector)
    #Peak_E_field = E_field_interp(decay_view_angle*180./np.pi, zenith_angle*180./np.pi, f_Lo, f_High, log10_tau_energy, dist_decay_to_detector)
    ####OLD: Peak_E_field = E_field_interp(decay_view_angle*180./np.pi, (np.pi/2. - exit_view_angle)*180./np.pi, f_Lo, f_High, log10_tau_energy, dist_decay_to_detector)
    #####DEBUG: Peak_E_field = np.zeros(len(GEOM_theta_exit)) + Noise_Voltage * 100.;
    for k in range(0,len(GEOM_theta_exit)):
        if( P_LUT[k] > 1.e-15):
            if( P_range[k] == 1.):
                Peak_Voltage = E_to_V_signal(Peak_E_field[k], Gain_dB)
                Peak_Voltage_SNR = Peak_Voltage / Noise_Voltage
                if(Peak_Voltage_SNR > threshold_voltage_snr):
                    #print len(triggered_events), 'Peak_Voltage %1.2e, Noise_Voltage %1.2e, SNR %1.2e, View_angle %1.2f'%(Peak_Voltage, Noise_Voltage, Peak_Voltage_SNR, exit_view_angle[k]*180./np.pi)
                    P_det[k] = 1.
                    triggered_events.append(np.array( [ log10_tau_energy[k], dist_exit_to_detector[k], X0_dist[k], dist_decay_to_detector[k], Peak_E_field[k], Peak_Voltage_SNR, exit_view_angle[k]*180./np.pi, P_LUT[k], GEOM_theta_exit[k], decay_view_angle[k]*180./np.pi,  zenith_angle[k]*180./np.pi ]))
                    #
		    #print LUT_E_tau[idx]
                    #print '%1.2f'%tau_energy
                #else:
                #    print "Didn't trigger on ", 'Peak_Voltage %1.2e, Noise_Voltage %1.2e, SNR %1.2e, View_angle %1.2f'%(Peak_Voltage, Noise_Voltage, Peak_Voltage_SNR, exit_view_angle[k]*180./np.pi)
	    #else:
	    	#print "\t", k, " didn't pass range cut", P_det[k]
        sum_P_exit                += P_LUT[k]
        sum_P_exit_P_range        += P_LUT[k] * P_range[k]
        sum_P_exit_P_range_P_det  += P_LUT[k] * P_range[k] * P_det[k]
	
        if(k%100000 == 0 and k>0):
            #if(k%1== 0 and k>0):
            print 'After %d events: %d events triggered: '%(len(GEOM_theta_exit[0:k]), len(triggered_events))
            print '\t %1.3e km^2 sr'%A_Omega
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit / float(k))
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range / float(k) )
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range_P_det / float(k) )
            print '\t ',np.array(triggered_events).shape
            print ''
        all_events.append([ log10_tau_energy[k], dist_exit_to_detector[k], X0_dist[k], dist_decay_to_detector[k], Peak_E_field[k], Peak_Voltage_SNR, exit_view_angle[k]*180./np.pi, P_LUT[k], GEOM_theta_exit[k], decay_view_angle[k]*180./np.pi,  zenith_angle[k]*180./np.pi ])

    print 'After all %d events, %d events triggered: '%(len(GEOM_theta_exit), len(triggered_events))
    print '\t %1.3e km^2 sr'%(A_Omega)
    print '\t %1.3e km^2 sr, %d events exited Earth'%(A_Omega*sum_P_exit/float(N_cut),N_tot*sum_P_exit/float(N_cut) )
    print '\t %1.3e km^2 sr, %d events decayed before detector'%(A_Omega*sum_P_exit_P_range/float(N_cut), N_tot*sum_P_exit_P_range/float(N_cut))
    print '\t %1.3e km^2 sr, %d events triggered'%(A_Omega*sum_P_exit_P_range_P_det/float(N_cut), N_tot*sum_P_exit_P_range_P_det/float(N_cut))

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
                            log10_tau_energy = log10_tau_energy,
                            triggered_events = np.array(triggered_events),
                            )
                            #all_events = np.array(all_events))

    print "Wrote ", outTag+'.npz and ', outTag+'_events.npz'

    exit()

