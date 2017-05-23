import numpy as np

tau_exit_num_events = 1.e6

# Define constants
Earth_radius   = 6371.      # km
tau_life_time  = 2.906e-13  # seconds
tau_mass       = 1776.82e6  # eV/c2
speed_of_light = 299792.458 # km / second

# Threshold voltage based on temperature:
T_ant  = 105.     # Kelvin (based on elog 655 entry by P. Gorham)
T_sys  = 140.     # System temperature at the lower part of the ANITA band from elog 119
kB     = 1.38e-23 # in Volts^2/Ohm/Hz/Kelvin
Z_load = 50.      # 50 Ohm system
BW     = 400.e6   # bandwidth for first two ANITA bands.
#threshold_voltage = 8.9e-06 # V, this is 8.9 microVolts for sky temp 290 Kelvin T_sys
#threshold_voltage =np.sqrt(kB*(T_ant+T_sys)*Z_load*BW) # V, this is 8.22 microVolts for the parameters given above

threshold_voltage = 72.e-6 # V, this is based on taking half the peak field for the weakest ANITA-1 EAS candidate. (0.466 mV/m convolved with antenna effective height at 300 MHz). 

print 'threshold_voltage', threshold_voltage
####################################################################################

def tau_lepton_decay_range(tau_lepton_energy):
    return speed_of_light * 10**(tau_lepton_energy)  / tau_mass * tau_life_time

####################################################################################

def E_field_model(view_angle_deg, log10_tau_energy, distance_km):
    # Lorentzian beam pattern based on Harm's results
    # Returns electric field peak in V/m
    
    # This is the older Lorentz fit to a lateral distribution function that was actually at 75 degree zenith exit angle rather than 65.
    #E_0 =  1.2e-4*10**(log10_tau_energy - 17.) # overall amplitude in ANITA 200 - 1200 MHz band
    #E_0 *= 85.2/distance_km # corrected for the distance the shower was simulated at.
    #arg =  (view_angle_deg - 1.1)/0.14 
    #return E_0 * 1. / (1+arg**2)

    # This is a more detailed fit to the lateral distribution function at 65 for a 200 - 1200 MHz band.
    E_0  = 0.183*1.e-3     # V/m at ground level for 10^17 eV tau lepton.
    peak_view_angle = 1.01 # degrees
    width = 0.160          # degrees
    E_1 = 0.00385*1.e-3    # V/m at ground level for 10^17 eV tau lepton.
    width2 = 1.14          # degrees
    v = (view_angle_deg-peak_view_angle)/width 
    frac_gauss = 0.825 # fraction of Gaussian peak ( 1-frac_gauss is lorentzian)
    # The E_field model is a linear combination of a Gaussian and a Lorentzian centered at the peak view angle. A Gaussian centered at 0 deg is added to account for the central lobe. 
    E_field = E_0*(frac_gauss*np.exp(-v**2/2.) + (1-frac_gauss)/(1+v**2)) + E_1*np.exp(-view_angle_deg**2/2./width2**2)
    E_field *= 86.4/distance_km   # distance to tau decay point correction
    E_field *= 10**(log10_tau_energy - 17.) # Energy scaling
    return E_field

####################################################################################

def E_to_V_signal(E_pk):
    return E_pk * (speed_of_light*1.e3)/300.e6 * np.sqrt(50./377. * 10./4/np.pi)

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
    x_pd = x_pos - x_det
    y_pd = y_pos - y_det
    z_pd = z_pos - z_det
    cos_view_angle = -( k_x * x_pd + k_y * y_pd + k_z * z_pd ) / np.sqrt(x_pd**2 + y_pd**2 + z_pd**2)
    return np.arccos(cos_view_angle)

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

def A_OMEGA_tau_exit(geom_file_name, LUT_file_name, cut_ang, outTag='test'):

    # 1. Load geometry file
    GEOM = np.load(geom_file_name)
    GEOM_inputs = parse_input_args(GEOM['input_args'])
    print '\nGEOMETRY FILE INPUT ARGUMENTS\n', 
    for key in GEOM_inputs.keys():
        print '\t',key.ljust(20), GEOM_inputs[key]
    cos_theta_exit  = GEOM['cos_theta_exit']
    exit_view_angle = GEOM['exit_view_angle']
    x_exit          = GEOM['x_exit']
    z_exit          = GEOM['z_exit']
    A_geom          = GEOM['A_geom']
    altitude = float(GEOM_inputs['altitude'])
    GEOM.close()
    N_tot = len(cos_theta_exit)
    print 'Number of Geometry Events', N_tot
    print ''

    # 2. Impose a viewing angle cut
    view_cut = exit_view_angle*180./np.pi<cut_ang
    N_cut = len(cos_theta_exit[view_cut])
    A_Omega = A_geom* float(N_cut) / float(N_tot)
    dist_to_detector = np.sqrt( x_exit**2 + (z_exit - altitude - Earth_radius)**2 ) [view_cut]
    GEOM_theta_exit = np.arccos(cos_theta_exit[view_cut])*180./np.pi
    exit_view_angle = exit_view_angle[view_cut]
    # 3. Load Energy Look-up Table
    LUT_th_exit, LUT_P_exit, LUT_log10_E_tau = load_tau_LUTs(LUT_file_name)

    # 4. Loop Through Geometry Events  
    sum_P_exit = 0.
    sum_P_exit_P_range = 0.
    sum_P_exit_P_range_P_det = 0.    
    triggered_events = []
    for k in range(0,len(GEOM_theta_exit)):
        # 4.1 Get LUT exit angle closest to geometry exit angle
        idx = np.argmin(np.abs(LUT_th_exit - GEOM_theta_exit[k]))

        # 4.2 Get tau lepton energy and decay position
        P_range = 0.                  # zero until it is proven to decay before it passes the detector
        P_det = 0.                    # zero until it is proven to be detectable
        if( LUT_P_exit[idx] > 1.e-15):  # make sure the probability of this event is non-zero
            log10_tau_energy = LUT_log10_E_tau[idx][np.random.randint(0,len(LUT_log10_E_tau[idx]))] # random tau energy
            decay_range = tau_lepton_decay_range(log10_tau_energy)                      # estimated decay range
            X0_dist = np.random.exponential(scale=decay_range)                    # sample exponentially distributed decay positions            
            if(X0_dist < dist_to_detector[k]): # If the event is contained the range probability is 1. 
                P_range = 1.
                #print exit_view_angle[k]*180./np.pi
                Peak_E_field = E_field_model(exit_view_angle[k]*180./np.pi, log10_tau_energy, dist_to_detector[k])
                Peak_Voltage = E_to_V_signal(Peak_E_field)
                if(Peak_Voltage > 3.*threshold_voltage):
                    #print 'Peak_Voltage %1.2e, view_angle %1.2f'%(Peak_Voltage, exit_view_angle[k]*180./np.pi)
                    P_det = 1.
                    triggered_events.append(np.array( [ log10_tau_energy, X0_dist, dist_to_detector[k], Peak_E_field, exit_view_angle[k]*180./np.pi, LUT_P_exit[idx], GEOM_theta_exit[k] ]))
            #print '%1.1f\t%1.1f\t%1.1f\t%1.1f\t%d\t%1.1f'%(log10_tau_energy, decay_range, X0_dist, dist_to_detector[k], X0_dist < dist_to_detector[k], Peak_Voltage/threshold_voltage)
        #print LUT_E_tau[idx]
        #print '%1.2f'%tau_energy
        sum_P_exit                += LUT_P_exit[idx] 
        sum_P_exit_P_range        += LUT_P_exit[idx] * P_range 
        sum_P_exit_P_range_P_det  += LUT_P_exit[idx] * P_range * P_det
        if(k%100000 == 0 and k>0):
            print '\t %1.3e km^2 sr'%A_Omega
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit / float(k))
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range / float(k))
            print '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range_P_det / float(k))
            print '\t ',np.array(triggered_events).shape
            print ''

    '\t %1.3e km^2 sr'%A_Omega
    '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit/float(N_cut))
    '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range/float(N_cut))
    '\t %1.3e km^2 sr'%(A_Omega*sum_P_exit_P_range_P_det/float(N_cut))

    np.savez(outTag+'.npz', A_Omega_start    = A_Omega,
                            A_Omega_exit     = A_Omega*sum_P_exit/float(N_cut),
                            A_Omega_range    = A_Omega*sum_P_exit_P_range/float(N_cut),
                            A_Omega_trig     = A_Omega*sum_P_exit_P_range_P_det/float(N_cut),
                            triggered_events = np.array(triggered_events))

    exit()

