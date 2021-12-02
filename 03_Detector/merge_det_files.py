import os
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

def read_det_acceptance_files(finame1, finame2, en1, en2):
    print("Reading ", finame1)
    fi1 = np.load(finame1)

    # acceptance calculated by the monte carlo at various levels
    # maximum geometric acceptance in km^2 sr
    A_Omega_start = fi1['A_Omega_start']
    # acceptance that includes the exit probability in km^2 sr
    A_Omega_exit  = fi1['A_Omega_exit']
    # acceptance that includes the tau decay range probabilitiy and tau exit probability in km^2 sr
    A_Omega_range = fi1['A_Omega_range']
    # acceptance that includes the tau exit, decay range probabilities and probability that the event triggers the detector in km^2 sr
    A_Omega_trig  = fi1['A_Omega_trig']
    # the energy in log( energy / eV)
    log_energy = np.log10(en1*pow(10., en2))
    # the number of events triggered
    N_triggered = fi1['N_triggered']
    # the number of events thrown in the Monte Carlo
    N_events_start = fi1['N_events_start']

    print("Reading ", finame2)
    fi2 = np.load(finame2)
    # the noise voltage from the sum of the noise powers from the galaxy, the amplitired, and thermal noise
    noise_voltage = fi2['noise_voltage']
    # an array of useful information about the events that triggered, see below
    triggered_events = fi2['triggered_events']
    # the binned emergence angles used to calculate the differential acceptance
    emerge_angle_bins = fi2['emerge_angle_bins']
    # the binned number of sims in the emergence angle bins
    nsim_emerge_angle = fi2['nsim_emerge_angle']

    return A_Omega_start, A_Omega_exit, A_Omega_range, A_Omega_trig, log_energy, N_triggered, N_events_start, noise_voltage, triggered_events, emerge_angle_bins, nsim_emerge_angle

def read_2det_acceptance_files(finame1, finame2, en1, en2):
    print("Reading ", finame1)
    fi1 = np.load(finame1)

    # acceptance calculated by the monte carlo at various levels
    # maximum geometric acceptance in km^2 sr
    A_Omega_start = fi1['A_Omega_start']
    # acceptance that includes the exit probability in km^2 sr
    A_Omega_exit  = fi1['A_Omega_exit']
    # acceptance that includes the tau decay range probabilitiy and tau exit probability in km^2 sr
    A_Omega_range = fi1['A_Omega_range']
    # acceptance that includes the tau exit, decay range probabilities and probability that the event triggers the detector in km^2 sr
    A_Omega_trig  = fi1['A_Omega_trig']
    # the energy in log( energy / eV)
    log_energy = np.log10(en1*pow(10., en2))
    # the number of events triggered
    N_triggered = fi1['N_triggered']
    # the number of events triggered the second detector
    N_triggered2 = fi1['N_triggered2']
    # the number of events thrown in the Monte Carlo
    N_events_start = fi1['N_events_start']
    # The number of events that triggered both detectors
    N_Det1Det2  = fi1['N_Det1Det2']

    print("Reading ", finame2)
    fi2 = np.load(finame2)
    # the noise voltage from the sum of the noise powers from the galaxy, the amplitired, and thermal noise
    noise_voltage = fi2['noise_voltage']
    # an array of useful information about the events that triggered, see below
    triggered_events = fi2['triggered_events']
    # an array of useful information about the events that triggered, see below
    triggered_events2 = fi2['triggered_events2']
    # the binned emergence angles used to calculate the differential acceptance
    emerge_angle_bins = fi2['emerge_angle_bins']
    # the binned number of sims in the emergence angle bins
    nsim_emerge_angle = fi2['nsim_emerge_angle']

    return A_Omega_start, A_Omega_exit, A_Omega_range, A_Omega_trig, log_energy, N_triggered, N_triggered2, N_Det1Det2, N_events_start, noise_voltage, triggered_events, triggered_events2, emerge_angle_bins, nsim_emerge_angle

def get_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain,
                        nphased, cross, eloss, threshold, start_event):
        finame1 = os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr_%d.npz"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold, start_event)
        finame2 = os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr_%d_events.npz"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold, start_event)
        finameout =  os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold)

        return finame1, finame2, finameout

def get_2det_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain,
                        nphased, cross, eloss, threshold, detsep, start_event):
        finame1 = os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr_detsep%1.1f_%d.npz"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold, detsep, start_event)
        finame2 = os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr_detsep%1.1f_%d_events.npz"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold, detsep, start_event)
        finameout =  os.environ['TAU_ACC_DET_DIR'] + "/detector_acceptance_altitude_%d_km_%1.1fkm_ice_midCS_stdEL_%de+%d_eV_%d-%dMHz_%2.1fdBi_%1.1fantennas_%2.1fsnr_detsep%1.1f"%(
                  altitude,icethick,energy1,energy2,f_Lo,f_High,gain,nphased, threshold, detsep)

        return finame1, finame2, finameout


def combine_det_files(nfiles, nstep, energy1, energy2, altitude, icethick, f_Lo, f_High, gain,
                        nphased, cross, eloss, threshold, remove_files=False):
    A_Omega_start = 0
    A_Omega_exit = 0
    A_Omega_range = 0
    A_Omega_trig = 0
    N_triggered = 0
    N_events_start = 0 
    triggered_events = np.array([])
    for start_event in range(0, (nfiles+1)*nstep, nstep):
        finame1, finame2, outTag = get_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain, nphased, cross, eloss, threshold, start_event)
        iA_Omega_start, iA_Omega_exit, iA_Omega_range, iA_Omega_trig, ilog_energy, iN_triggered, iN_events_start, inoise_voltage, itriggered_events, iemerge_angle_bins, insim_emerge_angle = read_det_acceptance_files(finame1, finame2, energy1, energy2)
        
        A_Omega_start += iA_Omega_start/nfiles
        A_Omega_exit += iA_Omega_exit/nfiles
        A_Omega_range += iA_Omega_range/nfiles
        A_Omega_trig += iA_Omega_trig/nfiles
        N_triggered += iN_triggered
        N_events_start += iN_events_start
	
        log_energy = ilog_energy
        noise_voltage = inoise_voltage
        triggered_events = np.array(list(triggered_events) + list(itriggered_events))
        emerge_angle_bins = iemerge_angle_bins
        nsim_emerge_angle = insim_emerge_angle

    print("Writing %s.npz and %s_events.npz"%( outTag, outTag))
    np.savez(outTag+'.npz', A_Omega_start    = A_Omega_start,
                            A_Omega_exit     = A_Omega_exit,
                            A_Omega_range    = A_Omega_range,
                            A_Omega_trig     = A_Omega_trig,
                            N_events_start   = N_events_start,
                            N_triggered      = N_triggered
			    )

    np.savez(outTag+'_events.npz', 
                            A_Omega_start    = A_Omega_start,
                            A_Omega_exit     = A_Omega_exit,
                            A_Omega_range    = A_Omega_range,
                            A_Omega_trig     = A_Omega_trig,
			    noise_voltage    = noise_voltage,
                            triggered_events = np.array(triggered_events),
			    emerge_angle_bins = emerge_angle_bins,
			    nsim_emerge_angle = nsim_emerge_angle,
			    )
    if( remove_files == True):
        for start_event in range(0, (nfiles+1)*nstep, nstep):
            finame1, finame2, outTag = get_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain, nphased, cross, eloss, threshold, start_event)
            print("rm "+finame1)
            print("rm "+finame2)
            os.remove(finame1)
            os.remove(finame2)

def combine_2det_files(nfiles, nstep, energy1, energy2, altitude, icethick, f_Lo, f_High, gain,
                        nphased, cross, eloss, threshold, detsep, remove_files=False):
    A_Omega_start = 0
    A_Omega_exit = 0
    A_Omega_range = 0
    A_Omega_trig = 0
    N_triggered = 0
    N_events_start = 0 
    N_triggered2 = 0
    N_triggered_det1det2 = 0
    triggered_events = np.array([])
    triggered_events2 = np.array([])
    for start_event in range(0, (nfiles+1)*nstep, nstep):
        finame1, finame2, outTag = get_2det_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain, nphased, cross, eloss, threshold, detsep, start_event)
        iA_Omega_start, iA_Omega_exit, iA_Omega_range, iA_Omega_trig, ilog_energy, iN_triggered, iN_triggered2,iN_triggered_det1det2, iN_events_start, inoise_voltage, itriggered_events, itriggered_events2, iemerge_angle_bins, insim_emerge_angle = read_2det_acceptance_files(finame1, finame2, energy1, energy2)
        
        A_Omega_start += iA_Omega_start/nfiles
        A_Omega_exit += iA_Omega_exit/nfiles
        A_Omega_range += iA_Omega_range/nfiles
        A_Omega_trig += iA_Omega_trig/nfiles
        N_triggered += iN_triggered
        N_triggered2 += iN_triggered2
        N_triggered_det1det2 += iN_triggered_det1det2
        N_events_start += iN_events_start
	
        log_energy = ilog_energy
        noise_voltage = inoise_voltage
        triggered_events = np.array(list(triggered_events) + list(itriggered_events))
        triggered_events2 = np.array(list(triggered_events2) + list(itriggered_events2))
        emerge_angle_bins = iemerge_angle_bins
        nsim_emerge_angle = insim_emerge_angle

    print("Writing %s.npz and %s_events.npz"%( outTag, outTag))
    np.savez(outTag+'.npz', A_Omega_start    = A_Omega_start,
                            A_Omega_exit     = A_Omega_exit,
                            A_Omega_range    = A_Omega_range,
                            A_Omega_trig     = A_Omega_trig,
                            N_events_start   = N_events_start,
                            N_triggered      = N_triggered,
			    N_triggered2     = N_triggered2,
			    N_Det1Det2 = N_triggered_det1det2
			    )

    np.savez(outTag+'_events.npz', 
                            A_Omega_start    = A_Omega_start,
                            A_Omega_exit     = A_Omega_exit,
                            A_Omega_range    = A_Omega_range,
                            A_Omega_trig     = A_Omega_trig,
			    noise_voltage    = noise_voltage,
                            triggered_events = np.array(triggered_events),
			    triggered_events2= np.array(triggered_events2),
			    emerge_angle_bins = emerge_angle_bins,
			    nsim_emerge_angle = nsim_emerge_angle,
			    )
    if( remove_files == True):
        for start_event in range(0, (nfiles+1)*nstep, nstep):
            finame1, finame2, outTag = get_2det_filenames(energy1, energy2, altitude, icethick, f_Lo, f_High, gain, nphased, cross, eloss, threshold, detsep, start_event)
            print("rm "+finame1)
            print("rm "+finame2)
            os.remove(finame1)
            os.remove(finame2)



def merge_files():
    nfiles = 10
    nstep = 1000000
    altitude = 37
    ice_thick = 2.0
    start_frequency = 200
    stop_frequency = 1200
    gain = 10.0 
    phased = 16.0
    threshold = 5.0
    cross = 'mid'
    eloss = 'std'
    remove_files = True

    for energy1 in [1,3]:
        for energy2 in range(15,22, 1):
             if energy1 == 3 and energy2 == 21:
                continue
             else:
                combine_det_files(nfiles, nstep, energy1, energy2, altitude, ice_thick, start_frequency,
                               stop_frequency, gain, phased, cross, eloss, threshold, remove_files=remove_files)

def merge_2det_files():
    nfiles = 10
    nstep = 1000000
    altitude = 3
    ice_thick = 0
    start_frequency = 200
    stop_frequency = 1200
    gain = 10.0 
    phased = 100.0
    threshold = 5.0
    cross = 'mid'
    eloss = 'std'
    detsep = 1.0
    remove_files = False

    for detsep in [0.5,1.0, 2.0,3.0,4.0,5.0,6.0,7.0]:
        for energy1 in [1]:
            for energy2 in range(18, 19, 1):
                if energy1 == 3 and energy2 == 21:
                    continue
                else:
                    combine_2det_files(nfiles, nstep, energy1, energy2, altitude, ice_thick, start_frequency,
                               stop_frequency, gain, phased, cross, eloss, threshold, detsep, remove_files=remove_files)




if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='Combines the output files from multiple sims')
    parser.add_argument("-e",    "--energy", default=1e19, help='energy in eV', type=float)
    parser.add_argument("-a",    "--altitude", default=3, help='altitude in km', type=float)
    parser.add_argument("-i",    "--ice_thick", default=0, help='ice_thick in km', type=float)
    parser.add_argument("-c",    "--cross", default='mid', help='cross section [low, mid, high]', type=str)
    parser.add_argument("-E",    "--eloss", default='std', help='energy loss [std, low]', type=str)
    parser.add_argument("-l",    "--start_frequency", default = 200., help="starting frequency (MHz)", type=float)
    parser.add_argument("-s",    "--stop_frequency", default = 1200., help="stopping frequency (MHz)", type=float)
    parser.add_argument("-g", "--gain", default = 10., help='gain of antenna in dBi; if phasing antennas, this is the gain of each antenna', type=float) 
    parser.add_argument("-p", "--phased", default=1, help='number of phased antennas', type=float) 
    parser.add_argument("-t", "--threshold", default=5.0, help='threshold SNR for trigger', type=float)
    parser.add_argument("-f",    "--nfiles", default=0, help='number of files', type=int)
    parser.add_argument("-n",    "--nstep", default=0, help='number of events in each sim', type=int)

    args = parser.parse_args()
    energy1 = float(str(args.energy)[0])
    energy2 = float(str(args.energy)[-2:])
    
    #combine_det_files(args.nfiles, args.nstep, energy1, energy2, args.altitude, args.ice_thick, args.start_frequency, args.stop_frequency, args.gain,
    #                    args.phased, args.cross, args.eloss, args.threshold)

    merge_files()
    #merge_2det_files()
