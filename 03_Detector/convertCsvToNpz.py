import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as pyp

def read_harm_sims(antenna_height, zenith, orientation, offangle, domain ):
    ''' Parameters:
	antenna_height: Height of detector in km. Sims are available for 1, 2, 3, 37 km. 
	orientation: x,y,z
	zenith: zenith angle. complement of the emergence angle
	offangle: angle from the shower axis. Sims are available in the range 0.1-4.0 deg in 0.1 deg bins. Input in decimal degrees.
	domain: frequency or time domain (freq or time)
    '''
    predirc = "/Users/wissels/Dropbox/MountainTop/harms_sims/"
    dirc = predirc + "%dkm/paper/DeltaPsi0.04/DecayID10/DecayHeight0/AntennaHeight%d/Zenith%d/Energy1e17/csv/"%(antenna_height, antenna_height*1000, zenith)
   
    if( domain == 'time'):
	finame = dirc + "antenna_%d_time_%s.csv"%(int(offangle/0.04)+1, orientation)
	print offangle, finame
   	return pd.read_csv(finame, names=['time_s','v_per_m'])
    elif( domain == 'freq'):
	finame = dirc + "antenna_%d_freq_%s.csv"%(int(offangle/0.04)+1, orientation)
   	return pd.read_csv(finame, names=['freq_hz','volts_per_m_hz'])

def read_efields(antenna_height, zenith):
    orientations = ['x','y', 'z']
    offangles = np.arange(0.0,0.04*80, 0.04)

    efield_td = {}
    efield_fd = {}
    for i_offang, offang in enumerate(offangles):
	eoffang_td = {}
	eoffang_fd = {}
	
	for pol in orientations:
		e_td = read_harm_sims(antenna_height, zenith, pol, offang, 'time')
		e_fd = read_harm_sims(antenna_height, zenith, pol, offang, 'freq')
		eoffang_td['time_s'] = e_td['time_s'].values # the time bases are the same for x,y,z
		eoffang_td["%s_v_per_m"%pol] = e_td['v_per_m'].values
		eoffang_fd['freq_hz'] = e_fd['freq_hz'].values # df's are the same for x,y,z
		eoffang_fd['%s_v_per_m_Hz'%pol] = e_fd['volts_per_m_hz'].values
	efield_td[i_offang] = eoffang_td
	efield_fd[i_offang] = eoffang_fd
	
    return efield_td, efield_fd

def write_npz_files(dirc, efield_td, efield_fd, antenna_height, zenith):
    finame = dirc + "/npz_files/altitude_%dkm_zenith_%d.npz"%(antenna_height, zenith)
    np.savez(finame, efield_td=efield_td, efield_fd=efield_fd)

def convert_csv_to_npz(antenna_height, zenith):
    efield_td, efield_fd = read_efields(antenna_height, zenith)
    write_npz_files( "/Users/wissels/Dropbox/MountainTop/harms_sims/", efield_td, efield_fd, antenna_height, zenith)

def read_npz_files(antenna_height, zenith):
    npzfile = np.load("/Users/wissels/Dropbox/MountainTop/harms_sims/npz_files/altitude_%dkm_zenith_%d.npz"%(antenna_height, zenith))
    efield_td = npzfile['efield_td'][()] ## for some reason just reading the arrays gives a 0-dim array, so skip this.
    efield_fd = npzfile['efield_fd'][()]
    return efield_td, efield_fd

def check_npz_files(antenna_height, zenith):
    efield_td, efield_fd = read_npz_files(h,z)
    off_angle_array = np.arange(0.0,80*0.04, 0.04)
    i_off_angle_array = range(0, len(off_angle_array))
    for i_offangle in i_off_angle_array:
        efield_td, efield_fd = read_npz_files(35,89)
        efield = efield_td[i_offangle]['y_v_per_m']
        time = efield_td[i_offangle]['time_s']
        pyp.plot(time, efield)
    pyp.show()

###################################
if __name__ == '__main__':
    antenna_heights = [37] #[1,2,3,35]
    zenith_angles = [55,60,65,70,75,80,85,87,89]
    #read_efields(1, 60)
    #convert_csv_to_npz(1, 60)
    for h in antenna_heights:
	for z in zenith_angles:
		convert_csv_to_npz(h, z)
		efield_td, efield_fd = read_npz_files(h,z)
		print "Height : ", h, " km ", " Zenith angle ", z
		print efield_td
		print efield_fd
		print "-------------------------------------------"	
