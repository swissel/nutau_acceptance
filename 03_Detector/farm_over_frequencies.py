import sys
import os
import re
import numpy as np
#/u/nobackup/bern/swissel/anaconda2/bin/python tau_Shower_Sim_TAU.py
#-G $TAU_ACC_GEOM_DIR/geom_3_km_5_deg_10M_ev.npz -L $TAU_ACC_LUT_DIR/0.0km_ice_lowCS_lowEL/LUT_1e+17_eV.npz -e $TAU_ACC_ZHAIRES_DIR/interpolator_efields_3.npz -l 100 -s 1000 -o $TAU_ACC_DET_DIR/detector_acceptance_altitude_3_km_0.0fkm_ice_lowCS_lowEL_1e+17_eV_10-1000MHz.npz
#setenv TAU_DIR="/u/project/bern/swissel/nutau"
#setenv TAU_ACC_GEOM_DIR="$TAU_DIR/geom_files"
#setenv TAU_ACC_ZHAIRES_DIR="$TAU_DIR/efield_files/"
#setenv TAU_ACC_DET_DIR="$TAU_DIR/det_files"
#setenv TAU_ACC_LUT_DIR="$TAU_DIR/tau_energy_lut/"

def edit_cluster_script(altitude, energy, f_Lo, f_High):
	dir = '/u/project/bern/swissel/nutau/nutau_acceptance/03_Detector/'
	os.chdir(dir)

        #edit the run command
        newlines = []
        prog1 = re.compile('  setenv ALTITUDE \d')
        prog2 = re.compile('  setenv FLO \d')
        prog3 = re.compile('  setenv FHIGH \d')
        prog4 = re.compile('  setenv ENERGY .*')
      
        with open("tau_Shower_Sim_TAU.py.cmd", 'r') as fi:
                for l in fi:
                        if( prog1.match(l)):
                                nl = '  setenv ALTITUDE %d\n'%altitude
                                print "Replacing :", nl,
                                newlines.append(nl)
                        elif(prog2.match(l)):
                                nl = '  setenv FLO %d\n'%f_Lo
                                print "Replacing :", nl,
                                newlines.append(nl)
                        elif(prog3.match(l)):
                                nl = '  setenv FHIGH %d\n'%f_High
                                print "Replacing :", nl,
                                newlines.append(nl)
                        elif(prog4.match(l)):
                                nl = '  setenv ENERGY "%de+%d"\n'%(energy[0], energy[1])
                                print "Replacing :", nl,
                                newlines.append(nl)
                        else:
                                newlines.append(l)

        with open("tau_Shower_Sim_TAU.py.cmd",'w') as fi:
                for l in newlines:
                        fi.write(l)

altitude = 3

for f_Lo in np.logspace(1, 3, 10):
	for BW in [30,100,300,1000]:
		for en1 in [1,3]:
			for en2 in range(15, 22):
				if( f_Lo + BW > 1660 or (en1 == 3 and en2 == 21)):
					continue
				else:
					energy = "%de+%d"%(en1, en2)
					f_High = f_Lo + BW
					edit_cluster_script(altitude, [en1,en2], f_Lo, f_High)
					cmd = "qsub tau_Shower_Sim_TAU.py.cmd"
					print cmd, "\n"
					os.system(cmd)
			
