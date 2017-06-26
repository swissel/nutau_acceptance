import os

altitude = 2
energies = ["1e+15", "3e+15", "1e+16", "3e+16", "1e+17", "3e+17", "1e+18", "3e+18", "1e+19", "3e+19", "1e+20", "3e+20", "1e+21"][::-1]
print energies

ice_thickness = 0.0
cross_section = "low"
elong_rate = "std"

for energy in energies:
	geom_file = os.environ['TAU_ACC_GEOM_DIR'] + '/geom_%d_km_5_deg_10M_ev.npz'%altitude
	lut_file  = os.environ['TAU_ACC_LUT_DIR'] + '/%1.1fkm_ice_%sCS_%sEL/LUT_%s_eV.npz'%(ice_thickness,cross_section, elong_rate,energy)
	efield_file = os.environ['TAU_ACC_ZHAIRES_DIR'] + '/interpolator_efields_%dkm.npz'%(altitude)

	for f_Lo in range(10, 1000, 30):
		for BW in [30,100,300,1000]:
			if( f_Lo + BW > 1660):
				continue
			else:
				out_tag = os.environ['TAU_ACC_DET_DIR'] + '/detector_acceptance_altitude_%d_km_%1.1fkm_ice_%sCS_%sEL_%s_eV_%d-%dMHz'%(altitude,ice_thickness,cross_section,elong_rate, energy, f_Lo, f_Lo+BW)
				cmd = "python tau_Shower_Sim_TAU.py -G %s -L %s -e %s -l %d -s %s -o %s"%(geom_file, lut_file, efield_file, f_Lo, f_Lo+BW, out_tag)
				print cmd
				os.system(cmd)
			
