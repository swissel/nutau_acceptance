import os

exe        = "/home/romerowo/nutau_acceptance/applications/ANITA_upper_bound/ANITA_upper_bound.py"
geom_file  = "/home/romerowo/nutau_acceptance/outputs/GEOM/geom_37_km_5_deg_10M_ev.npz"
LUT_dir    = "/home/romerowo/nutau_acceptance/data/tau_exit_LUTs"
out_dir    = "/aurora_nobackup/eva/romerowo/nutau_outputs/20170508"
inp_dir    = "/home/romerowo/nutau_acceptance/applications/ANITA_upper_bound/inputs"

if not os.path.isdir(inp_dir): os.makedirs(inp_dir)
if not os.path.isdir(out_dir): os.makedirs(out_dir)

count = 0
for ice_thickness in [0.0, 1.0, 2.0, 3.0, 4.0]:
    for CS_model in ['midCS', 'lowCS', 'uppCS']:
        for EL_model in ['stdEL', 'lowEL']:
            for energy in [1.e17, 3.e17, 1.e18, 3.e18, 1.e19, 3.e19, 1.e20, 3.e20, 1.e21]:
                LUT_file = "%s/%1.1fkm_ice_%s_%s/LUT_%1.0e_eV.npz"%(LUT_dir, ice_thickness, CS_model, EL_model, energy)
                out_tag  = "%s/ice_%1.1fkm_alt_37km_%1.0e_eV_%s_%s"%(out_dir, ice_thickness, energy, CS_model, EL_model)
                com = exe + " -G %s"%geom_file + " -L %s"%LUT_file + ' -c 5.' + " -out_tag %s"%out_tag + " > %s.log"%out_tag
                print com
                inp_fnm = "%s/input.%d"%(inp_dir,count)
                print inp_fnm
                fout = open(inp_fnm,"w")
                fout.write(com)
                fout.close()
                count+=1


