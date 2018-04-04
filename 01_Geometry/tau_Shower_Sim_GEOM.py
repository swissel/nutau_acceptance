#!/usr/bin/env python
'''
#!/u/nobackup/bern/swissel/anaconda2/bin/python
'''
import sys
sys.path.append("/home/romerowo/nutau_acceptance/02_TauDecay/")
sys.path.append("/home/romerowo/nutau_acceptance/src/")

import matplotlib
matplotlib.use('Agg')
from pylab import *
import argparse
#from tau_Shower_Sim_lib import *
import tau_Shower_Sim_lib as TSS
import os

if __name__ == "__main__":

  parser=argparse.ArgumentParser(description='Script the simulates tau neutrino induced tau lepton air shower detection from antennas at altitude')
  parser.add_argument("-nev",     "--num_events",        default = 10000000, help='number of events to simulate',         type=int)
  parser.add_argument("-tnev",    "--target_num_events", default = 1000000,  help='number of events to simulate',         type=int)
  parser.add_argument("-alt",     "--altitude",          default = 4.,       help='number of events to simulate',         type=float)
  parser.add_argument("-vc",      "--theta_view_cut",    default = 5.,       help='exit point view angle cut in degrees', type=float)
  parser.add_argument("-Er",      "--Earth_radius",      default = 6371.,    help='radius of Earth in km',                type=float)
  parser.add_argument("-out_tag", "--output_file_tag",   default = 'test',   help='tag for output files',                 type=str)

  #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
  for i, arg in enumerate(sys.argv):
    if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

  args=parser.parse_args()

  TSS.Earth_radius = args.Earth_radius # km
  print 'IN GEOM: TSS.Earth_radius', TSS.Earth_radius
  print ''
  print 'Simulating and observer at %1.2f km altitude'%(args.altitude)
  print 'Sampling %1.0e events per iteration'%args.num_events
  print 'Accepting only events with theta_view < %1.2f degrees'%args.theta_view_cut
  print ''
  # Calculate zeroth order acceptance
  # This is the acceptance you would get if you could see all the particles exiting from the ground.
  print ''
  print 'The acceptance for all isotropically distributed shower axes poking up the surface of the Earth is:'
  A0 = 2.*pi**2 * TSS.Earth_radius * args.altitude / (1. +  args.altitude / TSS.Earth_radius )
  print 'A0', A0, 'km^2 sr'

  ev_count = 0
  dat = TSS.GEOM(args.altitude, args.num_events, view_angle_cut_deg = args.theta_view_cut)
  dat = np.array(dat)
  ev_count += dat.shape[1]
  print ''
  print 'Iteration 0'
  print 'Hits %d of %1.0e'%(ev_count, args.num_events)


  iterations = 1 # just did it once before this loop to initizalize values, so start at 1 
  while(ev_count < args.target_num_events):
    #x_exit, z_exit, cos_th_E, sin_th_E, k_x, k_y, k_z, cos_theta_exit, view_angle 
    out = TSS.GEOM(args.altitude, args.num_events, view_angle_cut_deg = args.theta_view_cut)
    out = np.array(out)
    dat = np.concatenate([dat, out], axis=1)
    #print out.shape, dat.shape
    ev_count = dat.shape[1]
    print ''
    print 'Iteration %d'%(iterations)
    if(ev_count < 10000):
        print 'Hits %d of %1.0e'%(ev_count, (iterations+1)*args.num_events)
    else:
        print 'Hits %1.2e of %1.0e'%(ev_count, (iterations+1)*args.num_events)
    iterations += 1

    A =  A0 * float(ev_count) / float(iterations*args.num_events)
    print 'A %1.3e km^2 sr'%A
    # INTERMEDIATE SAVE IN THE WHILE LOOP
    np.savez(args.output_file_tag+'.npz', x_exit = dat[0], z_exit = dat[1],
                                          k_x = dat[2], k_y = dat[3], k_z = dat[4],
                                          cos_theta_exit = dat[5],
                                          exit_view_angle = dat[6],
                                          input_args = args,
                                          A0 = A0,
                                          A_geom = A)

  print 'status:', ev_count, iterations*args.num_events
  A =  A0 * float(ev_count) / float(iterations*args.num_events)
  print 'A', A
  # SAVE WHEN THE LOOP IS DONE
  np.savez(args.output_file_tag+'.npz', x_exit = dat[0], z_exit = dat[1],
                                        k_x = dat[2], k_y = dat[3], k_z = dat[4],
                                        cos_theta_exit = dat[5],
                                        exit_view_angle = dat[6],
                                        input_args = args,
                                        A0 = A0,
                                        A_geom = A)

