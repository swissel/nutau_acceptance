#!/usr/bin/env python

from pylab import *
rcParams['figure.facecolor']='white'
rcParams['font.size']=14
rcParams['legend.fontsize']=14
import matplotlib
#matplotlib.use('Agg')
import argparse
import os


if "TAU_ACCEPTANCE_DIR" not in os.environ:
    raise EnvironmentError('Environment variable TAU_ACCEPTANCE_DIR has not been set')

import sys
sys.path.append(os.environ['TAU_ACCEPTANCE_DIR']+'/src')
from tau_Shower_Sim_lib import *

if __name__ == "__main__":

  parser=argparse.ArgumentParser(description='Script the simulates tau neutrino induced tau lepton air shower detection from antennas at altitude')
  parser.add_argument("-G",    "--geom_file",     default = '/home/romerowo/nutau_sim/geom_outputs/geom_37_km_5_deg_10M_ev.npz', help='number of events to simulate', type=str)
  parser.add_argument("-L",    "--lut_file",      default = '/home/romerowo/nutau_sim/LUTs/2.0km_ice_lowCS_stdEL/LUT_1e+17_eV.npz',     help='number of events to simulate', type=str)
  parser.add_argument("-c",    "--cut_ang",       default = 5., help='number of events to simulate', type=float)
  parser.add_argument("-out_tag", "--output_tag", default = 'test', help='tag for output files', type=str)

  #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
  for i, arg in enumerate(sys.argv):
    if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

  args=parser.parse_args()

A_OMEGA_tau_exit(args.geom_file, args.lut_file, args.cut_ang, outTag = args.output_tag)


