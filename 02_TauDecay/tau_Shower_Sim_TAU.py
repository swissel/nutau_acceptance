#!/usr/bin/env python

from pylab import *
rcParams['figure.facecolor']='white'
rcParams['font.size']=14
rcParams['legend.fontsize']=14
import matplotlib
#matplotlib.use('Agg')
import argparse

from tau_Shower_Sim_lib import *
import os

if __name__ == "__main__":

  parser=argparse.ArgumentParser(description='Script the simulates tau neutrino induced tau lepton air shower detection from antennas at altitude')
  parser.add_argument("-G",    "--geom_file",     default = os.environ['TAU_ACC_GEOM_DIR'] + '/geom_37_km_5_deg_10M_ev.npz', help='number of events to simulate', type=str)
  parser.add_argument("-L",    "--lut_file",      default = os.environ['TAU_ACC_LUT_DIR'] + '/4.0km_ice_midCS_stdEL/LUT_3e+17_eV.npz',     help='number of events to simulate', type=str)
  parser.add_argument("-c",    "--cut_ang",       default = 5., help='number of events to simulate', type=float)
  parser.add_argument("-out_tag", "--output_tag", default = 'test', help='tag for output files', type=str)

  #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
  for i, arg in enumerate(sys.argv):
    if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

  args=parser.parse_args()

A_OMEGA_tau_exit(args.geom_file, args.lut_file, args.cut_ang, outTag = args.output_tag)


