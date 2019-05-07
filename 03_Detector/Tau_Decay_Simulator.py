#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 16:15:57 2019

@author: romerowo
"""

import numpy as np

class Tau_Decay_Simulator:
    def __init__(self):
        self.load_Tauola_Table()
        self.N = float(len(self.f_nu_tau))
        # The probabilities of muonic, electronic, and hadronic decay modes (PDG)
        self.shower_cut = np.logical_or(self.f_hadron>0., self.f_electron>0.)
        self.e_shower_cut = self.f_electron>0.
        self.h_shower_cut = self.f_hadron>0.
        self.P_shower     = float(np.sum(self.shower_cut))/self.N
        self.N_shower = np.sum(self.shower_cut)
        self.N_h_shower = np.sum(self.h_shower_cut)
        self.N_e_shower = np.sum(self.e_shower_cut)
        self.P_e_shower   = float(np.sum(self.e_shower_cut))/float(self.N_shower)
        self.P_h_shower   = float(np.sum(self.h_shower_cut))/float(self.N_shower)
        
        # sorted arrays for energy fraction sampling
        self.arr_h = np.sort(self.f_hadron[self.h_shower_cut])
        self.arr_e = np.sort(self.f_electron[self.e_shower_cut])
        
        self.tau_mass_eV = 1.77682e9 # eV/c^2

        
    def load_Tauola_Table(self):
        '''
        Tauola data from nutau_sim
        '''
        cc=0
        f_nu_tau   = []
        f_nu_mu    = []
        f_nu_e     = []
        f_hadron   = []
        f_muon     = []
        f_electron = []
        for line in file('tau_decay_tauola.data'):
            cc+=1
            if cc>1:
                f_nu_tau.append(float(line.split()[0]))
                f_nu_mu.append(float(line.split()[1]))
                f_nu_e.append(float(line.split()[2]))
                f_hadron.append(float(line.split()[3]))
                f_muon.append(float(line.split()[4]))
                f_electron.append(float(line.split()[5]))
                
        self.f_nu_tau   = np.array(f_nu_tau)
        self.f_nu_mu    = np.array(f_nu_mu)
        self.f_nu_e     = np.array(f_nu_e)
        self.f_hadron   = np.array(f_hadron)
        self.f_muon     = np.array(f_muon)
        self.f_electron = np.array(f_electron)
        
    def sample_shower_type(self):
        #0 is hadronic
        #1 is electron
        if np.random.uniform(0.,1.) > self.P_e_shower:
            return 0
        else:
            return 1

    def sample_energy_fraction(self, shower_type):
        #0 is hadronic
        #1 is electron
        if shower_type not in [0,1]:
            return 0.
        u = np.random.uniform(0., 1.)
        if shower_type==0:
            idx = u*float(len(self.arr_h)-1.) # so we don't hit the limit
            idx_lo = int(np.floor(idx))
            idx_hi = idx_lo+1
            val = self.lin_interp(idx, idx_lo, idx_hi, self.arr_h[idx_lo], self.arr_h[idx_hi])
            #return self.f_hadron[self.h_shower_cut][np.random.randint(0, self.N_h_shower)]
            return val
        if shower_type==1:
            idx = u*float(len(self.arr_e)-1.) # so we don't hit the limit
            idx_lo = int(np.floor(idx))
            idx_hi = idx_lo+1
            val = self.lin_interp(idx, idx_lo, idx_hi, self.arr_e[idx_lo], self.arr_e[idx_hi])
            #return self.f_electron[self.e_shower_cut][np.random.randint(0, self.N_e_shower)]
            return val
        
    def lin_interp(self, x, x1, x2, y1, y2):
        m = (y2-y1)/(x2-x1)
        b = y2 - m*x2
        return m*x + b 

    def sample_range(self, E_tau_eV, num_events):
        d_ref = 4.9 # km 
        E_ref = 1.e17 # eV
        d_val = E_tau_eV / E_ref * d_ref
        return np.random.exponential(d_val, num_events)
