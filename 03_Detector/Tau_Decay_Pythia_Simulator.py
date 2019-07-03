#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 16:15:57 2019

@author: romerowo
"""

import numpy as np
class DecayParticle:
    def __init__(self, name='', tag='', pythia_code=0, Aires_code=0,shower_type=None, frac_en=0.):
	self.name = name # a name you can use for plots and stuff
	self.tag  = tag # a more readable name
	self.pythia_code = pythia_code
	self.Aires_code = Aires_code
	self.shower_type = shower_type # does this particle produce a shower?
	self.frac_en = frac_en # fraction of the original particle's energy
	if( pythia_code != 0):
		self.codes(pythia_code)
		
    def codes(self,particle_type):
   	self.pythia_code = particle_type

	code_array = np.array([16, 13, -14, 111, -211, 11, -12, 211, 22, -311, -321, 321, 223, 311, 130, 310, 221, -323])
	code_array_Aires = [8,-3,-7,10,-11,-2,-6,11,1,-14,-14,14,15,14,13,12,15,12]
	name_array = [r'$\nu_{\tau}$',r'$\mu^{-}$',r'$\bar{\nu_{\mu}}$',r'$\pi^{0}$',r'$\pi^{-}$',r'$e^{-}$',r'$\bar{\nu_{e}}$',r'$\pi^{+}$',r'$\gamma$',r'$\bar{K^{0}}$',r'$K^{-}$',r'$K^{+}$',r'$\omega$',r'$K^{0}$',r'$K^{0}_{L}$',r'$K^{0}_S$',r'$\eta$',r'$\bar{K^{+}}$']
        tag_array = ['nu_tau', 'mu', 'anti_nu_mu', 'pi0', 'pi-', 'e-', 'anti_nu_e','pi+', 'gamma', 'anti-K0', 'K-', 'K+', 'omega', 'K0', 'K0_L', 'K0_S', 'eta', 'K+']
	shower_array = [None, None, None, 'hadron', 'hadron', 'em', None, 'hadron', 'em', 'hadron', 'hadron', 'hadron', 'hadron', 'hadron', 'hadron', 'hadron', 'hadron', 'hadron']

	idx = np.where( code_array == self.pythia_code)
	if( idx != None):
		idx = idx[0][0]
		self.Aires_code = code_array_Aires[idx]
		self.name = name_array[idx]
		self.tag = tag_array[idx]
		self.shower_type = shower_array[idx]
	else:
		self.Aires_code = 0
		self.name = ''
		self.tag = ''
		self.shower_type = None

class Tau_Decay_Simulator:
    def __init__(self, pol='negative'):
        self.load_Pythia_Table(pol=pol)   
        self.tau_mass_eV = 1.77682e9 # eV/c^2
     
    def load_Pythia_Table(self, pol='negative'):
        '''
        Pythia data run by Austin Cummings
        '''
        cc=0
        type = 0

	self.tau_decays = []
	part_types = []
	if( pol == 'negative'):
		f=open('negativepolnew.txt','r')
		print("Reading negativepolnew.txt")
	elif( pol == 'positive'):
		f=open('positivepolnew.txt','r')
		print("Reading positivepolnew.txt")
	else:
		print("Please specify either negative (default) or positive polarity for the Pythia runs")

	for line in f:
		all_parts = line.split(';')[:-1]
		single_tau_decay = []
		for i in range(0,len(all_parts)):
			b = all_parts[i].split(',')
			part_type = int(b[0])
			if part_type not in part_types:
				part_types.append(part_type)
			frac_en = float(b[1])

			# save the decay particle
			decay_particle = DecayParticle(pythia_code=part_type, frac_en=frac_en)
			single_tau_decay.append(decay_particle)
		self.tau_decays.append(single_tau_decay)
	print(("All possible daughters of the tau decay (within 10k decays): " + str(part_types)))	
	self.N = len(self.tau_decays)
	self.build_energy_distributions()

    def build_energy_distributions(self, nbins=10000): # use a a large number of bins to avoid numerical issues with the sampling
    	self.shower_energy_distribution(nbins=nbins, type='shower')
	self.shower_energy_distribution(nbins=nbins, type='hadron')
	self.shower_energy_distribution(nbins=nbins, type='em')

    def shower_energy_distribution(self, nbins=10000, type='shower'):
	#self.shower_energies = []
	energies = []
	for i, single_decay in enumerate(self.tau_decays):
		energy_fraction = 0
		choose_decay = False
		for j,particle in enumerate(single_decay):
			if( type == 'shower'):
				if(particle.shower_type != None):
					#print('got one', particle.shower_type, particle.tag, particle.frac_en)
					energy_fraction += particle.frac_en
					choose_decay = True
			else:
				if(particle.shower_type == type):
					#print('got one', particle.shower_type, particle.tag, particle.frac_en)
					energy_fraction += particle.frac_en
					choose_decay = True

		if choose_decay:
			energies.append(energy_fraction)

	if( type == 'shower'):
		self.shower_energies = energies
		self.shower_frac,self.shower_energybins = np.histogram(energies, bins=nbins, range=(0,1),weights = np.ones(len(energies))/float(len(energies)))
	elif(type == 'hadron'):
		self.hadron_energies = energies
		self.hadron_frac,self.hadron_energybins = np.histogram(energies, bins=nbins, range=(0,1),weights = np.ones(len(energies))/float(len(energies)))
	elif(type == 'em'):
		self.em_energies = energies
		self.em_frac,self.em_energybins = np.histogram(energies, bins=nbins, range=(0,1),weights = np.ones(len(energies))/float(len(energies)))


    def sample_energy_fraction(self, num_events=1, type='shower'):
    	# generate random samples from a weighted array using numpy.random.choice
	if(type == 'shower'):
		bins = (self.shower_energybins[:-1] + self.shower_energybins[1:])/2.
		weights = self.shower_frac
	elif(type == 'hadron'):
		bins = (self.hadron_energybins[:-1] + self.hadron_energybins[1:])/2.
		weights = self.hadron_frac
	elif(type == 'em'):
		bins = (self.em_energybins[:-1] + self.em_energybins[1:])/2.
		weights = self.em_frac

	samps = np.random.choice(bins,num_events, p=weights)
	return samps

    def sample_range(self, E_tau_eV, num_events):
        d_ref = 4.9 # km 
        E_ref = 1.e17 # eV
        d_val = E_tau_eV / E_ref * d_ref
        return np.random.exponential(d_val, num_events)

    def energy_fraction(self,type='shower', rebin_factor=200):
    	n_original_bins = len(self.shower_frac)
	n_new_bins = int(np.floor(n_original_bins/rebin_factor))
	
	if( type == 'shower'):
		original_bins = (self.shower_energybins[:-1] + self.shower_energybins[1:])/2.
		original_weights = self.shower_frac
	elif( type =='hadron'):
		original_bins = (self.hadron_energybins[:-1] + self.hadron_energybins[1:])/2.
		original_weights = self.hadron_frac
	elif( type =='em'):
		original_bins = (self.em_energybins[:-1] + self.em_energybins[1:])/2.
		original_weights = self.em_frac

	if( n_original_bins % rebin_factor != 0):
		print("I can only rebin the array of len(%d) by a number that is a integer divisor of it."%original_bins)
		return original_bins, original_weights

	new_bins = original_bins.reshape( n_original_bins/rebin_factor, rebin_factor).mean(1)
	new_weights = original_weights.reshape( n_original_bins/rebin_factor, rebin_factor).sum(1)

	return new_bins, new_weights

