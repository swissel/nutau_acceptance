# -*- coding: utf-8 -*-

"""
This module contains a set of convenience functions to parse the output 
of NEC2, modify the input (.nec) file, and rerun NEC as necessary.
NEC2 is the Numerical Electromagnetics Code, developed at LLNL.
The version of NEC2 this code currently assumes is `here <http://www.physics.otago.ac.nz/research/electronics/nec/index.html>`_.

Several NEC files are included with the LSL distribution for modeling the 
dipoles.  See the `README.NEC` file included in the LSL data directory for
more information about what is included.
"""

from numpy import *
from mathutil import regrid
import os
import re
import logging


__version__   = '0.1'
__revision__ = '$Rev: 1498 $'
__all__ = ['CloseTo', 'open_and_get_nec_freq', 'change_nec_freq', 'calcIME', 'NECImpedance', 'NECPattern', 'whichNEC4', '__version__', '__revision__', '__all__']
__author__    = 'P. S. Ray'
__maintainer__ = 'Jayce Dowell'

_NEC_UTIL_LOG = logging.getLogger('nec_util')
_NEC_UTIL_LOG.setLevel(logging.INFO)


def CloseTo( x, y, epsilon=0.005 ):
	"""
	Return True if two numbers are within a specified fractional, not 
	absolute, tolerance of each other.  Tolerance epsilon is a keyword 
	parameter with a default of 0.005.
	"""

	return ( 2.0*abs(x-y)/(x+y) < epsilon )


def open_and_get_nec_freq(fname):
	"""
	Open a NEC output file and return a tuple containing the open file 
	object and the first frequency found in the file (MHz).
	"""

	f = open(fname)

	# Start at beginning of file
	f.seek(0)

	# skip all lines until a line containing "STRUCTURE SPECIFICATION": this
	# effectively skips all comments (from CM cards) and text resulting from
	# reading Numerical Green's Function parts. Of course, if a user writes
	# "STRUCTURE SPECIFICATION" in his comment lines, this still fails...
	for line in f:
		if line.find('STRUCTURE SPECIFICATION') >= 0:
			break
	else:
		raise RuntimeError("STRUCTURE SPECIFICATION not found!")
		
	#  Now look for FREQUENCY and get the value
	for line in f:
		if line.find('FREQUENCY') >= 0:
			break
	for line in f:
		#print line
		if line.find('FREQUENCY=') >= 0:
			freq = float(line[line.find('=')+1:].split()[0])
			break
		if line.find('FREQUENCY :') >=0:
			freq = float(line[line.find(':')+1:].split()[0])
			break
	else:
		raise RuntimeError("Frequency value not found")
	
	_NEC_UTIL_LOG.debug("Found frequency %f MHz", freq)
	return (f, freq)


def change_nec_freq(necname, freq):
	"""
	Modify the FR card in a NEC input file to run at freq.
	"""

	f = open(necname,'r+')
	buf = f.read()
	lines = buf.splitlines(True)
	# Substitute the freq in the right field of the FR card
	for i in range(len(lines)):
		if lines[i][:2] == 'FR':
			_NEC_UTIL_LOG.debug("Found line : %s", lines[i])
			vals = re.split(',| +',lines[i])
			_NEC_UTIL_LOG.debug("Vals = %s", vals)
			vals[5] = "%.2f" % freq
			lines[i] = " ".join(vals)
			# Make sure this line ends in newline
			if lines[i][-1] != '\n':
				lines[i] += '\n'
	# Rewrite the file
	f.seek(0)
	f.writelines(lines)
	f.truncate()
	f.close()


def calcIME(necname, myfreqs = None, zpre = 100):
	"""
	Compute the impedance mismatch efficiency (IME), for a given NEC run 
	and write out the results in a file.  Assumes a default preamplifier input 
	impedance of 100 ohms, unless overridden by the zpre keyword argument.
	Returns the frequencies calculated by NEC, unless myfreqs is set, in 
	which case it interpolates onto that frequency grid.
	"""

	ant = NECImpedance(necname)
	gamma = (zpre - ant.z)/(zpre + ant.z)
	ime = (1 - abs(gamma)**2)
	if myfreqs is None:
		return (ant.freqs, ime)
	else:
		newime = regrid(ant.freqs, ime, myfreqs, method = 'linear')
		return (myfreqs, newime)


class NECImpedance:
	"""
	NECImpedance:
	Python class to read an array of impedance values from a NEC2 .out file The 
	.nec file should loop over a range of frequencies with an FR card like this:

		FR 0 91 0 0 10.0 1.0

	The RP card should be minimal to keep the runtime and output file size from 
	growing huge.  For example:

		RP 0,91,1,1000,0.,0.,1.0,1.0

	"""

	def __init__(self, necname):
		outname = os.path.splitext(necname)[0] + '.out'
		f = open(outname)

		# Start at beginning of file
		f.seek(0)
		
		# skip all lines until a line containing "STRUCTURE SPECIFICATION": this
		# effectively skips all comments (from CM cards) and text resulting from
		# reading Numerical Green's Function parts. Of course, if a user writes
		# "STRUCTURE SPECIFICATION" in his comment lines, this still fails...
		for line in f:
			if line.find('STRUCTURE SPECIFICATION') >= 0: break
		else:
			raise RuntimeError("STRUCTURE SPECIFICATION not found!")

		freqs = []
		impedances = []
		while (True):
			#  Now look for FREQUENCY and get the value
			for line in f:
				if line.find('FREQUENCY') >= 0:
					break
			for line in f:
				_NEC_UTIL_LOG.debug(line.strip())
				if line.find('FREQUENCY=') >= 0:
					freq = float(line[line.find('=')+1:].split()[0])
					break
				if line.find('FREQUENCY :') >=0:
					freq = float(line[line.find(':')+1:].split()[0])
					break
			else:
				_NEC_UTIL_LOG.debug("No more freqs...")
				break
			_NEC_UTIL_LOG.debug("Found frequency %f MHz",freq)
			for line in f:
				if line.find('ANTENNA INPUT PARAMETERS') >= 0:
					break
			gotimp = False
			for line in f:
				if line.find('IMPEDANCE') >= 0:
					gotimp = True
					break
			if not gotimp:
				raise RuntimeError("IMPEDANCE not found")
			for line in f:
				break
			for line in f:
				_NEC_UTIL_LOG.debug(line.strip())
				break

			# Here we need to add a space before - signs that
			# are not preceeded by an E, so it will parse
			line = re.sub(r'(\d)-',r'\1 -',line)
			re_z = float(line.split()[6])
			im_z = float(line.split()[7])
			freqs.append(freq)
			impedances.append(complex(re_z,im_z))
			
		self.freqs = array(freqs)
		self.z = array(impedances)


class NECPattern:
	"""
	NECPattern:
	Python class to read the pattern from a NEC2 .out file.  Note that the 
	.nec file should have an RP card to run over the full pattern, like this:

		RP 0,91,360,1000,0.,0.,1.0,1.0,0.

	The FR card should be a simple single frequency run:

		FR 0,1,0,0,74.0,1

	"""
        
	def __init__(self, necname, freq, rerun = True):
		# Modify NEC file to set FR card to use "freq"
		
		# Run NEC if necessary.
		# anntenna_pat_dB[az,alt] is the total gain or current in dB in the 
		# direction az, alt (integer degrees), where az (azimuth) runs from 
		# 0 to 359, where 0 is North and alt (altitude) runs from 0 to 89 , 
		# where 0 is the horizon The default pattern is all zeros (isotropic 
		# response)
		self.antenna_pat_dB=zeros(shape=(360,90),dtype=float_)

		outname = os.path.splitext(necname)[0] + '.out'
		try:
			f, filefreq = open_and_get_nec_freq(outname)
		except:
			print("NEC .out file not found! Running NEC")
			f = None
		
		if f is None or not CloseTo(filefreq,freq):
			if rerun:
				_NEC_UTIL_LOG.warning("NEC output file is at a different frequency \
					than the requested frequency: re-running")
				if f is not None:
					f.close()
				change_nec_freq(necname,freq)

				# Make sure we have NEC install
				if whichNEC4() is None:
					raise RuntimeError("NEC executable 'nec4d' not found in PATH")

				# Important NOTE:
				# This requires a modified version of NEC-4 that
				# takes 2 command line arguments instead of asking questions
				# interactively. See Paul Ray for info.
				cmdstr = "nec4d %s %s" % (necname, outname)
				ret=os.system(cmdstr)
				if ret != 0:
					raise RuntimeError("Bad return value from nec2++ call : %d" % ret)       
				f, filefreq = open_and_get_nec_freq(outname)
				if not CloseTo(filefreq,freq):
					raise ValueError("NEC failed to generate a file with the correct frequency.")
					
			else:
				raise ValueError("NEC output file is at a different frequency (%f) than the requested frequency (%f)." % \
					(filefreq, freq))
		
		#  Now look for RADIATION PATTERN or EXCITATION and read it
		radpat = None
		for line in f:
			if line.find('RADIATION PATTERN') >= 0:
				radpat = True
				break
			if line.find('EXCITATION') >= 0:
				radpat = False
				break
		else:
			raise RuntimeError("RADIATION PATTERN nor EXCITATION not found!")

		
		if radpat:
			self.__readRADIATION(f)
		else:
			self.__readEXCITATION(f)
		

	def __readRADIATION(self, f):
		"""
		Private function to read in a RADIATION PATTERN section of a NEC
		output file.
		"""

		# Some versions of NEC2 output extraneous data after "RADIATION PATTERNS" before the actual data
		# and column labels (e.g. RANGE = and EXP (-JKR) values ).  Discard until
		# the true bottom of the column labels (look for DB) */
		for line in f:
			if line.find('DB') >= 0:
				break

		n = 0
		for line in f:
			cols = line.split()
			if len(cols) < 4:
				break
			# Parse theta and phi into ints so we can use for indexing
			# Convert theta from zenith angle to altitude
			theta = 90-int(cols[0].split('.')[0])
			phi = int(cols[1].split('.')[0])
			if theta < 0 or theta > 89 or phi > 359:
				#print "Skipping ",phi,theta
				continue
			powgain= float(cols[4])
			#print phi, theta, powgain
			self.antenna_pat_dB[phi,theta] = powgain
			n += 1
			_NEC_UTIL_LOG.debug("theta %d phi %d gain %f", theta,phi,powgain)


	def __readEXCITATION(self, f):
		"""
		Private function to read in data stored in a collection of EXCITATION 
		sections in a NEC output file.
		"""

		n = 0
		# We have already read the first line of the first entry, so start there
		# The information we need is stored across the 15 lines following the 
		# EXCITATAION heading.  Read those lines into parts and then deal with the
		# results.  The keys lines are #2 (theta and phi) and #12 (induced currents)
		lineCount = 0
		for line in f:
			if lineCount % 16 == 0:
				parts = []
			parts.append( line )
			lineCount += 1
			if lineCount % 16 == 0:
				fieldsAngle = parts[1].split()
				fieldsCurrent = parts[11].split()

				# Direction of the incident radiation
				theta = 90 - int(float(fieldsAngle[3]))
				phi = int(float(fieldsAngle[6]))
				if theta < 0 or theta > 89 or phi > 359:
					pass
				else:
					# Get the absolute value and put it on a dB scale
					powcurr = float(fieldsCurrent[8])
					powcurr = 10.0*log10(powcurr)
					#print phi, theta, powcurr
					self.antenna_pat_dB[phi,theta] = powcurr
					n += 1
					_NEC_UTIL_LOG.debug("theta %d phi %d current %f", theta,phi,powcurr)


def whichNEC4():
	"""
	Return the path to the nec4d executable if it can be found in the 
	current path.  None otherwise.  This is useful for making sure that NEC
	is installed before trying to run something.
	"""

	def is_exe(fpath):
		return os.path.exists(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split('nec4d')
	if fpath:
		if is_exe('nec4d'):
			return 'nec4d'
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, 'nec4d')
			if is_exe(exe_file):
				return exe_file

	return None
