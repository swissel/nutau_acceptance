from numpy import *
from scipy.signal import butter, lfilter, bessel
from scipy.signal import blackman,bartlett, hamming, hanning, kaiser
from scipy import signal
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend.
mpl.use('Agg')
import matplotlib.pyplot as pyp
#######################################################################
# filters

def butter_bandstop(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    #print "Butter Filter: ", lowcut, highcut, fs, low, high, order 
    b, a = butter(order, [low, high], btype='bandstop')
    return b, a

def butter_bandstop_filter(data, lowcut, highcut, fs, order=3):
    #applies butterworth filter on timeseies and returns timeseries
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    #pyp.figure(12)
    #w, h = signal.freqz(b, a, worN=4096)
    #pyp.semilogy((fs * 0.5 / pi) * w, abs(h), label="order = %d" % order)
    y = lfilter(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    #print "Butter Filter: ", lowcut, highcut, fs, low, high, order 
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    #applies butterworth filter on timeseies and returns timeseries
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    #pyp.figure(12)
    #w, h = signal.freqz(b, a, worN=4096)
    #pyp.semilogy((fs * 0.5 / pi) * w, abs(h), label="order = %d" % order)
    y = lfilter(b, a, data)
    return y

def bessel_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 *fs
    low = lowcut / fs
    high = highcut /fs
    b, a = bessel(order, [low,high], btype='bandpass')
    return b,a

def bessel_bandpass_filter(data, lowcut, highcut, fs, order=3):
    #applies butterworth filter on timeseies and returns timeseries
    b, a = bessel_bandpass(lowcut, highcut, fs, order=order)
    #pyp.figure(12)
    #w, h = signal.freqz(b, a, worN=4096)
    #pyp.semilogy((fs * 0.5 / pi) * w, abs(h), label="order = %d" % order)
    y = lfilter(b, a, data)
    return y


def slac_butterworth_filter(x,fstart, fstop, fs):
	return butter_bandpass_filter(x, fstart, fstop, fs)
	#return bessel_bandpass_filter(x, fstart, fstop, fs)

'''def slac_butterworth_filter(x, fstart, fstop, fs):
    if( fstart == 200e6 and fstop == 1200e6 and fs == 5e9):
    	NPOLES = 10
    	GAIN = 4.55795945e+01
    	b = [-0.0112635125, 0.131958494, -0.696966183, 2.38373375, -5.82318248, 10.4051525, -13.903009, 14.0406871, -10.2444333, 4.71583406, 0]
    	a = [-1.0, 0.0, 5.0, 0.0, -10.0, 0.0, 10.0, 0.0, -5.0, 0.0, 1.0]
    	y = butterworth_filter(x,GAIN,a,b)
    	return y
    elif( fstart == 200e6 and fstop == 1200e6 and fs == 10e9):
	NPOLES = 10
	GAIN = 7.796705300e+02
	b = [-1.25430622e-01,1.43496775e+00,-7.51048919e+00,2.36821675e+01,-4.98200075e+01,7.30620523e+01,-7.56411835e+01,5.45748143e+01,-2.62446427e+01,7.58774853e+00,0]
	a = [-1.0, 0.0, 5.0, 0.0, -10.0, 0.0, 10.0, 0.0, -5.0, 0.0, 1.0]
    	y = butterworth_filter(x,GAIN,a,b)
    	return y
    elif( fstart == 0 and fstop == 1200e6 and fs == 10e9):
	NPOLES = 5
	GAIN = 3.614117003e+02
	b = [8.08708884e-02,	-6.00666748e-01, 1.85070748e+00,-2.99386863e+00, 2.57441532e+00, 0]
	a = [1,5,10,10,5,1]
	y = butterworth_filter(x,GAIN,a,b)
    	return y
    elif( fstart == 0 and fstop == 1200e6 and fs == 5e9 ):
	NPOLES = 5
	GAIN = 2.223750758e+01
	b = [1.76990902e-03,-5.77949601e-02,6.71216999e-02,-6.46993998e-01,1.96887159e-01]
	a =  [1,5,10,10,5,1]
	y = butterworth_filter(x,GAIN,a,b)
    	return y
    else:
	print "Error: Filter not defined for band ", fstart/1e6, " - ", fstop/1e6, "MHz amd sampling rate", fs/1e9, "GSPS"
	return x
'''
def butterworth_filter(x,GAIN,a,b,mode='td'):
    NSAMPLES = len(x)
    NCOEFF = len(a)-1
    xv=zeros(NCOEFF+1)
    yv=zeros(NCOEFF+1)
    
    y = zeros(len(x))
    for n in range(0,NSAMPLES):
        for i in range(NCOEFF):
            xv[i] = xv[i+1]
        xv[NCOEFF]=x[n] / GAIN
        for i in range(NCOEFF):
            yv[i] = yv[i+1]
        yv[NCOEFF] = sum(a*xv)
        yv[NCOEFF] += sum(b*yv)
        y[n] = yv[NCOEFF]
    if( mode == 'fd'):
        return fft.rfft(y)
    else:
        return y

def realRectFilter(fr, fft, fstart, fstop):
    ind1 = where( fr >= fstart)[0][0]
    ind2 = where( fr >= fstop )[0][0]
    
    fft[:ind1] = 1e-9
    fft[ind2:] = 1e-9 
    return fft 


def rectWindow(x, y, x1, x2):
    # x1 = start value 
    # x2 = stop value
    ind1 = where( x >= x1)[0][0]
    ind2 = where( x[ind1:] <  x2)[0][-1] + ind1
    
    ynew = zeros( len(y) )
    ynew[ind1:ind2] = y[ind1:ind2]
    return x, ynew

def rectWindow2(x, y, x1, x2):
    # x1 = negative distance from max(y)
    # x2 = positive distance from max(y)
    x = array(x, dtype=float64)
    y = array(y, dtype=complex64)
     
    if( abs( min(y) ) > abs( max(y) ) ):
        indm = where( y == min(y) )[0][0]
    else:
        indm = where( y == max(y) )[0][0]


    #indm = where( y == max(y))[0][0]
    
    dx = x[1]-x[0]
    ind1 = int(indm - x1/dx)
    ind2 = int(indm + x2/dx)
    
    ynew = zeros( len(y), dtype=complex64 )
    #ynew[ind1:ind2] = y[ind1:ind2].real + 1j*y[ind1:ind2].imag
    ynew[ind1:ind2] = y[ind1:ind2]
    return x, ynew

# useful for windowing

def gaussian_filter(length, t0, width, dt=1.):
    t = arange(0, length)*dt
    window = exp(-0.5*((t-t0)/width)**2)
    return window
