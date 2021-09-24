''' Compute the Lomb-Scargle periodograms of the lightcurve files 
    Fold on the best fit frequencies, generate multiple folded lightcurves
'''
import numpy as np
import glob, os

# Get the directory name with the lightcurve
lc_dir = os.getcwd() + '/../lightcurves'
lc_list = glob.glob(lc_dir + '/*txt')

# Construct a periodogram; find best fit freqs
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import math

def phase_fold(time, flux):
    f_nyq = 0.5 / (time[1] - time[0])
    ovr_sam = 3
    N = len(time)
    frequency = np.linspace(1/N, f_nyq, N*ovr_sam)
    power = LombScargle(time, flux).power(frequency=frequency)
    # Find the best fit frequencies
    ind = find_peaks(power, distance = ovr_sam,
                    height = np.max(power)/2,
                    prominence=max(power)/3)[0]
    ind = ind[np.argsort(power[ind])]
    plt.plot(frequency, power)
    plt.plot(frequency[ind], power[ind], 'rx')
    plt.savefig('figures/lombscargle.png')

    # Sorting by power does not give the best period usually; maybe also sort by period
    print(1/frequency[ind])
    print(ind, power[ind])
    phase = time * frequency[0] - math.floor(time * frequency[0])

class lightcurve:
    def __init__(self, **kwargs):
        self.flux = None
        self.time = None
        self.freq = None
        self.power = None
        self.period = None
        self.folded_flux = None
        self.phase = None
        if ('time' in kwargs): self.time = kwargs['time']
        if ('flux' in kwargs): self.time = kwargs['flux']
        if ('period' in kwargs): self.time = kwargs['period']
    # load the lightcurve from the provided lc dir
    def load_lc(self, lc_name, **kwargs):
        time, flux  = [], []
        if ('skiprows' in kwargs): skiprows = kwargs['skiprows']
        else: skiprows = 1
        with open(lc_name, "r") as f:
            for i, line in enumerate(f.readlines()):
                # Skip the first few lines
                if i < skiprows: pass
                else: 
                    data = line.split()
                    # Only load data points where third column is good (= 0)
                    if (float(data[3]) == 0):
                        time.append(float(data[0]))
                        flux.append(float(data[1]))
        self.flux = np.array(flux)
        self.time = np.array(time)
    # construct the periodogram for the lightcurve
    def lombscargle(self, **kwargs):
        if ((self.flux is not None) and (self.time is not None)):
            # oversampling factor for frequency grid
            if ('oversample' in kwargs): oversample = kwargs['oversample']
            else: oversample = 5
            # construct a frequency grid
            if ('freq' in kwargs): self.freq = kwargs['freq']
            else:
                # assuming a uniformly spaced data (TESS)
                f_low = 1 / (self.time[-1] - self.time[0])
                f_nyq = 0.5 /(self.time[1] - self.time[0])
                self.freq = np.linspace(f_low, f_nyq, oversample*len(self.time))
            # Now construct the periodogram
            self.power = LombScargle(self.time, self.flux).power(self.freq)
        else: raise ValueError('Flux and time not loaded in the lightcurve')
    # supplementary function to get the harmonics in an array with a certain tolerance
    def get_harmonics(arr, **kwargs):
        # sort the array by frequency
        arr = np.sort(arr)
        tot_counts = [] # total number of harmonics present in an array
        if ('tol' in kwargs): tol = kwargs['tol']
        else: tol = 0.02
        for i, a in enumerate(arr[:-1]):
            # get the fractional remainder from the lowest frequency
            _ = (arr[i+1:] % a) / a
            ctr = 0 # count number of harmonics in the array
            if ((f < tol) for f in _): ctr += 1
            tot_counts.append(ctr)
        n_harmonics = max(tot_counts)
        if n_harmonics == 0: 
            print("Found no harmonics")
            return -1.
        else: 
            # return the frequency with the most harmonics present
            f_best = arr[np.argmax(tot_counts)]
            print("Found %i harmonics with best fit frequency as %d" %n_harmonics, f_best)
            return f_best
    # estimate the period of the signal if none is provided
    ''' if 'sufficiently' strong harmonics are present then the best fit frequency is
    the lowest of the harmoics otherwise the period with max power is selected'''
    def guess_period(self, **kwargs):
        if not(self.power): self.lombscargle()
        # narrow down on the top few frequencies from LombScargle
        # some of these values are artificially set to make sense
        peak_distance = 5 * int(len(self.freq) / len(self.time))    # dist between peaks
        peak_height = max(self.power) / 2                           # min peak height 
        prominence = max(self.power) / 3                            # min peak prominence from its neighbors
        indices = find_peaks(self.power, distance = peak_distance, 
                            height = peak_height, prominence=prominence)[0]
        # sort peak indices by power
        indices = indices[np.argsort(self.power[indices])]
        max_peaks = 10                                              # maximum frequencies to look at
        indices = indices[-1:-max_peaks:-1]
        test_frequencies = self.freq[indices]
        # check if mutiple harmonics are returned; otherwise use the best frequency
        self.period = 1 / self.get_harmonics(test_frequencies, **kwargs)
        if self.period == -1: self.period = 1 / test_frequencies[-1]
        print("Guessed period of the lightcurve: %d" %self.period)
    # phase fold on the best guessed frequency
    def phase_fold(self, **kwargs):
        if not(self.power): self.lombscargle()
        if not(self.period): self.guess_period(**kwargs)
        self.phase = (self.time % self.period) / self.period
        _ = np.argsort(self.phase)
        self.phase = self.phase[_]
        self.folded_flux = self.flux[_]
    # binning on the phase folded lightcurve
    def bin_lc(self, **kwargs):
        if ('N' in kwargs): N = kwargs['N']
        else: N = len(self.flux) / 100
        if not(self.folded_flux): self.phase_fold(self)
        # throw away points that are more than 4 std from binned flux
        bins = np.split(self.phase, 100)
        
lc = lightcurve()
lc.load_lc(lc_name = lc_list[0])
lc.lombscargle()