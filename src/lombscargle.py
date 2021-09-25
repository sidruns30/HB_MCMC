''' Compute the Lomb-Scargle periodograms of the lightcurve files 
    Fold on the best fit frequencies, generate multiple folded lightcurves
'''
from typing_extensions import final
import numpy as np
import glob, os

from numpy.core.numeric import indices

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
        # required
        self.flux = None
        self.time = None
        # remove initial points / white data
        self.clean_flag = False
        # for period finding
        self.freq = None
        self.power = None
        self.period = None
        self.folded_flux = None
        self.phase = None
        if ('time' in kwargs): self.time = kwargs['time']
        if ('flux' in kwargs): self.flux = kwargs['flux']
        if ('period' in kwargs): self.period = kwargs['period']
    # load the lightcurve from the provided lc dir
    def load_lc(self, lc_name, **kwargs):
        self.id = lc_name[len(lc_dir):-4]
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
    # clean the lightcurve flux
    def clean_lc(self):
        initial_points = len(self.time)
        if (self.flux is None) or (self.time is None): raise ValueError ("Lightcuve not loaded properly")
        # Recenter the time array
        self.time -= self.time[0]
        # Throw the first 6 hours of data
        tmp = self.time
        self.time = self.time[tmp>0.25]
        self.flux = self.flux[tmp>0.25]
        dt = 24*60*(self.time[1] - self.time[0])
        # Search for gaps in the lightcurve. If present (> 1 day), throw away the next 6 hours
        gap_times = []
        for i, t in enumerate(self.time[:-1]):
            if self.time[i+1] - self.time[i] >= 1: gap_times.append(t)
        for gap in gap_times:
            tmp = self.time
            self.time = self.time[(tmp < gap ) | (tmp - gap > 0.25)]
            self.flux = self.flux[(tmp < gap ) | (tmp - gap > 0.25)]
        final_points = len(self.time)
        self.clean_flag = True
    # construct the periodogram for the lightcurve
    def lombscargle(self, **kwargs):
        if (not self.clean_flag): self.clean_lc()
        if ((self.flux is not None) and (self.time is not None)):
            # oversampling factor for frequency grid
            if ('oversample' in kwargs): oversample = kwargs['oversample']
            else: oversample = 10
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
    def get_harmonics(self, arr, **kwargs):
        # sort the array by frequency
        arr = np.sort(arr)
        #print("Test frequencies are: ", arr)
        tot_counts = [] # total number of harmonics present in an array
        if ('tol' in kwargs): tol = kwargs['tol']
        else: tol = 0.07
        for i, a in enumerate(arr[:-1]):
            # get the fractional remainder from the lowest frequency
            _ = (arr[i+1:] % a) / a
            ctr = 0 # count number of harmonics in the array
            for f in _:
                if (f < tol): ctr += 1
            tot_counts.append(ctr)
        n_harmonics = max(tot_counts)
        if n_harmonics < 1: 
            #print("Found no harmonics")
            return -1.
        else: 
            # return the frequency with the most harmonics present
            f_best = arr[np.argmax(tot_counts)]
            #print("Found %i harmonics with best fit frequency as %f" %(n_harmonics, f_best))
            return f_best
    # estimate the period of the signal if none is provided
    ''' if 'sufficiently' strong harmonics are present then the best fit frequency is
    the lowest of the harmoics otherwise the period with max power is selected'''
    def guess_period(self, **kwargs):
        if self.power is None: self.lombscargle()
        if self.period is not None: return
        # narrow down on the top few frequencies from LombScargle
        # some of these values are artificially set to make sense
        peak_distance = 5 * int(len(self.freq) / len(self.time))    # dist between peaks
        peak_height = max(self.power) / 10                           # min peak height 
        prominence = max(self.power) / 5                            # min peak prominence from its neighbors
        indices = find_peaks(self.power, distance = peak_distance, 
                            height = peak_height, prominence=prominence)[0]
        # sort peak indices by power
        indices = indices[np.argsort(self.power[indices])]
        max_peaks = 10                                              # maximum frequencies to look at
        indices = indices[-1:-max_peaks:-1]                         # indices sorted from highest to lowest in power
        test_frequencies = self.freq[indices]
        # check if mutiple harmonics are returned; otherwise use the best frequency
        self.period = 1 / self.get_harmonics(np.array(test_frequencies), **kwargs)
        if self.period == -1: self.period = 1 / test_frequencies[-1]
        print("Guessed period of the lightcurve after LombScargle: %f" %self.period)
        # Now perform a second iteration of lombscargle, based on the guessed frequency
        if ('refine' in kwargs): refine = kwargs['refine']
        # refine 1: use lombscargle on a more resolved frequency grid
        if refine == 1:
            freq_grid = np.linspace(0.98 / self.period, 1.02 / self.period, 2000)
            power = LombScargle(self.time, self.flux).power(freq_grid)
            self.period = 1 / freq_grid[np.argmax(power)]
            print("Period after refining using LombScargle is: %f" %self.period)
        # refine 2: minimize the sum(|dy|) of the phase folded lightcurve
        if refine == 2:
            freq_grid = np.linspace(0.98 / self.period, 1.02 / self.period, 2000)
            sum_abs_dflux = []
            # phase fold over each frequency
            for freq in freq_grid:
                phase = (self.time % (1 / freq)) * freq
                _ = np.argsort(phase)
                phase = phase[_]
                flux = self.flux[_]
                # Take the sum of absolute differences between neighboring elements
                # should be minimum for the correctly folded lightcurve? 
                sum_abs_dflux.append(np.sum(np.abs(flux[1:] - flux[:-1])))
            self.period = 1 / freq_grid[np.argmin(sum_abs_dflux)]
            print("Period after refining using dispersion min is: %f" %self.period)
    # phase fold on the best guessed frequency
    def phase_fold(self, **kwargs):
        if self.power is None: self.lombscargle()
        if self.period is None: self.guess_period(**kwargs)
        self.phase = (self.time % self.period) / self.period
        _ = np.argsort(self.phase)
        self.phase = self.phase[_]
        self.folded_flux = self.flux[_]
    # binning on the phase folded lightcurve
    def bin_lc(self, **kwargs):
        if ('N' in kwargs): N = kwargs['N']
        else: N = len(self.flux) / 100
        if self.folded_flux is None: self.phase_fold(self)
        # throw away points that are more than 4 std from binned flux
        bins = np.split(self.phase, 100)
    # Plot the orignal and phase folded lightcurve
    def plot_lc(self):
        if self.phase is not None:
            fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7,10))
            axes[0].plot(self.time, self.flux, "x")
            axes[0].set_title("Time and flux")
            axes[1].plot(self.phase, self.folded_flux, '.', markersize=1)
            axes[1].set_title("Folded flux folded at period %f" %self.period)
            axes[2].plot(self.freq, self.power)
            axes[2].set_title("LombScargle periodogram for the lightcurve")
        else:
            fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(7,7))
            axes[0].plot(self.time, self.flux, "x")
            axes[0].set_title("Time and flux")
            axes[1].plot(self.freq, self.power)
            axes[1].set_title("LombScargle periodogram for the lightcurve")
        plt.savefig("../figures/%s_fig.png" % self.id)

for lk_name in lc_list:
    lc = lightcurve()
    print("Loading lightcurve: ", lk_name)
    try:    lc.load_lc(lc_name = lk_name)#lc_dir + '/110602878.txt')#lc_list[0])
    except: pass
    lc.lombscargle()
    lc.guess_period(refine=2)
    lc.phase_fold()
    lc.plot_lc()

'''
class helper_functions:
    # I don't calculate with 
    def get_chi_square(arr):
        arr = 
'''