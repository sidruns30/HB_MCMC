''' Compute the Lomb-Scargle periodograms of the lightcurve files 
    Fold on the best fit frequencies, generate multiple folded lightcurves
'''
from typing_extensions import final
import numpy as np
import glob, os, sys
from numpy.core.numeric import indices

# Get the directory name with the lightcurve
lc_dir = os.getcwd() + '/../lightcurves/original'
lc_list = glob.glob(lc_dir + '/*txt')

# Construct a periodogram; find best fit freqs
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib
import math
matplotlib.use('agg')

class functions:
    def phase_fold(period, time, flux):
        phase = (time % period) / period
        _ = np.argsort(phase)
        return phase[_], flux[_]
    # Average array with Npts per bin; returned shape is len(arr)//Npts
    def avg_bin(arr, Npts):
        return np.array([np.average(arr[i:i+Npts]) for i in range(0, len(arr) - len(arr)%Npts, Npts)])
    # Take std dev of array with Npts per bin
    def std_bin(arr, Npts):
        return np.array([np.nanstd(arr[i:i+Npts]) for i in range(0, len(arr) - len(arr)%Npts, Npts)])
    # Subtract a staright line fit from every n points; returned shape is len(arr) [last few points are not subtracted]
    def remove_line(xarr, yarr, Npts):
        new_arr = np.array([])
        for i in range(0, len(xarr) - len(xarr)%Npts, Npts):
            xpts = xarr[i:i+Npts]
            ypts = yarr[i:i+Npts]
            X = np.average(xpts)
            Y = np.average(ypts)
            m = np.average((xpts - X)*(ypts - Y))/np.average(np.square(xpts - X))
            b = Y - m*X
            new_arr = np.concatenate((new_arr, np.array(ypts - m*xpts  - b)), axis=-1)
        return new_arr.ravel()
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
        self.id = lc_name[len(lc_dir)+1:-4]
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
                f_nyq = 0.5 /(30 / 60 / 24)         # sampling time is 30 min
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
        if len(tot_counts) > 0: n_harmonics = max(tot_counts)
        else: n_harmonics = -1
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
        peak_height = max(self.power) / 20                          # min peak height 
        prominence = max(self.power) / 5                            # min peak prominence from its neighbors
        rel_height = 0.2                                            # relative prominence from which to calculate widths
        indices, prop = find_peaks(self.power, distance = peak_distance, 
                            height = peak_height, prominence=prominence, width=rel_height)
        # integrate by power (?)
        # sort peak indices by areas
        '''IMPORTANT: FREQUENCIES ARE NOT MONOTONIC IF SORTED BY AREA! (Change that to avoid errors in harmonic finding)'''
        #peak_areas = prop['prominences'] * prop['width_heights']
        #indices = indices[np.argsort(peak_areas)]
        #max_peaks = 10                                              # maximum frequencies to look at
        #indices = indices[-1:-max_peaks:-1]                         # indices sorted from highest to lowest in area
        test_frequencies = self.freq[indices]
        # setting reasonable limits on test frequenies
        f_low_limit = 0.1   # period of 10 days
        f_high_limit = 6    # period of 4 hours
        test_frequencies = test_frequencies[(test_frequencies > f_low_limit) & (test_frequencies < f_high_limit)]
        # expend limits if no frequencies are found
        if len(test_frequencies) < 1: test_frequencies = self.freq[indices]
        # check if mutiple harmonics are returned; otherwise use the best frequency
        self.period = 1 / self.get_harmonics(np.array(test_frequencies), **kwargs)
        if self.period == -1: self.period = 1 / test_frequencies[-1]
        print("Guessed period of the lightcurve after LombScargle: %f" %self.period)
        # Now perform a second iteration of lombscargle, based on the guessed frequency
        if ('refine' in kwargs): 
            refine = kwargs['refine']
            freq_grid = np.linspace(0.95 / self.period, 1.05 / self.period, 2000)
        # refine 1: use lombscargle on a more resolved frequency grid
        if refine == 1:
            power = LombScargle(self.time, self.flux).power(freq_grid)
            self.period = 1 / freq_grid[np.argmax(power)]
            print("Period after refining using LombScargle is: %f" %self.period)
        # refine 2: minimize the sum(|dy|) of the phase folded lightcurve
        if refine == 2:
            sum_abs_dflux = []
            # phase fold over each frequency
            for freq in freq_grid:
                phase, flux = functions.phase_fold(1/freq, self.time, self.flux)
                # Take the sum of absolute differences between neighboring elements
                sum_abs_dflux.append(np.sum(np.abs(flux[1:] - flux[:-1])))
            self.period = 1 / freq_grid[np.argmin(sum_abs_dflux)]
            print("Period after refining using L1 min is: %f" %self.period)
        # refine 3: minimize standard deviation after selecting a few periods
        if refine == 3:
            sum_std = []
            # phase fold over each frequency
            for freq in freq_grid:
                phase, flux = functions.phase_fold(1/freq, self.time, self.flux)
                # now bin this flux with 10 pts per bin and find std
                sum_std.append(np.sum(functions.std_bin(flux - functions.remove_line(phase, flux, 10))))
            self.period = 1 / freq_grid[np.argmin(sum_std)]
            print("Period after minimizing std in binned flux is: %f" %self.period)
    # refine an input guessed period (basically copies the methods from top but uses better schemes)
    def refine_guess(self, refine, **kwargs):
        print("Enter guessed period for lightcurve: %s" %self.id)
        period_guess = float(input())
        # first construct a frequency grid over which to sample periods
        fguess = 1/ period_guess
        fmin = fguess * 0.67    # Tmax = 1.5 guessed period
        fmax = fguess * 1.5     # Tmin = 0.67 the guessed period
        freq_grid = np.linspace(fmin, fmax, 10000)
        # refine 1: use lombscargle on a more resolved frequency grid
        if refine == 1:
            power = LombScargle(self.time, self.flux).power(freq_grid)
            self.period = 1 / freq_grid[np.argmax(power)]
            print("Period after refining using LombScargle is: %f" %self.period)
        # refine 2: minimize the sum(|dy|) of the phase folded lightcurve
        if refine == 2:
            sum_abs_dflux = []
            # phase fold over each frequency
            for freq in freq_grid:
                phase, flux = functions.phase_fold(1/freq, self.time, self.flux)
                # Take the sum of absolute differences between neighboring elements
                # should be minimum for the correctly folded lightcurve? 
                sum_abs_dflux.append(np.sum(np.abs(flux[1:] - flux[:-1])))
            self.period = 1 / freq_grid[np.argmin(sum_abs_dflux)]
            print("Period after refining using dispersion min is: %f" %self.period)
        # refine 3: minimize standard deviation after selecting a few periods
        if refine == 3:
            sum_std = []
            freq_grid = np.linspace(fguess*.95, fguess*1.05, 500)
            # phase fold over each frequency
            for freq in freq_grid:
                phase, flux = functions.phase_fold(1/freq, self.time, self.flux)
                # now bin this flux with 10 pts per bin and find std
                sum_std.append(np.sum(functions.std_bin(flux[:(10*(len(self.folded_flux)//10))] - functions.remove_line(phase, flux, 10), 10)))
                #sum_std.append(np.sum([np.nanstd([element for element in flux[i:i+10]]) for i in range(0, len(flux)//10, 10)]))
            self.period = 1 / freq_grid[np.argmin(sum_std)]
            print("Period after minimizing std in binned flux is: %f" %self.period)
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
        else: N = 10       # Number of points per bin
        if self.folded_flux is None: self.phase_fold(self)
        # Future: throw away points that are more than 4 std from binned flux
        arr, ph = self.folded_flux, self.phase
        self.binned_phase = functions.avg_bin(ph, N)
        self.binned_folded_flux = functions.avg_bin(arr, N)
        print(len(functions.remove_line(self.phase, self.folded_flux, N)))
        self.binned_folded_flux_err = functions.std_bin(self.folded_flux[:(10*(len(self.folded_flux)//10))] - functions.remove_line(self.phase, self.folded_flux, N), N)
    # Plot the orignal and phase folded lightcurve
    def plot_lc(self, savedir=None):
        if self.phase is None: self.phase_fold()
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7,15))
        axes[0].plot(self.time[(self.flux > 0.95) & (self.flux < 1.05)], self.flux[(self.flux > 0.95) & (self.flux < 1.05)], "x")
        axes[0].set_title("Time and flux", family="serif")
        axes[0].set_xlabel("Time [JD]", family="serif")
        axes[0].set_ylabel("Flux", family="serif")
        axes[0].set_xlim((0,20))
        axes[1].plot(self.phase[(self.folded_flux > 0.95) & (self.folded_flux < 1.05)], self.folded_flux[(self.folded_flux > 0.95) & (self.folded_flux < 1.05)],
                                 '.', color='orange', markersize=1)
        axes[1].set_title("Folded flux folded at period %f" %self.period, family="serif")
        axes[1].errorbar(x=self.binned_phase , 
                        y=self.binned_folded_flux, 
                        yerr=self.binned_folded_flux_err, marker='o',ecolor='red', 
                        color='k', label='binned folded flux', ls='none')
        axes[1].set_xlabel("Phase", family="serif")
        axes[1].set_ylabel("Flux", family="serif")
        axes[2].plot(1/self.freq, self.power, color='pink')
        axes[2].set_title("LombScargle periodogram for the lightcurve", family="serif")
        axes[2].set_xlabel("Period [day]", family="serif")
        axes[2].set_ylabel("Power", family="serif")
        axes[2].set_xlim((0, 30))
        for ax in axes: ax.legend()
        if not savedir: savedir = "../figures/%s_fig.png" % self.id
        plt.savefig(savedir)
    # write the phase folded lightcurve to a file + append period to the period file
    def write_lc(self, **kwargs):
        if self.binned_folded_flux is None: self.bin_lc(**kwargs)
        fname = "../lightcurves/folded_lightcurves/%s" % self.id
        with open(fname, "w") as f:
            # Number of points in LC
            N = len(self.binned_folded_flux)
            f.write("%d\n" % N)
            # Write for two periods
            for iter in range(0, 2, 1):
                for phase, flux, flux_err in zip(self.binned_phase, self.binned_folded_flux, self.binned_folded_flux_err):
                    phase_to_time = iter * self.period + phase * self.period
                    f.write("%f\t%f\t%f\n" % (phase_to_time, flux, flux_err))
        fname = "../lightcurves/periods.txt"
        with open(fname, "a") as f:
            # Entries are lc id, period and period flag
            good_period = int(input("Period found? [0 and 1]"))
            f.write("%d\t%f\t%d\n" %(int(self.id), self.period, int(good_period)))
    # Wrapper to call all of the above functions with just one call
    def generate_lc(self, lc_name):
        self.load_lc(lc_name)
        self.lombscargle()
        self.guess_period(refine=3)
        self.phase_fold()
        self.bin_lc()
        self.plot_lc()
        self.write_lc()
    # Wrapper function to refine period
    def user_period(self, lc_name):
        self.load_lc(lc_name)
        print("LC ID IS: ", lc.id)
        self.lombscargle()
        self.guess_period(refine=2)
        self.phase_fold()
        self.bin_lc()
        self.plot_lc(savedir="../figures/tst/%s_fig.png" % self.id)
        cont = input("Press enter to pass or any other key to try different period")
        while(cont):
            refine = int(input("Choice of period refinement [1,2,3]"))
            self.refine_guess(refine)
            self.phase_fold()
            self.bin_lc()
            self.plot_lc(savedir="../figures/tst/%s_fig.png" % self.id)
            cont = input("Try different period? ")
            
        self.write_lc()
        print("\n*****************************************\n")
start = int(sys.argv[1])

for lk_name in lc_list[start:]:
    print("lc_number: %i" %start)
    lc = lightcurve()
    lc.user_period(lk_name)
    start += 1

'''
class helper_functions:
    # I don't calculate with 
    def get_chi_square(arr):
        arr = 
'''