from subprocess import call
import glob, os, sys, numpy as np

periods = {}

# Open the text file and scan the period
def fetch_period(lc_id):
    data = np.loadtxt("../data/lightcurves/periods.txt")
    found = False
    for name, period, flag in data:
        if int(name) == int(lc_id):
            found = True
            break
        else: pass
    if found:
        period = np.log10(period)
        print("Success, period is %f log 10 days" %period)
        return period
    else: 
        print ("LC not found")
        return None  

'''
Reads the lightcurve given the path
'''
def read_lc(lc_id, lc_path):
    lc_time = []
    lc_flux = []
    print("Loading TIC: %s" % lc_id)
    with open(lc_path, "r") as fobj:
        for i, line in enumerate(fobj.readlines()):
            if i == 0: pass
            else:
                data = line.split()
                # only read points where 4th column is 0
                if (float(data[3]) == 0.):
                    lc_time.append(float(data[0]))
                    lc_flux.append(float(data[1]))
    return lc_id, np.array(lc_time), np.array(lc_flux)

'''
Folds lightcurve onto its period and repeats the array twice
'''
def fold_lc(lc_id, lc_time, lc_flux):
    period = (periods[lc_id])
    if period is not None:
        period = np.power(10., period)
    else:
        period = (lc_time[-1] - lc_time[0])
        
    print("Folding TIC %s on %f days" % (lc_id, period))
    phases = (lc_time % period)
    indices = np.argsort(phases)
    
    phases = phases[indices]
    folded_flux = lc_flux[indices]

    return phases, folded_flux

'''
Auxiliary function: Fit a line to a function
'''
def fit_line(xarr, yarr):
    xarr = np.array(xarr)
    yarr = np.array(yarr)
    X = np.average(xarr)
    Y = np.average(yarr)
    m = np.average((xarr - X)*(yarr - Y))/np.average(np.square(xarr - X))
    b = Y - m*X

    return m*xarr + b
'''
Auxiliary function: Remove a line from bins
'''
def remove_line(xarr, yarr, Npts):
    xarr = np.array(xarr)
    yarr = np.array(yarr)
    new_arr = np.array([])
    stddev = []
    for i in range(0, len(xarr) - len(xarr)%Npts, Npts):
        xpts = xarr[i:i+Npts]
        ypts = yarr[i:i+Npts]
        new_arr = np.concatenate((new_arr, 
                                  np.array(ypts - fit_line(xpts, ypts))),
                                 axis=-1)
        stddev.append(np.nanstd(np.array(ypts - fit_line(xpts, ypts))))
        
    return new_arr.ravel(), np.repeat(np.array(stddev), Npts)

'''
Cleans the folded lightcurve and throws away the outliers
'''
def clean_lc(lc_phases, lc_flux):
    bin_len = 20

    #Throw away obvious outliers
    flux_mean = np.mean(lc_flux)
    flux_std = np.nanstd(lc_flux - flux_mean)
    
    new_flux = []
    new_phase = []

    for phase, flux in zip(lc_phases, lc_flux):
        if abs(flux - flux_mean) < 7 * flux_std:
            new_flux.append(flux)
            new_phase.append(phase)
    
    # Now start binning and removing the outliers
    N = len(new_phase)
    line_flux = np.ones(N)
    good_index = []

    for i in range(0, N - N%bin_len, bin_len):
        y = fit_line(new_phase[i:i+bin_len], new_flux[i:i+bin_len])
        y.reshape(-1)
        y_std = np.nanstd(y)
        line_flux[i:i+bin_len] = y
        # Throw points more than 3 sigma from line fits
        for j in range(0, bin_len, 1):
            if abs(y[j] - new_flux[i + j]) < 7*y_std:
                good_index.append(i + j)
    
    new_phase = np.array(new_phase)[good_index]
    new_flux = np.array(new_flux)[good_index]
    
    N = len(new_phase)
    double_phases = np.zeros(2*N)
    double_flux = np.zeros(2*N)
    double_phases[:N] = new_phase
    double_phases[N:] = new_phase + new_phase[-1]
    double_flux[:N] = new_flux
    double_flux[N:] = new_flux
    return double_phases, double_flux

'''
Generates error by binning the lightcurves and finding the stddev
'''
def bin_lc(phase, flux):
    bin_len = 5
    N = len(phase)
    binned_phase = []
    binned_flux = []
    binned_err = []
    for i in range(0, N - N%bin_len, bin_len):
        x = phase[i:i+bin_len]
        y = flux[i:i+bin_len]
        binned_err.append(np.nanstd(y - fit_line(x, y)))
        binned_phase.append(np.average(x))
        binned_flux.append(np.average(y))
    map (np.array, (binned_phase, binned_flux, binned_err))
    return binned_phase, binned_flux, binned_err

'''
Writes to a file
'''
def write_lc(lc_id, phase, flux, err):
    fname = "../data/lightcurves/folded_lightcurves/" + lc_id + "_new.txt"
    with open(fname, "w") as f:
        f.write("%d\n" %len(phase))
        for p, fl, e in zip(phase, flux, err):
            f.write(f"{p}\t{fl}\t{e}\n")

'''
Combines all of the above functions to produce the final lightcurve
'''
def generate_lc(lc_id, lc_path):
    lc_phases, lc_folded_flux = fold_lc(*read_lc(lc_id, lc_path))
    c_phases, c_fflux = clean_lc(lc_phases, lc_folded_flux)
    b_phase, b_flux, b_err = bin_lc(c_phases, c_fflux)
    write_lc(lc_id, b_phase, b_flux, b_err)


lc_list = ["102289966", "110106107", "110602878", "116066534" , "118305806", "120684604", "127079833", "145974291",
           "147380831", "19971063", "365831185", "431305729"]

periods = {}

for lc in lc_list:
    period = fetch_period(lc)
    periods[lc] = period

for lc in lc_list:
    lc_path = "../data/lightcurves/original/" + lc + ".txt"
    generate_lc(lc_id=lc, lc_path=lc_path) 


##################################################
# MCMC RUN STARTS HERE
##################################################
'''

lc_dir = os.getcwd() + '/../lightcurves/original'
lc_list = glob.glob(lc_dir + '/*txt')

iter = int(sys.argv[1]) # selects the lightcurve #
NCHAINS = 2000000
tol = 0.01
burn_in  = 1

# check if the period is good or not, change the burn in accordingly
def check_period(ID):
    with open("../lightcurves/periods.txt") as f:
        lines = f.read().splitlines()
        for line in lines:
            line = line.split("\t")
            if (int(line[0]) == int(ID) and (int(line[2]) == 1)): return 2, float(line[1])
        # Set burn in to 1 otherwise
        else: return 1, 0

for lc in lc_list[iter:iter+2]:
    lc_id = lc[len(lc_dir)+1:-4]
    burn_in, period = check_period(lc_id)
    print("id, burn in and period are (%s, %d, %f)" %(lc_id, burn_in, period))
    # convert period to log10 days
    if burn_in == 2: period = np.log10(period)
    call(["./mcmc", "%d" %NCHAINS, lc_id, "%f" %tol, "%d" %burn_in, "%f" %period])
    '''