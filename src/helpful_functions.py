import numpy as np
import warnings
import os
import matplotlib.pyplot as plt
import sys
import time
import re
import json
import pprint
import requests
import glob
from urllib.parse import quote as urlencode
from scipy.optimize import curve_fit
from astropy.timeseries import BoxLeastSquares
from scipy.signal import find_peaks

## Sub-functions for lightcurve handling
'''
Reads the lightcurve given the path
'''
def read_lc(lc_path):
    lc_time = []
    lc_flux = []
    print("Loading lc %s" % lc_path)
    with open(lc_path, "r") as fobj:
        for i, line in enumerate(fobj.readlines()):
            if i == 0: pass
            else:
                data = line.split()
                # only read points where 4th column is 0
                if (float(data[3]) == 0.):
                    lc_time.append(float(data[0]))
                    lc_flux.append(float(data[1]))
    lc_time = np.array(lc_time)
    lc_flux = np.array(lc_flux)
    # Skips points with gaps
    ind = np.where(lc_time[1:] - lc_time[:-1] > 0.3)[0][0]
    lc_time = np.concatenate((lc_time[:ind-10], lc_time[ind+10:]))
    lc_flux = np.concatenate((lc_flux[:ind-10], lc_flux[ind+10:]))
    return lc_time, lc_flux

'''
Find the period using Box Least Squares
'''
def find_period(time, flux):
    model = BoxLeastSquares(time, flux, err)
    periodogram = model.autopower(0.1, oversample=2)
    return periodogram.period[np.argmax(periodogram.power)]

'''
Folds lightcurve onto its period and repeats the array twice
'''
def fold_lc(lc_time, lc_flux, period):
    if period is not None:
        period = np.power(10., np.log10(period))
    else:
        period = find_period(lc_time, lc_flux)
        
    print("Folding lightcurve on %f days" % (period))
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
Auxiliary Function: Average array with Npts per bin; returned shape is len(arr)//Npts
'''
def avg_bin(arr, Npts):
    return np.array([np.average(arr[i:i+Npts]) for i in range(0, len(arr) - len(arr)%Npts, Npts)])
'''
Auxiliary Function: Take std dev of array with Npts per bin
'''
def std_bin(arr, Npts):
    return np.array([np.nanstd(arr[i:i+Npts]) for i in range(0, len(arr) - len(arr)%Npts, Npts)])

'''
Auxiliary Function: Remove a parabola from bins
'''
def remove_par(xarr, yarr, Npts):
    new_arr = np.array([])
    for i in range(0, len(xarr) - len(xarr)%Npts, Npts):
        xpts = xarr[i:i+Npts]
        ypts = yarr[i:i+Npts]
        model = lambda x, a, b, c: a*x**2 + b*x + c
        pars, _ = curve_fit(model, xpts, ypts)
        new_arr = np.concatenate((new_arr, np.array(ypts - model(xpts, *pars))), axis=-1)
    return new_arr.ravel()


'''
Cleans the folded lightcurve and throws away the outliers
lc_period should be in days
'''
def clean_lc(lc_phases, lc_flux, bin_len = 20):

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

    good_index = []

    for i in range(0, N - N%bin_len, bin_len):
        y = fit_line(new_phase[i:i+bin_len], new_flux[i:i+bin_len])
        y.reshape(-1)
        y_std = np.nanstd(y)
        # Throw points more than 7 sigma from line fits
        for j in range(0, bin_len, 1):
            if abs(y[j] - new_flux[i + j]) < 7*y_std:
                good_index.append(i + j)
    
    new_phase = np.array(new_phase)[good_index]
    new_flux = np.array(new_flux)[good_index]

    return new_phase, new_flux

'''
Generates error by binning the lightcurves and finding the stddev
'''
def bin_lc(phase, flux, lc_period, method="linear", bin_len=4, uniform_err=False):
    binned_phase = avg_bin(phase, bin_len)
    binned_folded_flux = avg_bin(flux, bin_len)
    
    # For brevity
    N = bin_len
    if method == "linear":
        binned_err = std_bin(flux[:(N*(len(flux)//N))] - 
                                                    remove_line(phase, flux, N), N)
    elif method == "quadratic":
        binned_err = std_bin(flux[:(N*(len(flux)//N))] - 
                                                    remove_par(phase, flux, N), N)
    
    print(binned_err.shape)
    if uniform_err: binned_err = np.ones(flux.shape) * uniform_err
    
    # double the phases anf fluxes
    # notic we change N again  
    N = len(binned_phase)
    double_phases = np.zeros(2*N)
    double_flux = np.zeros(2*N)
    double_err = np.zeros(2*N)
    double_phases[:N] = binned_phase
    double_phases[N:] = binned_phase + lc_period
    double_flux[:N] = binned_folded_flux
    double_flux[N:] = binned_folded_flux
    double_err[:N] = binned_err
    double_err[N:] = binned_err
    
    return double_phases, double_flux, double_err

'''
Writes to a file
'''
def write_lc(lc_id, phase, flux, err):
    lcdir = os.getcwd() + "/../data/lightcurves/folded_lightcurves/"
    fname = lcdir + lc_id + "_new.txt"
    with open(fname, "w") as f:
        f.write("%d\n" %len(phase))
        for p, fl, e in zip(phase, flux, err):
            f.write(f"{p}\t{fl}\t{e}\n")

'''
 Lookup the TIC in the data file - check if the correct period was found
 if not then raise a warning and use the original lightcurve instead. 
 Otherwise load the old lightcurve, clean it using the found period
 write the lightcurve to a file. Returns the burn in too

 Parameters:
             TIC (str): TIC id of the lightcurve
 Returns:
             burn_in:   boolean to check if correct period was found
                        or not. If it was then a new lightcurve file
                        is also created in the data directory
'''
def process_lc(TIC, plot=False):

    # Read the lightcurve data file
    lc_file_path = os.getcwd() + "/../data/lightcurves/original/" + TIC + ".txt"
    
    if (not os.path.isfile(lc_file_path)): 
        raise ValueError ("File %s not found " % lc_file_path)
    
    lc_time, lc_flux = read_lc(lc_file_path)
     
    # lc period lookup
    lc_periods = np.loadtxt("../data/lightcurves/periods.txt")
    for i, lc_id in enumerate(lc_periods[:, 0]):
        success = 0
        if int(lc_id) == int(TIC):
            success = 1
            period = lc_periods[i, 1]
            break
 
    if success:
        print("TIC %d period in file is %f days" % (int(TIC), period))
    
    if (not success):
        warnings.warn("TIC not found in the lightcurve file; estimating period from Box Least squares")
        period = find_period(lc_time, lc_flux)
        print("Guessed period is %f", period)
    

    phases, folded_flux = fold_lc(lc_time, lc_flux, period) 
    clean_phases, clean_flux = clean_lc(phases, folded_flux, 
                                        bin_len = 20)
    binned_phase, binned_flux, binned_err = bin_lc(clean_phases,
                                                clean_flux, 
                                                period,
                                                bin_len=4)
    write_lc(TIC, binned_phase, binned_flux, binned_err)
    
    print("Successfully wrote binned lc to file")
    
    if plot:
        plt.errorbar(binned_phase, binned_flux, binned_err, ecolor="r")
        plt.savefig("test_lc.png")
    return 1, period

########################MAST QUERY FUNCTIONS######################
def mast_query(request):
    """Perform a MAST query.
    
        Parameters
        ----------
        request (dictionary): The MAST request json object
        
        Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    # Base API url
    request_url='https://mast.stsci.edu/api/v0/invoke'    
    
    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    
    # Perform the HTTP request
    resp = requests.post(request_url, data="request="+req_string, headers=headers)
    
    # Pull out the headers and response content
    head = resp.headers
    content = resp.content.decode('utf-8')

    return head, content

def resolve_name(TIC): 
        request = {'service':'Mast.Name.Lookup',
                   'params':{'input':'TIC ' + TIC,
                             'format':'json'},
        }
    
        headers, out_string = mast_query(request)
    
        out_data = json.loads(out_string)
        return out_data

'''
    Now construct the color and the magnitude information
    Search MAST for the given TIC, extract color and distance
    information, and write the info to a file
    If the any of the errors are not known then set them
    to be a large number and raise a warning
'''
def get_magnitudes(TIC):
    print("Resolving object RA and Dec...")
    # Get the RA and Dec data
    resolved_tic = resolve_name(TIC)
    try:
        ra = resolved_tic['resolvedCoordinate'][0]['ra']
        dec = resolved_tic['resolvedCoordinate'][0]['decl']
    except:
        raise ValueError ("An exception occured while resolving the name ")
        
    print("Fetching magnitude data from MAST...")
    # Create MAST request object
    mast_request = { 'service':'Mast.Catalogs.Tic.Cone',
                        'params': {'ra':ra, 'dec':dec, 'radius':0.006},
                        'format':'json',
                        'pagesize':2000,
                        'removenullcolumns':True,
                        'timeout':30,
                        'clearcache':True,
                        'removecache':True
                        }
    head, mast_data = mast_query(mast_request)
    data = json.loads(mast_data)
    print("Number of objects found = %d" % len(data["data"]))
    
    # The star of interest is assumed to have the highest Gaia magnitude
    data = data["data"]
    
    Gmags = []
    for obj in data:
        Gmags.append(obj["GAIAmag"])
    
    star_index = np.argmin(Gmags)
    
    Gmags = np.sort(Gmags)
    # estimate blending
    if (len(Gmags > 1)):
        f1_over_f2 = pow(10., (Gmags[0] - Gmags[1]) / (-2.5)) 
        blending = 1 / (f1_over_f2 + 1)
    else:
        blending = 0
    
    print("Estimated blending is: %f" %blending)
    
    # Now extract the magnitude data from the list
    for i, obj in enumerate(data):
        if i == star_index:
            Bmag = obj["Bmag"]
            Bmag_e = obj["e_Bmag"]
            Vmag = obj["Vmag"]
            Vmag_e = obj["e_Vmag"]
            Gmag = obj["GAIAmag"]
            Gmag_e = obj["e_GAIAmag"]
            Tmag = obj["Tmag"]
            Tmag_e = obj["e_Tmag"]
            Distance = obj["d"]
            Distance_e = obj["e_d"]
            break
            
    # Compute colors and errors
    if Gmag_e and Distance and Distance_e:
        Gmag_e = Gmag_e + (5 * Distance_e / Distance / np.log(10))
        
    fdir = os.getcwd() + "/../data/magnitudes/"
    fname = fdir + "%s.txt" % TIC
    
    with open(fname, "w") as file:
        if not Distance:
            warnings.warn("Distance measurement of the object does not exist")
            file.write(f"{1000.}\n")
            
        else:  file.write(f"{Distance}\n")
        
        if not Gmag: 
            warnings.warn("GAIA magnitude of the object does not exist")
            file.write(f"{10.}\t{10000.}\n")
            
        else:
            if not Gmag_e:
                    file.write(f"{Gmag}\t{1.}\n")
            else:   file.write(f"{Gmag}\t{Gmag_e}\n")

        if (not Bmag) or (not Vmag):
            warnings.warn("B/V magnitude of the object does not exist")
            file.write(f"{0.}\t{1000.}\n")

        else:  
            if (not Bmag_e) or (not Vmag_e):
                file.write(f"{Bmag-Vmag}\t{1.}\n")
            else: file.write(f"{Bmag-Vmag}\t{abs(Bmag_e) + abs(Vmag_e)}\n")
            
        if (not Vmag) or (not Gmag):
            warnings.warn("V/G magnitude of the object does not exist")
            file.write(f"{0.}\t{1000.}\n")
        
        else:  
            if (not Vmag_e) or (not Gmag_e):
                file.write(f"{Vmag-Gmag}\t{1.}\n")
            else:  file.write(f"{Vmag-Gmag}\t{abs(Vmag_e) + abs(Gmag_e)}\n")
                
        if (not Gmag) or (not Tmag):
            warnings.warn("G/T magnitude of the object does not exist")
            file.write(f"{0.}\t{1000.}\n")
            
        else:    
            if (not Gmag_e) or (not Tmag_e):
                file.write(f"{Gmag-Tmag}\t{1.}\n")
            else:  file.write(f"{Gmag-Tmag}\t{abs(Gmag_e) + abs(Tmag_e)}\n")
    
        print("Magnitudes written to file %s" %fname)

        return blending

