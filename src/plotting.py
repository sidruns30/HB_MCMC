import matplotlib.pyplot as plt
import matplotlib
import glob, os, sys
import numpy

matplotlib.use('agg')

id = sys.argv[1]
# Plot the lightcurves from the first 10000 chains of mcmc
data_dir = os.getcwd() + "/../data/lightcurves/mcmc_lightcurves/"
lc =  glob.glob(data_dir + id + ".out")[0]


try:
    print("Loading lc: %s" %id)
    data = []
    print(lc)
    # Read lightcurve data
    with open(lc, "r") as f:
        print(lc)
        for i,line in enumerate(f.readlines()):
            if (i==0): N = line
            else: data.append(line.split('  ', 2))
    print("Lines: ", N)
    time, lc, model = [], [], [] 
    for d in data:
        time.append(float(d[0]))
        lc.append(float(d[1]))
        model.append(float(d[2]))
    
    # Read model parameters 
    '''
    par_file = os.getcwd() + "/../chains/chain.%s.dat" %id
    N_iter = 199
    params = []
    with open(par_file, "r") as f:
        for i, line in enumerate(f.readlines()):
            if (i != N_iter): pass
            else: 
                params.append(line.split(' '))
    param_list = ["$M_1$", "$M_2$", "$P$", "$e$", "$i$", "$\\omega$", "$\\Omega_0$", "$T_0$", "$F_0$", "$r_{R1}$", "$r_{R2}$"]
    some_string = ""
    for i, val in enumerate(params[0][2:-1]):
        if i == 2: val = str(10**(float(val)))
        some_string += " | " + str(param_list[i]) + " = " + val 

    '''
    # stop plotting after the first gap of a day
    time, lc, model = map(numpy.array, (time, lc, model))
    gaps = time[1:] - time[:-1]
    try: index = numpy.where(gaps > 10)[0][0]
    except: index = len(time)
    time, lc, model = time[:index], lc[:index], model[:index]

    # Plot
    plt.style.use("dark_background")
    plt.figure(figsize=(10, 5))
    plt.title("Computed vs model lightcurve (%s)" %id, fontsize=30, family="serif")
    #plt.gcf().text(0.2, 0.78, some_string[:140], fontsize=8, color='lime')
    #plt.gcf().text(0.25, 0.73, some_string[140:], fontsize=8, color='lime')
    plt.plot(time, lc, 'r.', linewidth=0.5, label="Real Lightcurve")
    plt.plot(time, model, '.-', color='yellow', linewidth=1, label="Computed lightcurve")
    plt.xlabel("Time", fontsize=15, family="serif")
    plt.ylabel("Flux", fontsize=15, family="serif")
    #plt.xlim((2225, 2255))
    #plt.ylim((0.95, 1.05))
    plt.grid(alpha=0.5)
    plt.legend()
    plt.savefig(os.getcwd() + '/../data/figures/mcmc/model_lightcurve_%s.png' %id)
    plt.close()
    print("Generated plot: %s" % id)
except:pass
