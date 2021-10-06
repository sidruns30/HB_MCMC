import matplotlib.pyplot as plt
import matplotlib
import glob, os
import numpy

matplotlib.use('agg')

# Plot the lightcurves from the first 10000 chains of mcmc
data_dir = os.getcwd() + "/../lightcurves/folded_lightcurves"
lc_data =  glob.glob(data_dir + "/*")
lc_ids = [lc[len(data_dir)+1:] for lc in lc_data]
#lc_ids = [lc[len(data_dir)+1:-4] for lc in lc_data]

for lc, id in zip(lc_data, lc_ids):
    try:
        print("Loading lc: %s" %id)
        data = []
        # Read lightcurve data
        with open(lc, "r") as f:
            for i,line in enumerate(f.readlines()):
                if (i==0): N = line
                else: data.append(line.split('\t', 2))
        print("Lines: ", N)
        time, lc, model = [], [], [] 
        for d in data:
            time.append(float(d[0]))
            lc.append(float(d[1]))
            model.append(float(d[2]))
        print(len(time), len(lc), len(model))
        # Read model parameters 
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

        # stop plotting after the first gap of a day
        time, lc, model = map(numpy.array, (time, lc, model))
        gaps = time[1:] - time[:-1]
        try: index = numpy.where(gaps > 10)[0][0]
        except: index = len(time)
        time, lc, model = time[:index], lc[:index], model[:index]

        # Plot
        plt.style.use("dark_background")
        plt.figure(figsize=(10, 5))
        plt.title("Binned and folded lightcurve (%s)" %id, fontsize=30, family="serif")
        #plt.gcf().text(0.2, 0.78, some_string[:140], fontsize=8, color='lime')
        #plt.gcf().text(0.25, 0.73, some_string[140:], fontsize=8, color='lime')
        plt.plot(time,lc, color='r', marker='o', linewidth=0.5, label="Real Lightcurve")
        #plt.plot(time, model, 'o', color='yellow', linewidth=1, label="Computed lightcurve")
        plt.xlabel("Time", fontsize=15, family="serif")
        plt.ylabel("Flux", fontsize=15, family="serif")
        #plt.xlim((2225, 2255))
        plt.ylim((0.95, 1.05))
        plt.grid(alpha=0.5)
        plt.savefig(os.getcwd() + '/../figures/mcmc/model_lightcurve_%s.png' %id)
        plt.close()
        print("Generated plot: ", id)
    except: pass
