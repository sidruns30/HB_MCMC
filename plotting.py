import matplotlib.pyplot as plt

filename = "110602878.out"
data = []

# Read lightcurve data
with open(filename, "r") as f:
    for i,line in enumerate(f.readlines()):
        if (i==0): pass
        else: data.append(line.split('  ', 2))

time, lc, model = [], [], [] 
for d in data:
    time.append(float(d[0]))
    lc.append(float(d[1]))
    model.append(float(d[2]))

# Read model parameters 
par_file = "chain.110602878.dat"
N_iter = 19887
params = []
with open(par_file, "r") as f:
    for i, line in enumerate(f.readlines()):
        if (i != N_iter): pass
        else: 
            params.append(line.split(' '))
param_list = ["$M_1$", "$M_2$", "$P$", "$e$", "$i$", "$\\omega$", "$\\Omega_0$", "$T_0$", "$F_0$", "$r_{R1}$", "$r_{R2}$"]
some_string = ""


for i, val in enumerate(params[0][2:-1]): 
    some_string += " | " + str(param_list[i]) + " = " + val 

# Plot
plt.style.use("dark_background")
plt.figure(figsize=(10, 5))
plt.title("Computed vs model lightcurve", fontsize=30, font="serif")
plt.gcf().text(0.2, 0.8, some_string[:140], fontsize=7)
plt.gcf().text(0.25, 0.75, some_string[140:], fontsize=7)
plt.plot(time, lc, 'r-', linewidth=0.5, label="Real Lightcurve")
plt.plot(time, model, 'w--', linewidth=0.5, label="Computed lightcurve")
plt.xlabel("Time", fontsize=15, font="serif")
plt.ylabel("Flux", fontsize=15, font="serif")
plt.xlim((2225, 2255))
plt.ylim((0.8, 1.2))
plt.savefig('Generated_Lightcurve.png')
plt.show()