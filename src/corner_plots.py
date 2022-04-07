import corner
import numpy as np
import sys, os
import matplotlib.pyplot as plt

lc_id = sys.argv[1]

# Load the mcmc chain
def load_chain(lc_id, remove_period=True):
    cwd = os.getcwd()
    chain_dir = cwd + "/../data/chains/chain.%s.txt" % lc_id
    chain = np.loadtxt(chain_dir, max_rows=110000)
    # remove the last column
    chain = chain[:-1,:]
    if chain.shape[-1] > 11: alpha_flag = True
    logL = chain[:,1]
    if remove_period:
        # The 5th element is the period 
        samples = np.hstack((chain[:,2:4], chain[:,5:]), )
        if not alpha_flag: labels = ["$\\log_{10}M_1$", "$log_{10}M_2$", "$e$", "i", "$\\Omega$", "$\\omega$", 
                                    "$T_0$", "$\\log_{10}rr_1$", "$\\log_{10}rr_2$", "$$"]
        else:
            labels = ["$\\log_{10}M_1$", "$log_{10}M_2$", "$e$", "i", "$\\Omega$", "$\\omega$", 
            "$T_0$", "$\\log_{10}rr_1$", "$\\log_{10}rr_2$", "$\\mu_1$", "$\\tau_1$", 
            "$\\mu_2$", "$\\tau_2$", "$\\alpha_{ref1}$", "$\\alpha_{ref2}$"]
    else: 
        samples = chain[:,2:]
        if not alpha_flag: labels = ["$\\log_{10}M_1$", "$log_{10}M_2$", "P", "$e$", "i", "$\\Omega$", 
                                    "$\\omega$", "$T_0$", "$\\log_{10}rr_1$", "$\\log_{10}rr_2$"]
        else:
            labels = ["$\\log_{10}M_1$", "$log_{10}M_2$", "P", "$e$", "i", "$\\Omega$", 
                    "$\\omega$", "$T_0$", "$\\log_{10}rr_1$", "$\\log_{10}rr_2$", "$\\mu_1$", "$\\tau_1$",
                          "$\\mu_2$", "$\\tau_2$", "$\\alpha_{ref1}$", "$\\alpha_{ref2}$"]
    return logL, samples, labels

logL, samples, labels = load_chain(lc_id)
samples = samples[20000:,:]
logL_0 = logL[20000:]
samples[:,4] = np.random.rand(samples.shape[0])
print(np.min(samples), np.max(samples), np.argmin(samples), np.argmax(samples))
iternum = 10*np.arange(len(logL_0))

#plt.plot(iternum, logL_0)
#print(samples.shape)
plt.hist(np.power(10, samples[1000:,0]) / np.power(10, samples[1000:,1]))
#fig = corner.corner(samples, labels=labels,quantiles=[0.16, 0.5, 0.84],
#                       show_titles=True, title_kwargs={"fontsize": 12})
#plt.savefig("../data/figures/corner/%s.png" %lc_id)
#plt.close()
#plt.plot(logL_0)
#plt.savefig("test.png")