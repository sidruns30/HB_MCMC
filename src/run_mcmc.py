from subprocess import call
import glob, os, sys, numpy

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
    if burn_in == 2: period = numpy.log10(period)
    call(["./mcmc", "%d" %NCHAINS, lc_id, "%f" %tol, "%d" %burn_in, "%f" %period])