from subprocess import call
import glob
import os
import sys

lc_dir = os.getcwd() + '/../lightcurves/original'
lc_list = glob.glob(lc_dir + '/*txt')

iter = int(sys.argv[1])

NCHAINS = 20000
tol = 0.01
burn_in  = 1

for lc in lc_list[iter:iter+1]:
    lc_id = lc[len(lc_dir)+1:-4]
    print(lc_id)
    call(["./mcmc", "%d" %NCHAINS, lc_id, "%f" %tol, "%d" %burn_in])