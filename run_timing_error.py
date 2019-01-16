from __future__ import print_function
import numpy as np
import subprocess
import os

# params : L spin B N
names = ['L', 's', 'B', 'N']
L = 1024
for spin in [0, 2]:
    for B in [2]:
        for N in [5]:
            params = [L, spin, B, N]
            outfilename = '_'.join(['timings', 'errors'] + [nm + str(x) for nm, x in zip(names, params)]) + '.csv'
            command = ' '.join(['bin/s2let_test_csv'] + [str(x) for x in params])
            print('Executing', command)
            print('Writing to', outfilename)
            outfile = open(outfilename, 'w+')
            p = subprocess.Popen(command, stdout=outfile, shell=True)
