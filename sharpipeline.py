import sys, string, os
import numpy as np

import sharpener as sharpy

file = sys.argv[1]
cfg = open(file)

spar=sharpy.sharpener(file)

run = spar.go(spar.cfg_par)

if run == 0: 
    print '\t+------+\n\t Done \n\t+------+'