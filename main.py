# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 14:49:03 2020

@author: Admin
"""


import numpy as np
import sys
import skimage
from skimage import io
import optargv
import optfun

# read data from command line
if len(sys.argv) > 1 and (sys.argv[1] == "--help"or sys.argv[1] == "-h"):
    options = optargv.default_options()
    print("Usage:")
    print("      python main.py opt1=value1 opt2=value2...")
    print("      e.g python main.py --image_path=image/peppers256.png\n")
    print("Available options:default values")
    for k, v in options.items():
        print("", k, ":", v)
    sys.exit()
options, globalnames = optargv.setoptions(argv=sys.argv[1:], kw=None)
globals().update(globalnames)

#read image and convert to [0,1]
f = skimage.io.imread(image_path)/256
f = np.asarray(f, dtype=float)



