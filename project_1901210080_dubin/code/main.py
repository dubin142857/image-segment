# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 16:15:00 2020

@author: Admin
"""

import  sys
import numpy as np
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

#initial image
image = np.ones((N,N))
exp_ham=optfun.metropolis(image,J,h,q,max_iter,k_B,init_temp,proposal_type)

# if not output_image == "None":
#     import matplotlib.pyplot as plt
#     fig = plt.figure()
#     ax = fig.add_subplot(1,3,1)
#     ax.imshow(f, cmap='gray')
#     ax.set_title("$u_{true}$", fontsize=20); ax.set_xticks([]); ax.set_yticks([])
#     ax = fig.add_subplot(1,3,2)
#     ax.imshow(f_noise, cmap='gray')
#     ax.set_title("f", fontsize=20); ax.set_xticks([]); ax.set_yticks([])
#     ax = fig.add_subplot(1,3,3)
#     ax.imshow(u1_heat, cmap='gray')
#     ax.set_title("$u_{compute}$", fontsize=20); ax.set_xticks([]); ax.set_yticks([])
#     fig.set_size_inches([15,5])
#     fig.savefig('heat-'+output_image, dpi=200)
