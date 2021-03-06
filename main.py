# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 14:49:03 2020

@author: Admin
"""


import numpy as np
import sys
import skimage
from skimage import  color
from skimage import io
import optargv
import optfun
import matplotlib.pyplot as plt

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
f_raw = skimage.io.imread(image_path)
f     = color.rgb2gray(f_raw)
f     = np.asarray(f, dtype=float)
#观察图像不加噪声模糊
image=f.copy()
grad_image_x,grad_image_y=optfun.central_diff(image)
grad_image_norm_2=grad_image_x**2+grad_image_y**2

#level set method
u=optfun.initalizePhi(image,d)
g_force_image=optfun.g_force(grad_image_norm_2,g_force_lambda)
grad_x_g,grad_y_g=optfun.central_diff(g_force_image)
max_grad_xg=np.maximum(grad_x_g,0)
min_grad_xg=np.minimum(grad_x_g,0)
max_grad_yg=np.maximum(grad_y_g,0)
min_grad_yg=np.minimum(grad_y_g,0)

for i in range(max_iter):
    cur_matrix,sqrt_grad_u=optfun.gauss_cur(u)
    upwind_plus               =optfun.non_osc_upwind_plus(u)
    upwind_minus               =optfun.non_osc_upwind_minus(u)
    diff_x_fw,diff_y_fw       =optfun.forward_diff(u)
    diff_x_bc,diff_y_bc       =optfun.back_diff(u)
    u_next=u+delta_t*(g_force_image*cur_matrix*sqrt_grad_u\
                      +level_alpha*g_force_image*upwind_minus+\
                max_grad_xg*diff_x_fw+min_grad_xg*diff_x_bc\
                + max_grad_yg*diff_y_fw+min_grad_yg*diff_y_bc)
    #if not output_image =="None":
    print(i)
    if np.linalg.norm(u_next-u)<tol:
        u=u_next.copy()
        break
    u=u_next.copy()
    if i%restart_num==0:
        u=optfun.reinitial2d(u,level_steps)
    
    if i%iprint==0:
        fig=plt.figure()
        ax=fig.add_subplot()
        ax.imshow(f,cmap='gray')
        ax.set_title("image and segment line",fontsize=20);ax.set_xticks([]);ax.set_yticks([])
        ax.contour(u,[0],colors='r',linewidth=2)
        #ax.show(block=False);ax.pause(0.01)
        fig.savefig('i'+output_image,dpi=200)
        
    
    


