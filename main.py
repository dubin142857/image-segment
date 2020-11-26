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
f_raw = skimage.io.imread(image_path)/256
f     = color.rgb2gray(f_raw)
f     = np.asarray(f, dtype=float)
#观察图像不加噪声模糊
image=f.copy()
grad_image_x,grad_image_y=optfun.gradient_matrix(image)
grad_image_norm_2=grad_image_x**2+grad_image_y**2

#level set method
u=optfun.initalizePhi(image,d)
g_force_image=optfun.g_force(grad_image_norm_2)
max_g=np.maximum(g_force_image,0)
min_g=np.minimum(g_force_image,0)
grad_x_g,grad_y_g=optfun.gradient_matrix(g_force_image)
max_grad_xg=np.maximum(grad_x_g,0)
min_grad_xg=np.minimum(grad_x_g,0)
max_grad_yg=np.maximum(grad_y_g,0)
min_grad_yg=np.minimum(grad_y_g,0)
# row=image.shape[0]
# col=image.shape[1]
# h=1/(max(row,col)-1)
# delta_t=h*cfl
for i in range(max_iter):
    cur_matrix,cur_multiy_grad=optfun.gauss_cur(u)
    upwind_plus               =optfun.non_osc_upwind_plus(u)
    upwind_minus              =optfun.non_osc_upwind_minus(u)
    diff_x_fw,diff_y_fw       =optfun.forward_diff(u)
    diff_x_bc,diff_y_bc       =optfun.back_diff(u)
    u+=delta_t*(g_force_image*cur_multiy_grad+level_alpha*max_g*upwind_plus+\
                min_g*upwind_minus+max_grad_xg*diff_x_bc+\
                min_grad_xg*diff_x_fw + max_grad_yg*diff_y_bc+\
        min_grad_yg*diff_y_fw)
    #if not output_image =="None":
    print(i)
    if i%iprint==0:
        fig=plt.figure()
        ax=fig.add_subplot()
        ax.imshow(f,cmap='gray')
        ax.set_title("image and segment line",fontsize=20);ax.set_xticks([]);ax.set_yticks([])
        ax.contour(u,[0],colors='r',linewidth=2)
        #ax.show(block=False);ax.pause(0.01)
        fig.savefig('i'+output_image,dpi=200)
        
    
    


