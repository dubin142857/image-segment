# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 14:54:38 2020

@author: Admin
"""


import numpy as np
import matplotlib.pyplot as plt

def initalizePhi(image,d):
    """
    Usage: use to initalize the curve
    Parameters:
        image:the image that is used to be exam
        d    :level set function,outside is +d ,inside is -d
    """
    u=np.ones_like(image)
    u[0,:]=0;u[:,0]=0;u[-1,:]=0;u[:,-1]=0
    u=-d*u
    return u

def forward_diff(image):
    """
    Usage:use to the right difference of each point in image
          and the right difference of the first col is all zero 
    """
    diff_x_fw=np.zeros_like(image)
    diff_x_fw[:,1:]=np.diff(image,1,1)
    
    diff_y_fw=np.zeros_like(image)
    diff_y_fw[1:,:]=np.diff(image,0,1)
    
    return diff_x_fw,diff_y_fw

def back_diff(image):
    """
    Usage:use to the left difference of each point in image
          and the left difference of the last col is all zero 
    """
    diff_x_bc=np.zeros_like(image)
    diff_x_bc[:,0:-1]=-np.diff(image,1,1)
    
    diff_y_bc=np.zeros_like(image)
    diff_y_bc[0:-1,:]=-np.diff(image)
    
    return diff_x_bc,diff_y_bc

def g_force(image):
    """
    Usage:define force g function,which form is 1/(1+s),s is square of gradient image 
    """
    return 1./(np.ones_like(image)+image)

def boundry_pad(image,mid,boundry_type="periodic"):
    """
    Usage: use to pad image boundry by the way of boundry_type.

    Parameters
    ----------
    image : s*s image matrix.
    mid : the number which each edge needs to pad.
    boundry_type : the boundry pad type. The default is "periodic".

    Returns
    -------
    (s+2*mid)*(s+2*mid) pading image matrix.

    """
    if boundry_type="periodic":
        #填补周期边界
        observe_image=np.zeros((image.shape[0]+2*n,image.shape[1]+2*n))
        observe_image[0:mid,0:mid]=image[-mid:,-mid:].copy()
        observe_image[0:mid,mid:-mid]=image[-mid:,:].copy()
        observe_image[0:mid,-mid:]=image[-mid:,0:mid].copy()
        observe_image[mid:-mid,0:mid]=image[:,-mid:].copy()
        observe_image[mid:-mid,mid:-mid]=image.copy()
        observe_image[mid:-mid,-mid:]=image[:,0:mid].copy()
        observe_image[-mid:,0:mid]=image[0:mid,-mid:].copy()  
        observe_image[-mid:,mid:-mid]=image[0:mid,:].copy()
        observe_image[-mid:,-mid:]=image[0:mid,0:mid].copy()
    return observe_image
       
def central_diff(image):
    """
    Usage: use to compute the central differences of each point in image
    """
    image_xx=np.zeros_like(image)
    image_yy=np.zeros_like(image)
    image_xy=np.zeros_like(image)
    image_x =np.zeros_like(image)
    image_y =np.zeros_like(image)
    pad_image=boundry_pad(image, mid=1)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            image_xx[i,j]=pad_image[i+2,j+1]+pad_image[i,j+1]-2*pad_image[i+1,j+1]
            image_yy[i,j]=pad_image[i+1,j+2]+pad_image[i+1,j]-2*pad_image[i+1,j+1]
            image_xy[i,j]=pad_image[i+2,j+2]+pad_image[i,j]-pad_image[i,j+2]-\
                          pad_image[i+2,j]
            image_xy[i,j]=0.25*image_xy[i,j]
            image_x      =0.5*(pad_image[i+2,j+1]-pad_image[i,j+1])
            image_y      =0.5*(pad_image[i+1,j+2]-pad_image[i+1,j])
    return image_xx,image_yy,image_xy,image_x,image_y

def gauss_cur(image):
    """
    Usage: use to compute the gauss curvature of each point in image

    Parameters
    ----------
    image : matrix.

    Returns
    -------
    curvature matrix and matrix which record curvature*the norm of gradent of image

    """
    cur_matrix=np.zeros_like(image)
    image_xx,image_yy,image_xy,image_x,image_y=central_diff(image)
    image_x_2=image_x**2
    image_y_2=image_y**2
    grad_matrix=image_x_2+image_y_2
    cur_matrix=image_xx*image_y_2-2*image_y*image_x*image_xy+\
               image_yy*image_x_2
    cur_matrix=cur_matrix/(grad_matrix**(1.5))
    cur_multiy_grad=cur_matrix/grad_matrix
    return cur_matrix,cur_multiy_grad

def non_osc_upwind_plus(image):
    """
    Uaage:use to compute non-oscillotary upwind discretization while c>0

    Parameters
    ----------
    image : matrix.

    Returns
    -------
    upwind matrix when c>0.

    """
    image_x_bc,image_y_bc=back_diff(image)
    image_x_fw,image_y_fw=forward_diff(image)
    upwind_plus=np.zeros_like(image)
    upwind_plus=np.maximum(image_x_bc,0)**2+np.minimum(image_x_fw,0)**2+\
                np.maximum(image_y_bc,0)**2+np.minimum(image_y_fw,0)**2
    upwind_plus=np.sqrt(upwind_plus)
    return upwind_plus

def non_osc_upwind_minus(image):
    """
    Uaage:use to compute non-oscillotary upwind discretization while c<=0

    Parameters
    ----------
    image : matrix.

    Returns
    -------
    upwind matrix when c>0.

    """
    image_x_bc,image_y_bc=back_diff(image)
    image_x_fw,image_y_fw=forward_diff(image)
    upwind_minus=np.zeros_like(image)
    upwind_minus=np.maximum(image_x_fw,0)**2+np.minimum(image_x_bc,0)**2+\
                np.maximum(image_y_fw,0)**2+np.minimum(image_y_bc,0)**2
    upwind_minus=np.sqrt(upwind_minus)
    return upwind_minus
    
    