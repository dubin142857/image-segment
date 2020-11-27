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
    # u=np.ones_like(image)
    # u=-d*u
    # u[0:3,:]=d;u[:,0:3]=d;u[-3:,:]=d;u[:,-3:]=d
    
    row=image.shape[0]
    col=image.shape[1]
    center_x=0.5*row
    center_y=0.5*col
    u=np.zeros_like(image)
    for i in range(row):
        for j in range(col):
            u[i,j]=np.sqrt((i-center_x)**2+(j-center_y)**2)-200
    return u

def reinitial2d(phi,steps):
    """
    reinitialization
    """
    row=phi.shape[0]
    col=phi.shape[1]
    h=1/(max(row,col)-1)
    eps=1e-9
    cfl=0.5
    dt=h*cfl
    dx=1/(row-1)
    dy=1/(col-1)
    for k in range(steps):
        xp,yp=forward_diff(phi)
        xn,yn=back_diff(phi)
        xp/=dx;yp/=dy;
        xn/=dx;yn/=dy;
        phi_p=np.int64(phi>0)
        phi_n=1-phi_p
        godnov_p=np.sqrt(np.maximum((np.maximum(xn,0))**2,(np.minimum(xp,0))**2)+np.maximum((np.maximum(yn,0))**2,\
                 np.minimum(yp,0)**2))
        godnov_n=np.sqrt(np.maximum((np.minimum(xn,0))**2,(np.maximum(xp,0))**2)+np.maximum((np.minimum(yn,0))**2,\
                 np.maximum(yp,0)**2))
        phi = phi-dt*phi_p*(godnov_p-1)*phi/(np.sqrt(phi**2+godnov_p**2*dx*dy+eps))\
            -dt*phi_n*(godnov_n-1)*phi/(np.sqrt(phi**2+godnov_n**2*dx*dy+eps))
    return phi
            

def gradient_matrix(image):
    """
    compute gradient

    Parameters
    ----------
    image : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """  
    image_x_fw,image_y_fw=forward_diff(image)
    image_x_bc,image_y_bc=back_diff(image)
    gradient_x=image_x_fw+image_x_bc
    gradient_x[:,1:-1]=0.5*(gradient_x[:,1:-1])
    gradient_y=image_y_fw+image_y_bc
    gradient_y[1:-1,:]=0.5*(gradient_y[1:-1,:])
    return gradient_x,gradient_y          
def forward_diff(image):
    """
    Usage:use to the right difference of each point in image
          and the right difference of the first col is all zero 
    """
    dx=1
    dy=1
    # row=image.shape[0]
    # col=image.shape[1]
    # dx=1/(row-1)
    # dy=1/(col-1)
    diff_x_fw=np.zeros_like(image)
    diff_x_fw[:,0:-1]=np.diff(image,1,1)/dx
    diff_x_fw[:,-1]=-image[:,-1]
    diff_y_fw=np.zeros_like(image)
    diff_y_fw[0:-1,:]=np.diff(image,1,0) /dy  
    diff_y_fw[0:-1,:]=-image[-1,:]
    return diff_x_fw,diff_y_fw

def back_diff(image):
    """
    Usage:use to the left difference of each point in image
          and the left difference of the last col is all zero 
    """
    dx=1
    dy=1
    # row=image.shape[0]
    # col=image.shape[1]
    # dx=1/(row-1)
    # dy=1/(col-1)
    diff_x_bc=np.zeros_like(image)
    diff_x_bc[:,1:]=np.diff(image,1,1)/dx
    diff_x_bc[:,0]=image[:,0]
    
    diff_y_bc=np.zeros_like(image)
    diff_y_bc[1:,:]=np.diff(image,1,0)/dy
    diff_y_bc[0,:]=image[0,:]
    return diff_x_bc,diff_y_bc

def g_force(image,g_force_lambda):
    """
    Usage:define force g function,which form is 1/(1+s),s is square of gradient image 
    """
    return 1./(1+g_force_lambda*image)

def boundry_pad(image,mid,boundry_type="zero_pad"):
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
    if boundry_type=="periodic":
        #填补周期边界
        observe_image=np.zeros((image.shape[0]+2*mid,image.shape[1]+2*mid))
        observe_image[0:mid,0:mid]=image[-mid:,-mid:].copy()
        observe_image[0:mid,mid:-mid]=image[-mid:,:].copy()
        observe_image[0:mid,-mid:]=image[-mid:,0:mid].copy()
        observe_image[mid:-mid,0:mid]=image[:,-mid:].copy()
        observe_image[mid:-mid,mid:-mid]=image.copy()
        observe_image[mid:-mid,-mid:]=image[:,0:mid].copy()
        observe_image[-mid:,0:mid]=image[0:mid,-mid:].copy()  
        observe_image[-mid:,mid:-mid]=image[0:mid,:].copy()
        observe_image[-mid:,-mid:]=image[0:mid,0:mid].copy()
    if boundry_type=="mirror":
        observe_image=np.zeros((image.shape[0]+2*mid,image.shape[1]+2*mid))
        observe_image[0,0]=image[0,0].copy()
        observe_image[0,1:-1]=image[0,:].copy()
        observe_image[0,-1]=image[0,-1].copy()
        observe_image[1:-1,0]=image[:,0].copy()
        observe_image[1:-1,1:-1]=image.copy()
        observe_image[1:-1,-1]=image[:,-1].copy()
        observe_image[-1,0]=image[-1,0].copy()  
        observe_image[-1,1:-1]=image[-1,:].copy()
        observe_image[-1,-1]=image[-1,-1].copy()
    if boundry_type=="zero_pad":
        observe_image=np.zeros((image.shape[0]+2*mid,image.shape[1]+2*mid))
        observe_image[1:-1,1:-1]=image.copy()
    return observe_image
def central_diff(image):
    """
    Usage: use to compute the central differences of each point in image
    """
    dx=1
    dy=1
    # row=image.shape[0]
    # col=image.shape[1]
    # dx=1/(row-1)
    # dy=1/(col-1)
    # image_xx=np.zeros_like(image)
    # image_yy=np.zeros_like(image)
    # image_xy=np.zeros_like(image)
    # image_x =np.zeros_like(image)
    # image_y =np.zeros_like(image)
    # pad_image=boundry_pad(image, mid=1)
    # for i in range(image.shape[0]):
    #     for j in range(image.shape[1]):
    #         image_xx[i,j]=pad_image[i+2,j+1]+pad_image[i,j+1]-2*pad_image[i+1,j+1]
    #         image_yy[i,j]=pad_image[i+1,j+2]+pad_image[i+1,j]-2*pad_image[i+1,j+1]
    #         image_xy[i,j]=pad_image[i+2,j+2]+pad_image[i,j]-pad_image[i,j+2]-\
    #                       pad_image[i+2,j]
    #         image_xy[i,j]=0.25*image_xy[i,j]
    #         image_x [i,j]=0.5*(pad_image[i+2,j+1]-pad_image[i,j+1])
    #         image_y [i,j]=0.5*(pad_image[i+1,j+2]-pad_image[i+1,j])
    # return image_xx/(dx**2),image_yy/(dy**2),image_xy/(dx*dy),image_x/dx,image_y/dy
    image_x=np.zeros_like(image)
    image_y=np.zeros_like(image)
    pad_image=boundry_pad(image,mid=1)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            image_x[i,j]=0.5*(pad_image[i+2,j+1]-pad_image[i,j+1])
            image_y[i,j]=0.5*(pad_image[i+1,j+2]-pad_image[i+1,j])
    return image_x,image_y
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
    # cur_matrix=np.zeros_like(image)
    # image_xx,image_yy,image_xy,image_x,image_y=central_diff(image)
    # image_x_2=image_x**2
    # image_y_2=image_y**2
    # grad_matrix=image_x_2+image_y_2
    # grad_matrix1=grad_matrix+1e-6*np.ones_like(image)
    # cur_matrix=image_xx*image_y_2- 2*image_y*image_x*image_xy+\
    #            image_yy*image_x_2
    # cur_matrix=cur_matrix/(grad_matrix1**(1.5))
    # cur_multiy_grad=cur_matrix/grad_matrix1
    # return cur_matrix,cur_multiy_grad,np.sqrt(grad_matrix)
    cur_matrix=np.zeros_like(image)
    eps=1e-6
    image_center_x,image_center_y  =central_diff(image)
    norm_g=np.sqrt(image_center_x**2+image_center_y**2)
    norm_g1=norm_g+eps
    image_back_x,image_back_y      =back_diff(image)
    image_back_x/=norm_g1
    image_back_y/=norm_g1
    k1,nomatter1=forward_diff(image_back_x)
    nomatter2,k2=forward_diff(image_back_y)
    return k1+k2,norm_g
    
    
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
    
    