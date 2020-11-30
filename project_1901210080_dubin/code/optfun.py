# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 16:20:24 2020

optional funtions which may be used.

@author: Admin
"""
import numpy as np
from random import  choice
def signle_flip(image,q):
    """
    Usage: this function use single flip proposal to decide matrix Q
    Parameters:
        image:N*N matrix,which record the state of  Potts model.
    
        q    : each point in matrix image have q states.
    """
    N=int(image.shape[0])
    random_x=np.random.randint(N)
    random_y=np.random.randint(N)
    q_xy=image[random_x,random_y]
    state_list=list(range(1,q+1))
    state_list.remove(q_xy)
    # old_state=image[random_x,random_y]
    image[random_x,random_y]=choice(state_list)
    # delta_state=image[random_x,random_y]-old_state
    # new_image=image.copy()
    # new_image[random_x,random_y]=delta_state
    # #delta_hamilton=hamiltonian(new_image, J, h)
    return image

def periodic_pad(image):
    """
    Usage:use to periodic pad n*n image to (n+2)*(n+2) new_image 
    """
    n=int(image.shape[0])
    new_image=np.zeros((n+2,n+2))
    new_image[0,0]      =image[-1,-1]
    new_image[0,1:-1]   =image[-1,:]
    new_image[0,-1]     =image[-1,0]
    new_image[1:-1,0]   =image[:,-1]
    new_image[1:-1,1:-1]=image
    new_image[1:-1,-1]  =image[:,0]
    new_image[-1,0]     =image[0,-1]
    new_image[-1,1:-1]  =image[0,:]
    new_image[-1,-1]    =image[0,0]
    return new_image
def kronecker(a,b):
    """
    kronecker delta funtion
    """
    if a==b:
        return 1
    else:
        return 0
    
def window33(image_part):
    """
    Usage:use to compute the hamiltonian of x which is at the middle of 3*3 matrix
    image_part
    """
    result=0
    result=kronecker(image_part[1,1],image_part[1,0])+\
           kronecker(image_part[1,1],image_part[1,2])+\
           kronecker(image_part[1,1],image_part[2,1])+\
           kronecker(image_part[1,1],image_part[0,1])
    return result

def hamiltonian(image,J,h):
    """"
    Usage: this function is used to compute the hamiltonian of system image
    """
    N=int(image.shape[0])
    new_image=periodic_pad(image)
    hamilton=0
    for i in range(1,N+1):
        for j in range(1,N+1):
            hamilton+=window33(new_image[i-1:i+2,j-1:j+2])
    hamilton=-J*hamilton
    if h!=0:
        hamilton+=-h*(np.sum(image))
    return hamilton

def metropolis(image,J=1,h=0,q=3,max_iter=300,k_B=1,init_temp=1,proposal_type="signle_flip"):
    """
    Usage:use to metropolis algorithm to solve MC intergal.
    proposal_maxtix is default "signle_flip".
    Parameters:
        image        :the state of system
        J            :hamiltonian parameters
        h            :hamiltonian parameters
        q            :the number of particle states
        k_B          :pdf parameters
        init_temp    :initial temperature
        max_iter     ：given temperature,MC method get max_iter states to compute integral.
        proposal_type:proposal matrix type
    """
    beta=1/(k_B*init_temp)
    N=int(image.shape[0])
    new_image=np.zeros_like(image)
    hamilton_list=[]
    hamilton_list.append(hamiltonian(image, J, h))
    if proposal_type =="signle_flip":
        for i in range(max_iter):
            new_image=signle_flip(image,q)
            hamilton_list.append(hamiltonian(new_image,J,h))
            delta_hamilton=hamilton_list[-1]-hamilton_list[-2]
            if delta_hamilton<=0:
                continue
            else:
                r=np.random.rand()
                if r<=np.exp(-beta*delta_hamilton):
                    continue
                else:
                    hamilton_list[-1]=hamilton_list[-2]
    expectation_hamilton=sum(hamilton_list)/max_iter
    return expectation_hamilton
        
    
            

