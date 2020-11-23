# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:28:37 2020

@author: Admin
"""


import numpy as np

u=np.ones((10,10))
u[0,:]=0;u[:,0]=0;
a=np.exp(np.array([1,2]))