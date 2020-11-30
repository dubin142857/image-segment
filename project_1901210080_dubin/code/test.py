# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 17:07:22 2020

@author: Admin
"""

from random import choice
import numpy as np
r=np.random.rand()
str=list(range(1,9))
a=choice(str)
b=[1,2,3]
b[-1]=b[-2]
b[-2]=9
str