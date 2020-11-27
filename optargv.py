# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 09:00:28 2020
@author: Admin
该程序是用来接收命令行参数
"""


import sys, getopt

def  default_options():
    options = {
        "--image_path":"image/strawberry.jpg",
        "--level_steps":10,
        "--d":2,
        "--cfl":0.5,
        "--delta_t":0.01,
        "--max_iter":500,
        "--level_alpha":30,
        "--g_force_lambda":-1,
        "--restart_num":10,
        "--tol":1e-6,
        "--output_image":"None",
        "--iprint":10
        }
    return options

def _setopt(options):
    options.pop('-f',1)
    default =default_options()
    for k in default:
        if isinstance(default[k],(float,int)): 
            if isinstance(options[k],str):
                options[k]=eval(options[k])
            if isinstance(default[k], float):
                options[k]=float(options[k])
            else:
                options[k]=int(options[k])
        else:
            options[k]=type(default[k])(options[k])
    return None
def setoptions(*,argv=None,kw=None):
    options=default_options()
    longopts=list(k[2:]+'=' for k in options)
    argv =({} if argv is None else
           dict(getopt.getopt(argv,shortopts='f',longopts=longopts)[0]))
    kw=({} if kw is None else kw)
    
    options.update(kw)
    options.update(argv)
    
    _setopt(options)
    
    globalnames={}
    for k,v in options.items():
        globalnames[k[2:]]=v
    return options,globalnames