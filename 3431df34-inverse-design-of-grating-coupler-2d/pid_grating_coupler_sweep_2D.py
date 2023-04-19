"""
    Copyright (c) 2020 Ansys Inc. """
    
######## IMPORTS ########
# General purpose imports
import os,sys
sys.path.append("C:\\Program Files\\Lumerical\\v222\\api\\python\\")
import scipy as sp
import numpy as np
import json
from lumjson import LumEncoder, LumDecoder

import lumapi

######## OPTIMIZABLE GEOMETRY ########
lambda_c = 1.55e-6 
bandwidth_in_nm = 0     #< Only optimize for center frequency of 1550nm
F0 = 0.95
height = 220e-9
etch_depth = 80e-9
y0 = 0
x_begin = -5.1e-6
x_end = 22e-6
n_grates = 25

indexSi = 3.47668
indexSiO2 = 1.44401

data_file = "pid_grating_coupler_initial_params.json"
base_file = "pid_grating_coupler_2D_TE_base.fsp"


def grating_params_pos(params, n_grates):
    y0      = 0

    y3 = y0+height
    y1 = y3-etch_depth

    x_start = params[0]*1e-6  #< First parameter is the starting position
    R  = params[1]*1e6        #< second parameter (unit is 1/um)
    a  = params[2]            #< Third parameter (dim-less)
    b  = params[3]            #< Fourth parameter (dim-less)

    x0 = x_start
  
    verts = np.array( [[x_begin,y0],[x_begin,y3],[x0,y3],[x0,y1]] )
    
    ## Iterate over all but the last tooth
    for i in range(n_grates-1):
        F = F0-R*(x0-x_start)
        Lambda = lambda_c / (a+F*b)
        x1 = x0 + (1-F)*Lambda    #< Width of the etched region
        x2 = x0 + Lambda          #< Rest of cell
        verts = np.concatenate((verts,np.array([[x1,y1],[x1,y3],[x2,y3],[x2,y1]])),axis=0)
        x0 = x2

    ## Last tooth is special
    F = F0-R*(x0-x_start)
    Lambda = lambda_c / (a+F*b)
    x1 = x0 + (1-F)*Lambda        #< Width of the etched region
    verts = np.concatenate((verts,np.array([[x1,y1],[x1,y3],[x_end,y3],[x_end,y0]])),axis=0) 

    return verts


if __name__ == "__main__":
    with open(data_file) as fh:
        initial_params = json.load(fh, cls=LumDecoder)["initial_params"][0]
                
    # Alternate starting point
    # initial_params = [ -2.5, 0.03, 2.4, 0.5369]

    vtx = grating_params_pos(initial_params, 25)

    if os.path.exists(base_file):
        with lumapi.FDTD(filename=base_file) as fdtd:
            fdtd.addpoly()
            fdtd.set("vertices", vtx)
            fdtd.set("x", 0)
            fdtd.set("y", 0)
            fdtd.set("index", indexSi)
            fdtd.setglobalsource("center wavelength", lambda_c)
            fdtd.setglobalsource("wavelength span", 0)
            fdtd.save()
            fdtd.runsweep("sweep source position")
            sweep_pos = fdtd.getsweepdata("sweep source position", "x")
            sweep_T = fdtd.getsweepdata("sweep source position", "T")
            
            Tmax = np.amax(sweep_T)
            Tpos = sweep_pos[np.where(sweep_T == np.amax(sweep_T))[0][0]][0]
            
            print("Max transmission:", Tmax*100, "%")
            print("Position", Tpos*1e6, "um")
            
            fdtd.setnamed("source", "x", Tpos)
            
            fdtd.select("polygon")
            fdtd.delete()
            fdtd.save()            
    else:
        print("base file doesn't exist...")
