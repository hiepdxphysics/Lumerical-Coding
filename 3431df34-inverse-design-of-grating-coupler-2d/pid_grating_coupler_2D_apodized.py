##"""
 ##   Copyright (c) 2020 Ansys Inc. """

######## IMPORTS ########
# General purpose imports
import os,sys
import numpy as np
import scipy as sp
import json
sys.path.append("C:\\Program Files\\Lumerical\\v222\\api\\python\\")
from lumjson import LumEncoder, LumDecoder


# Optimization specific imports
from lumopt.utilities.load_lumerical_scripts import load_from_lsf
from lumopt.utilities.wavelengths import Wavelengths
from lumopt.geometries.parameterized_geometry import ParameterizedGeometry
from lumopt.geometries.polygon import FunctionDefinedPolygon
from lumopt.figures_of_merit.modematch import ModeMatch
from lumopt.optimizers.generic_optimizers import ScipyOptimizers
from lumopt.optimization import Optimization
from lumopt.utilities.materials import Material

import lumapi

cur_path = os.path.dirname(os.path.realpath(__file__))


# Optimization global parameters
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

params_file = "pid_grating_coupler_initial_params.json"
base_sim_2d = "pid_grating_coupler_2D_TE_base.fsp"
base_sim_apodized_2d = "pid_grating_coupler_2D_TE_base_apodized.fsp"
base_script_2d = 'pid_grating_coupler_2D_TE_base.lsf'
params_file_apod = "pid_optim_1.json"



def grating_params_pos(params):
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
    
    
def get_vertices_from_distances(self, x_start, distances):
        #x_start = params[0]*1e-6  #< First parameter is the starting position
        #R  = params[1]*1e6        #< second parameter (unit is 1/um)
        #a  = params[2]            #< Third parameter (dim-less)
        #b  = params[3]            #< Fourth parameter (dim-less)
        
        x_begin = self.x_min-2e-6
        x = np.cumsum(np.concatenate(([x_start],distances)))        

        y0 = 0
        y3 = y0+self.wg_height
        y1 = y3-self.etch_depth
    
        verts = np.array( [[x_begin,y0],[x_begin,y3],[x[0],y3]] )       
    
        ## Iterate over all but the last tooth
        for i in range(1,len(x),2):
            verts = np.concatenate((verts,np.array([[x[i-1],y1],[x[i],y1],[x[i],y3],[x[i+1],y3]])),axis=0)
    
        ## Close off the polygon
        verts = np.concatenate((verts,np.array([[x[-1],y0]])),axis=0) 
        return verts


if __name__ == "__main__":   
    with open(os.path.join(cur_path, params_file)) as fh:
        initial_params = json.load(fh, cls=LumDecoder)["initial_params"][0]
        
    with lumapi.FDTD(filename = os.path.join(cur_path, base_sim_2d), hide = True) as fdtd:
        vtx = grating_params_pos(initial_params)
        fdtd.addpoly()
        fdtd.set("vertices", vtx)
        fdtd.set("x", 0)
        fdtd.set("y", 0)
        fdtd.set("index", indexSi)
        fdtd.save(os.path.join(cur_path, base_sim_apodized_2d))
        # extracting the distsnces from vertices
        vx=np.array(vtx)
        vertsx=vx[2:(len(vx)-2):2,0]
        vert_x = np.concatenate(([np.array(vx[2,0])],(vertsx[1:]-vertsx[0:len(vertsx)-1]))) 
        vert_x = vert_x *1e6
               
    with open(os.path.join(cur_path, params_file_apod), "w") as fh:
        json.dump({ "initial_params": vert_x }, fh, cls=LumEncoder, indent = 4)
        
        
        
        
        
        
        
        
         
        
       
