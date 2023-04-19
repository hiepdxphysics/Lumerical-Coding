"""
    Copyright (c) 2020 Ansys Inc. """

######## IMPORTS ########
# General purpose imports
import os,sys
sys.path.append("C:\\Program Files\\Lumerical\\v212\\api\\python\\")
import numpy as np
import json
from lumjson import LumEncoder, LumDecoder

import lumapi

# Optimization global parameters
height = 220e-9
etch_depth = 80e-9
y0 = 0
n_grates = 25

params_size = 2*n_grates

indexSi = 3.47668
indexSiO2 = 1.44401

base_sim_3D = "pid_grating_coupler_3D_TE_base.fsp"
params_file = "pid_optim_final.json"
sim_3D_file = "pid_grating_coupler_3D.fsp"
gds_file = "pid_grating_coupler_3d.gds"

cur_path = os.path.dirname(os.path.realpath(__file__))


def focusing_grating(params, fdtd, update_only = False, rotation_angle = 0.0, in_taper_len = 25.0e-6, out_taper_len = 5.0e-6):
    assert params.size == params_size
    group_name = 'grating'
    Si_layer_name = 'silicon'
    theta = np.arctan2(12.0e-6, in_taper_len+1e-6*params.sum())*(180.0/np.pi)
    if not update_only:
        fdtd.addstructuregroup()
        fdtd.set('name', group_name)
        fdtd.setnamed(group_name, 'x', 0.0)
        fdtd.setnamed(group_name, 'y', 0.0)
        fdtd.setnamed(group_name, 'z', 0.0)
        fdtd.setnamed(group_name, 'construction group', False)
        fdtd.addring()
        fdtd.set('name', Si_layer_name)
        fdtd.setnamed(Si_layer_name, 'index', indexSi)
        fdtd.setnamed(Si_layer_name, 'override mesh order from material database', True)
        fdtd.setnamed(Si_layer_name, 'mesh order', 2)
        fdtd.setnamed(Si_layer_name, 'x', 0.0)
        fdtd.setnamed(Si_layer_name, 'y', 0.0)
        fdtd.setnamed(Si_layer_name, 'z', 0.0)
        fdtd.setnamed(Si_layer_name, 'z span', height)
        fdtd.setnamed(Si_layer_name, 'inner radius', 0.0)
        fdtd.setnamed(Si_layer_name, 'theta start', -theta)
        fdtd.setnamed(Si_layer_name, 'theta stop', theta)
        fdtd.addtogroup(group_name)
    Si_layer_name = group_name + '::' + Si_layer_name
    inner_radius = float(in_taper_len) + params[0]*1e-6
    for idx in range(n_grates):
        etch_len = params[2*idx+1]*1e-6
        try:
            tooth_len = params[2*idx+2]*1e-6
        except:
            tooth_len = out_taper_len
        # etch_len = (params[1]+params[2]*idx/(num_grates-1))*1e-6
        # tooth_len = params[3]*1e-6 - etch_len
        etch_name = 'etch_{}'.format(idx)
        if not update_only:
            fdtd.addring()
            fdtd.set('name', etch_name)
            fdtd.setnamed(etch_name, 'material', '<Object defined dielectric>')
            fdtd.setnamed(etch_name, 'index', indexSiO2)
            fdtd.setnamed(etch_name, 'override mesh order from material database', True)
            fdtd.setnamed(etch_name, 'mesh order', 1)
            fdtd.setnamed(etch_name, 'x', 0.0)
            fdtd.setnamed(etch_name, 'y', 0.0)
            fdtd.setnamed(etch_name, 'z max', height/2.0)
            fdtd.setnamed(etch_name, 'z min', height/2.0 - etch_depth)
            fdtd.setnamed(etch_name, 'theta start', -theta)
            fdtd.setnamed(etch_name, 'theta stop', theta)
            fdtd.addtogroup(group_name)
        etch_name = group_name + '::' + etch_name
        fdtd.setnamed(etch_name, 'inner radius', inner_radius)
        inner_radius += etch_len
        fdtd.setnamed(etch_name, 'outer radius', inner_radius)
        inner_radius += tooth_len
    inner_radius += out_taper_len
    fdtd.setnamed(Si_layer_name, 'outer radius', inner_radius)
    if not update_only:
        x0, y0, z0 = -in_taper_len, 0.0, height/2.0
        rotation_angle_in_rad = (np.pi/180.0)*rotation_angle
        y0_rot = y0*np.cos(rotation_angle_in_rad)-z0*np.sin(rotation_angle_in_rad)
        z0_rot = y0*np.sin(rotation_angle_in_rad)+z0*np.cos(rotation_angle_in_rad)
        fdtd.setnamed(group_name, 'x', x0)
        fdtd.setnamed(group_name, 'y', y0_rot)
        fdtd.setnamed(group_name, 'z', z0_rot)
        fdtd.setnamed(group_name, 'first axis', 'x')
        fdtd.setnamed(group_name, 'rotation 1', rotation_angle)

def gds_export_script(fdtd, gdsfile):
    fdtd.eval("gds_filename = '{0}';".format(gdsfile))
    fdtd.eval("top_cell = 'model';")
    fdtd.eval("layer_def = [1, 0.0, {0}; 2 , {1}, {0}];".format(height, height-etch_depth))
    fdtd.eval('n_circle = 64;')
    fdtd.eval('n_ring = 64;')
    fdtd.eval('n_custom = 64;')
    fdtd.eval('n_wg = 64;')
    fdtd.eval('round_to_nm = 1;')
    fdtd.eval('grid = 1e-9;')
    fdtd.eval('max_objects = 10000;')
    fdtd.eval('Lumerical_GDS_auto_export;')

if __name__ == "__main__":    
    with open(params_file) as fh:
        initial_params = json.load(fh, cls=LumDecoder)["initial_params"] 
   
    with lumapi.FDTD(filename = os.path.join(cur_path, base_sim_3D), hide = False) as fdtd:
        focusing_grating(initial_params, fdtd)
        fdtd.save(os.path.join(cur_path, sim_3D_file))
        gds_export_script(fdtd, gds_file)