import os
import numpy as np
import lumapi

np.random.seed(seed = 98765)

from lumopt.utilities.wavelengths import Wavelengths
from lumopt.utilities.materials import Material
from lumopt.geometries.polygon import FunctionDefinedPolygon
from lumopt.figures_of_merit.modematch import ModeMatch
from lumopt.optimizers.generic_optimizers import ScipyOptimizers
from lumopt.optimization import Optimization

######## SPECTRAL RANGE #########
wavelengths = Wavelengths(start = 1530.0e-9, stop = 1570.0e-9, points = 25)

######## OPTIMIZABLE GEOMETRY ########
n_grates = 20
wg_height = 220.0e-9
wg_length = 30.0e-6
etch_depth = 0.8
x0 = -6.0e-6
y0 = 0.0
def grate_function(params):
    y2 = y0+wg_height
    y1 = y2-etch_depth*wg_height
    x_begin = x0-wg_length
    verts = np.array( [ [x_begin,y0], [x_begin,y2], [x0,y2], [x0,y1] ] )
    xp = float(x0)
    for idx in range(n_grates):
        a = params[2*idx]*1e-6
        b = params[2*idx+1]*1e-6
        pitch = a+b
        verts = np.concatenate((verts, [[xp+a, y1], [xp+a, y2], [xp+pitch, y2], [xp+pitch, y1]] ), axis = 0)
        xp += pitch
    verts = np.concatenate((verts, [[xp, y1], [xp, y0]]), axis = 0)
    return verts
initial_params = np.zeros(2*n_grates)
for idx in range(n_grates):
    initial_params[2*idx] = 0.1+0.2*np.sin(np.pi/2.0*idx/n_grates)*np.random.random()
    initial_params[2*idx+1] = 0.7-initial_params[2*idx]
bounds = [(0.1, 0.9)] * (2*n_grates)
Si = Material(base_epsilon = 3.47668**2, mesh_order = 2)
geometry = FunctionDefinedPolygon(func = grate_function, 
                                  initial_params = initial_params,
                                  bounds = bounds,
                                  z = 0.0,
                                  depth = wg_height,
                                  eps_out = 1.0 ** 2,
                                  eps_in = Si,
                                  edge_precision = 5,
                                  dx = 1.0e-5)

######## FIGURE OF MERIT ########
fom = ModeMatch(monitor_name = 'fom',
                mode_number = 'fundamental TM mode',
                direction = 'Backward',
                target_T_fwd = lambda wl: 0.5*np.ones(wl.size),
                norm_p = 1,
                target_fom = 0.0)

######## OPTIMIZATION ALGORITHM ########
optimizer = ScipyOptimizers(max_iter = 200,
                            method = 'L-BFGS-B',
                            scaling_factor = 1.0,
                            pgtol = 1.0e-5,
                            ftol = 1.0e-5,
                            scale_initial_gradient_to = 0.0)

######## BASE 2-D SIMULATION ########
proj_path = os.path.join(os.path.dirname(__file__))
base_sim = os.path.join(proj_path, 'grating_coupler_base.fsp')

######## PUT 2-D OPTIMIZATION TOGETHER ########
opt = Optimization(base_script = base_sim,
                   wavelengths = wavelengths,
                   fom = fom,
                   geometry = geometry,
                   optimizer = optimizer,
                   use_var_fdtd = False,
                   hide_fdtd_cad = False,
                   use_deps = True,
                   plot_history = True,
                   store_all_simulations = True,
                   save_global_index = False,
                   label = None)

######## RUN THE 2-D OPTIMIZATION ########
res = opt.run()

######## 2-D SIMULATION WITH OPTIMIZED STRUCTURE ########
final_sim = os.path.join(proj_path, 'optimized_grating_coupler_etch_0p8_bandwidth_40nm.fsp')
with lumapi.FDTD() as fdtd:
    fdtd.load(base_sim)
    fdtd.addpoly()
    fdtd.setnamed('polygon', 'vertices', grate_function(res[1]))
    fdtd.setnamed('polygon', 'material', Si.name)
    fdtd.setnamed('polygon', 'index', np.sqrt(Si.base_epsilon))
    fdtd.setnamed('polygon', 'x', x0)
    fdtd.setnamed('polygon', 'y', y0)
    fdtd.setnamed('back', 'enabled', True)