

import autograd.numpy as np 

param = {}

# time
param['t0'] = 0.
param['tf'] = 20
param['dt'] = 0.05
T = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['T'] = np.reshape(T, (len(T),-1))
param['nt'] = len(param.get('T'))

# agents
param['na'] = 1
param['nb'] = 1
param['ni'] = param.get('na') + param.get('nb')

# dimension
param['nd'] = 2
param['n'] = 2 * param.get('nd') * param.get('ni')
param['m'] = param.get('nd') * param.get('nb')

# dynamics
param['model'] = 'reynolds'
param['kx'] = 1.
param['kv'] = 1.
param['R_des'] = 1.
param['gamma'] = 3
param['lambda_a'] = 0.8
param['R_comm'] = 5.

# controller
param['controllers'] = ['clf'] #[ 'fdbk','clf','scp','mpc']
param['mpc_horizon'] = 30
param['mpc_update'] = 5
param['control_max'] = 1
param['k_fdbk'] = 10

# objective
param['objective'] = 'tracking'
param['case_xd'] = 2

# solver
param['max_iters'] = 2000
param['n_scp_iter'] = 5
param['scp_cost_tol'] = 0.01
param['scp_state_tol'] = 0.5
param['p_u'] = 1000 
param['p_v'] = 1000
param['p_tau'] = 1000
param['tau_x'] = 0.05
param['tau_u'] = 0.05
param['nl_correction_on'] = True

# init
param['plim'] = 2
param['vlim'] = 1 
param['min_dist'] = 0.1

# plotting
param['ControlAgentColor'] = 'r'
param['FreeAgentColor'] = 'b'
param['DesiredTrajectoryColor'] = 'g'
param['CentroidColor'] = 'k'
param['start_marker'] = 's'
param['stop_marker'] = 'x'
param['plot_buffer'] = '0.2'
param['gif_on'] = False
param['nframes_gif'] = '50'
param['fn_gif'] = 'gif.gif'
param['fn_plots'] = 'plots.pdf'