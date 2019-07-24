import numpy as np 

param = {}

# time
param['t0'] = 0.
param['tf'] = 20.
param['dt'] = 0.05
T = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['T'] = np.reshape(T, (len(T),-1))
param['nt'] = len(param.get('T'))

# agents
param['na'] = 2
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

# controller
# [ 'empty', 'fdbk', 'clf', 'scp-clf']
param['controllers'] = [ 'fdbk','clf','mpc']
param['mpc_horizon'] = 30 # in timesteps
param['mpc_update'] = 5 
param['control_max'] = 1
param['k_fdbk'] = 100

# desired trajectory
param['case_xd'] = 2

# solver
param['max_iters'] = 2000
param['n_scp_iter'] = 5
param['scp_tol'] = 0.1
param['p_u'] = 0 
param['p_v'] = 1000
param['tau_trust'] = 0.1

# init
param['plim'] = 1 
param['vlim'] = 1 
param['min_dist'] = 1 

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