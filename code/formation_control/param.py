import numpy as np 

param = {}

# time
param['t0'] = 0.
param['tf'] = 20.
param['dt'] = 0.1
param['T'] = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['nt'] = len(param.get('T'))

# number of agents
param['na'] = 3
param['nb'] = 4
param['ni'] = param.get('na') + param.get('nb')

# number of dimensions
param['nd'] = 2 # spatial dimension (ie x,y in R^2)
param['n'] = 2*param.get('nd')*param.get('ni')
param['m'] = 2*param.get('nd')*param.get('nb')

# dynamics
param['kx'] = 1.
param['kv'] = 1.
param['R_des'] = 1.5
param['lambda_a'] = 1.
param['R_comm'] = 2.

# controller
param['mpc_horizon'] = 10
param['mpc_update'] = 5
param['controller'] = 'mpc'
param['control_max'] = 1.

# solver
param['tau_x'] = 3
param['tau_u'] = 3
param['scp_cost_tol'] = 0.1
param['n_iter_scp'] = 5

# desired formation
param['case_xd'] = '3_node_triangle'
#: '3_node_triangle'
#: 'single_node'

# init
param['plim'] = 2
param['vlim'] = 1 
param['min_dist'] = 0.1

# plotting
param['fn_plots'] = 'plots.pdf'
param['start_marker'] = 's'
param['stop_marker'] = 'x'
param['ControlAgentColor'] = 'r'
param['FreeAgentColor'] = 'b'