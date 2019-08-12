
import autograd.numpy as np 
import dynamics
import plotter
import utilities as util 
from param import param
import cvxpy as cp

def get_u(x,t):
	# input

	print(param.get('controller'))

	if param.get('controller') is 'empty':
		u = np.zeros( [param.get('m'),1])

	elif param.get('controller') is 'mpc':
		start_idx = np.where( param.get('T') == t)[0][0]
		end_idx = np.min( (start_idx + param.get('mpc_horizon'), param.get('nt')-1))
		T = param.get('T')[start_idx:end_idx]
		u = get_mpc_controller( x, T)

	return u 

def get_mpc_controller(x0,T):

	xbar, ubar = get_initial_trajectory(x0,T)
	Vs = util.augment(util.get_Vs_a(param.get('xd')))
	VsT_T_pi_a = np.dot( np.dot( np.dot( 
		np.transpose(Vs), util.get_T_a()), util.get_pi_a()), util.get_pv_a() )

	i_iter = 0 
	cost_curr = np.inf
	cost_diff = np.inf
	while i_iter < param.get('n_iter_scp') and \
		np.abs(cost_diff) > param.get('scp_cost_tol'):

		print('Iteration: ' + str(i_iter))

		u = cp.Variable(param.get('m'),len(T))
		x = cp.Variable(param.get('n'),len(T))
		delta = cp.Variable(2*param.get('nd')*(param.get('na')-1),len(T)) 

		constraints = []
		constraints.append(x[:,0] == x0)

		for k,t in enumerate(T):

			# synchronization
			constraints.append(
				VsT_T_pi_a * x[:,k] == delta[:,k])

			# control bounds
			constraints.append( 
				cp.abs( u[:,k]) <= param.get('control_max'))

			# trust region
			# constraints.append( 
			# 	cp.abs( x[:,k] - xbar[:,k]) <= param.get('tau_x'))			
			# constraints.append( 
			# 	cp.abs( u[:,k] - ubar[:,k]) <= param.get('tau_u'))

			# dynamics
			if k < len(T)-1:

				F_k, B_k, d_k = dynamics.get_linear_dynamics( xbar[:,k], ubar[:,k], T[k])
				constraints.append(
					x[:,k+1] == F_k*x[:,k] + B_k*u[:,k] + d_k)

		obj = cp.Minimize( 100*cp.sum_squares(delta) + cp.sum_squares(u))
		prob = cp.Problem( obj, constraints)

		# prob.solve( verbose = True, solver = cp.GUROBI) 
		prob.solve( verbose = True) 

		plotter.plot_VsTx( np.asarray(delta.value), T, title = 'VsTx')
		fig, ax = plotter.make_fig()
		plotter.plot_SS(xbar, T, 
			title = 't: ' + str(T[0]) + '\n iter: ' + str(i_iter), fig = fig, ax = ax)
		plotter.plot_SS(np.asarray(x.value), T, 
			title = 't: ' + str(T[0]) + '\n iter: ' + str(i_iter), fig = fig, ax = ax)

		ubar = np.asarray(u.value)
		xbar = np.asarray(x.value)

		cost_next = prob.value
		cost_diff = cost_next - cost_curr
		print('Current Cost:' + str(cost_next))
		print('Cost Difference: ' + str(cost_diff))

		i_iter += 1
		cost_curr = cost_next

	param['xbar'] = xbar
	param['ubar'] = ubar
	param['delta'] = np.asarray(delta.value)
	return ubar


def get_initial_trajectory(x0,T):


	if param.get('ubar') is not None:
		U = np.hstack((
			param.get('ubar')[:,param.get('mpc_update'):],
			np.zeros((param.get('m'), param.get('mpc_update'))) ))

	else:
		U = np.zeros( (param.get('m'), param.get('nt')))
	
	# forward propagate
	x_curr = x0
	for k,t in enumerate(T):
		u_curr = U[:,k]
		x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, T[k])*param.get('dt')
		
		try:
			X = np.hstack( (X, x_curr))
		except:
			X = x_curr

	return X,U
