

import autograd.numpy as np 
import dynamics
import utilities as util 
from param import param
import cvxpy as cp
import plotter

def get_u(x,t):
	# input

	print(param.get('controller'))

	if param.get('controller') is 'empty':
		u = np.zeros( [param.get('m'),1])

	elif param.get('controller') is 'fdbk': 
		u = get_fdbk_controller( x,t)

	elif param.get('controller') is 'clf':
		u = get_clf_controller( x,t)

	elif param.get('controller') is 'scp':
		start_idx = np.where( param.get('T') == t)[0][0]
		end_idx = param.get('nt')-1
		T = param.get('T')[start_idx:end_idx]
		u = get_scp_clf_controller( x, T)

	elif param.get('controller') is 'mpc':
		start_idx = np.where( param.get('T') == t)[0][0]
		end_idx = np.min( (start_idx + param.get('mpc_horizon'), param.get('nt')-1))
		T = param.get('T')[start_idx:end_idx]
		u = get_scp_clf_controller( x, T)

	return u 

def get_fdbk_controller( x,t):

	k = np.where(param.get('T') == t)[0][0]
	my_1 = util.get_my_1()
	eta = dynamics.get_eta( x,t)
	dvdota_dvb = dynamics.get_dvdota_dvb(x,t)
	xtilde_dot = dynamics.get_xtildedot( x,t)
	dvdota_dxtilde = dynamics.get_dvdota_dxtilde( x,t)

	A = np.matmul( np.transpose( my_1), \
		np.matmul( dvdota_dxtilde, xtilde_dot)) - \
		util.list_to_vec( param.get('ad')[k,:])
	B = np.matmul( np.transpose( my_1), dvdota_dvb)
	K = param.get('k_fdbk')*np.kron( np.eye(param.get('nd')), \
		np.ones((1,param.get('gamma'))))
	u = np.matmul(np.linalg.pinv(B), - A - np.matmul(K,eta))

	u = np.clip( u, -param.get('control_max'), param.get('control_max'))
	return u

def get_clf_controller( x,t):

	LgV = dynamics.get_LgV( x,t) 
	LfV = dynamics.get_LfV( x,t)
	V   = dynamics.get_V( x,t)
	lambda_v = util.get_stabilization_rate()

	# cvxpy
	u = cp.Variable(param.get('m'))
	delta_v = cp.Variable()
	constraints = [ 
		LfV + LgV*u + lambda_v * V <= delta_v,
		cp.abs( u)  <= param.get('control_max'),
		delta_v >= 0
	]
	obj = cp.Minimize( \
		cp.sum_squares(u) + \
		param.get('p_v') * delta_v \
		)
	
	prob = cp.Problem(obj, constraints)
	prob.solve(verbose = False, solver = cp.GUROBI)

	u = np.array(u.value)
	return u

def get_scp_clf_controller( x0,T):

	lambda_v = util.get_stabilization_rate()
	xbar, ubar = get_scp_initial_trajectory(x0, T)
	xbar = np.transpose( xbar, (1,0,2))
	ubar = np.transpose( ubar, (1,0,2))

	i_iter = 0
	cost_curr = np.inf
	cost_diff = np.inf
	tx_curr = xbar
	state_diff = np.inf
	C = []

	while i_iter < param.get('n_scp_iter') and \
		np.abs(cost_diff) > param.get('scp_cost_tol'):

		if not (i_iter == 0):
			if param.get('nl_correction_on'):
				xbar = []
				x_curr = x0
				for k,t in enumerate(T):
					u_curr = np.reshape( u_next[:,k], (param.get('m'),1))
					x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')
					xbar.append(x_curr)
					x_curr = x_next
				xbar = np.transpose( np.asarray(xbar), (1,0,2))
				ubar = u_next
			else:
				ubar = u_next
				xbar = tx_next


		print('Iteration: ' + str(i_iter) + '/' + str(param.get('n_scp_iter')))
	
		u  = cp.Variable( param.get('m'), len(T) )
		tx = cp.Variable( param.get('n'), len(T) )
		tV = cp.Variable( len(T), 1)
		delta_v = cp.Variable( len(T)-1, 1)

		constraints = []
		constraints.append( tx[:,0] == xbar[:,0])

		for k in range( len(T)):

			F_k, B_k, d_k = dynamics.get_linear_dynamics( xbar[:,k], ubar[:,k], T[k])
			R_k, w_k = dynamics.get_linear_lyapunov( xbar[:,k], ubar[:,k], T[k])

			constraints.append( 
				tV[k] == R_k * tx[:,k] + w_k)
			constraints.append( 
				cp.abs( tx[:,k] - xbar[:,k]) <= param.get('tau_x') * np.power( 0.95, i_iter))
			constraints.append( 
				cp.abs( u[:,k] - ubar[:,k]) <= param.get('tau_u') * np.power( 0.95, i_iter))
			constraints.append( 
				cp.abs( u[:,k]) <= param.get('control_max'))
		
			if k < len(T)-1:
				constraints.append( 
					tV[k+1] <= (1-lambda_v*param.get('dt')) * tV[k] + delta_v[k]*param.get('dt'))
				constraints.append(
					tx[:,k+1] == F_k * tx[:,k] + B_k * u[:,k] + d_k)
				constraints.append(
					delta_v[k] >= 0)


		# obj = cp.Minimize( cp.sum_squares(u) + param.get('p_v')*sum(delta_v)) 
		# obj = cp.Minimize( param.get('p_v')*sum(delta_v)) 
		obj = cp.Minimize( sum(tV) + cp.sum_squares(u)) 
		prob = cp.Problem( obj, constraints)

		prob.solve(verbose = False, solver = cp.GUROBI, \
			BarQCPConvTol = 1e-6, BarHomogeneous = True)
		# prob.solve(verbose = False)

		tx_next = np.asarray(tx.value).reshape(param.get('n'), len(T), -1)
		u_next = np.asarray(u.value).reshape(param.get('m'), len(T), -1)
		tV_next = np.asarray(tV.value).reshape(len(T), -1)

		cost_next = prob.value
		cost_diff = cost_next - cost_curr
		state_diff = sum( np.linalg.norm( tx_next - tx_curr, axis = 1))

		i_iter += 1
		C.append( cost_curr)
		cost_curr = cost_next
		tx_curr = tx_next
		print('Curr Cost - Prev Cost: ' + str(cost_diff))
		print('State Change: ' + str(state_diff))
		print('Timestep: ' + str(T[0]) + '/' + str(param.get('T')[-1]) + '\n')

	plotter.debug_scp_iteration_plot( tx_next, u_next, xbar, ubar, x0, T, i_iter)
	plotter.plot_cost_iterations(C)
	U = np.transpose( ubar, (1,0,2))
	return U


def get_scp_initial_trajectory(x0, T):
	# use clf solution

	X = []
	U = []
	x_curr = x0

	for t in T:
		u_curr = get_clf_controller( x_curr, t) 
		x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')

		X.append(x_curr)
		U.append(u_curr)
		x_curr = x_next

	X = np.asarray(X)
	U = np.asarray(U)
	return X,U


def calc_cost( U):

	lambda_v = util.get_stabilization_rate()
	cost = 0
	x_curr = param.get('x0')
	for k,t in enumerate( param.get('T')):
		if k < param.get('nt')-1:
			u_curr = np.reshape( U[k,:], (param.get('m'),1))
			x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')

			v_curr = dynamics.get_V(x_curr, t)
			v_next = dynamics.get_V(x_next, param.get('T')[k+1])
			delta_v = np.max((0, \
				(v_next - v_curr)/param.get('dt') + lambda_v * v_curr))

			cost += np.linalg.norm( u_curr) + param.get('p_v')*delta_v
			x_curr = x_next
	return cost
