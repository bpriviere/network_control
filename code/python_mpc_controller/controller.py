

import numpy as np 
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
		# u = get_mpc_scp_clf_controller( param.get('x0'), t)
		u = get_scp_clf_controller()

	elif param.get('controller') is 'mpc':
		u = get_mpc_scp_clf_controller( x, t)

	return u 

def get_fdbk_controller( x,t):
	k 			   = np.where(param.get('T') == t)[0][0]
	my_1 		   = util.get_my_1()
	eta 		   = dynamics.get_eta( x,t)
	dvdota_dvb 	   = dynamics.get_dvdota_dvb(x,t)
	xtilde_dot	   = dynamics.get_xtildedot( x,t)
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

	u = util.list_to_vec([x for x in u.value])
	return u

def get_scp_clf_controller():

	lambda_v = util.get_stabilization_rate()
	T = param.get('T')
	p_v = param.get('p_v')
	control_max = param.get('control_max')
	tau_trust = param.get('tau_trust')
	dt = param.get('dt')
	max_iters = param.get('max_iters')

	i_iter = 0
	cost_curr = np.inf
	cost_diff = np.inf
	while i_iter < param.get('n_scp_iter') and cost_diff > param.get('scp_tol'):
	
		print('SCP Iteration: ' + str(i_iter) + '/' + str(param.get('n_scp_iter')))

		if i_iter == 0:
			xbar, ubar = get_scp_initial_trajectory(param.get('x0'), param.get('T'))
		else:
			xbar = np.transpose( np.asarray([y for y in tx.value]), (0,2,1))
			ubar = np.transpose( np.asarray([y for y in u.value]),  (0,2,1))

		if np.mod( i_iter, 1) == 0:
			plotter.plot_SS( xbar, param.get('T'), title = 'SCP Iteration: ' + str(i_iter))

		u  = cp.Variable( len(T), param.get('m') )
		tx = cp.Variable( len(T), param.get('n') )
		tV = cp.Variable( len(T),1)
		delta_v = cp.Variable( len(T)-1,1 )

		constraints = []
		constraints.append( tx[0,:].T == param.get('x0'))

		for k in range( len(T)):
			if k < len(T)-1:
				F_k, B_k, d_k = dynamics.get_linear_dynamics( xbar[k], ubar[k], T[k])

				constraints.append( 
					tV[k+1] <= (1-lambda_v*dt) * tV[k] + delta_v[k]*dt)

				constraints.append(
					tx[k+1,:].T == F_k * tx[k,:].T + B_k * u[k,:].T + d_k)

				constraints.append(
					delta_v[k] >= 0)

			R_k, w_k = dynamics.get_linear_lyapunov( xbar[k], ubar[k], T[k])

			constraints.append( 
				tV[k] == R_k * tx[k,:].T + w_k)

			constraints.append( 
				cp.abs(tx[k,:].T-xbar[k]) <= tau_trust)

			constraints.append( 
				cp.abs(u[k,:]) <= control_max)

		obj = cp.Minimize( p_v*sum(delta_v) + cp.sum_squares(u))
		prob = cp.Problem( obj, constraints)

		prob.solve(verbose = False, solver = cp.GUROBI,  max_iters = max_iters)

		i_iter += 1
		cost_next = prob.value
		cost_diff = cost_next - cost_curr
		cost_curr = cost_next

		print(prob.value)
		print(cost_diff)

	X = np.squeeze(np.asarray([y for y in tx.value]))
	U = np.asarray([y for y in u.value])

	# print(np.shape(U))
	return U

def get_mpc_scp_clf_controller( x0,t):
	
	lambda_v = util.get_stabilization_rate()
	p_v = param.get('p_v')
	control_max = param.get('control_max')
	tau_trust = param.get('tau_trust')
	dt = param.get('dt')
	max_iters = param.get('max_iters')

	start_idx = np.where( param.get('T') == t)[0][0]
	end_idx = np.min( (start_idx + param.get('mpc_horizon'), param.get('nt')-1))
	T = param.get('T')[start_idx:end_idx]

	i_iter = 0
	cost_curr = np.inf
	cost_diff = np.inf
	if len(T) > 1:
		while i_iter < param.get('n_scp_iter') and np.abs(cost_diff) > param.get('scp_tol'):
		
			print('SCP Iteration: ' + str(i_iter) + '/' + str(param.get('n_scp_iter')))

			if i_iter == 0:
				xbar, ubar = get_scp_initial_trajectory(x0, T)
			else:
				xbar = np.transpose( np.asarray([y for y in tx.value]), (0,2,1))
				ubar = np.transpose( np.asarray([y for y in u.value]),  (0,2,1))

			# if np.mod( i_iter, 1) == 0:
			# 	plotter.plot_SS( xbar, T, title = 'SCP Iteration: ' + str(i_iter))

			u  = cp.Variable( len(T), param.get('m') )
			tx = cp.Variable( len(T), param.get('n') )
			tV = cp.Variable( len(T), 1)
			delta_v = cp.Variable( len(T)-1, 1)

			constraints = []
			constraints.append( tx[0,:].T == x0)

			for k in range( len(T)):
				if k < len(T)-1:
					F_k, B_k, d_k = dynamics.get_linear_dynamics( xbar[k], ubar[k], T[k])

					constraints.append( 
						tV[k+1] <= (1-lambda_v*dt) * tV[k] + delta_v[k]*dt)

					constraints.append(
						tx[k+1,:].T == F_k * tx[k,:].T + B_k * u[k,:].T + d_k)

					constraints.append(
						delta_v[k] >= 0)

				R_k, w_k = dynamics.get_linear_lyapunov( xbar[k], ubar[k], T[k])

				constraints.append( 
					tV[k] == R_k * tx[k,:].T + w_k)

				constraints.append( 
					cp.abs(tx[k,:].T-xbar[k]) <= tau_trust)

				constraints.append( 
					cp.abs(u[k,:]) <= control_max)

			obj = cp.Minimize( p_v*sum(delta_v) + cp.sum_squares(u))
			prob = cp.Problem( obj, constraints)

			prob.solve(verbose = True, solver = cp.GUROBI,  max_iters = 1000, BarQCPConvTol = 1e-6)
			
			i_iter += 1
			cost_next = prob.value
			cost_diff = cost_next - cost_curr
			cost_curr = cost_next
			
			print('Objective: ' + str(prob.value))
			print('Objective Change: ' + str(cost_diff))


		U = np.transpose( np.asarray([y for y in u.value]), (0,2,1))
	else:
		U = np.zeros( (param.get('m'), len(T) ))
	# print(np.shape(U))
	return U


def get_scp_initial_trajectory(x0, T):
	# use feedback linearizing solution

	X = []
	U = []
	x_curr = x0

	for t in T:
		u_curr = get_fdbk_controller( x_curr, t) 
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