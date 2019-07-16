

import numpy as np 
import utilities as util
from param import param

def get_dxdt(x,u,t):

    f = get_f(x)
    g = get_g()
    dxdt = f + np.matmul(g,u);
    return dxdt

def get_f(x):
	# drift dynamics
	f = []
	if param.get('model') is 'reynolds':

		# free agents
		for i in range(param.get('na')):
			v_i = util.get_v(x,i)
			a_i = reynolds( x, i)
			f.append( np.vstack( (v_i, a_i)))

		# control agents
		for i in range(param.get('nb')):
			v_i = util.get_v(x,i+param.get('na'))
			a_i = np.zeros([param.get('nd'),1])
			f.append( np.vstack( (v_i, a_i)))

	# flattens list of list into a single list, then converts to numpy array
	return util.list_of_list_to_vec(f)

def get_g():
	# control matrix
	g = np.zeros( [param.get('n'), param.get('m')])

	for i in range(param.get('nb')):
		row_idx = 2*param.get('nd')*i + \
			2*param.get('na')*param.get('nd') + param.get('nd')
		col_idx = i*param.get('nd') 
		g[row_idx:row_idx + param.get('nd'), col_idx:col_idx + param.get('nd')] = \
			np.eye(param.get('nd')) 
	return g

def get_LgV( x,t):
	dVdeta = get_dVdeta(x,t)
	detadx = get_detadx(x,t)
	g 	   = get_g()
	return np.matmul( dVdeta, np.matmul( detadx, g))

def get_LfV( x,t):
	dVdeta = get_dVdeta(x,t)
	detadx = get_detadx(x,t)
	detadt = get_detadt(x,t)
	f 	   = get_f(x)
	return np.matmul( dVdeta, np.matmul( detadx, f) + detadt)

def get_V( x,t):
	# Lyapunov (not velocity)
	P = util.get_lyapunov_metric()
	eta = get_eta( x,t)
	return np.matmul( np.transpose(eta), np.matmul(P, eta))

def get_dVdeta(x,t):
	return 2 * np.matmul( np.transpose(get_eta(x,t)), \
		util.get_lyapunov_metric())

def get_detadx(x,t):
	
	my_1 = util.get_my_1()
	O_a = np.zeros( (param.get('nd'), param.get('nd')*param.get('na')))
	O_b = np.zeros( (param.get('nd'), param.get('nd')*param.get('nb')))
	dvdota_dpa = get_dvdota_dpa(x,t)
	dvdota_dva = get_dvdota_dva(x,t)
	dvdota_dpb = get_dvdota_dpb(x,t)
	dvdota_dvb = get_dvdota_dvb(x,t)

	row_1 = np.hstack( (np.transpose(my_1), O_a, O_b, O_b))
	row_2 = np.hstack( (O_a, np.transpose(my_1), O_b, O_b))
	row_3 = np.hstack(( \
		np.matmul(np.transpose(my_1), dvdota_dpa), \
		np.matmul(np.transpose(my_1), dvdota_dva), \
		np.matmul(np.transpose(my_1), dvdota_dpb), \
		np.matmul(np.transpose(my_1), dvdota_dvb)))
	return util.permute_rows(np.vstack( (row_1, row_2, row_3)))

def get_detadt( x,t):
	k = np.where(param.get('T') == t)[0][0]
	yp  = - util.list_to_vec( param.get('vd')[k,:])
	ypp = - util.list_to_vec( param.get('ad')[k,:]) 
	yppp= - util.list_to_vec( param.get('jd')[k,:]) 
	detadt = np.vstack( (yp, ypp, yppp));
	return util.permute_rows( detadt);

def get_eta( x,t):
	# tracking output dynamics

	k = np.where(param.get('T') == t)[0][0]
	my_1 = util.get_my_1()
	pose_a = util.get_p_a( x);
	vel_a  = util.get_v_a( x);
	vdot_a = get_vdot_a( x,t);

	y   = np.matmul( np.transpose(my_1), pose_a) \
		- util.list_to_vec( param.get('pd')[k,:])
	yp  = np.matmul( np.transpose(my_1), vel_a) \
		- util.list_to_vec( param.get('vd')[k,:])
	ypp = np.matmul( np.transpose(my_1), vdot_a) \
		- util.list_to_vec( param.get('ad')[k,:])

	eta = np.vstack( (y, yp, ypp));
	eta = util.permute_rows( eta);

	return eta

def get_vdot_a( x,t):
	vdot_a = []
	for i in range(param.get('na')):
		vdot_a.append( reynolds( x,i))
	return util.list_to_vec(np.asarray(vdot_a).flatten())

def get_dvdota_dpa( x,t):

	dvdota_dpa = np.zeros( (param.get('na')*param.get('nd'),\
		param.get('na')*param.get('nd')))

	for i in range(param.get('na')):
		for j in range(param.get('na')):
			if i is not j:
				row_idx_0 = i*param.get('nd')
				row_idx_1 = row_idx_0 + param.get('nd')
				col_idx_0 = j*param.get('nd') 
				col_idx_1 = col_idx_0 + param.get('nd')
				dvdota_dpa[ row_idx_0:row_idx_1, col_idx_0:col_idx_1] = \
					get_reynolds_position_derivative_ij( x, i, j)
	return dvdota_dpa

def get_dvdota_dpb( x,t):

	dvdota_dpb = np.zeros( (param.get('na')*param.get('nd'),\
		param.get('nb')*param.get('nd')))

	for i in range(param.get('na')):
		for j in range(param.get('nb')):
			row_idx_0 = i*param.get('nd')
			row_idx_1 = row_idx_0 + param.get('nd')
			col_idx_0 = j*param.get('nd') 
			col_idx_1 = col_idx_0 + param.get('nd')
			dvdota_dpb[ row_idx_0:row_idx_1, col_idx_0:col_idx_1] = \
				get_reynolds_position_derivative_ij( x, i, j + param.get('na'))
	return dvdota_dpb

def get_xtildedot( x,t):
	v_a 	= util.get_v_a(x)
	v_b 	= util.get_v_b(x)
	vdot_a  = get_vdot_a( x,t)
	return np.vstack( (v_a, vdot_a, v_b))

def get_dvdota_dxtilde( x,t):
	dvdotadpa = get_dvdota_dpa( x, t)
	dvdotadva = get_dvdota_dva( x, t)
	dvdotadpb = get_dvdota_dpb( x, t)
	return np.hstack( (dvdotadpa, dvdotadva, dvdotadpb))

def get_dvdota_dva( x, t):
	dvdota_dva = np.zeros( (param.get('na')*param.get('nd'),\
		param.get('na')*param.get('nd')))

	for i in range(param.get('na')):
		for j in range(param.get('na')):
			if i is not j:
				row_idx_0 = i*param.get('nd')
				row_idx_1 = row_idx_0 + param.get('nd')
				col_idx_0 = j*param.get('nd') 
				col_idx_1 = col_idx_0 + param.get('nd')
				dvdota_dva[ row_idx_0:row_idx_1, col_idx_0:col_idx_1] = \
					get_reynolds_velocity_derivative_ij( x, i, j)
	return dvdota_dva

def get_dvdota_dvb( x, t):
	dvdota_dvb = np.zeros( (param.get('na')*param.get('nd'),\
		param.get('nb')*param.get('nd')))

	for i in range(param.get('na')):
		for j in range(param.get('nb')):
			row_idx_0 = i*param.get('nd')
			row_idx_1 = row_idx_0 + param.get('nd')
			col_idx_0 = j*param.get('nd') 
			col_idx_1 = col_idx_0 + param.get('nd')
			dvdota_dvb[ row_idx_0:row_idx_1, col_idx_0:col_idx_1] = \
				get_reynolds_velocity_derivative_ij( x, i, j + param.get('na'))
	return dvdota_dvb

def reynolds( x, i):
	A   = util.get_A(x)
	p_i = util.get_p(x,i)
	v_i = util.get_v(x,i)	
	a_i = np.zeros([param.get('nd'),1])
	for j in range(param.get('ni')):
		if i is not j:
			p_j = util.get_p(x,j)
			v_j = util.get_v(x,j)

			r_ij = p_j - p_i
			dist = np.linalg.norm(r_ij)

			a_i = a_i + A[i,j]*( \
                param.get('kv')*(v_j - v_i) + \
                param.get('kx')*r_ij*(1 - param.get('R_des')/dist) 
                )
	return a_i 

def get_reynolds_position_derivative_ij( x, i, j):

	A = util.get_A(x)
	p_i = util.get_p( x,i)
	p_j = util.get_p( x,j)
	r_ij = p_j - p_i
	dist = np.linalg.norm(r_ij)
	
	return A[i,j]*param.get('kx')* \
		( 1 - param.get('R_des')/dist*( np.eye(param.get('nd')) - \
		(np.matmul( r_ij, np.transpose( r_ij)))/dist))

def get_reynolds_velocity_derivative_ij( x, i, j):
	A = util.get_A(x)
	return A[i,j]*param.get('kv')*np.eye(param.get('nd'))

def get_linearized_dynamics( xbar, ubar, t):

	#  xdot = h(x,u) = f(x) + g(x)u
	dhdx = get_dfdx(xbar)
	dhdu = get_g()

	# continuous linearized dynamics
	A  = dhdx 
	B  = dhdu
	d  = np.matmul( -dhdx, xbar) + np.matmul( -dhdu, ubar) + \
		get_f(xbar) + np.matmul( get_g(), ubar)

	# discretize
	Ad = np.eye( param.get('n')) + A*param.get('dt')
	Bd = B*param.get('dt')
	dd = d*param.get('dt')
	return Ad, Bd, dd

def get_dfdx(x):

	O_a = np.zeros((param.get('na')*param.get('nd'), \
		param.get('na')*param.get('nd')))
	O_b = np.zeros((param.get('nb')*param.get('nd'), \
		param.get('nb')*param.get('nd')))
	I_a = np.eye(param.get('na')*param.get('nd'))
	I_b = np.eye(param.get('nb')*param.get('nd'))
	
	dvdotadpa = get_dvdota_dpa( x, t)
	dvdotadva = get_dvdota_dva( x, t)
	dvdotadpb = get_dvdota_dpb( x, t)
	dvdotadvb = get_dvdota_dvb( x, t)

	row1 = np.hstack( (O_a, I_a, O_a, O_a))
	row2 = np.hstack( (dvdotadpa, dvdotadva, dvdotadpb, dvdotadvb))
	row3 = np.hstack( (O_b, O_b, I_b, O_b))
	row4 = np.hstack( (O_b, O_b, O_b, O_b))

	return util.permute_rows(np.vstack( (row1, row2, row3, row4)))