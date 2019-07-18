
from param import param
from scipy.linalg import solve_lyapunov, block_diag
import numpy as np 


def get_min_dist( x):
	min_dist = np.inf 
	for i in range(param.get('ni')):
		pose_i = get_p(x,i)
		for j in range(param.get('ni')):
			if i is not j:
				pose_j = get_p(x,j)
				dist = np.linalg.norm( pose_i - pose_j)
				if dist < min_dist:
					min_dist = dist
	return min_dist

def get_A( x):
	# TODO
	A = np.ones([param.get('ni'), param.get('ni')])
	return A

def get_v( x, i):
	# get velocity of agent i
	state_idx = 2*param.get('nd')*i + param.get('nd')
	return x[state_idx:state_idx+param.get('nd')]

def get_p( x, i):
	# get position of agent i
	state_idx = 2*param.get('nd')*i
	return x[state_idx:state_idx+param.get('nd')]

def get_p_a( x):
	# get stacked positions of all free agents
	p_a = []
	for i in range(param.get('na')):
		p_a.append( get_p( x, i))
	return list_to_vec(np.asarray(p_a).flatten())

def get_v_a( x):
	# get stacked velocity of all free agents
	v_a = []
	for i in range(param.get('na')):
		v_a.append( get_v( x, i))
	return list_to_vec(np.asarray(v_a).flatten())

def get_v_b( x):
	v_b = []
	for i in range(param.get('nb')):
		v_b.append( get_v( x, i + param.get('na')))
	return list_to_vec(np.asarray(v_b).flatten())

def get_P( X, i):
	# get position of agent i for all time
	state_idx = 2*param.get('nd')*i
	return X[:,state_idx:state_idx+param.get('nd')]

def get_P_bar( X,T):
	# get centroid position of free agents for all time 
	P_bar = []
	my_1 = get_my_1()
	for i_t in range(len(T)):
		p_a = get_p_a( X[i_t])
		P_bar.append( np.matmul( np.transpose(my_1), p_a))
	return np.asarray(P_bar)

def get_V( X, i):
	# get velocity of agent i for all time 
	state_idx = 2*param.get('nd')*i + param.get('nd')
	return X[:,state_idx:state_idx+param.get('nd')]

def get_lyapunov_metric():

	# closed loop output dynamics for single output
	a_cl = np.eye(param.get('gamma') - 1)
	a_cl = np.hstack(( \
		np.zeros( (param.get('gamma')-1,1)), a_cl \
		))
	a_cl = np.vstack(( \
		a_cl, -param.get('k_fdbk')*np.ones((1, param.get('gamma'))) \
		))

	# full closed loop output dynamics
	A_cl = a_cl
	for _ in range(param.get('nd')-1):
		A_cl = block_diag( A_cl, a_cl)

	# solution to lyapunov equation
	P = solve_lyapunov( np.transpose(A_cl), \
		-np.eye(param.get('nd')*param.get('gamma')))
	return P

def get_stabilization_rate():
	P = get_lyapunov_metric()
	l, v = np.linalg.eig(P)
	return 1/np.max(l)

def get_my_1():
	return np.kron(np.ones((param.get('na'),1)),np.eye(param.get('nd')))/ \
		param.get('na')

def list_to_vec(x):
	return np.reshape(x, (len(x),-1))

def list_of_list_to_vec(list_of_list):
	return list_to_vec( [x for sublist in list_of_list for x in sublist])

def get_plot_lim(X,T):
	xmin = np.inf
	ymin = np.inf
	xmax = -np.inf
	ymax = -np.inf

	# check agents
	for t in range(len(T)):
		x = X[t]
		for i in range(param.get('ni')):
			pose_i = get_p(x,i)
			if xmin > pose_i[0]:
				xmin = pose_i[0]
			if xmax < pose_i[0]:
				xmax = pose_i[0]
			if ymin > pose_i[1]:
				ymin = pose_i[1]
			if ymax < pose_i[1]:
				ymax = pose_i[1]

	# check desired
	for t in range(len(T)):
		pose_d = param.get('pd')[t,:]
		if xmin > pose_d[0]:
			xmin = pose_d[0]
		if xmax < pose_d[0]:
			xmax = pose_d[0]
		if ymin > pose_d[1]:
			ymin = pose_d[1]
		if ymax < pose_d[1]:
			ymax = pose_d[1]

	# add buffer
	if xmin < 0:
		xmin = xmin*(1+param.get('plot_buffer'))
	else:
		xmin = xmin*(1-param.get('plot_buffer'))
	if ymin < 0:
		ymin = ymin*(1+param.get('plot_buffer'))
	else:
		ymin = ymin*(1-param.get('plot_buffer'))
	if xmax < 0:
		xmax = xmax*(1-param.get('plot_buffer'))
	else:
		xmax = xmax*(1+param.get('plot_buffer'))
	if ymax < 0:
		ymax = ymax*(1-param.get('plot_buffer'))
	else:
		ymax = ymax*(1+param.get('plot_buffer'))

	return xmin, xmax, ymin, ymax

def permute_rows( x):
	perm_mat = np.zeros( (len(x), len(x)))
	gamma = param.get('gamma') # relative degree

	for dim_idx in range( param.get('nd')):
		row_idx = dim_idx* gamma
		for gamma_idx in range( gamma):
			col_idx = gamma_idx* param.get('nd') + dim_idx
			perm_mat[row_idx, col_idx] = 1
			row_idx += 1
	return np.matmul( perm_mat, x)
	# return x

def permute_rows_2( x):

	perm_mat = np.zeros( (param.get('n'), param.get('n')))
	for i_a in range(param.get('na')):
		old_p_idx = i_a*param.get('nd')
		old_v_idx = i_a*param.get('nd') + param.get('na')*param.get('nd')
		new_p_idx = i_a*2*param.get('nd')
		new_v_idx = i_a*2*param.get('nd') + param.get('nd')
		for i_d in range(param.get('nd')):
			perm_mat[ old_p_idx+i_d, new_p_idx+i_d] = 1
			perm_mat[ old_v_idx+i_d, new_v_idx+i_d] = 1

	for i_b in range(param.get('nb')):
		old_p_idx = i_b*param.get('nd') + 2*param.get('na')*param.get('nd')
		old_v_idx = old_p_idx + param.get('nb')*param.get('nd')
		new_p_idx = i_b*2*param.get('nd') + 2*param.get('na')*param.get('nd')
		new_v_idx = new_p_idx + param.get('nd')
		for i_d in range(param.get('nd')):
			perm_mat[ old_p_idx+i_d, new_p_idx+i_d] = 1
			perm_mat[ old_v_idx+i_d, new_v_idx+i_d] = 1
		
	return np.matmul( perm_mat,x)


def permute_cols( x):

	perm_mat = np.zeros( (param.get('n'), param.get('n')))
	for i_a in range(param.get('na')):
		old_p_idx = i_a*param.get('nd')
		old_v_idx = i_a*param.get('nd') + param.get('na')*param.get('nd')
		new_p_idx = i_a*2*param.get('nd')
		new_v_idx = i_a*2*param.get('nd') + param.get('nd')
		for i_d in range(param.get('nd')):
			perm_mat[ old_p_idx+i_d, new_p_idx+i_d] = 1
			perm_mat[ old_v_idx+i_d, new_v_idx+i_d] = 1

	for i_b in range(param.get('nb')):
		old_p_idx = i_b*param.get('nd') + 2*param.get('na')*param.get('nd')
		old_v_idx = old_p_idx + param.get('nb')*param.get('nd')
		new_p_idx = i_b*2*param.get('nd') + 2*param.get('na')*param.get('nd')
		new_v_idx = new_p_idx + param.get('nd')
		for i_d in range(param.get('nd')):
			perm_mat[ old_p_idx+i_d, new_p_idx+i_d] = 1
			perm_mat[ old_v_idx+i_d, new_v_idx+i_d] = 1
		
	return np.matmul( x, perm_mat)

def get_x0():

	min_dist = np.inf
	count = 0 
	while min_dist > param.get('min_dist'):
		
		x0 = []
		for i in range(param.get('ni')):
			p_i = param.get('plim')*np.random.rand(param.get('nd'),1) \
				- param.get('plim')/2.
			v_i = param.get('vlim')*np.random.rand(param.get('nd'),1) \
				- param.get('vlim')/2.
			x0.append( np.vstack((p_i, v_i)))

		count += 1
		min_dist = get_min_dist( list_of_list_to_vec(x0))

		if count > 100:
			print('Error: Incompatible Initial Conditions')
			return 

	param['x0'] = list_of_list_to_vec(x0)

def get_xd():

	param['pd'] = np.zeros((len(param.get('T')),2))
	# point, x = y = 1 for all t
	if param.get('case_xd') == 0:
		param['pd'] = np.ones((len(param.get('T')),2))
	elif param.get('case_xd') == 1:
		param['pd'][:,0] = np.squeeze(param.get('T'))
		param['pd'][:,1] = np.squeeze(param.get('T'))
	elif param.get('case_xd') == 2:
		param['pd'][:,0] = np.cos( 2*np.pi*np.squeeze(param.get('T'))/param.get('T')[-1])
		param['pd'][:,1] = np.sin( 2*np.pi*np.squeeze(param.get('T'))/param.get('T')[-1])

	param['vd'] = np.gradient(param.get('pd'), param.get('dt'), axis = 0)
	param['ad'] = np.gradient(param.get('vd'), param.get('dt'), axis = 0)
	param['jd'] = np.gradient(param.get('ad'), param.get('dt'), axis = 0)
