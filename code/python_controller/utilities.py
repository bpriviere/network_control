
from param import param
from scipy.linalg import solve_lyapunov, block_diag
import autograd.numpy as np 

def get_min_dist( x):
	min_dist = np.inf 
	for i in range(param.get('ni')):
		pose_i = np.dot(get_p_i(i), x)
		for j in range(param.get('ni')):
			if i is not j:
				pose_j = np.dot(get_p_i(j), x)
				dist = np.linalg.norm( pose_i - pose_j)
				if dist < min_dist:
					min_dist = dist
	return min_dist

def get_A( x):
	
	A = np.ones([param.get('ni'), param.get('ni')])
	# A = np.zeros([param.get('ni'), param.get('ni')])
	# for i in range(param.get('ni')):
	# 	p_i = np.dot( get_p_i(i), x)
	# 	for j in range(param.get('ni')):
	# 		p_j = np.dot( get_p_i(j), x)
	# 		dist = np.linalg.norm( p_i - p_j)
	# 		if dist < param.get('R_comm'):

	# 			print('Adj')
	# 			print(np.shape(x))
	# 			print(np.shape(p_i))
	# 			print(np.shape(p_j))
	# 			print(np.shape(dist))
	# 			print(param.get('lambda_a'))
	# 			A[i,j] = np.exp( -param.get('lambda_a')*dist)
	# print(A)
	return A

def get_p_a():
	# get matrix to extract stacked positions of all free agents
	pi = np.zeros( (param.get('na'), 2*param.get('ni') ))
	for i in range(param.get('na')):
		pi[i,2*param.get('nd')*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_v_a():
	# get matrix to extract stacked velocities of all free agents
	pi = np.zeros( (param.get('na'), 2*param.get('ni') ))
	for i in range(param.get('na')):
		pi[i,2*i + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_v_b():
	# get matrix to extract stacked velocities of all control agents
	pi = np.zeros( (param.get('nb'), 2*param.get('ni')))
	for i in range(param.get('nb')):
		pi[i,2*(i+param.get('ni')) + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_v_i(i):
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_p_i(i):
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_P_bar( X,T):

	my_1 = get_my_1()
	P_bar = []
	for k,t in enumerate(T):
		P_bar.append( \
			np.dot( np.transpose( my_1), \
			np.dot( get_p_a(), X[k,:])))
	return np.asarray(P_bar)

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

def permute_eta_rows():
	perm_mat = np.zeros( (param.get('nd')*param.get('gamma'), \
		param.get('nd')*param.get('gamma')))
	gamma = param.get('gamma') # relative degree

	for dim_idx in range( param.get('nd')):
		row_idx = dim_idx* gamma
		for gamma_idx in range( gamma):
			col_idx = gamma_idx* param.get('nd') + dim_idx
			perm_mat[row_idx, col_idx] = 1
			row_idx += 1
	return perm_mat

def permute_x_rows( x):

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
	
	return np.dot( perm_mat, x)


def permute_x_cols():

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
		
	return perm_mat

def get_x0():

	min_dist = -np.inf
	count = 0 
	while min_dist < param.get('min_dist'):
		
		for i in range(param.get('ni')):
			p_i = param.get('plim')*np.random.rand(param.get('nd'),1) \
				- param.get('plim')/2.
			v_i = param.get('vlim')*np.random.rand(param.get('nd'),1) \
				- param.get('vlim')/2.
			try:
				x0 = np.vstack((x0, p_i, v_i))
			except:
				x0 = np.vstack((p_i, v_i))

		count += 1
		min_dist = get_min_dist( x0)

		if count > 100:
			print('Error: Incompatible Initial Conditions')
			return 

	param['x0'] = x0

def get_xd():

	pd = []
	vd = []
	ad = []
	for t in param.get('T'):
		pd.append( get_pd_t(t))
		vd.append( get_vd_t(t))
		ad.append( get_ad_t(t))
	param['pd'] = np.asarray(pd)
	param['vd'] = np.asarray(vd)
	param['ad'] = np.asarray(ad)

def get_pd_t(t):

	if param.get('case_xd') == 0:
		return np.ones((2,1))
	elif param.get('case_xd') == 1:
		return t*np.ones((2,1))
	elif param.get('case_xd') == 2:
		return np.array([ 
				np.cos( 2*np.pi/param.get('tf')*t), 
				np.sin( 2*np.pi/param.get('tf')*t)
			])
	return 

def get_vd_t(t):

	if param.get('case_xd') == 0:
		return np.zeros((2,1))
	elif param.get('case_xd') == 1:
		return np.ones((2,1))
	elif param.get('case_xd') == 2:
		return (2*np.pi/param.get('tf'))*np.array( [
				-np.sin( 2*np.pi/param.get('tf')*t),
				 np.cos( 2*np.pi/param.get('tf')*t)
			])
	return 

def get_ad_t(t):

	if param.get('case_xd') == 0:
		return np.zeros((2,1))
	elif param.get('case_xd') == 1:
		return np.zeros((2,1))
	elif param.get('case_xd') == 2:
		return np.power( 2*np.pi/param.get('tf'),2.)*np.array( [
				-np.cos( 2*np.pi/param.get('tf')*t),
				-np.sin( 2*np.pi/param.get('tf')*t)
			])
	return 
