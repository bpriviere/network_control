
from param import param
import autograd.numpy as np 
from autograd.extend import primitive, defvjp  # defvjp is now a function
# from autograd import primitive
from scipy.linalg import block_diag as block_diag


def get_A(x):
	
	# pose
	P = np.dot(get_p(),x)
	
	X = np.dot( 
		np.ones((param.get('ni'),1)),
		to_vec(P[np.mod(np.arange(0,len(P)),2)==False]).T)
	Y = np.dot( 
		np.ones((param.get('ni'),1)),
		to_vec(P[np.mod(np.arange(0,len(P)),2)==True]).T)
	
	D = np.sqrt( 
		np.power( X - X.T + 1e-9,2) + 
		np.power( Y - Y.T + 1e-9,2))

	# print(D)
	
	A = np.exp( -param.get('lambda_a') * D)

	A = A / np.linalg.norm(A, ord=2, axis=1)

	return A


def get_A_a(x):

	A = get_A(x)[0:param.get('na'),0:param.get('na')]
	return A

def get_L(x):
	# Laplacian Matrix
	A = get_A(x)
	D = np.zeros( np.shape(A))
	for i in range( np.shape(A)[0]):
		D[i,i] = sum( A[i,:])
	return D - A

def get_L_a(x):
	A = get_A_a(x)
	D = np.zeros( np.shape(A))
	for i in range( np.shape(A)[0]):
		D[i,i] = sum( A[i,:])
	return D - A

def get_Vs_a(x):
	L = get_L_a(x)
	l, v = np.linalg.eig(L)
	l = np.round(l, 5)
	return v[:, l!=0]

def get_Vs(x):
	L = get_L(x)
	l, v = np.linalg.eig(L)
	l = np.round(l, 5)
	return v[:, l!=0]

def augment(x):
	return np.kron( x, np.eye(2*param.get('nd')))

def to_vec(x):
	return np.array(np.reshape(x, (len(x), 1)))

def get_p_a():
	# get matrix to extract stacked positions of all free agents
	pi = np.zeros( (param.get('na'), 2*param.get('ni') ))
	for i in range(param.get('na')):
		pi[i,2*param.get('nd')*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_p():
	# get matrix to extract stacked positions of all agents
	pi = np.zeros( (param.get('ni'), 2*param.get('ni') ))
	for i in range(param.get('ni')):
		pi[i,2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_pv_a():
	# get matrix to extract stacked positions and velocity of all free agents
	pi = np.eye( param.get('na')*2*param.get('nd'))
	pi = np.hstack(( pi, 
		np.zeros(( param.get('na')*2*param.get('nd'), param.get('nb')*2*param.get('nd'))) ))
	return pi

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
		pi[i,2*(i + param.get('na')) + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_v_i(i):
	# get row vector to extract velocity of a single agent
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_p_i(i):
	# get row vector to extract position of a single agent
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))

def get_p_i_idx(i):
	return i*2*param.get('nd') + np.arange( 0, param.get('nd'))

def get_v_i_idx(i):
	return (i*2*param.get('nd') + param.get('nd')) + np.arange( 0, param.get('nd'))	

def get_R(phase):
	# rotation matrix, 2D
	return np.vstack(( 
		np.hstack(( np.cos(phase), -np.sin(phase) )),
		np.hstack(( np.sin(phase),  np.cos(phase) )) ))

def get_pi():
	# permutation matrix
	# return np.kron( np.eye(param.get('ni')), np.eye(2*param.get('nd')))
	pi = np.asarray([
		[0,0,1,0,0],
		[1,0,0,0,0],
		[0,1,0,0,0],
		[0,0,0,1,0],
		[0,0,0,0,1]
		])
	return np.kron( pi, np.eye(2*param.get('nd')))
	# return np.eye(param.get('n'))

def get_pi_a():
	# pi = np.asarray([
	# 	[0,0,1],
	# 	[0,1,0],
	# 	[1,0,0]
	# 	])
	pi = np.eye(param.get('na'))
	return np.kron( pi, np.eye(2*param.get('nd')))

def get_T():
	
	T = block_diag( \
			param.get('radius_d')[0]*get_R(param.get('phase_d')[0]),\
			np.eye(param.get('nd')))
	for i in range(1,param.get('na')): 
		t = block_diag( \
			param.get('radius_d')[i]*get_R(param.get('phase_d')[i]),\
			np.eye(param.get('nd')))
		T = block_diag(T, t)
		
	for i in range(param.get('nb')):
		T = block_diag(T, np.zeros((2*param.get('nd'),2*param.get('nd'))))

	return T

def get_T_a():
	
	T = block_diag( \
			param.get('radius_d')[0]*get_R(param.get('phase_d')[0]),\
			np.eye(param.get('nd')))
	for i in range(1,param.get('na')): 
		t = block_diag( \
			param.get('radius_d')[i]*get_R(param.get('phase_d')[i]),\
			np.eye(param.get('nd')))
		T = block_diag(T, t)
		
	return T	

def get_xd():

	if param.get('case_xd') is "3_node_triangle":

		param['phase_d'] = [
			np.pi*2/3,
			np.pi*4/3,
			np.pi*6/3,
		]

		param['radius_d'] = [
			1,
			1,
			1
		]

		param['ni_d'] = 3

	elif param.get('case_xd') is "single_node":
		param['phase_d'] = 0
		param['radius_d'] = 1
		param['ni_d'] = 1

	if param.get('na') is not param.get('ni_d'):
		print('Error: Number of Agents Mismatch')
		exit()

	if param.get('ni_d') > 1:
		xd = np.zeros( param.get('n'))
		for i in range( param.get('na')):
			xd[get_p_i_idx(i)] = param.get('radius_d')[i] * np.asarray( [
				np.cos( param.get('phase_d')[i]), np.sin( param.get('phase_d')[i]) ])
		param['xd'] = xd
	else:
		param['xd'] = np.zeros(param.get('n')) 


def get_x0():
	min_dist = -np.inf
	count = 0 
	while min_dist < param.get('min_dist'):
		
		x0 = np.zeros(param.get('n'))
		for i in range(param.get('ni')):
			p_i = param.get('plim')*np.random.rand(param.get('nd')) \
				- param.get('plim')/2.
			v_i = param.get('vlim')*np.random.rand(param.get('nd')) \
				- param.get('vlim')/2.
			
			x0[get_p_i_idx(i)] = p_i
			x0[get_v_i_idx(i)] = v_i

		count += 1
		min_dist = get_min_dist( x0)

		if count > 100:
			print('Error: Incompatible Initial Conditions')
			return 

	param['x0'] = to_vec(x0)

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



# @primitive
# def get_A(x):
# 	# Adjacency Matrix

# 	A = np.zeros((param.get('ni'), param.get('ni')))
# 	for i in range(param.get('ni')):
# 		p_i = np.dot(get_p_i(i), x)
# 		for j in range(param.get('ni')):
# 			p_j = np.dot(get_p_i(j), x)
# 			dist = np.linalg.norm( p_i - p_j)
# 			# if dist < param.get('R_comm'):
# 			A[i,j] = np.exp( -param.get('lambda_a') * dist)
	
# 	# A = np.ones((param.get('ni'), param.get('ni')))
# 	return A

# def get_A_vjp(ans, x):

# 	try:
# 		x = x._value
# 	except:		
# 		pass

# 	def grad_A(g):

# 		grad_A = np.zeros( (param.get('ni'), param.get('ni'), param.get('n')))
# 		for i in range(param.get('ni')):
# 			p_i = np.dot(get_p_i(i), x)
# 			p_i_idx = get_p_i_idx(i)
# 			for j in range(param.get('ni')):
# 				if j != i:
# 					p_j = np.dot(get_p_i(j), x)
# 					p_j_idx = get_p_i_idx(j)
# 					dist = np.linalg.norm( p_i - p_j)
# 					# if dist < param.get('R_comm'):
# 					C = np.squeeze(
# 						np.exp( -param.get('lambda_a')*dist)*
# 						-param.get('lambda_a')/dist*
# 						(p_i - p_j))

# 					grad_A[i,j,p_i_idx] = C
# 					grad_A[i,j,p_j_idx] = -C

# 		return np.tensordot( g, grad_A)

# 	return grad_A

# def my_check_grad( func, dfunc, x):

# 	fx = func(x)
# 	eps = 1e-11
# 	err = np.zeros((*np.shape(fx)[:], len(x)))
# 	for i in range( len(x)):
# 		print('x_idx: ' + str(i))
# 		x1 = x 
# 		x1[i] = x[i] + eps
# 		fx1 = func(x1)
# 		df = (fx1 - fx) / eps
# 		err[:,:,i] = df - dfunc(x)[:,:,i]
	
# 	print(err)
# 	print(err.shape)