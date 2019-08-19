
import numpy as np 
import cvxpy as cp
import dynamics 
import utilities as util
from param import param


def get_centralized_ta(x,t):

	if param.get('ta_on'):
		c = get_cij(x,t)
		my_1 = np.ones((param.get('ni'),1))
		pi = cp.Variable((param.get('ni'),param.get('ni')))
		
		constraints = []

		constraints.append(pi >= 0)
		constraints.append(pi <= 1)
		for i in range(param.get('ni')):
			constraints.append(
				sum(pi[i,:]) == 1 )
			constraints.append(
				sum(pi[:,i]) == 1 )

		obj = cp.Minimize(cp.sum( cp.multiply(c,pi)))
		prob = cp.Problem(obj,constraints)
		# prob.solve(solver = 'ECOS', verbose = True)
		prob.solve()
		pi = pi.value
	else:
		pi = np.eye(param.get('ni'))

	return np.round(pi)

def get_cij(x,t):
	# euclidean distance heuristic 
	pl = dynamics.get_xl(x,t)[0:param.get('nd')]
	pb = dynamics.get_xb(x,t)[0:param.get('nd')]
	T = dynamics.get_T(x,t)
	c = np.zeros((param.get('ni'),param.get('ni')))
	for i in range(param.get('ni')):
		# agent positions
		pi = np.dot(util.get_p_i(i),x)

		for j in range(param.get('ni')):
			# extract configuration
			my_1 = np.zeros((param.get('ni'),1))
			my_1[j] = 1
			my_1 = np.kron(my_1, np.eye(param.get('nd')))
			Tj = np.dot(np.dot(
				my_1.T, T), my_1)

			# desired positions
			pj = pl + \
				np.dot(Tj,pb)

			c[j,i] = np.linalg.norm(
				pi - pj
				)

	return c