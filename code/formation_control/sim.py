

from param import param
import dynamics
import controller
import autograd.numpy as np 
from autograd.extend import primitive, defvjp 
from autograd import primitive, jacobian
import utilities as util
import plotter
import os, timeit, subprocess

def main():

	np.random.seed(88)

	# 
	util.get_x0()
	util.get_xd()
	dynamics.compute_jacobians()
	defvjp(util.get_A_a, util.get_A_a_vjp)
	defvjp(util.get_A, util.get_A_vjp)	

	x_curr = param.get('x0')

	Vs = util.augment(util.get_Vs_a(param.get('xd')))		
	VsT_T_pi_a = np.dot( np.dot( np.dot( 
		np.transpose(Vs), util.get_T_a()), util.get_pi_a()), util.get_pv_a() )
	delta = np.dot( VsT_T_pi_a, param.get('xd'))
	print(delta)
 

	X = np.zeros( (param.get('n'), param.get('nt')))
	U = np.zeros( (param.get('m'), param.get('nt')))
	D = np.zeros( (2*param.get('nd')*(param.get('na')-1), param.get('nt')))
	for k,t in enumerate( param.get('T')):

		if param.get('controller') is 'mpc':
			if np.mod( k, param.get('mpc_update')) == 0:
				U_mpc = controller.get_u( x_curr, t)
			try:
				u_curr = util.to_vec(U_mpc[:,np.mod( k, param.get('mpc_update'))])
			except:
				u_curr = np.zeros( (param.get('m'),1))
		else:
			u_curr = controller.get_u( x_curr, t)

		x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')
		
		Vs = util.augment(util.get_Vs_a(param.get('xd')))		
		VsT_T_pi_a = np.dot( np.dot( np.dot( 
			np.transpose(Vs), util.get_T_a()), util.get_pi_a()), util.get_pv_a() )
		delta = np.dot( VsT_T_pi_a, x_curr)

		X[:,k] = np.squeeze(x_curr)
		U[:,k] = np.squeeze(u_curr)
		D[:,k] = np.squeeze(delta)

		x_curr = x_next
		print( '\t' + str(t) + '/' + str(param.get('T')[-1]))

	plotter.plot_SS( X, param.get('T'))
	plotter.plot_VsTx( D, param.get('T'), title = 'VsTx')
	plotter.plot( param.get('T'), np.linalg.norm(D, axis = 0))
	plotter.save_figs()

if __name__ == '__main__':

	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))
	if os.path.exists(pdf_path):
		os.remove(pdf_path)
	
	start = timeit.default_timer()	
	main()
	stop = timeit.default_timer()
	print('Total Time: ', stop - start)

	subprocess.call(["xdg-open", pdf_path])





	# x_curr_nl = param.get('x0')
	# x_curr_l = param.get('x0') + .05 * np.ones( (param.get('n'), 1))
	# for k,t in enumerate( param.get('T')):

	# 	u_curr = np.sin( 2*np.pi* t / param.get('T')[-1])*np.ones( (param.get('m'),1))

	# 	x_next_nl = x_curr_nl + dynamics.get_dxdt( x_curr_nl, u_curr, t) * param.get('dt')

	# 	F_k, B_k, d_k = dynamics.get_linear_dynamics( x_curr_l, u_curr, t)
	# 	x_next_l = np.dot(F_k, x_curr_l) + np.dot(B_k, u_curr) + d_k 

	# 	# x_next_l = x_curr_l + dynamics.get_dxdt( x_curr_l, u_curr, t) * param.get('dt')

	# 	try:
	# 		X_nl = np.hstack( (X_nl, x_curr_nl))
	# 		X_l = np.hstack( (X_l, x_curr_l))
	# 	except:
	# 		X_nl = x_curr_nl
	# 		X_l = x_curr_l

	# 	x_curr_nl = x_next_nl
	# 	x_curr_l = x_next_l
	# 	print( '\t' + str(t) + '/' + str(param.get('T')[-1]))

	# plotter.plot( param.get('T'), np.linalg.norm( X_nl - X_l, axis = 0))