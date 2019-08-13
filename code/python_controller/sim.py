
from param import param
import numpy as np 
import dynamics
import controller
import plotter
import utilities as util
import subprocess, os, timeit 

def main():

	np.random.seed(88)

	util.get_x0()
	util.get_xd()
	dynamics.compute_jacobians()

	C = []
	print('Sim...')
	for u_type in param.get('controllers'):

		param['controller'] = u_type

		print('Controller: ' + str(u_type))
		
		X = []
		U = []
		V = []
		pbar = []
		vbar = []
		abar = []
		x_curr = param.get('x0')
		count = 0 
		for k,t in enumerate(param.get('T')):

			# try:
			
			if param.get('controller') is 'scp':
				if count == 0:
					U_scp = controller.get_u( x_curr, t)
				try:
					u_curr = U_scp[count]
				except:
					u_curr = np.zeros((param.get('m'),1))

			elif param.get('controller') is 'mpc':
				if np.mod( count, param.get('mpc_update')) == 0:
					count = 0
					U_mpc = controller.get_u( x_curr, t)
				try:
					u_curr = U_mpc[count]
				except:
					u_curr = np.zeros((param.get('m'),1))

			else:
				u_curr = controller.get_u( x_curr,t)

			# except:
			# 	break

			x_next = x_curr + dynamics.get_dxdt(x_curr,u_curr,t) * param.get('dt')

			X.append(x_curr)
			U.append(u_curr)
			V.append(dynamics.get_V( x_curr,t))

			# for error plot 
			pbar.append( np.dot( \
				np.transpose( util.get_my_1()), np.dot( util.get_p_a(), x_curr)))
			vbar.append( np.dot( \
				np.transpose( util.get_my_1()), np.dot( util.get_v_a(), x_curr)))
			abar.append( np.dot( \
				np.transpose( util.get_my_1()), dynamics.get_vdot_a(x_curr)))

			x_curr = x_next
			count += 1 
			print( '\t' + str(t) + '/' + str(param.get('T')[-1]))

			# if k > 2:
			# 	break


		X = np.squeeze(np.asarray(X))
		U = np.squeeze(np.asarray(U))
		V = np.asarray(V)
		pbar = np.asarray(pbar)
		vbar = np.asarray(vbar)
		abar = np.asarray(abar)

		C.append( controller.calc_cost(U))

		print('Plots...')
		plotter.plot_SS(X, param.get('T'), title = str(u_type) + ' State Space' )
		plotter.plot_V( V, param.get('T'), title = str(u_type) + ' Lyapunov' )
		plotter.plot_U( U, param.get('T'), title = str(u_type) + ' Controller')
		plotter.plot_test2( \
			np.squeeze(pbar) - np.squeeze(param.get('pd')), \
			np.squeeze(vbar) - np.squeeze(param.get('vd')), \
			np.squeeze(abar) - np.squeeze(param.get('ad')), \
			param.get('T'), title = str(u_type) + ' Errors')
		print('Plots Complete')


	# plotter.plot_cost( C)
	plotter.save_figs()


	if param.get('gif_on'):
		print('Gif...')
		plotter.make_gif(X,param.get('T'))
		print('Gif Complete')

if __name__ == "__main__":

	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))
	gif_path = os.path.join( os.getcwd(), param.get('fn_gif'))

	# remove if exists 
	if os.path.exists(gif_path): 
		os.remove(gif_path)
	if os.path.exists(pdf_path):
		os.remove(pdf_path)

	start = timeit.default_timer()	
	main()
	stop = timeit.default_timer()
	print('Total Time: ', stop - start)

	subprocess.call(["xdg-open", pdf_path])

	if param.get('gif_on'):
		proc = subprocess.Popen(['animate', gif_path])
		proc.communicate()