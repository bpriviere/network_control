
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
	C = []
	print('Sim...')
	for u in param.get('controllers'):

		try:
			param['controller'] = u

			print('Controller: ' + str(u))
			
			X = []
			U = []
			V = []
			pbar = []
			vbar = []
			abar = []
			x_curr = param.get('x0')
			cost = 0
			count = 0 
			for t in param.get('T'):
				
				if param.get('controller') is 'scp':
					if count == 0:
						U_scp, C_scp = controller.get_u( x_curr, t)
					u_curr = np.transpose(U_scp[count])

				else:
					u_curr, cost_t = controller.get_u( x_curr,t)

				x_next = x_curr + dynamics.get_dxdt(x_curr,u_curr,t) * param.get('dt')

				X.append(x_curr)
				U.append(u_curr)
				V.append(dynamics.get_V( x_curr,t))

				# temp 
				pbar.append( np.matmul( \
					np.transpose( util.get_my_1()), util.get_p_a(x_curr)))
				vbar.append( np.matmul( \
					np.transpose( util.get_my_1()), util.get_v_a(x_curr)))
				abar.append( np.matmul( \
					np.transpose( util.get_my_1()), dynamics.get_vdot_a(x_curr,t)))

				x_curr = x_next
				count += 1 
				print( '\t' + str(t) + '/' + str(param.get('T')[-1]))

			X = np.squeeze(np.asarray(X))
			U = np.squeeze(np.asarray(U))
			V = np.asarray(V)
			pbar = np.asarray(pbar)
			vbar = np.asarray(vbar)
			abar = np.asarray(abar)

			C.append( dynamics.calc_cost(U))

			print('Plots...')
			plotter.plot_SS(X, param.get('T'), title = str(u) + ' State Space' )
			plotter.plot_V( V, param.get('T'), title = str(u) + ' Lyapunov' )
			plotter.plot_U( U, param.get('T'), title = str(u) + ' Controller')
			plotter.plot_test2( \
				np.squeeze(pbar) - param.get('pd'), \
				np.squeeze(vbar) - param.get('vd'), \
				np.squeeze(abar) - param.get('ad'), \
				param.get('T'), title = str(u) + ' Errors')
			print('Plots Complete')
		except:
			pass

	plotter.plot_cost( C, C_scp)
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