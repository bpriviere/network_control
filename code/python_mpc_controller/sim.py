
from param import param
import numpy as np 
import dynamics
import controller
import plotter
import utilities as util
import subprocess, os, timeit 

def main():

	util.get_xd()

	X = []
	U = []
	V = []
	LgV = []
	LfV = []

	x_curr = util.get_x0()
	print('Sim...')
	for t in param.get('T'):
		
		u_curr = controller.get_u( x_curr,t)
		x_next = x_curr + dynamics.get_dxdt(x_curr,u_curr,t) * param.get('dt')

		X.append(x_curr)
		U.append(u_curr)
		V.append(dynamics.get_V( x_curr,t))

		# temp 
		LgV.append( dynamics.get_LgV(x_curr, t))
		LfV.append( dynamics.get_LfV(x_curr, t))

		x_curr = x_next
		print( '\t' + str(t) + '/' + str(param.get('T')[-1]))

		return

	X = np.squeeze(np.asarray(X))
	U = np.squeeze(np.asarray(U))
	V = np.asarray(V)
	LgV = np.asarray(LgV)
	LfV = np.asarray(LfV)

	print('Plots...')
	plotter.plot_SS(X,param.get('T'))
	plotter.plot_V( V, param.get('T'))
	plotter.plot_U( U,param.get('T'))
	plotter.plot_test( V, LgV, LfV, U, param.get('T'))
	# plotter.show()
	plotter.save_figs()
	print('Plots Complete')

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