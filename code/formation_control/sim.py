

from param import param
import dynamics
import utilities as util
import plotter
import os, timeit, subprocess
import numpy as np


def main():

	# random seed
	np.random.seed(88)

	# init
	util.get_x0()
	util.get_xd()

	# preallocate
	X = np.zeros( (param.get('n'), param.get('nt')))

	# go sim 
	x_curr = param.get('x0')
	for k,t in enumerate( param.get('T')):

		u_curr = util.to_vec( np.zeros(param.get('m')))
		x_next = x_curr + \
			dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')		
		X[:,k] = np.squeeze(x_curr)

		x_curr = x_next

	plotter.plot_SS( X, param.get('T'))
	plotter.plot_ss( X[:,-1], param.get('T')[-1])
	plotter.save_figs()

if __name__ == '__main__':

	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))
	if os.path.exists(pdf_path):
		os.remove(pdf_path)
	
	main()

	if os.path.exists(pdf_path):
		subprocess.call(["xdg-open", pdf_path])                    