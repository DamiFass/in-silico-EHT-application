import sys
import numpy as np
import subprocess


def interp(initial,final):

	print('Loading initial fibers...')
	v_fib = np.loadtxt(initial + '.lon', skiprows=1)

	C = {}

	for i in range(v_fib.shape[1]):
		tmp = v_fib[:,i]
		np.savetxt('col_{}.dat'.format(i), tmp, fmt='%.6f')
		cm = 'meshtool interpolate elemdata -omsh={} -imsh={} -idat=col_{}.dat -odat=new_col_{}'.format(final,initial,i,i)
		subprocess.check_call(cm, shell = True, universal_newlines = True)
		C[str(i)]=np.loadtxt('new_col_{}'.format(i))
		if i == 0:
			new_fib = np.zeros((C[str(i)].shape[0],6))
		new_fib[:,i] = np.copy(C[str(i)])

	# Set to 0 the rows whose elements have tag 0:
	# Merged_elems = np.loadtxt('P100-final_i.elem', skiprows=1, usecols=[1,2,3,4,5])

	# new_fib[np.where(Merged_elems[:,4]==0)[0],:] = 0

	print('Saving final fibers...')
	np.savetxt(final+'.lon', new_fib, fmt='%.6f')

	cm = "sed  -i '1i 2' {}.lon".format(final)
	subprocess.check_call(cm, shell = True, universal_newlines = True)


def main():

	initial = sys.argv[1]  #P100
	final = sys.argv[2]    #P100_epilayer_refined_myoc.

	interp(initial,final)

if __name__ == '__main__':
	main()
