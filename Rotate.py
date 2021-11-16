import numpy as np
import array
import sys
import subprocess


def rotation_Z(meshname, theta):

	R = np.array([[np.cos(theta),-np.sin(theta),0], [np.sin(theta),np.cos(theta),0], [0,0,1]])

	nodes = np.loadtxt(str(meshname) + '.pts', skiprows=1)
	elems = np.loadtxt(str(meshname) + '.elem', skiprows=1, usecols=[1,2,3,4])
	fib = np.loadtxt(str(meshname) + '.lon', skiprows=1)

	nodes_rotated = np.zeros((nodes.shape[0], nodes.shape[1]))
	for i in range(nodes.shape[0]):
		nodes_rotated[i,:] = np.dot( R,  np.atleast_2d(nodes[i,:]).transpose() ).transpose()

	nodes_nuovi = np.zeros((nodes.shape[0], nodes.shape[1]))
	nodes_nuovi[:,0] = nodes_rotated[:,0] + (np.abs(np.mean(nodes[:,0])-np.mean(nodes_rotated[:,0])))
	nodes_nuovi[:,1] = nodes_rotated[:,1] - (np.abs(np.mean(nodes[:,1])-np.mean(nodes_rotated[:,1])))
	nodes_nuovi[:,2] = nodes_rotated[:,2] 

	np.savetxt('bbox.pts', nodes_nuovi, fmt='%.6f' )

	com5 = "sed  -i '1i {}' bbox.pts".format(nodes_nuovi.shape[0])
	subprocess.check_call(com5, shell = True, universal_newlines = True)

	com6 = 'meshtool convert -imsh=bbox -ifmt=carp_txt -omsh=bbox_rotated'+ str(theta) +' -ofmt=vtk'
	subprocess.check_call(com6, shell = True, universal_newlines = True) 

	return 

def main():

	meshname = sys.argv[1]
	theta=73

	rotation_Z(meshname,theta)

if __name__ == '__main__':
	main()

