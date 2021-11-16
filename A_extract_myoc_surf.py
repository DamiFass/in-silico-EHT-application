import numpy as np
import os.path, time
import timeit
import sys

def extract(meshname):

	print('Loading mesh...')
	t1 = timeit.default_timer()

	nodes = np.loadtxt(str(meshname) + '.pts', skiprows=1)
	elems = np.loadtxt(str(meshname) + '.elem', skiprows=1, usecols=[1,2,3,4,5],dtype=int)
	fib = np.loadtxt(str(meshname) + '.lon', skiprows=1)

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t1)/60 ))

	# Myocardium tags: 1,4.  Bath tags: 0,3   ---> We want the scar (3) to have fibers:

	# copy the fibres:
	new_fib = np.copy(fib)
	# Put default fibers (1,0,0) where the tag is 3:
	new_fib[np.where(elems[:,4]==3)[0],:] = np.array([1,0,0,1,0,0])

	nmeshname = meshname+'_scar_with_fibers'
	# Copy pts and elem files:
	print('Copying .pts and .elem files... ')
	t2 = timeit.default_timer()

	cmd = "ln -s %s %s"%(meshname+".pts",nmeshname+".pts")
	os.system(cmd)
	cmd = "ln -s %s %s"%(meshname+".elem",nmeshname+".elem")
	os.system(cmd)

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t2)/60 ))

	# Create new fibres file:
	print('Writing new fibers file...')
	t3 = timeit.default_timer()

	f = open(nmeshname+".lon","w")
	f.write("2\n")
	for i in range(new_fib.shape[0]):
		f.write("%f %f %f %f %f %f\n"%(float(new_fib[i,0]),float(new_fib[i,1]),float(new_fib[i,2]),float(new_fib[i,3]),float(new_fib[i,4]),float(new_fib[i,5])))
	f.close()

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	# Extract myocardium from the new mesh:
	print('Extracting myocardium...')
	t3 = timeit.default_timer()

	myoc_meshname = 'myocardium'
	cmd = "meshtool extract myocard -msh={} -ifmt=carp_txt -ofmt=carp_txt -submsh={}".format(nmeshname,myoc_meshname)
	os.system(cmd)

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	# Extract surface from the myocardial mesh:
	print('Extracting myocardial surface...')
	t3 = timeit.default_timer()

	surfname = 'myocardial_surface'
	cmd = "meshtool extract surface -msh={} -ofmt=carp_txt -surf={}".format(myoc_meshname,surfname)
	os.system(cmd)

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	# Compute the elements centres: 
	print('Computing myocardial surface elements center...')
	t3 = timeit.default_timer()

	ctrname = "elemCentres"
	cmd = "GlElemCenters -m {}.surfmesh -o {}".format(surfname,ctrname)
	os.system(cmd)

	print('Done in {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	# Convert myocardial_surface mesh in vtk to visualize:
	cmd = "meshtool convert -ifmt=carp_txt -ofmt=vtk -imsh={}.surfmesh -omsh={}".format(surfname,surfname)
	os.system(cmd)

	cmd = "rm myocardium.*"
	os.system(cmd)

	cmd = "rm {}_scar_with_fibers.*".format(meshname)
	os.system(cmd)

	return


def main():

	t_tot = timeit.default_timer()

	meshname = sys.argv[1]

	extract(meshname)

	print('Everything done in {:.1f} min'.format( (timeit.default_timer() - t_tot)/60 ))
	print('\n')
	print('OPEN THE MESH SURFACE ("myocardial_surface.vtk") IN PARAVIEW AND CHOOSE A PLANE SLIGHTLY BELOW THE BASE (norm vector oriented towards the apex)')

if __name__ == '__main__':
	main()