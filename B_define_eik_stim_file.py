import numpy as np
import os.path, time
import timeit
from numba import jit
from scipy.spatial import distance
import sys


def create_stim_file(origin, normal, meshname):

	print('Loading centers...')
	ctrname = "elemCentres"
	elemcentr = np.loadtxt(str(ctrname), skiprows=1)
	# Compute dot product:
	dotproducts = []
	for i in range(elemcentr.shape[0]):
		tmp = np.copy(elemcentr[i,:])
		v = tmp - origin
		dotproducts.append(np.dot(v,normal))

	print('Loading myocardial surface...')
	surfname = 'myocardial_surface'
	# Tag differently elements above the surface (i.e. negative dot product):
	surf_elem = np.loadtxt(str(surfname) + '.surfmesh.elem', skiprows=1, usecols=[1,2,3,4])
	surf_elem[np.where(np.array(dotproducts)<0),3] = 8

	print('Saving the new elem file...')
	# Save the elem file:
	f = open(surfname+".surfmesh.elem","w")
	f.write("%d\n" % surf_elem.shape[0])
	for i in range(surf_elem.shape[0]):
		f.write("Tr %d %d %d %d\n"%(int(surf_elem[i,0]),int(surf_elem[i,1]),int(surf_elem[i,2]),int(surf_elem[i,3])))
	f.close()

	# Eliminate the base:
	no_base_name = 'No_base'
	cmd = "meshtool extract mesh -msh={}.surfmesh -ofmt=carp_txt -tags=1,3,4 -submsh={}".format(surfname,no_base_name)
	os.system(cmd)

	# Create the base:
	base_name = 'BASE'
	cmd = "meshtool extract mesh -msh={}.surfmesh -ofmt=carp_txt -tags=8 -submsh={}".format(surfname,base_name)
	os.system(cmd)

	# Divide mesh to get endo and epicardial surface:
	layers = "layers"
	cmd = "meshtool extract unreachable -msh={} -ifmt=carp_txt -submsh={} -ofmt=carp_txt".format(no_base_name, layers)
	os.system(cmd)

	print('Loading epicardial mesh points...')
	# Load points of the epicardial mesh:
	nodes_epi = np.loadtxt(str(layers)+'.part1.pts', skiprows=1)

	print('Loading endocardial mesh points...')
	# Load points of the epicardial mesh:
	nodes_endo = np.loadtxt(str(layers)+'.part0.pts', skiprows=1)

	print('Loading BASE mesh points...')
	# Load points of the epicardial mesh:
	nodes_base = np.loadtxt(str(base_name)+'.pts', skiprows=1)

	@jit(nopython=True, parallel=True)
	def loop(nodes, xcoord, ycoord, zcoord):
		giusta = []
		not_found = []
		for i in range(nodes.shape[0]):
			tmp_array = np.intersect1d( np.intersect1d( np.where(xcoord==nodes[i,0])[0], np.where(ycoord==nodes[i,1])[0] ), np.where(zcoord==nodes[i,2])[0] )
			if tmp_array.shape[0] > 0:
				giusta.append(tmp_array[0])
			else:
				not_found.append(i)
		return giusta, not_found
		
	print('Loading refined mesh points...')
	nodes_big = np.loadtxt(meshname+'.pts', skiprows=1)

	xcoord = np.copy(nodes_big[:,0])
	ycoord = np.copy(nodes_big[:,1])
	zcoord = np.copy(nodes_big[:,2])

	t3 = timeit.default_timer()
	print('Numba for loop for epicardial nodes...')
	gg, fail = loop(nodes_epi, xcoord, ycoord, zcoord)

	print('Loop done in: {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	t3 = timeit.default_timer()
	print('Numba for loop for endocardial nodes...')
	gg_endo, fail_endo = loop(nodes_endo, xcoord, ycoord, zcoord)

	print('Loop done in: {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	t3 = timeit.default_timer()
	print('Numba for loop for base nodes...')
	gg_base, fail_base = loop(nodes_base, xcoord, ycoord, zcoord)

	print('Loop done in: {:.1f} min'.format( (timeit.default_timer() - t3)/60 ))

	new = []
	for i in fail:
		new.append( np.intersect1d( np.where(xcoord==nodes_epi[i,0])[0], np.where(ycoord==nodes_epi[i,1])[0] )[0] )

	new_endo = []
	for i in fail_endo:
		new_endo.append( np.intersect1d( np.where(xcoord==nodes_endo[i,0])[0], np.where(ycoord==nodes_endo[i,1])[0] )[0] )

	new_base = []
	for i in fail_base:
		new_base.append( np.intersect1d( np.where(xcoord==nodes_base[i,0])[0], np.where(ycoord==nodes_base[i,1])[0] )[0] )

	# gg_giusto = gg + new 

	# gg_endo_giusto = gg_endo + new_endo

	# gg_base_giusto = gg_base + new_base

	# Save the epicardial nodes as vtx file:
	stim_file_name= meshname+'_EPI'
	with open(stim_file_name+'.vtx', 'w') as f:
		f.write('{}\n'.format(len(gg)))
		f.write('extra\n')
		np.savetxt(f, gg, fmt='%.0f')

	# Save the endocardial nodes as vtx file:
	stim_file_name_endo= meshname+'_ENDO'
	with open(stim_file_name_endo+'.vtx', 'w') as f:
		f.write('{}\n'.format(len(gg_endo)))
		f.write('extra\n')
		np.savetxt(f, gg_endo, fmt='%.0f')

	# Save the base nodes as vtx file:
	stim_file_name_base= meshname+'_BASE'
	with open(stim_file_name_base+'.vtx', 'w') as f:
		f.write('{}\n'.format(len(gg_base)))
		f.write('extra\n')
		np.savetxt(f, gg_base, fmt='%.0f')

	eik_mesh = meshname+'_eik'
	# Copy pts and elem files:
	cmd = "ln -s %s %s"%(meshname+".pts",eik_mesh+".pts")
	os.system(cmd)
	cmd = "ln -s %s %s"%(meshname+".elem",eik_mesh+".elem")
	os.system(cmd)

	# Put the fibres in ALL elements, using meshtool generate fibres.
	cmd = "meshtool generate fibres -msh={} -outmsh={}".format(meshname,eik_mesh)
	os.system(cmd)

	### Generate healthy mesh:

	print('Loading mesh elements to retag and generate healthy mesh...')
	elems = np.loadtxt(str(meshname) + '.elem', skiprows=1, usecols=[1,2,3,4,5],dtype=int)
	elemtags = np.copy(elems[:,4])
	# Put tag 1 where the tag was 3 (scar):
	elems[np.where(elemtags==3)[0],4] = 1
	# Put tag 1 where the tag was 4 (BZ):
	elems[np.where(elemtags==4)[0],4] = 1

	print('Saving elements with new tags...')
	# Save elements with new tag:
	outputname = meshname + "_HEALTHY"

	f = open(outputname+".elem","w")
	f.write("%d\n" % elems.shape[0])
	for i in range(elems.shape[0]):
		f.write("Tt %d %d %d %d %d\n"%(int(elems[i,0]),int(elems[i,1]),int(elems[i,2]),int(elems[i,3]),int(elems[i,4])))
	f.close()

	# Copy pts and elem files:
	cmd = "cp %s %s"%(meshname+".pts",outputname+".pts")
	os.system(cmd)
	cmd = "cp %s %s"%(meshname+".lon",outputname+".lon")
	os.system(cmd)

	# Clean unsued files from previous script: 
	cmd = "rm layers.part*"
	os.system(cmd)
	cmd = "rm myocardi*"
	os.system(cmd)
	cmd = "rm No_base.*"
	os.system(cmd)

	return


def main():

	t_tot = timeit.default_timer()

	a = float(sys.argv[1])
	b = float(sys.argv[2])
	c = float(sys.argv[3])
	d = float(sys.argv[4])
	e = float(sys.argv[5])
	f = float(sys.argv[6])

	meshname = sys.argv[7]

	origin = np.array([a,b,c]) 
	normal = np.array([d,e,f])

	create_stim_file(origin, normal, meshname)

	eik_mesh = meshname+'_eik'
	stim_file_name= meshname+'EPI'
	print('Everything done in {:.1f} min'.format( (timeit.default_timer() - t_tot)/60 ))
	print('\n')
	print('USE THE FILES {}_eik.pts, {}_eik.elem, {}_eik.lon AND THE STIMULUS FILE {}.vtx TO RUN THE EIKONAL SIMULATION!'.format(eik_mesh, eik_mesh, eik_mesh, stim_file_name))

	# scp -r $meshname'_eik.*' dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/
	# scp -r $stim_file_name'.vtx' dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/
	## AFTER THE SIMULATION:
	# scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/Pxyz_eik/vm_act_seq.dat /home/dfa18/Desktop/D/CREATE_PATCHED_MESHES
	
	## TO VISUALIZE AND CHECK THE EIKONAL ACTIVATION:
	# meshtool collect -imsh=Pxyx_eik -ifmt=carp_txt -omsh=CHECK_FIRS_EIKONAL -nod=vm_act_seq.dat -ofmt=vtk

if __name__ == '__main__':
	main()