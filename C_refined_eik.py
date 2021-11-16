import numpy as np
import os.path, time
import timeit
from scipy.spatial import distance
import sys

def second_eik(meshname, eik_ATs, pthick):

	# Interpolate the eikonal ATs to the elements: 
	eik_ATs_elem = eik_ATs+'.elems.dat'
	cmd = "meshtool interpolate node2elem -omsh={} -idat={} -odat={}".format(meshname,eik_ATs,eik_ATs_elem)
	os.system(cmd)

	print('Loading eikonal solution on elements...')
	# Load EK solution interpolated on elements
	epibath = np.loadtxt(eik_ATs_elem)

	print('Loading big mesh elements...')
	elems = np.loadtxt(str(meshname) + '.elem', skiprows=1, usecols=[1,2,3,4,5],dtype=int)

	# Separate nodes indexes and element tags:
	elementi = np.copy(elems[:,0:4])
	elemtags = np.copy(elems[:,4])
	nelems = elementi.shape[0]

	# Put a different tag (7) on the elements for which the eikonal solution is >-1 (i.e. towards outer bath) and <= a certain thicnkess., 
	# AND the tag is 0! (Otherwise it will retag also the myocardial nodes!)
	elemtags[np.where(elemtags==7)[0]] = 0  
	epilayer = np.logical_and(np.logical_and(epibath>-1,epibath<=pthick),elemtags==0)  
	elemtags[epilayer] = 7

	print('Saving elements with new tag...')
	# Save elements with new tag:
	outputname = meshname + "_epilayer"

	f = open(outputname+".elem","w")
	f.write("%d\n" % nelems)
	for i in range(nelems):
		f.write("Tt %d %d %d %d %d\n"%(int(elems[i,0]),int(elems[i,1]),int(elems[i,2]),int(elems[i,3]),int(elemtags[i])))
	f.close()

	# Copy points and lon files with the new meshname:
	cmd = "cp {}.pts {}.pts".format(meshname,outputname)
	os.system(cmd)
	cmd = "cp {}.lon {}.lon".format(meshname,outputname)
	os.system(cmd)

	# refine epi layer
	ref_nmeshname = outputname + "_refined"
	cmd = "meshtool resample mesh -msh={} -tags=7 -ifmt=carp_txt -ofmt=carp_txt -avrg=300 -postsmth=0 -outmsh={}".format(outputname,ref_nmeshname)
	os.system(cmd)

	# Now we have a different tag just outised the epicardium, so we can extract the surface easily, without splitting the plane, etc.. like done before.
	# extract new epi surface, from the mesh: PXYZ-epilayer-refined.*
	epi_name = ref_nmeshname+'_stimulus'
	cmd = "meshtool extract surface -msh={} -surf={} -op=1,3,4,5:7 -ofmt=carp_txt".format(ref_nmeshname,epi_name)
	os.system(cmd)

	# generate fibers everywhere for the eikonal simulation:
	cmd = "meshtool generate fibres -msh={} -outmsh={}".format(ref_nmeshname, ref_nmeshname)
	os.system(cmd)

	cmd = "rm {}_epilayer.*".format(meshname)
	os.system(cmd)

	return ref_nmeshname, epi_name

def main():

	t_tot = timeit.default_timer()

	meshname = sys.argv[1]
	eik_ATs = 'vm_act_seq.dat'
	# Define how thick is the layer we want to refine (from epicadium, outward)
	pthick = 3. # mm 

	ref_nmeshname, epi_name = second_eik(meshname, eik_ATs, pthick)

	print('Everything done in {:.1f} min'.format( (timeit.default_timer() - t_tot)/60 ))
	print('\n')
	print('USE THE FILES {}.pts/.elem/.lon AND THE STIMULUS FILE {}.vtx TO RUN THE SECOND EIKONAL SIMULATION!'.format(ref_nmeshname, ref_nmeshname, ref_nmeshname, epi_name))

	# scp -r P***_epilayer_refined.elem P***_epilayer_refined.pts P***_epilayer_refined.lon  dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/
	# scp -r P***_epilayer_refined_stimulus.surf.vtx  dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/
	## AFTER THE SIMULATION:
	# scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/P***_epilayer_refined/vm_act_seq.dat /home/dfa18/Desktop/D/CREATE_PATCHED_MESHES/


if __name__ == '__main__':
	main()
