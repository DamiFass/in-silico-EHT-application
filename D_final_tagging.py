import numpy as np
import os.path, time
import timeit
from scipy.spatial import distance
import sys

def tag_EHT_CP(meshname,eik_ATs,rscaling,pthick,ehtthick):
	
	ref_nmeshname = meshname+'_epilayer_refined'
	# Interpolate again, but remember to use the epilayer-refined mesh: 
	eik_ATs_elem = eik_ATs+'.elems.dat'
	cmd = "meshtool interpolate node2elem -omsh={} -idat={} -odat={}".format(ref_nmeshname,eik_ATs,eik_ATs_elem)
	os.system(cmd)

	print('Loading eikonal solution on elements...')
	# Load EK interpolated solution again:
	epibath = np.loadtxt(eik_ATs_elem)

	# compute element centres of the epilayer-refined mesh:
	ctrname = ref_nmeshname + ".centres.vpts"
	cmd = "GlElemCenters -m {} -o {}".format(ref_nmeshname,ctrname)
	os.system(cmd)

	print('Loading elements centers...')
	# Load element centres 
	elemctr = np.loadtxt(ctrname,skiprows=1)

	print('Loading refined mesh elements...')
	# Load elements of the epilayer-refined mesh:
	elems = np.loadtxt(str(ref_nmeshname) + '.elem', skiprows=1, usecols=[1,2,3,4,5],dtype=int)

	# Separate nodes indexes and element tags:
	elementi = np.copy(elems[:,0:4])
	elemtags = np.copy(elems[:,4])
	nelems = elementi.shape[0]

	print('Getting epicardial scar elements...')
	# GET EPI SCAR ELEMENTS:
	# nodes index of the epi surface
	episurf = ref_nmeshname + "_stimulus.surf"
	epivtx = np.loadtxt(episurf+".vtx",skiprows=2,dtype=int)   
	# True only on the rows where tag is 3 or 4:
	scar = np.logical_or(elemtags==3,elemtags==4)  
	# True only on the rows of elements that have at least a node belonging to the epicardial surface
	epi = np.any(np.isin(elems,epivtx),axis=1)   
	# True only on the rows of elements that have an epicardial node AND a scar node
	episcar = np.logical_and(scar,epi)  

	print('Computing elements distances from scar center...')
	# GET ELEMENTS DISTANCE FROM SCAR CENTRE:
	# Mean of the coordinates of the centers of the elements that have an epicardial AND a scar node:
	scentroid = np.mean(elemctr[episcar,:],axis=0)
	# Store the distancese of each element centre from the scar centroid:
	elemdist = np.sqrt(np.sum((elemctr-scentroid)**2, axis=1))

	print('Computing patch radius...')
	# COMPUTE PATCH RADIUS:
	# Get the most right scar element centroid:
	stopright = np.array([np.amax(elemctr[episcar,:],axis=0)])
	# Get the most left scar element centroid:
	sbottomleft = np.array([np.amin(elemctr[episcar,:],axis=0)])
	# Compute the distances between the farthest elements, and divide to get the radius.
	sradius = np.sqrt(np.sum((stopright-sbottomleft)**2, axis=1))/2.
	# Scale the radius based on how much I want the patch to be covering the scar. 
	sradius = sradius*rscaling

	print('Labelling EHT and CP...')
	#LABEL EHT AND CP ELEMENTS:
	# True only on the rows of elements whose eikonal solution is smaller than thickness, distance is less than radius, and were labelled as epilayer from the first eikonal simulation:
	cp = np.logical_and(np.logical_and(np.logical_and(epibath>-1,epibath<=pthick),elemdist<=sradius),elemtags==7)
	# True only on the rows of elements whose eikonal solution is smaller than eht thickness, distance is less than radius, and were labelled as epilayer from the first eikonal simulation:
	eht = np.logical_and(np.logical_and(np.logical_and(epibath>-1,epibath<=ehtthick),elemdist<=sradius*0.95),elemtags==7) # 0.95: CP must envelop EHT
	# cp contains True values also for EHT elements, but then I am assigning different tags to the EHT elements later, so it does not matter:
	elemtags[cp] = 6 # CP
	elemtags[eht] = 5 # EHT

	print('Saving elements with new tags...')
	# Save elements with new tag:
	cmd = "rm {}.elem".format(ref_nmeshname)
	os.system(cmd)
	f = open(ref_nmeshname+".elem","w")
	f.write("%d\n" % nelems)
	for i in range(nelems):
		f.write("Tt %d %d %d %d %d\n"%(int(elems[i,0]),int(elems[i,1]),int(elems[i,2]),int(elems[i,3]),int(elemtags[i])))
	f.close()

	# Now the fibers have [1,0,0] in all the elements! So I need to get back the correct fibers for the myocardium (tags 1,3,4), set [1,0,0] for EHT (tag 5) and 0 everywhere else.

	# Extract submesh with tags 1,3,4 from the final mesh (because extract myocardium does not work, because the fibers are all [1,0,0])
	submeshname = ref_nmeshname+"_myoc"
	cmd = "meshtool extract mesh -msh={} -ifmt=carp_txt -ofmt=carp_txt -tags=1,3,4 -submsh={}".format(ref_nmeshname, submeshname)
	os.system(cmd)
	# Run Interp_fibers_final.py
	cmd = "python3 Interp_fibers_final.py {} {}".format(meshname, submeshname)
	os.system(cmd)
	# meshtool insert submesh, con il nuovo .lon file
	finalname = meshname+'-final-{}'.format(rscaling)
	cmd = "meshtool insert submesh -submsh={} -msh={} -ofmt=carp_txt -outmsh={}".format(submeshname, ref_nmeshname, finalname)
	os.system(cmd)

	# Now the fibers are correct in the myocardium, and [1,0,0] everywhere else: 
	print('Loading last mesh fibers...')
	fib_ref = np.loadtxt(str(finalname) + '.lon', skiprows=1)
	nelems = fib_ref.shape[0]
	print('Setting bath elements fibers to 0...')
	# Set fibers (1,0,0) in the EHT elements, so the elements tagged with 5.
	fib_ref[np.where(elemtags==5)[0],0] = 1
	# Set fibers to 0 in the CP:
	fib_ref[np.where(elemtags==6)[0],:] = 0
	# Set fibers to 0 in the scar:
	fib_ref[np.where(elemtags==3)[0],:] = 0
	# Set fibers to 0 in the refined layer:
	fib_ref[np.where(elemtags==7)[0],:] = 0
	# Set fibers to 0 in the bath:
	fib_ref[np.where(elemtags==0)[0],:] = 0
	# # Set fibers to 0 in the border zone: 
	# fib_ref[np.where(elemtags==4)[0],:] = 0 # Border zones stay in the tissue!

	print('Saving new fibers...')
	lonname = finalname + ".lon"
	fs = open(lonname,"w")
	fs.write("2\n")
	for i in range(nelems):
		fs.write("%f %f %f %f %f %f\n"%(fib_ref[i,0],fib_ref[i,1],fib_ref[i,2],fib_ref[i,3],fib_ref[i,4],fib_ref[i,5]))
	fs.close()

	return

def main():

	t_tot = timeit.default_timer()

	meshname = sys.argv[1]
	eik_ATs = 'vm_act_seq.dat'
	# Define th extension of the patch: 0 to 1 (1=on all the scar):
	rscaling = 0.9
	# pthick is the thickness of both! CP and EHT summed up together! CAN'T BE MORE THAN 3mm!
	pthick = 2. # mm --> (i.e. 2mm, could be 1 + 1 or 0.5mm for EHT and 1.5 for CP) 
	# Thickness of the EHT: 
	ehtthick = 1.3 # mm --> thickness of CP will be pthick - ehtthick 

	tag_EHT_CP(meshname, eik_ATs,rscaling,pthick,ehtthick)

	print('Everything done in {:.1f} min'.format( (timeit.default_timer() - t_tot)/60 ))
	print('\n')
	print('DONE! THE FINAL MESH SHOULD BE NAMED: {}-final-{}.pts/elem/lon'. format(meshname,rscaling))

	# Convert in vtk to check! meshtool convert -imsh=P***-final-** -ifmt=carp_txt -omsh=P***-final-** -ofmt=vtk

	# Maybe increase the quality!!
	# meshtool clean quality -msh=P100-final-1_lowq -thr=0.5 -surf=P100_epilayer_refined_stimulus  -outmsh=P100-final-1_preserved_surf

if __name__ == '__main__':
	main()


############ Now I can use the mesh to run simulation, but I have to define the stimulus area first!: 

# Send the mesh on tom2-login:

################ ON TOM2 ##################

###Extract myocardium to define the stimulus points and visualize vm later:
#meshtool extract myocard -msh=$meshname -ifmt=carp_txt -submsh=$meshname'_i' -ofmt=carp_txt

#Convert extracted myocardial mesh in vtk:
#meshtool convert -imsh=$meshname -ifmt=carp_txt -omsh=$meshname -ofmt=vtk

########### BACK ON DESKTOP ############

#scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/P***_i.vtk /home/dfa18/Desktop/

#### Open with Paraview to select the points. Save them as csv file. 

#Create the vtx file for the stimulus:
#python3 Write_stimolo_vtx.py stimolo.csv

#Send the vtx file to TOM2:
#scp -r Stim.vtx dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/

#Map the stim.vtx file back onto the big mesh:
#meshtool map -submsh=$meshname'_i' -files=Stim.vtx -outdir=/scratch/dfa18/ -mode=s2m

##### RUN THE SIMULATION ON TOM2 ######

#Copy on Desktop the myocardial mesh in carp format (or convert the vtk one already on Desktop)
#scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/P166_i.elem (and .pts and .lon) /home/dfa18/Desktop/

#Copy on Desktop vm.igb file to visualize with meshalyzer
#scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/P166/vm.igb /home/dfa18/Desktop/



# It works with:
# - with stimulus defined with volume and not with the .vtx file
