#!/bin/bash

ventricle=$1

meshname=${ventricle:0:4}
 
#Create bbox triangular mesh:
meshtool generate bboxmesh -msh=$ventricle -outmsh=bbox_triangle -ofmt=carp_txt -scale=1.25

#Create bbox tetrahedral mesh:
meshtool generate mesh -surf=bbox_triangle -outmsh=bbox_surface -ifmt=carp_txt -ofmt=carp_txt

#Resample the bbox to 2mm edge average:
meshtool resample mesh -msh=bbox_surface -ifmt=carp_txt -avrg=2000 -ofmt=carp_txt -outmsh=bbox_surface_2mm

#Copy the elements file and the lon file:
cp bbox_surface_2mm.elem bbox.elem
cp bbox_surface_2mm.lon bbox.lon

#Rotate the bbox tetrahedral mesh:
python3 Rotate.py bbox_surface_2mm

#Extract surface from the rotated bbox:
meshtool extract surface -msh=bbox_rotated73.vtk -ifmt=vtk -ofmt=vtk -surf=bbox

#Extract surfaces from the ventricle mesh:
meshtool extract surface -msh=$ventricle -ifmt=vtk -op=1 -surf=one -ofmt=vtk
meshtool extract surface -msh=$ventricle -ifmt=vtk -op=3 -surf=three -ofmt=vtk
meshtool extract surface -msh=$ventricle -ifmt=vtk -op=4 -surf=four -ofmt=vtk

#Generate mesh using all the surfaces: 
meshtool generate mesh -surf=bbox.surfmesh,one.surfmesh,three.surfmesh,four.surfmesh -ifmt=vtk -ofmt=carp_txt -ins_tag=0,1,3,4 -outmsh=Merged

#Convert ventricle in carp_txt format:
meshtool convert -imsh=$ventricle -ifmt=vtk -omsh=v -ofmt=carp_txt

#Add fibers: (modified to set the fibers to 0 where the taf is 0)
python3 Interp_fibers.py v 

#Convert merged mesh with new fibers in vtk:
meshtool convert -imsh=Merged -ifmt=carp_txt -omsh=$meshname'_coarse' -ofmt=vtk

#Send to TOM2:
scp -r $meshname'_coarse.vtk' dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/

# Check that myocardium=1,4 and bath=0,3
meshtool query tags -msh=$meshname'_coarse.vtk' -ifmt=vtk

# Remove files we don't need:
rm bbox* 
rm col_* 
rm four.* 
rm new_col_*
rm one.*
rm three.*
rm v.*
rm Merged.*
rm $meshname'_coarse.vtk'
rm $meshname'-0.8mm-vol-final-labeled.fcon'

################ ON TOM2 ##################

###Resample the myocardium (tags 1,4) to average edge lenght 300um:
#meshtool resample mesh -msh=$meshname'_coarse' -ifmt=vtk -tags=1,4 -avrg=300 -ofmt=carp_txt -outmsh=/scratch/dfa18/$meshname

################ BACK ON DESKTOP ##################

### Bring the resampled mesh (in carp_txt format) back on DESKTOP:
#scp -r dfa18@tom2-login.hpc.isd.kcl.ac.uk:/scratch/dfa18/Pxyz_.* /home/dfa18/Desktop/D/CREATE_PATCHED_MESHES/


######## Now it's ready to be modified for the eikonal simulation and add the EHT patch ########### ---> Mod_mesh_for_eik.py P***



