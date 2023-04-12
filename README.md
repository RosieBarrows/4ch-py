# 4ch-py

-------------------------------------------------------------------------------------
This section of the pipeline takes the mesh as an input and returns a folder ready to be copied to a machine for simulation. It extracts all the necessary surfaces, calculates UVCs on the ventricles and atria, adds ventricular and atrial fibres and prepares the mesh for simulation. 
-------------------------------------------------------------------------------------

Set up a conda environment with the packages listed in README_conda_setup.md 

Deactivate your base environment and activate your conda environment. 

-------------------------------------------------------------------------------------
You should already have a heart folder (i.e. example_mesh) with the following structure:

example_mesh:
		- meshing
				- myocardium_OUT
						- myocardium.elem
						- myocardium.pts
						- myocardium.lon
		- segmentations (not needed but might be present from previous steps)


If the folder containing your mesh is not called myocardium_OUT (or your mesh is not called myocardium), you will need to change this parameter within the file main_surf.py. 
-------------------------------------------------------------------------------------

Within the parfiles folder (within this scripts folder) there is a parfiles folder. Change the name in heartFolder.txt to the name of your heart folder i.e. example heart. 

Check that the tags in tags_vent_fibres.json match those of your mesh and change if necessary. 

Run bash 0_extract_surfs.sh in terminal

Use the la.vtk (in surfaces_uvc_LA/la/) to select a point for the LA apex and a point for the LA septum. Save the Point IDs in la.lvapex.vtx and la.rvspet_pt.vtx, respectively (blank files provided in surfaces_uvc_LA/la/).

Repeat this step for the RA using ra.vtk. 

Within the ./parfiles/etags/ folder:
	- Open the etags.sh file. Change the value of T_LV to be the LV label in your mesh and change the value of T_RV to be the RV label in your mesh. 

	- Open the etags_la.sh file. Change the value of T_LV to be the !! LA !! label in your mesh.

	- Open the etags_ra.sh file. Change the value of T_LV to be the !! RA !! label in your mesh.

-------------------------------------------------------------------------------------
MAKE SURE YOU HAVE A VALID CARP LICENSE (used for mguvc command)
-------------------------------------------------------------------------------------

Run bash 1_calculate_UVCs.sh

In the bash 2_add_vent_fibres.sh file, check that the four_chamber_name parameter matches the name of your four chamber mesh.
Also check that the path to your CARP_FOLDER is correct. 




