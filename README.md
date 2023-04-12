# 4ch-py

This section of the pipeline takes the mesh as an input and returns a folder ready to be copied to a machine for simulation. It extracts all the necessary surfaces, calculates UVCs on the ventricles and atria, adds ventricular and atrial fibres and prepares the mesh for simulation. 

Set up a conda environment with the packages listed in README_conda_setup.md 

Deactivate your base environment and activate your conda environment. 

You should already have a heart folder (i.e. example_mesh) with the following structure:

example_mesh:
		- meshing
				- myocardium_OUT
						- myocardium.elem
						- myocardium.pts
						- myocardium.lon
		- segmentations (not needed but might be present from previous steps)

If the folder containing your mesh is not called myocardium_OUT (or your mesh is not called myocardium), you will need to change this parameter within the file main_surf.py. 

Within the parfiles folder (within this scripts folder) there is a parfiles folder. Change the name in heartFolder.txt to the name of your heart folder i.e. example heart. 

Check that the tags in tags_vent_fibres.json match those of your mesh and change if necessary. 


Run bash 0_extract_surfs.sh in terminal

Use the la.vtk to select a 


