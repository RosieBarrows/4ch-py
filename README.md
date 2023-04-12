# 4ch-py


*This section of the pipeline takes the mesh as an input and returns a folder ready to be copied to a machine for simulation. It extracts all the necessary surfaces, calculates UVCs on the ventricles and atria, adds ventricular and atrial fibres and prepares the mesh for simulation.*


Set up a conda environment with the packages listed in README_conda_setup.md 

## Quick start 
Set up your anaconda environment with `conda create -n 4ch python=3.9 -y`.
Then install the requirements using `pip`: 
```sh
pip install -r requirements.txt
```
Confirm everything has been installed. First, activate the environment `conda activate 4ch`, then type: 
```sh
python docker/check_installs.py 
```

Deactivate your base environment and activate your conda environment. 


*You should already have a heart folder (i.e. example_mesh) with the following structure:*

example_mesh:
		- meshing
				- myocardium_OUT
						- myocardium.elem
						- myocardium.pts
						- myocardium.lon
		- segmentations (not needed but might be present from previous steps)


*If the folder containing your mesh is not called myocardium_OUT (or your mesh is not called myocardium), you will need to change this parameter within the file main_surf.py.*


Within the parfiles folder (within this scripts folder) there is a parfiles folder. Change the name within heartFolder.txt to the path of your heart folder i.e. /data/Dropbox/example heart. 

**Check that the tags in tags_vent_fibres.json match those of your mesh and change if necessary**

Run:  
   `bash 0_extract_surfs.sh in terminal`

Use the la.vtk (in surfaces_uvc_LA/la/) to select a point for the LA apex and a point for the LA septum. Save the Point IDs in la.lvapex.vtx and la.rvspet_pt.vtx, respectively (blank files provided in surfaces_uvc_LA/la/).

Repeat this step for the RA using ra.vtk.

Finally, select the apex of the RAA (right atrial appendage) and save the 3 coordinates of this point (with spaces in between each coordinate) in ${heart_folder}/raa_apex.vtx

Within the ./parfiles/etags/ folder:

1. Open the etags.sh file. Change the value of T_LV to be the LV label in your mesh and change the value of T_RV to be the RV label in your mesh. 

2. Open the etags_la.sh file. Change the value of T_LV to be the !! LA !! label in your mesh.

3. Open the etags_ra.sh file. Change the value of T_LV to be the !! RA !! label in your mesh.

**MAKE SURE YOU HAVE A VALID CARP LICENSE (used for mguvc command)**
**NOTE - YOUR VERSION OF CARP MUST BE THE VERSION THAT INCLUDES A --CUSTOM-APEX FLAG FOR THE UVCS**

Run:  
   `bash 1_calculate_UVCs.sh`

In the bash 2_add_vent_fibres.sh file, **check that the four_chamber_name parameter matches the name of your four chamber mesh.**
Also **check that the path to your CARP_FOLDER is correct.** 

Run:  
   `bash 2_add_vent_fibres.sh`

**Check that the labels in tags_atrial_fibres.json match the labels in your mesh.**

Run:  
   `bash 3_create_surfaces_endo_landmarks.sh`  
   `bash 4_la_4ch_endo.sh`  
   `bash 5_ra_4ch_endo.sh`  
   `bash 6_map_2d_to_3d.sh`

Use the vtk (myocardium_fibres_l.vtk in the atrial_fibres folder) to check that the ventricular and atrial fibres appear as expected. 

**Check that the labels in tags_presim.json match the labels in your mesh**

*If you would like to change the location of the Bachmann bundle, this can be done in ./parfiles/bachmann_bundle_fec_settings.json*

Run:  
   `bash 7_define_tags.sh`  
   `bash 8_extract_surfs.sh`

If you need different labels for the fast endocardial conduction zones in the LV and RV (because you need to assign distinct material/conduction parameters in these regions), run:  
   `bash 9_split_fec.sh`

To complete the simulation folder, use the vtk (myocardium_AV_FEC_BB.vtk or myocardium_AV_FEC_BB_lvrv.vtk) to select a point for the apex and the sinoatrial node. Save these within the simulation folder as myocardium_AV_FEC_BB_apex.vtx and myocardium_AV_FEC_BB_SA.vtx, respectively. 
