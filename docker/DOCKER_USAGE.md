# Usage of docker container in conjuction with CARP 
The docker container `cermg/4ch:TAG` cannot be used to run the entire pipeline. 
The user **must** have CARP installed and be willing to use it. 

This document is a translation of how the container can be used in conjunction with CemrgApp
or via the command line. The following subsections correspond to the steps in order that need to happen
to run the 4ch pipeline. 

## General notes
+ `BASE_DIR` is a folder where the entire project is saved. This has a pre-defined structure and naming conventions.
+ We assume `BASE_DIR` will be passed to docker as a `--volume=$BASE_DIR:/data`, so all other paths need to be _relative_ to it.
+ The container assumes that its folder `/data` is the base folder. 
+ The container parses the input arguments, then passes the right arguments to a function inside `common_4ch.process_handler`, which are modelled after the `main_*.py` scripts in this repository

## Step 0. Extract Surfaces
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surfs --par-folder parfiles --input-tags-setup file.json --apex-septum-setup apex_septum_setup
```
The container will modify the inputs and run function `common_4ch.extract_surfs`, which is modelled after script `main_surfs.py`.

## Step 1. Calculate UVCs 
> This cannot be run inside container as it calls CARP. 
This step is carried out by `main_UVCs.py`.

## Step 2. Add Ventricular Fibres 
> This cannot be run inside container as it calls CARP
This step is carried out by `main_fibres.py`

The function `correct_fibres.py` **can** be run inside the container: 
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch correcfibres --mesh-path surfaces_uvc/BiV 
```

## Step 3. Create Surfaces Endo Landmarks 
This has to be run twice, once per atrium:
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surf2vol --atrium la --fibres-endo Fibre_endo_l --fibres-epi Fibre_epi_l --mesh-path atrial_fibres/UAC
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surf2vol --atrium ra --fibres-endo Fibre_endo_l --fibres-epi Fibre_epi_l --mesh-path atrial_fibres/UAC
```

Important Assumptions: 
+ UAC files are stored inside the LA_endo and RA_endo folders.    
+ mesh is stored in `{base_dir}/{mesh_path}/{atrium}/{meshname}`, e.g `/data/atrial_fibres/UAC/la/la`

## Step 4: LA fibres / Step 5: RA fibres
> This cannot be run on the container because it depends on another container 

## Step 6: Map 2D to 3D
> This can be run partly in the container. 

The user will run the `laplace_prep` mode of operation, then a call to CARP outside the container, and 
finally call the mode `surf`.

Important Assumptions: 
+ UAC files are stored inside the LA_endo and RA_endo folders.    
+ surf files are stored in f"{base_dir}/{mesh_path}/{atrium}/{meshname}", e.g /data/atrial_fibres/UAC/la/la

#### `laplace_prep`
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch laplace_prep --atrium la --surf-endo LA_endo --surf-epi LA_epi --mesh-path atrial_fibres/UAC  
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch laplace_prep --atrium ra --surf-endo RA_endo --surf-epi RA_epi --mesh-path atrial_fibres/UAC 
```

#### Call to CARP outside container 

#### `surf2vol`
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surf2vol --atrium la --fibres-endo Fibre_endo_l --fibres-epi Fibre_epi_l --mesh-path atrial_fibres/UAC 
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surf2vol --atrium ra --fibres-endo Fibre_endo_l --fibres-epi Fibre_epi_l --mesh-path atrial_fibres/UAC 
```

## Step 7: Define tags
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch tags --data-subdir surfaces_uvc --mesh-path atrial_fibres --meshname myocardium_fibres_l --par-folder parfiles --input-tags-setup tags_presim.json --input-bb-settings bachmann_bundle_fec_settings.json
```
Important Assumptions: 
+ We assume there's a BiV, la, and ra folders inside the directory/data-subdir.

## Step 8: Extract Surfaces (Presim)
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch surfs_presim --data-subdir surfaces_uvc --par-folder parfiles --input-tags-setup tags_presim.json --map-settings map_settings.json --fch-apex fch_apex.txt --fch-sa fch_sa.txt
```
Assumptions:
+ We assume there are other subfolders with name `{subfolder}_LA` and `{subfolder}_RA` inside the directory

## Step 9: FEC
```sh
docker run --rm --volume=$BASE_DIR:/data cemrg/4ch fec --meshname myocardium --mesh-path meshing/myocardium_OUT --par-folder parfiles --input-tags-setup tags_presim.json --lvrv-tags lvrv_tags.json
```