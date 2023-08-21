#!/bin/bash

clear

INPUT_heartFolder=$(cat /data/Dropbox/4ch-py/parfiles/heartFolder.txt)
mesh="${INPUT_heartFolder}/atrial_fibres/myocardium_AV_FEC_BB"
BiV_folder="${INPUT_heartFolder}/surfaces_uvc/BiV"
input_tags="./parfiles/tags_presim.json"
atria_map_settings="./parfiles/atria_map_settings.json"
fch_apex="./parfiles/myocardium_AV_FEC_BB_apex.vtx"
fch_sa="./parfiles/myocardium_AV_FEC_BB_SA.vtx"

LA_folder="${INPUT_heartFolder}/surfaces_uvc_LA/la/"
RA_folder="${INPUT_heartFolder}/surfaces_uvc_RA/ra/"

CMD="python main_surfs_presim.py --heartFolder ${INPUT_heartFolder}
						 	 --meshPath ${mesh}
						 	 --bivFolder ${BiV_folder}
						 	 --laFolder ${LA_folder}
						 	 --raFolder ${RA_folder}
						 	 --input_tags_setup ${input_tags}
						 	 --map_settings ${atria_map_settings}
						 	 --fch_apex ${fch_apex}
						 	 --fch_sa ${fch_sa}"

eval $CMD