#!/bin/bash

clear

INPUT_heartFolder=$(cat /data/Dropbox/scripts_cemrgapp/pipeline/vent_fibres/parfiles/heartFolder.txt)
mesh_path="${INPUT_heartFolder}/surfaces_uvc/BiV/"
four_chamber_path="${INPUT_heartFolder}/surfaces_uvc/"
four_chamber_name="myocardium"

CARP_FOLDER="/home/rb21/Software/CARPentry_KCL_latest/bin/"
alphaENDO=60	
alphaEPI=-60
betaENDO=-65
betaEPI=25

CMD="python main_fibres.py --heartFolder ${INPUT_heartFolder}
						   --mshPath ${mesh_path}
						   --fchPath ${four_chamber_path}
						   --fchName ${four_chamber_name}
						   --CARPFOLDER ${CARP_FOLDER}
						   --alpha_endo ${alphaENDO}
						   --alpha_epi ${alphaEPI}
						   --beta_endo ${betaENDO}
						   --beta_epi ${betaEPI}"
eval $CMD