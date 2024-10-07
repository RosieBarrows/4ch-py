#!/bin/bash
set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <carp_root>'
    exit 1
fi

clear

INPUT_heartFolder=$1
mesh_path="${INPUT_heartFolder}/surfaces_uvc/BiV/"
four_chamber_path="${INPUT_heartFolder}/surfaces_uvc/"
four_chamber_name="myocardium"

carp_root=$2
CARP_FOLDER="${carp_root}/bin/"
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