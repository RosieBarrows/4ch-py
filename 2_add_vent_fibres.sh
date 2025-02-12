#!/bin/bash
set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: "2_add_vent_fibres".sh <INPUT_heartFolder> <carp_root> [MESHPATH]'
    exit 1
fi

# clear

INPUT_heartFolder=$1
mesh_path="${INPUT_heartFolder}/surfaces_uvc/BiV/"
four_chamber_path="${INPUT_heartFolder}/surfaces_uvc/"
four_chamber_name="myocardium"

carp_root=$2
mesh_input=${3:-meshing/myocardium_OUT/myocardium}
CARP_FOLDER="${carp_root}/bin/"
alphaENDO=60	
alphaEPI=-60
betaENDO=-65
betaEPI=25

echo "MESH INPUT: ${mesh_input}"
echo "MESH PATH FOR FIBRES: ${mesh_path}"

CMD="python main_fibres_refact.py --heartFolder ${INPUT_heartFolder}
						   --mshPath ${mesh_path}
						   --fchPath ${four_chamber_path}
						   --fchName ${four_chamber_name}
						   --CARPFOLDER ${CARP_FOLDER}
						   --alpha_endo ${alphaENDO}
						   --alpha_epi ${alphaEPI}
						   --beta_endo ${betaENDO}
						   --beta_epi ${betaEPI}
						   --msh-path-in ${mesh_input}"
echo $CMD >> second_half.log
eval $CMD
