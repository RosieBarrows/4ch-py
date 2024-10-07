#!/bin/bash

set -euo pipefail

if [ $# -lt 1 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 7_define_tags.sh <heart_folder>'
    exit 1
fi

clear

INPUT_heartFolder=$1
mesh="${INPUT_heartFolder}/atrial_fibres/myocardium_fibres_l"
BiV_folder="${INPUT_heartFolder}/surfaces_uvc/BiV"
input_tags="./parfiles/tags_presim.json"
BB_settings="./parfiles/bachmann_bundle_fec_settings.json"

LA_folder="${INPUT_heartFolder}/surfaces_uvc_LA/la/"
RA_folder="${INPUT_heartFolder}/surfaces_uvc_RA/ra/"

CMD="python main_tags.py --heartFolder ${INPUT_heartFolder}
						 	 --meshPath ${mesh}
						 	 --bivFolder ${BiV_folder}
						 	 --laFolder ${LA_folder}
						 	 --raFolder ${RA_folder}
						 	 --input_tags_setup ${input_tags}
						 	 --input_BB_settings ${BB_settings}"

eval $CMD