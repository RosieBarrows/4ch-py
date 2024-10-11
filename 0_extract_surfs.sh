#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
	>&2 echo 'Insufficient arguments supplied'
	>&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <files_folder>'
	exit 1
fi

clear

INPUT_heartFolder=$1
input_tags="$2/tags_vent_fibres.json"
apex_septum="$2/apex_septum_templates/"

CMD="python main_surfs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --apex_septum_setup ${apex_septum}"
eval $CMD

echo ;
echo ;
echo " ### !! USER ACTION REQUIRED !! ### "
echo " ### You must now select a point for the apex and a point for the septum on both the LA and RA ### "
echo " ### Then you must select a point for apex of the right atrial appendage ###"
echo ;
echo ;


