#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
	>&2 echo 'Insufficient arguments supplied'
	>&2 echo 'Usage: 0_extract_surfs.sh <INPUT_heartFolder> <files_folder> [meshing_path]'
    >&2 echo ''
    >&2 echo 'meching_path = meshing/myocardium_OUT/myocardium (default)' 
	exit 1
fi

#clear

INPUT_heartFolder=$1
input_tags="$2/tags_vent_fibres.json"
apex_septum="$2/apex_septum_templates/"
meshname="${3:-meshing/myocardium_OUT/myocardium}" 

echo "MESHNAME: ${meshname}"


CMD="python main_surfs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --apex_septum_setup ${apex_septum}
                             -msh ${meshname}"

echo $CMD > ${INPUT_heartFolder}/0_extract_surfs.log
eval $CMD

echo ;
echo ;
echo " ### !! USER ACTION REQUIRED !! ### "
echo " ### You must now select a point for the apex and a point for the septum on both the LA and RA ### "
echo " ### Then you must select a point for apex of the right atrial appendage ###"
echo ;
echo ;


