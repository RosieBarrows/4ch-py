#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
	>&2 echo 'Insufficient arguments supplied'
	>&2 echo 'Usage: 0_1_extract_surfs_and_UVCs.sh <INPUT_heartFolder> <files_folder> [meshing_path]'
    >&2 echo ''
    >&2 echo 'meching_path = meshing/myocardium_OUT/myocardium (default)' 
	exit 1
fi

#clear

INPUT_heartFolder=$1
input_tags="$2/tags_vent_fibres.json"
etags="$2/etags/"
apex_septum="$2/apex_septum_templates/"
meshname="${3:-meshing/myocardium_OUT/myocardium}" 

echo "MESHNAME: ${meshname}"


CMD="python main_surfs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --apex_septum_setup ${apex_septum}
                             -msh ${meshname}"

echo $CMD > ${INPUT_heartFolder}/second_half.log
eval $CMD

echo ;

CMD="python main_UVCs_refact.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --etags_setup ${etags}
						 	 --apex_septum_setup ${apex_septum}
                             --meshing-folder ${meshname}"
                             
echo $CMD >> ${INPUT_heartFolder}/second_half.log
eval $CMD

