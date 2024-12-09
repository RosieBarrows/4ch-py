#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <files_folder> [mesh_path]'
    >&2 echo ''
    >&2 echo 'mesh_path = meshing/myocardium_OUT/myocardium (default)'
    exit 1
fi

# clear

INPUT_heartFolder=$1
input_tags="$2/tags_vent_fibres.json"
etags="$2/etags/"
apex_septum="$2/apex_septum_templates/"
meshing_folder=${3:-meshing/myocardium_OUT/myocardium}

echo "MESHING FOLDER: ${meshing_folder}"

CMD="python main_UVCs_refact.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --etags_setup ${etags}
						 	 --apex_septum_setup ${apex_septum}
                             --meshing-folder ${meshing_folder}"
                             
echo $CMD
eval $CMD