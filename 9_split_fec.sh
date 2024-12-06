#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
	>&2 echo 'Insufficient arguments supplied'
	>&2 echo 'Usage: 9_split_fec.sh <INPUT_heartFolder> <files_folder> [mesh_subfolder] '
    >&2 echo ''
    >&2 echo 'mesh_subfolder=myocardium_OUT (default)'
	exit 1
fi


INPUT_heartFolder=$1
files_folder=$2
mesh_subfolder=${3:-myocardium_OUT}
original_mesh="${INPUT_heartFolder}/meshing/${mesh_subfolder}/myocardium"
input_tags="$files_folder/tags_presim.json"
lvrv_tags="$files_folder/tags_lvrv.json"

echo "MESH: $original_mesh"

CMD="python main_fec.py --heartFolder ${INPUT_heartFolder}
						 	 --original_mesh ${original_mesh}
						 	 --input_tags_setup ${input_tags}
						 	 --lvrv_tags ${lvrv_tags}"

eval $CMD
