#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
	>&2 echo 'Insufficient arguments supplied'
	>&2 echo 'Usage: 9_split_fec.sh <INPUT_heartFolder> <files_folder>'
	exit 1
fi

clear

INPUT_heartFolder=$1
files_folder=$2
original_mesh="${INPUT_heartFolder}/meshing/myocardium_OUT/myocardium"
input_tags="$files_folder/tags_presim.json"
lvrv_tags="$files_folder/tags_lvrv.json"

CMD="python main_fec.py --heartFolder ${INPUT_heartFolder}
						 	 --original_mesh ${original_mesh}
						 	 --input_tags_setup ${input_tags}
						 	 --lvrv_tags ${lvrv_tags}"

eval $CMD