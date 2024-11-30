#!/bin/bash

clear

INPUT_heartFolder=$1
original_mesh="${INPUT_heartFolder}/meshing/myocardium_OUT/myocardium"
input_tags="./parfiles/tags_presim.json"
lvrv_tags="./parfiles/tags_lvrv.json"

CMD="python main_fec.py --heartFolder ${INPUT_heartFolder}
						 	 --original_mesh ${original_mesh}
						 	 --input_tags_setup ${input_tags}
						 	 --lvrv_tags ${lvrv_tags}"

eval $CMD