#!/bin/bash

clear

INPUT_heartFolder=$(cat /data/Dropbox/scripts_cemrgapp/pipeline/pre_sim/parfiles/heartFolder.txt)
original_mesh="${INPUT_heartFolder}/meshing/myocardium_OUT/myocardium"
input_tags="./parfiles/input_tags_setup.json"
lvrv_tags="./parfiles/tags_lvrv.json"

CMD="python main_fec.py --heartFolder ${INPUT_heartFolder}
						 	 --original_mesh ${original_mesh}
						 	 --input_tags_setup ${input_tags}
						 	 --lvrv_tags ${lvrv_tags}"

eval $CMD