#!/bin/bash

INPUT_heartFolder=$1

scripts_path="/home/croderog/Desktop/IC_projects/barrows_preprocessing/"
heart_name=${INPUT_heartFolder}
echo ${heart_name}

cmd="mkdir -p ${heart_name}/atrial_fibres"
eval $cmd

MESHNAME="${heart_name}/surfaces_uvc/myocardium_bayer_60_-60"
UACFOLDER="${heart_name}/atrial_fibres/UAC/"

# CMD="python main_mesh.py --meshname ${MESHNAME}
# 						 --outdir ${UACFOLDER}
# 						 --raa_apex_file ${heart_name}/raa_apex.txt
# 						 --input_tags_setup ${scripts_path}/parfiles/tags_atrial_fibres.json
# 						 --surface endo"

CMD="python main_mesh_v3.py --meshname ${MESHNAME}
						 --outdir ${UACFOLDER}
						 --raa_apex_file ${heart_name}/raa_apex.txt
						 --input_tags_setup ${scripts_path}/parfiles/tags_atrial_fibres.json
						 --surface endo"

eval $CMD

CMD="cp ${heart_name}/raa_apex.txt ${UACFOLDER}/ra/raa_apex.txt"
eval $CMD
