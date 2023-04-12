#!/bin/bash

INPUT_heartFolder=$(cat /data/Dropbox/scripts_cemrgapp/pipeline/vent_fibres/parfiles/heartFolder.txt)

scripts_path="/data/Dropbox/4ch-py"
heart_name=${INPUT_heartFolder}

cmd="mkdir ${heart_name}/atrial_fibres"
eval $cmd

MESHNAME="${heart_name}/surfaces_uvc/myocardium_bayer_60_-60"
UACFOLDER="${heart_name}/atrial_fibres/UAC/"

CMD="python main_mesh.py --meshname ${MESHNAME}
						 --outdir ${UACFOLDER}
						 --raa_apex_file ${heart_name}/raa_apex.txt
						 --input_tags_setup ${scripts_path}/parfiles/tags_atrial_fibres.json
						 --surface endo"

eval $CMD

CMD="mv ${heart_name}/raa_apex.txt ${UACFOLDER}/ra/raa_apex.txt"
eval $CMD
