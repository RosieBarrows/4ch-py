#!/bin/bash

scripts_path="/data/Dropbox/cemrg_docker_scripts/4ch-py"
heart_name=$1

cmd="mkdir ${heart_name}/atrial_fibres"

MESHNAME="${heart_name}/surfaces_uvc/myocardium_bayer_60_-60"
UACFOLDER="${heart_name}/atrial_fibres/UAC/"

CMD="python main_mesh.py --meshname ${MESHNAME}
						 --outdir ${UACFOLDER}
						 --raa_apex_file ${heart_name}/atrial_fibres/raa_apex.txt
						 --input_tags_setup ${scripts_path}/parfiles/tags_atrial_fibres.json
						 --surface endo"

eval $CMD
