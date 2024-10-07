#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <parfiles>'
    exit 1
fi

INPUT_heartFolder=$1

parfiles=$2
heart_name=${INPUT_heartFolder}
echo ${heart_name}

cmd="mkdir -p ${heart_name}/atrial_fibres"
eval $cmd

MESHNAME="${heart_name}/surfaces_uvc/myocardium_bayer_60_-60"
UACFOLDER="${heart_name}/atrial_fibres/UAC/"

CMD="python main_mesh.py --meshname ${MESHNAME}
						 --outdir ${UACFOLDER}
						 --raa_apex_file ${heart_name}/raa_apex.txt
						 --input_tags_setup ${parfiles}/tags_atrial_fibres.json
						 --surface endo"

eval $CMD

CMD="cp ${heart_name}/raa_apex.txt ${UACFOLDER}/ra/raa_apex.txt"
eval $CMD
