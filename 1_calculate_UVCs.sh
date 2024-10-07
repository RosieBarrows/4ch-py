#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <files_folder>'
    exit 1
fi

clear

INPUT_heartFolder=$1
input_tags="$2/tags_vent_fibres.json"
etags="$2/etags/"
apex_septum="$2/apex_septum_templates/"

CMD="python main_UVCs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --etags_setup ${etags}
						 	 --apex_septum_setup ${apex_septum}"
eval $CMD