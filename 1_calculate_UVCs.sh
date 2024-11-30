#!/bin/bash

clear

INPUT_heartFolder=$1
input_tags="./parfiles/tags_vent_fibres.json"
etags="./parfiles/etags/"
apex_septum="./parfiles/apex_septum_templates/"

# CMD="python main_UVCs.py --heartFolder ${INPUT_heartFolder}
# 						 	 --input_tags_setup ${input_tags}
# 						 	 --etags_setup ${etags}
# 						 	 --apex_septum_setup ${apex_septum}"

CMD="python main_UVCs_automatic.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --etags_setup ${etags}
						 	 --apex_septum_setup ${apex_septum}"
eval $CMD