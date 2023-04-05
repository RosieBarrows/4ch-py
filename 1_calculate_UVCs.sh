#!/bin/bash

clear

INPUT_heartFolder=$(cat /data/Dropbox/scripts_cemrgapp/pipeline/vent_fibres/parfiles/heartFolder.txt)
input_tags="./parfiles/input_tags_setup"
etags="./parfiles/etags/"
apex_septum="./parfiles/apex_septum_templates/"

CMD="python main_UVCs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --etags_setup ${etags}
						 	 --apex_septum_setup ${apex_septum}"
eval $CMD