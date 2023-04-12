#!/bin/bash

clear

INPUT_heartFolder=$(cat /data/Dropbox/scripts_cemrgapp/pipeline/vent_fibres/parfiles/heartFolder.txt)
input_tags="./parfiles/tags_vent_fibres.json"
apex_septum="./parfiles/apex_septum_templates/"

CMD="python main_surfs.py --heartFolder ${INPUT_heartFolder}
						 	 --input_tags_setup ${input_tags}
						 	 --apex_septum_setup ${apex_septum}"
eval $CMD

echo ;
echo ;
echo " ### !! USER ACTION REQUIRED !! ### "
echo " ### You must now select a point for the apex and a point for the septum on both the LA and RA ### "
echo " ### Then you must select a point for apex of the right atrial appendage ###"
echo ;
echo ;