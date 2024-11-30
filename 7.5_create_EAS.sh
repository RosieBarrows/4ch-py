#!/bin/bash

clear

INPUT_heartFolder=$1
fascicles_settings=$INPUT_heartFolder/parfiles/fascicles_settings.json

cp ./parfiles/fascicles_settings.json $fascicles_settings

CMD="python main_electrodes.py --heartFolder ${INPUT_heartFolder}
						--fascicles_settings ${fascicles_settings}
                        "

eval $CMD