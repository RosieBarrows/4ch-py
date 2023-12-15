#!/bin/bash

clear

INPUT_heartFolder=$(cat /home/croderog/Desktop/IC_projects/barrows_preprocessing/parfiles/heartFolder.txt)
fascicles_settings=$INPUT_heartFolder/fascicles_settings.json

CMD="python main_electrodes.py --heartFolder ${INPUT_heartFolder}
						--fascicles_settings ${fascicles_settings}
                        "

eval $CMD