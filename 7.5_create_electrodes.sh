#!/bin/bash

clear

INPUT_heartFolder=$1
fascicles_settings=$INPUT_heartFolder/fascicles_settings.json

CMD="python main_electrodes.py --heartFolder ${INPUT_heartFolder}
						--fascicles_settings ${fascicles_settings}
                        "

eval $CMD