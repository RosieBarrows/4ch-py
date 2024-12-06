#!/bin/bash

set -euo pipefail

if [ $# -eq 0 ] ; then
    >&2 echo 'No arguments supplied'
    >&2 echo '    EXAMPLE_FOLDER'
    exit 1
fi

INPUT_heartFolder=$1
fascicles_settings=$INPUT_heartFolder/fascicles_settings.json

cp ${INPUT_heartFolder}/parfiles/fascicles_settings.json ${fascicles_settings} 

CMD="python main_electrodes.py --heartFolder ${INPUT_heartFolder}
						--fascicles_settings ${fascicles_settings}
                        "

eval $CMD
