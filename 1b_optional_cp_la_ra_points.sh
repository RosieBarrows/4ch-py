#!/bin/bash

set -euo pipefail

if [ $# -lt 2 ] ; then
    >&2 echo 'Insufficient arguments supplied'
    >&2 echo 'Usage: 1_calculate_UVCs.sh <INPUT_heartFolder> <GT_folder>'
    exit 1
fi

INPUT_heartFolder=$1
gt_heart_folder=$2

echo "Making directories"
mkdir -p $INPUT_heartFolder/surfaces_uvc_RA/ra/uvc
mkdir -p $INPUT_heartFolder/surfaces_uvc_LA/la/uvc

# only for debugging 
# if third argument is given, then copy files to loactions 
echo "Copying RA/LA apex and septum files to from GT folder to correct to locations"
cp $gt_heart_folder/surfaces_uvc_RA/ra/ra.rvsept_pt.vtx $INPUT_heartFolder/surfaces_uvc_RA/ra/ra.rvsept_pt.vtx
cp $gt_heart_folder/surfaces_uvc_RA/ra/ra.lvapex.vtx $INPUT_heartFolder/surfaces_uvc_RA/ra/ra.lvapex.vtx

cp $gt_heart_folder/surfaces_uvc_LA/la/uvc/la.lvapex.vtx $INPUT_heartFolder/surfaces_uvc_LA/la/uvc/la.lvapex.vtx
cp $gt_heart_folder/surfaces_uvc_LA/la/uvc/la.rvsept_pt.vtx $INPUT_heartFolder/surfaces_uvc_LA/la/uvc/la.rvsept_pt.vtx

cp $gt_heart_folder/surfaces_uvc_LA/la/la.rvsept_pt.vtx $INPUT_heartFolder/surfaces_uvc_LA/la/la.rvsept_pt.vtx
cp $gt_heart_folder/surfaces_uvc_LA/la/la.lvapex.vtx $INPUT_heartFolder/surfaces_uvc_LA/la/la.lvapex.vtx

mkdir -p $INPUT_heartFolder/atrial_fibres/UAC/ra
cp $gt_heart_folder/atrial_fibres/UAC/ra/raa_apex.txt $INPUT_heartFolder/atrial_fibres/UAC/ra/raa_apex.txt
cp $gt_heart_folder/atrial_fibres/UAC/ra/raa_apex.txt $INPUT_heartFolder/raa_apex.txt

cp $gt_heart_folder/fascicles_settings.json $INPUT_heartFolder/fascicles_settings.json
