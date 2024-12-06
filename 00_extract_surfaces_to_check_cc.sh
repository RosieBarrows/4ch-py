#!/usr/bin/env bash
set -euo pipefail

if [ $# -eq 0 ] ; then
    >&2 echo 'No arguments supplied'
    exit 1
fi

INPUT_heartfolder=$1
mesh_location=$2
meshname="${INPUT_heartfolder}/meshing/${mesh_location}/myocardium" 
outfolder="${INPUT_heartfolder}/meshing/${mesh_location}/surface_check"

cmd="mkdir -p ${outfolder}"
echo $cmd
eval $cmd 

cmd="meshtool extract surface -msh=${meshname} -surf=${INPUT_heartfolder}/meshing/${mesh_location}/surface_check/surface_heart -ofmt=carp_txt"
echo $cmd 
eval $cmd 

cmd="meshtool extract surface -msh=${meshname} -surf=${INPUT_heartfolder}/meshing/${mesh_location}/surface_check/surface_heart -ofmt=vtk"
echo $cmd
eval $cmd 

cmd="meshtool extract unreachable -msh=${INPUT_heartfolder}/meshing/${mesh_location}/surface_check/surface_heart.surfmesh -submsh=${INPUT_heartfolder}/meshing/${mesh_location}/surface_check/surface_heart_CC -ofmt=vtk"
echo $cmd 
eval $cmd 

echo ""

file_count=$(find "$outfolder" -type f -name "surface_heart_CC.part*" | wc -l)
echo "NUMBER OF PARTS found: ${file_count}"
echo "$(basename "${INPUT_heartfolder}"): ${file_count}" >> counts.txt 
