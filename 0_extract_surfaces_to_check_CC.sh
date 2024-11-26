#!/bin/bash

meshfolder=$1

meshname="${meshfolder}/meshing/myocardium_OUT/myocardium"

cmd="mkdir -p ${meshfolder}/meshing/myocardium_OUT/surface_check/"
eval $cmd

cmd="meshtool extract surface -msh=${meshname} -surf=${meshfolder}/meshing/myocardium_OUT/surface_check/surface_heart -ofmt=carp_txt"
eval $cmd

cmd="meshtool extract surface -msh=${meshname} -surf=${meshfolder}/meshing/myocardium_OUT/surface_check/surface_heart -ofmt=vtk"
eval $cmd

cmd="meshtool extract unreachable -msh=${meshfolder}/meshing/myocardium_OUT/surface_check/surface_heart.surfmesh -submsh=${meshfolder}/meshing/myocardium_OUT/surface_check/surface_heart_CC -ofmt=vtk"
eval $cmd
