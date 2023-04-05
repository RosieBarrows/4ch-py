import os
import sys

import argparse
import warnings

from py_pre_sim.json_utils import *
from py_pre_sim.mesh_utils import *
from py_pre_sim.meshtools_utils import *

def main(args):

	os.system("clear")

	warnings.warn("MAKE SURE INPUT TAGS ARE CORRECT")

	input_tags = load_json(args.input_tags_setup)

	heartFolder = args.heartFolder
	bivFolder = args.bivFolder
	laFolder = args.laFolder
	raFolder = args.raFolder
	mesh = args.meshPath
	map_settings = args.map_settings 
	fch_apex = args.fch_apex
	fch_sa = args.fch_sa

	presimFolder = args.heartFolder+"/pre_simulation/"
	meshName = "myocardium_AV_FEC_BB"
	mesh = presimFolder+meshName

	# # ----------------------------------------------------------------------------------------------
	# # Extracting the surface for the perciardium boundary condition
	# # ----------------------------------------------------------------------------------------------
	# meshtool_extract_peri(mesh,presimFolder,input_tags)
	# connected_component_to_surface(presimFolder+"/peri_surface_CC.part1",
	# 							   presimFolder+"/peri_surface.surf",
	# 							   presimFolder+"/epicardium_for_sim")

	# surf2vtk(mesh,presimFolder+"/epicardium_for_sim"+".surf",presimFolder+"/epicardium_for_sim"+".vtk")

	# os.system("rm "+presimFolder+"/*CC*")

	# # ----------------------------------------------------------------------------------------------
	# # Extracting the epi and LV/RV/LA/RA endo surfaces
	# # ----------------------------------------------------------------------------------------------
	# os.system("mkdir "+presimFolder+"/surfaces_simulation")
	# meshtool_extract_epi_endo_surfs(mesh,presimFolder,input_tags)

	# os.system("rm "+presimFolder+"/surfaces_simulation/surface_heart_CC.*")

	# # ----------------------------------------------------------------------------------------------
	# # Extracting the surfaces of the rings
	# # ----------------------------------------------------------------------------------------------
	# os.system("mkdir "+presimFolder+"/surfaces_simulation/surfaces_rings")
	# meshtool_extract_rings(mesh,presimFolder,input_tags)

	# # ----------------------------------------------------------------------------------------------
	# # Setting up the pericardium scale
	# # ----------------------------------------------------------------------------------------------
	# print("Setting up the pericardium scale...")
	# set_pericardium(mesh,presimFolder,heartFolder)


	# # ----------------------------------------------------------------------------------------------
	# # Setting up the pericardium scale for the atria
	# # ----------------------------------------------------------------------------------------------
	# print("Setting up the pericardium for the atria...")
	# la_mesh=heartFolder+"/surfaces_uvc_LA/la"
	# uvcs=la_mesh+"/uvc/"
	# os.system("python3 motion_atria_BCs.py --mesh "+la_mesh+" --uvcs "+uvcs+" --chamber la --map_settings "+map_settings)

	# ra_mesh=heartFolder+"/surfaces_uvc_RA/ra"
	# uvcs=ra_mesh+"/uvc/"
	# os.system("python3 motion_atria_BCs.py --mesh "+ra_mesh+" --uvcs "+uvcs+" --chamber ra --map_settings "+map_settings)


	# # ----------------------------------------------------------------------------------------------
	# # Combining the pericardium scaling maps for the ventricles and atria
	# # ----------------------------------------------------------------------------------------------
	# print("Combining the pericardium scaling maps for the ventricles and the atria...")
	# combine_elem_dats(heartFolder,presimFolder)

	# # ----------------------------------------------------------------------------------------------
	# # Setting up a folder with the simulation-ready mesh
	# # ----------------------------------------------------------------------------------------------
	# print("Setting up a folder with the simulation-ready mesh")
	# setup_sim(heartFolder,presimFolder,fch_apex,fch_sa)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--meshPath', type=str, default=None,
                        help='Provide path to mesh with atrial fibres (if applicable)')

    parser.add_argument('--bivFolder', type=str, default=None,
                        help='Provide path to folder containing BiV mesh (and uvc folder)')

    parser.add_argument('--laFolder', type=str, default=None,
                        help='Provide path to folder containing LA mesh (and uvc folder)')

    parser.add_argument('--raFolder', type=str, default=None,
                        help='Provide path to folder containing RA mesh (and uvc folder)')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--map_settings', type=str, default="./parfiles/atria_map_settings.json",
                        help='Provide json file with settings for the atrial pericardium map')

    parser.add_argument('--fch_apex', type=str, default=None,
    					help='Provide vtx file with ID of four chamber apex')

    parser.add_argument('--fch_sa', type=str, default=None,
    					help='Provide vtx file with ID of four chamber SA nodes')

    args = parser.parse_args()

    main(args)