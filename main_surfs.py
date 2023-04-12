import os
import sys

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *

import argparse
import warnings

def main(args):

	os.system("clear")

	warnings.warn("MAKE SURE INPUT TAGS ARE CORRECT")

	input_tags = load_json(args.input_tags_setup)
	apexFolder = args.apex_septum_setup

	heartFolder = args.heartFolder
	mesh=heartFolder+"/meshing/myocardium_OUT/myocardium"
	surf_folder=heartFolder+"/surfaces_uvc/"
	surf_folder_la=heartFolder+"/surfaces_uvc_LA/"
	surf_folder_ra=heartFolder+"/surfaces_uvc_RA/"

	os.system("mkdir "+surf_folder)
	os.system("mkdir "+surf_folder+"/tmp")
	os.system("mkdir "+surf_folder+"/BiV")
	os.system("mkdir "+surf_folder_la)
	os.system("mkdir "+surf_folder_la+"/tmp")
	os.system("mkdir "+surf_folder_la+"/la")
	os.system("mkdir "+surf_folder_ra)
	os.system("mkdir "+surf_folder_ra+"/tmp")
	os.system("mkdir "+surf_folder_ra+"/ra")

	# ----------------------------------------------------------------------------------------------
	# Extract the base
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the base ## ")
	meshtool_extract_base(mesh,surf_folder,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Extract the epicardium, LV endo and RV endo
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting epi, LV endo and RV endo ## ")
	meshtool_extract_surfaces_lv_rv_epi(mesh,surf_folder,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Extract the RV septum
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting septum ## ")
	meshtool_extract_septum(mesh,surf_folder,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Mapping surfaces
	# ----------------------------------------------------------------------------------------------
	print(" ## Mapping surfaces ## ")
	mapping_surfaces(mesh,surf_folder,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Removing the septum
	# ----------------------------------------------------------------------------------------------
	print(" ## Removing the septum ## ")
	remove_sept(mesh,surf_folder)

	# ----------------------------------------------------------------------------------------------
	# Prepare vtx for UVCs
	# ----------------------------------------------------------------------------------------------
	print(" ## Preparing vtx files for UVCs")
	prepare_vtx_for_uvc(surf_folder)

	# ----------------------------------------------------------------------------------------------
	# Extracting surfaces in the LA
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the LA base ## ")
	meshtool_extract_la_base(mesh,surf_folder_la,input_tags)

	print(" ## Extracting the LA epi and LA endo")
	meshtool_extract_la_surfaces(mesh,surf_folder_la,input_tags)

	print(" ## Mapping surfaces LA")
	mapping_surfaces_la(mesh,surf_folder_la,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Extracting surfaces in the RA
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the RA base ## ")
	meshtool_extract_ra_base(mesh,surf_folder_ra,input_tags)

	print(" ## Extracting the RA epi and RA endo")
	meshtool_extract_ra_surfaces(mesh,surf_folder_ra,input_tags)

	print(" ## Mapping surfaces RA")
	mapping_surfaces_ra(mesh,surf_folder_ra,input_tags)

	# ----------------------------------------------------------------------------------------------
	# Extracting the BiV mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the biventricular mesh ## ")
	meshtool_extract_biv(mesh,surf_folder,input_tags)

	print(" ## Mapping vtx files from four-chamber mesh to BiV mesh ## ")
	meshtool_map_vtx(surf_folder)

	print(" ## Renaming files ## ")
	renaming_myo_files(surf_folder)

	# ----------------------------------------------------------------------------------------------
	# Extracting the la mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the left atrial mesh ##")
	meshtool_extract_la_for_UVCs(mesh,surf_folder_la,input_tags)

	print(" ## Mapping vtx files from four-chamber mesh to left atrial mesh ## ")
	meshtool_map_vtx_la(surf_folder_la)

	print(" ## Copying blank files for LA apex and septum ## ")
	os.system("cp "+apexFolder+"/la.lvapex.vtx "+surf_folder_la+"/la/la.lvapex.vtx")
	os.system("cp "+apexFolder+"/la.rvsept_pt.vtx "+surf_folder_la+"/la/la.rvsept_pt.vtx")


	# ----------------------------------------------------------------------------------------------
	# Extracting the ra mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the right atrial mesh ##")
	meshtool_extract_ra_for_UVCs(mesh,surf_folder_ra,input_tags)

	print(" ## Mapping vtx files from four-chamber mesh to right atrial mesh ## ")
	meshtool_map_vtx_ra(surf_folder_ra)

	print(" ## Copying blank files for RA apex and septum ## ")
	os.system("cp "+apexFolder+"/ra.lvapex.vtx "+surf_folder_ra+"/ra/ra.lvapex.vtx")
	os.system("cp "+apexFolder+"/ra.rvsept_pt.vtx "+surf_folder_ra+"/ra/ra.rvsept_pt.vtx")

	print(" ## Copying blank file for RAA apex ## ")
	os.system("cp "+apexFolder+"/raa_apex.txt "+heartFolder+"/raa_apex.txt")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--apex_septum_setup', type=str, default="./parfiles/apex_septum_templates",
                        help='Provide folder with templates for LA/RA apex and septum vtx files')

    args = parser.parse_args()

    main(args)