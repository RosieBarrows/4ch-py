import os
import sys

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *
from common_4ch.file_utils import mycp, mymkdir

import argparse
import warnings

def main(args):
	
	os.system("clear")

	warnings.warn("MAKE SURE INPUT TAGS ARE CORRECT")

	input_tags = load_json(args.input_tags_setup)
	apexFolder = args.apex_septum_setup

	heartFolder = args.heartFolder
	mesh=os.path.join(heartFolder, args.mesh) # default = "meshing/myocardium_OUT/myocardium"
	surf_folder=heartFolder+"/surfaces_uvc/"
	surf_folder_la=heartFolder+"/surfaces_uvc_LA/"
	surf_folder_ra=heartFolder+"/surfaces_uvc_RA/"
	
	list_of_folders = [os.path.join(surf_folder, subf) for subf in ["", "tmp", "BiV"] ]
	list_of_folders += [os.path.join(surf_folder_la, subf) for subf in ["", "tmp", "la"] ]
	list_of_folders += [os.path.join(surf_folder_ra, subf) for subf in ["", "tmp", "ra"] ]
      
	for folder in list_of_folders:
		mymkdir(folder)

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
	mycp(f"{apexFolder}/la.lvapex.vtx", f"{surf_folder_la}/la/la.lvapex.vtx")
	mycp(f"{apexFolder}/la.rvsept_pt.vtx", f"{surf_folder_la}/la/la.rvsept_pt.vtx")

	# ----------------------------------------------------------------------------------------------
	# Extracting the ra mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Extracting the right atrial mesh ##")
	meshtool_extract_ra_for_UVCs(mesh,surf_folder_ra,input_tags)

	print(" ## Mapping vtx files from four-chamber mesh to right atrial mesh ## ")
	meshtool_map_vtx_ra(surf_folder_ra)

	print(" ## Copying blank files for RA apex and septum ## ")
	mycp(f"{apexFolder}/ra.lvapex.vtx", f"{surf_folder_ra}/ra/ra.lvapex.vtx")
	mycp(f"{apexFolder}/ra.rvsept_pt.vtx", f"{surf_folder_ra}/ra/ra.rvsept_pt.vtx")

	print(" ## Copying blank file for RAA apex ## ")
	mycp(f"{apexFolder}/raa_apex.txt", f"{heartFolder}/raa_apex.txt")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--apex_septum_setup', type=str, default="./parfiles/apex_septum_templates",
                        help='Provide folder with templates for LA/RA apex and septum vtx files')
	
    parser.add_argument('-msh', '--mesh', type=str, default='meshing/myocardium_OUT/myocardium',
					 	help="Path to the mesh file (relative to heartFolder)")

    args = parser.parse_args()

    main(args)