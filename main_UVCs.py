import os
import sys

from py_vent_fibres.json_utils import *
from py_vent_fibres.meshtools_utils import *

import argparse
import warnings

def main(args):

	os.system("clear")

	warnings.warn("MAKE SURE YOU HAVE SELECTED APEX/SEPTUM IN THE LA/RA")
	warnings.warn("MAKE SURE TAGS IN ETAGS FILES ARE CORRECT")

	etagsFolder = args.etags_setup

	heartFolder = args.heartFolder
	mesh=heartFolder+"/meshing/myocardium_OUT/myocardium"
	surf_folder=heartFolder+"/surfaces_uvc/"
	surf_folder_la=heartFolder+"/surfaces_uvc_LA/"
	surf_folder_ra=heartFolder+"/surfaces_uvc_RA/"


	# ----------------------------------------------------------------------------------------------
	# Calculating UVCs for the BiV mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Calculating UVCs for the BiV mesh ##")
	os.system("cp "+etagsFolder+"/etags.sh "+surf_folder+"/BiV/etags.sh")
	os.system("mguvc --model-name "+surf_folder+"BiV/BiV --input-model biv --output-model biv --np 20 --tags-file "+surf_folder+"BiV/etags.sh --output-dir "+surf_folder+"BiV/uvc/ --laplace-solution")

	os.system("GlVTKConvert -m "+surf_folder+"BiV/BiV -n "+surf_folder+"BiV/uvc/BiV.uvc_phi.dat -n "+surf_folder+"BiV/uvc/BiV.uvc_z.dat -n "+surf_folder+"BiV/uvc/BiV.uvc_ven.dat -n "+surf_folder+"BiV/uvc/BiV.uvc_rho.dat -o "+surf_folder+"BiV/uvc/uvc --trim-names")
	os.system("GlVTKConvert -m "+surf_folder+"BiV/BiV -n "+surf_folder+"BiV/uvc/BiV.sol_apba_lap.dat -n "+surf_folder+"BiV/uvc/BiV.sol_rvendo_lap.dat -n "+surf_folder+"BiV/uvc/BiV.sol_endoepi_lap.dat -n "+surf_folder+"/BiV/uvc/BiV.sol_lvendo_lap.dat -o "+surf_folder+"BiV/uvc/laplace --trim-names")

	# ----------------------------------------------------------------------------------------------
	# Calculating UVCs for the la mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Calculating UVCs for the la mesh ##")
	os.system("cp "+etagsFolder+"/etags_la.sh "+surf_folder_la+"/la/etags.sh")
	os.system("mguvc --ID=$MESH/UVC_ek --model-name "+surf_folder_la+"la/la --input-model lv --output-model lv --np 20 --tags-file "+surf_folder_la+"la/etags.sh --output-dir "+surf_folder_la+"la/uvc/ --laplace-solution --custom-apex")

	os.system("GlVTKConvert -m "+surf_folder_la+"la/la -n "+surf_folder_la+"la/uvc/la.uvc_phi.dat -n "+surf_folder_la+"la/uvc/la.uvc_z.dat -n "+surf_folder_la+"la/uvc/la.uvc_ven.dat -n "+surf_folder_la+"la/uvc/la.uvc_rho.dat -o "+surf_folder_la+"la/uvc/uvc --trim-names")

	# ----------------------------------------------------------------------------------------------
	# Calculating UVCs for the ra mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Calculating UVCs for the ra mesh ##")
	os.system("cp "+etagsFolder+"/etags_ra.sh "+surf_folder_ra+"/ra/etags.sh")
	os.system("mguvc --ID=$MESH/UVC_ek --model-name "+surf_folder_ra+"ra/ra --input-model lv --output-model lv --np 20 --tags-file "+surf_folder_ra+"ra/etags.sh --output-dir "+surf_folder_ra+"ra/uvc/ --laplace-solution --custom-apex")

	os.system("GlVTKConvert -m "+surf_folder_ra+"ra/ra -n "+surf_folder_ra+"ra/uvc/ra.uvc_phi.dat -n "+surf_folder_ra+"ra/uvc/ra.uvc_z.dat -n "+surf_folder_ra+"ra/uvc/ra.uvc_ven.dat -n "+surf_folder_ra+"ra/uvc/ra.uvc_rho.dat -o "+surf_folder_ra+"ra/uvc/uvc --trim-names")



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--etags_setup', type=str, default="./parfiles/etags",
                        help='Provide folder with etags files for BiV, LA and RA')

    parser.add_argument('--apex_septum_setup', type=str, default="./parfiles/apex_septum_templates",
                        help='Provide folder with templates for LA/RA apex and septum vtx files')

    args = parser.parse_args()

    main(args)