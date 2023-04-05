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

	heartFolder = args.heartFolder
	input_tags = args.input_tags_setup
	lvrv_tags = args.lvrv_tags
	original_mesh = args.original_mesh

	# ----------------------------------------------------------------------------------------------
	# Splitting the FEC
	# ----------------------------------------------------------------------------------------------
	print(" ## Splitting the FEC (fast endocardial conduction) zone ## ")
	print(" ## This allows material/fibre properties to be assigned to the LV and RV independently ## ")

	f_input = open(input_tags,"r")
	original_tags = json.load(f_input)
	f_input.close()

	f_input = open(lvrv_tags,"r")
	new_tags = json.load(f_input)
	f_input.close()


	separate_FEC_lvrv(original_mesh+".elem",
                           heartFolder+'/sims_folder/myocardium_AV_FEC_BB.elem',
                           heartFolder+'/sims_folder/LV_endo.surf',
                           heartFolder+'/sims_folder/RV_endo.surf',
                           heartFolder+'/sims_folder/myocardium_AV_FEC_BB_lvrv.elem',
                           original_tags,
                           new_tags)                                                                       

	os.system("cp "+heartFolder+"/sims_folder/myocardium_AV_FEC_BB.lon "+heartFolder+"/sims_folder/myocardium_AV_FEC_BB_lvrv.lon")
	os.system("cp "+heartFolder+"/sims_folder/myocardium_AV_FEC_BB.pts "+heartFolder+"/sims_folder/myocardium_AV_FEC_BB_lvrv.pts")
	os.system("meshtool convert -imsh="+heartFolder+"/sims_folder/myocardium_AV_FEC_BB_lvrv -omsh="+heartFolder+"/sims_folder/myocardium_AV_FEC_BB_lvrv -ofmt=vtk_bin")



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--original_mesh', type=str, default=None,
                        help='Provide path to original mesh (before tag reassignment')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--lvrv_tags', type=str, default=None,
                        help='Provide json file with input tags settings for split FEC')


    args = parser.parse_args()

    main(args)