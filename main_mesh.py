import os
import sys

from py_atrial_fibres.json_utils import *
from py_atrial_fibres.meshtools_utils import *

import argparse
import warnings

def main(args):

	os.system("clear")

	warnings.warn("The following executables should be in your PATH: carp.pt, igbextract, GlVTKConvert.")

	input_tags = load_json(args.input_tags_setup)
	output_tags = load_json(args.output_tags_setup)

	meshname = args.meshname
	output_folder = args.outdir

	choice = input('Do you want only to recompute the RAA base? [y/n]')
	while choice not in ['y','n']:
		choice = input('Please type y for recomputing only the RAA base and n for recomputing everything.')

	if choice == 'n':

		biatrial_folder = output_folder+'/biatrial/'
		la_folder = output_folder+'/la/'
		ra_folder = output_folder+'/ra/'	

		os.system('mkdir -p '+biatrial_folder)
		os.system('mkdir -p '+la_folder)
		os.system('mkdir -p '+ra_folder)	

		meshtool_extract_biatrial(meshname,
								  output_folder,
								  input_tags)	

		meshtool_extract_LA(output_folder,
							input_tags)	

		meshtool_extract_RA(output_folder,
							input_tags)	

		meshtool_extract_surfaces(meshname,
			 					  output_folder,
			 					  input_tags,
			 					  export_sup_inf=True,
			 					  rm_vv_from_aa=True,
			 					  surface=args.surface)	

		export_vtk_meshes_caroline(output_folder,raa_apex_file=None,surface=args.surface)

		os.system("mkdir "+output_folder+"/Landmarks/")
		os.system("mkdir "+output_folder+"/Landmarks/LA/")
		os.system("mkdir "+output_folder+"/Landmarks/RA/")	

		os.system("cp "+output_folder+"/LA_epi/prodLaLandmarks.txt "+output_folder+"/Landmarks/LA/prodRaLandmarks.txt")
		os.system("cp "+output_folder+"/LA_epi/prodLaRegion.txt "+output_folder+"/Landmarks/LA/prodRaRegion.txt")	

		os.system("cp "+output_folder+"/RA_epi/prodRaLandmarks.txt "+output_folder+"/Landmarks/RA/")
		os.system("cp "+output_folder+"/RA_epi/prodRaRegion.txt "+output_folder+"/Landmarks/RA/")	

		os.system("rm -r "+output_folder+"/tmp/")

		print('=========================================================================================')
		print('- Select the right atrial appendage tip')
		print('- Save the coordinate in a .txt file')
		print('- Re-run this code by setting the raa_apex_file argument pointing to this file')
		print('- When asked if you want to recompute only the RAA base, say "y"')
		print('=========================================================================================')

	else:
		recompute_raa_base(output_folder,
					   args.raa_apex_file,
					   output_folder+"/RA_epi/prodRaRegion.txt",
					   scale=0.001,
					   surface=args.surface)
		scale_landmarks(output_folder+"/RA_epi/prodRaLandmarks.txt",
					    scale=0.001)

		os.system("mkdir "+output_folder+"/Landmarks/")
		os.system("mkdir "+output_folder+"/Landmarks/LA/")
		os.system("mkdir "+output_folder+"/Landmarks/RA/")	

		scale_landmarks(output_folder+"/LA_epi/prodLaLandmarks.txt",
					    scale=0.001)
		scale_landmarks(output_folder+"/LA_epi/prodLaRegion.txt",
					    scale=0.001)

		os.system("cp "+output_folder+"/LA_epi/prodLaLandmarks.txt "+output_folder+"/LA_endo/prodLaLandmarks.txt")
		os.system("cp "+output_folder+"/LA_epi/prodLaRegion.txt "+output_folder+"/LA_endo/prodLaRegion.txt")	
		os.system("cp "+output_folder+"/LA_epi/prodLaLandmarks.txt "+output_folder+"/Landmarks/LA/prodRaLandmarks.txt")
		os.system("cp "+output_folder+"/LA_epi/prodLaRegion.txt "+output_folder+"/Landmarks/LA/prodRaRegion.txt")	

		os.system("cp "+output_folder+"/RA_epi/prodRaLandmarks.txt "+output_folder+"/RA_endo/prodRaLandmarks.txt")
		os.system("cp "+output_folder+"/RA_epi/prodRaRegion.txt "+output_folder+"/RA_endo/prodRaRegion.txt")
		os.system("cp "+output_folder+"/RA_epi/prodRaLandmarks.txt "+output_folder+"/Landmarks/RA/")
		os.system("cp "+output_folder+"/RA_epi/prodRaRegion.txt "+output_folder+"/Landmarks/RA/")

		os.system("meshtool convert -imsh="+output_folder+"/LA_endo/LA_endo.vtk -ofmt=carp_txt -omsh="+output_folder+"/LA_endo/LA_only")
		os.system("meshtool convert -imsh="+output_folder+"/LA_epi/LA_epi.vtk -ofmt=carp_txt -omsh="+output_folder+"/LA_epi/LA_only")
		os.system("meshtool convert -imsh="+output_folder+"/RA_endo/RA_endo.vtk -ofmt=carp_txt -omsh="+output_folder+"/RA_endo/RA_only")
		os.system("meshtool convert -imsh="+output_folder+"/RA_epi/RA_epi.vtk -ofmt=carp_txt -omsh="+output_folder+"/RA_epi/RA_only")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--meshname', type=str, default=None,
                        help='Provide path to the carp four-chamber mesh')

    parser.add_argument('--surface', type=str, default="epi",
                        help='Surface to compute the landmarks (endo or epi)')

    parser.add_argument('--input_tags_setup', type=str, default="./parfiles/input_tags_setup.json",
                        help='Provide json file with input tags settings')

    parser.add_argument('--output_tags_setup', type=str, default="./parfiles/output_tags_setup.json",
                        help='Provide json file with output tags settings (for Milan fibres)')

    parser.add_argument('--raa_apex_file', type=str, default=None,
                        help='Provide txt file with coordinate of the right atrial appendage tip')

    parser.add_argument('--outdir', type=str, default=None,
                        help='Provide output directory')

    args = parser.parse_args()

    main(args)