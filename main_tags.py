import os
import sys

import argparse
import warnings

from common_4ch.json_utils import *
from common_4ch.mesh_utils import *

def main(args):

	os.system("clear")

	warnings.warn("MAKE SURE INPUT TAGS ARE CORRECT")

	input_tags = load_json(args.input_tags_setup)
	bb_settings = load_json(args.input_BB_settings)

	heartFolder = args.heartFolder
	bivFolder = args.bivFolder
	laFolder = args.laFolder
	raFolder = args.raFolder
	mesh=args.meshPath

	# ----------------------------------------------------------------------------------------------
	# Setting up folder structure
	# ----------------------------------------------------------------------------------------------

	presimFolder = heartFolder+"/pre_simulation/"
	os.system("mkdir "+presimFolder)

	os.system("cp "+mesh+".elem "+presimFolder+"/myocardium.elem")
	os.system("cp "+mesh+".pts "+presimFolder+"/myocardium.pts")
	os.system("cp "+mesh+".lon "+presimFolder+"/myocardium.lon")

	MESHNAME="myocardium"
	MESHNAME_AV="myocardium_AV"
	MESHNAME_FEC="myocardium_AV_FEC"
	MESHNAME_BB="myocardium_AV_FEC_BB"

	# ----------------------------------------------------------------------------------------------
	# Defining the AV separating plane
	# ----------------------------------------------------------------------------------------------
	define_AV_separation(presimFolder+"/"+MESHNAME+".elem",
				  		 			input_tags,
				  		 			input_tags["AV_plane"],
				  		 			presimFolder+"/"+MESHNAME_AV+".elem")

	# ----------------------------------------------------------------------------------------------
	# Defining the FEC layer
	# ----------------------------------------------------------------------------------------------

	biv_mesh=bivFolder+"/BiV"
	Zbiv_file=bivFolder+"/uvc/BiV.uvc_z.dat"
	RHObiv_file=bivFolder+"/uvc/BiV.uvc_rho.dat"

	define_FEC(presimFolder+"/"+MESHNAME_AV+".elem",
						  biv_mesh,
						  Zbiv_file,
						  RHObiv_file,
						  presimFolder+"/"+MESHNAME_FEC+".elem",
						  input_tags["FEC"],
						  include_septum=bivFolder+"/BiV.rvsept.surf",
						  FEC_height=bb_settings["FEC_height"])

	# ----------------------------------------------------------------------------------------------
	# Defining the BB area
	# ----------------------------------------------------------------------------------------------

	laMesh=laFolder+"/la"
	raMesh=raFolder+"/ra"
	Zla_file=laFolder+"uvc/la.uvc_z.dat"
	Zra_file=raFolder+"uvc/ra.uvc_z.dat"
	PHIla_file=laFolder+"uvc/la.uvc_phi.dat"
	PHIra_file=raFolder+"uvc/ra.uvc_phi.dat"

	define_BB(presimFolder+"/"+MESHNAME_FEC+".elem",
						 laMesh,
						 raMesh,
						 Zla_file,
						 Zra_file,
						 PHIla_file,
						 PHIra_file,
						 bb_settings,
						 input_tags,
						 presimFolder+"/"+MESHNAME_BB+".elem")

	os.system("cp  "+presimFolder+"/"+MESHNAME+".pts "+presimFolder+"/"+MESHNAME_BB+".pts")
	os.system("cp  "+presimFolder+"/"+MESHNAME+".lon "+presimFolder+"/"+MESHNAME_BB+".lon")
	os.system("meshtool convert -imsh="+presimFolder+"/"+MESHNAME_BB+" -omsh="+presimFolder+"/"+MESHNAME_BB+".vtk")


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

    parser.add_argument('--input_BB_settings', type=str, default="./parfiles/bachmann_bundle_fec_settings.json",
                        help='Provide json file with settings for the Bachmann bundle')

    args = parser.parse_args()

    main(args)