import os
import sys

from common_4ch.file_utils import *
from common_4ch.mesh_utils import *
from common_4ch.linalg_utils import *

import argparse
import warnings

def main(args):

	meshname_3d = args.meshname
	meshname_2d = args.meshname_uac

	print('================================================================')
	print('Mapping fibres from 2D to 3D on '+meshname_3d+'...')
	print('================================================================')

	warnings.warn("The following executables should be in your PATH: meshtool.")

	print('================================================================')
	print('Mapping endo-epi Laplace solution from nodes to elements...')
	print('================================================================')

	laplace_endo2elem(meshname_3d,
					  args.endo_epi_laplace)

	print('================================================================')
	print('Computing element centers on both meshes 3D and 2D meshes...')
	print('================================================================')

	compute_elemCenters(meshname_3d,el_type="Tt")
	compute_elemCenters(meshname_2d,el_type="Tr")

	print('================================================================')
	print('Mapping fibres...')
	print('================================================================')

	map_fibres_3d(args.endo_epi_laplace[:-4]+"_el.dat",
				  meshname_3d+"_elemCenters.pts",
				  meshname_2d+"_elemCenters.pts",
				  args.endo_fibres,
				  args.epi_fibres,
				  args.outmeshname+".lon",
				  map_tr_tt_file=meshname_3d+"_endo_to_3d.map")

	print('================================================================')
	print('Finding sheet direction...')
	print('================================================================')

	find_transmural_direction(meshname_3d,
							  meshname_2d,
							  meshname_3d+"_elemCenters.pts",
				  			  meshname_2d+"_elemCenters.pts",
							  meshname_3d+"_transmural.lon")

	find_rotation_axes(args.outmeshname+".lon",
					   meshname_3d+"_transmural.lon",
					   meshname_3d+"_rotation_axes.lon")

	make_sheet_orthogonal(args.outmeshname+".lon",
						  meshname_3d+"_transmural.lon",
						  meshname_3d+"_rotation_axes.lon",
						  args.outmeshname+"_sheet.lon")

	print('================================================================')
	print('Converting to VTK for visualization...')
	print('================================================================')	

	os.system("cp "+meshname_3d+".elem "+args.outmeshname+"_sheet.elem")
	os.system("cp "+meshname_3d+".pts "+args.outmeshname+"_sheet.pts")

	cmd = "meshtool convert -imsh="+args.outmeshname+"_sheet -omsh="+args.outmeshname+"_sheet.vtk"
	os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--meshname', type=str, default=None,
                        help='Provide path to the carp mesh')

    parser.add_argument('--meshname_uac', type=str, default=None,
                        help='Provide path to the carp mesh with the fibres you want to map')

    parser.add_argument('--endo_fibres', type=str, default=None,
                        help='Provide path to the lon file with the endo fibres')

    parser.add_argument('--epi_fibres', type=str, default=None,
                        help='Provide path to the lon file with the epi fibres')

    parser.add_argument('--endo_epi_laplace', type=str, default=None,
                        help='Provide dat file with the endo-epi Laplace solution (defined on nodes)')

    parser.add_argument('--outmeshname', type=str, default=None,
                        help='Provide path to the output carp mesh')

    args = parser.parse_args()

    main(args)