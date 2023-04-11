import os
import sys

from common_4ch.file_utils import *
from common_4ch.mesh_utils import *

import argparse
import warnings

def main(args):

	meshname_3d = args.meshname
	meshname_2d = args.meshname_uac

	warnings.warn("The following executables should be in your PATH: meshtool.")

	print('===================================================================')
	print('Mapping tags from nodes to triangles in  '+meshname_2d+'...')
	print('===================================================================')

	tagfile_list = args.tagfiles
	tagfile_list_el_cc = []

	for tag_file in tagfile_list:

			tag_node2elem(meshname_2d,
						  tag_file,
						  tag_file[:-4]+"_el.dat")	

			tag_connected_component_2d(meshname_2d,
								   	   tag_file[:-4]+"_el.dat",
								   	   tag_file[:-4]+"_el_cc.dat")

			tagfile_list_el_cc.append(tag_file[:-4]+"_el_cc.dat")

	print('===================================================================')
	print('Unifying tags in one file using '+args.tag_settings_file+' settings...')
	print('===================================================================')

	tag_settings = load_json(args.tag_settings_file)
	labels_list = args.tag_labels
	tags_list = [int(tag_settings[l]) for l in labels_list]

	combine_tags(meshname_2d,
				 tagfile_list_el_cc,
				 tags_list,
				 tag_settings[args.default_tag_label],
				 args.outmeshname+"_2d.dat")

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

	if args.tag_labels_endo is not None:
		labels_list_endo = args.tag_labels_endo
		tags_list_endo = [int(tag_settings[l]) for l in labels_list_endo]
	else:
		tags_list_endo = []

	print('================================================================')
	print('Mapping tags from 2D to 3D...')
	print('================================================================')

	map_tags_3d(args.endo_epi_laplace[:-4]+"_el.dat",
				meshname_3d+"_elemCenters.pts",
				meshname_2d+"_elemCenters.pts",
				args.outmeshname+"_2d.dat",
				args.outmeshname+".dat",
				endo_tags=tags_list_endo,
				map_tr_tt_file=meshname_3d+"_endo_to_3d.map")

	elem = read_elem(meshname_3d+".elem",el_type='Tt',tags=False)
	region_tags = np.loadtxt(args.outmeshname+".dat",dtype=int)
	write_elem_caroline(elem,
			   region_tags,
			   args.outmeshname+".elem",
			   el_type='Tt')

	print('================================================================')
	print('Converting to VTK for visualization...')
	print('================================================================')	

	os.system("cp "+meshname_3d+".lon "+args.outmeshname+".lon")
	os.system("cp "+meshname_3d+".pts "+args.outmeshname+".pts")

	cmd = "meshtool convert -imsh="+args.outmeshname+" -omsh="+args.outmeshname+".vtk"
	os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--meshname', type=str, default=None,
                        help='Provide path to the carp mesh you want to map the tags to')

    parser.add_argument('--meshname_uac', type=str, default=None,
                        help='Provide path to the carp mesh you want to map the tags from')

    parser.add_argument('--tag_settings_file', type=str, default="./parfiles/region_tags.json",
                        help='Provide the json file with the tag settings')

    parser.add_argument('--tag_labels', nargs='+', 
                        help='Provide the list of labels in tag_settings_file corresponding to the tagfiles.')

    parser.add_argument('--tag_labels_endo', nargs='+', 
                        help='Provide the list of labels in tag_settings_file corresponding to the tagfiles that you want to map only to the endo.')

    parser.add_argument('--default_tag_label', type=str, default=None,
                        help='Provide the label containing the default tag.')

    parser.add_argument('--tagfiles', nargs='+', 
    	                help='Provide the list of files containing the tags you want to map (space separated)')

    parser.add_argument('--endo_epi_laplace', type=str, default=None,
                        help='Provide dat file with the endo-epi Laplace solution (defined on nodes)')

    parser.add_argument('--outmeshname', type=str, default=None,
                        help='Provide path to the carp mesh with tagged outputs')

    args = parser.parse_args()

    main(args)