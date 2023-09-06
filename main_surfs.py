import os
import argparse

from common_4ch.process_handler import extract_surfs
from common_4ch.config import configure_logging
milog = configure_logging(log_name=__name__)

def main(args):
	
	os.system("clear")

	milog.warning("MAKE SURE INPUT TAGS ARE CORRECT")

	extract_surfs(args.heartFolder, args.input_tags_setup, args.apex_septum_setup, args.mesh)


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