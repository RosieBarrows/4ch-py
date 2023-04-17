import os
import argparse
import logging
import pathlib

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *

from common_4ch.custom_commands import my_cp, my_mkdir

def _surfs(base_dir, input_tags_setup, apex_septum_setup, meshname="myocardium", debug=False, help=False): 
    """
    Mode of operation: surfs
    Extract surfaces from the mesh using meshtool

    Parameters:
    --input-tags-setup: relative path to the json file containing the input tags
    --apex-septum-setup: relative path to the folder containing the apex and septum tags
    --meshname: name of the mesh file (default: myocardium)

    Usage:
    surfs --input-tags-setup sub/dir/to/file.json --apex-septum-setup apex_septum_setup
    """

    if help :
        print(_surfs.__doc__)
        return

    if debug :  logging.debug(f"Importing extract_surfs function")
    from common_4ch.custom_commands import extract_surfs

    extract_surfs(base_dir, input_tags_setup, apex_septum_setup, meshname, debug=False)


    

def _correct_fibres(base_dir, mesh_path, debug, help=False) :
    """
    Mode of operation: correctfibres
    Correct the fibres

    Parameters:
    --mesh-path: full path to the folder containing the BiV mesh

    Usage:
    correctfibres --mesh-path sub/dir/to/folder

    """
    if help :
        print(_correct_fibres.__doc__)
        return
    
    if debug :  logging.debug(f"Importing correct_fibres function")
    from common_4ch.custom_commands import correct_fibres 

    correct_fibres(f"{base_dir}/{mesh_path}/BiV") 

def _surf2vol(meshname, meshname_uac, endo_fibres, epi_fibres, endo_epi_laplace, output, debug, help=False) :
    """
    Mode of operation: surf2vol
    Convert surfaces to volume, including fibres orientations using an endo-to-epi laplace field

    Parameters:
    

    """


def main(args):

    logging.info("Starting 4ch docker container...")

    mode=args.operation
    myhelp=args.help
    output=args.output_dir
    debug=args.debug
    base_dir=args.dev_base_dir
    codes_d=args.dev_code_dir
    local=args.dev_run_local

    if mode == "surfs":
        input_tags=args.input_tags_setup
        apex_septum=args.apex_septum_setup
        meshname=args.meshname
        _surfs(base_dir, input_tags, apex_septum, meshname, debug, help=myhelp)
    
    elif mode == "correctfibres": 
        mesh_path=args.mesh_path
        _correct_fibres(base_dir, mesh_path, debug, help=myhelp)

    elif mode == "surf2vol":
        lon_file = args.sv2_fibres
        lon_file += ".lon" if not lon_file.endswith(".lon") else ""
        lon_location = lon_file[0:lon_file.find(".lon")]
        atrium = args.s2v_atrium
        layer = args.s2v_layer

        meshname = f"{base_dir}/{args.mesh_path}/{atrium}/{args.meshname}" # meshname = la 
        meshname_uac = f"{base_dir}/{args.mesh_path}/{atrium.upper()}_endo/{lon_location}" 
        endo_fibres = f"{base_dir}/{args.mesh_path}/{atrium.upper()}_endo/{lon_location}"
        


   
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser(prog="docker run --rm --volume=/path/to/data:/data cemrg/uac:TAG",
                                           description="4ch docker container entrypoint",
                                           usage="%(prog)s [surfs|correctfibres|surf2vol|latfield|stim|labels|vis|getparfile] [options]",
                                           epilog="$ docker run cemrg/4ch:TAG MODE help \n# for specific help about the operation mode")
    input_parser.add_argument("operation",
                              metavar="mode_of_operation",
                              choices=["surfs", "correctfibres", "surf2vol", "latfield",
                                       "labels", "stim", "vis", "getparfile"],
                              type=str, help="Modes of operation [surfs|correctfibres|surf2vol|latfield|stim|labels|vis|getparfile]")

    input_parser.add_argument("help", nargs='?', type=bool, default=False, help="Help page specific to each mode")

    input_parser.add_argument("--input-tags-setup", metavar="input_tags_setup", nargs='?', type=str)
    input_parser.add_argument("--apex-septum-setup", metavar="apex_septum_setup", nargs='?', type=str)
    input_parser.add_argument("--mesh-path", metavar="meshPath", nargs='?', type=str) 
    input_parser.add_argument("--meshname", metavar="meshName", nargs='?', type=str, default="myocardium")
    input_parser.add_argument("--fch-path", metavar="fchPath", nargs='?', type=str)

    input_parser.add_argument("--s2v-atrium", metavar="atrium", choices=['la', 'ra'], nargs='?', type=str)
    input_parser.add_argument("--s2v-layer", metavar="Layer", choices=['endo', 'epi'],  nargs='?', type=str)
    input_parser.add_argument("--s2v-fibres", metavar="fibre_file.lon", nargs='?', type=str, help="lon file with fibres")

    input_parser.add_argument("--output-dir", metavar="OutputName", nargs='?', type=str, default="")
    input_parser.add_argument("--debug", action='store_true', help="Only show command to run")

    input_parser.add_argument("--dev-base-dir", "-bdir", metavar="dir", nargs='?', default='/data', type=str, help="(only DEVs) Data path")
    input_parser.add_argument("--dev-code-dir", "-code",  metavar="dir", nargs='?', default='/code', type=str, help="(only DEVs) Code path")
    input_parser.add_argument("--dev-run-local", "-local",  action='store_true', help="(only DEVs) Run locally toggle")

    args = input_parser.parse_args()
    
    main(args)