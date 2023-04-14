import os
import argparse
import logging

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *

FORMAT = '[%(funcName)s:%(levelname)s] %(message)s'
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(FORMAT)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

def my_cp(src, dst, debug=False) : 
    if debug :
        logger.debug(f"cp {src} {dst}")
    os.system(f"cp {src} {dst}")

def my_mkdir(folder_path, create=True, debug=False) :
    if debug :
        logger.debug(f"mkdir -p {folder_path}") 
    res = False
    if os.path.exists(folder_path) : 
        res = True
    elif create : 
        os.makedirs(folder_path)
        res = True
    return res

def _surfs(base_dir, input_tags_setup, apex_septum_setup, meshname="myocardium", debug, help=False): 
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

    mesh=f"{base_dir}/meshing/myocardium_OUT/{meshname}"
    surf_folder=f"{base_dir}/surfaces_uvc/"
    surf_folder_la=f"{base_dir}/surfaces_uvc_LA/"
    surf_folder_ra=f"{base_dir}/surfaces_uvc_RA/" 
    
    my_mkdir(surf_folder, debug)
    my_mkdir(f"{surf_folder}/tmp", debug)
    my_mkdir(f"{surf_folder}/BiV", debug)
    my_mkdir(surf_folder_la, debug)
    my_mkdir(f"{surf_folder_la}/tmp", debug)
    my_mkdir(f"{surf_folder_la}/la", debug)
    my_mkdir(surf_folder_ra, debug)
    my_mkdir(f"{surf_folder_ra}/tmp", debug)
    my_mkdir(f"{surf_folder_ra}/ra", debug)

    logger.info("Extracting the base")
    meshtool_extract_base(mesh, surf_folder, input_tags_setup)

    logger.info("Extracting epi, LV endo and RV endo")
    meshtool_extract_surfaces_lv_rv_epi(mesh, surf_folder, input_tags_setup)

    logger.info("Extracting septum")
    meshtool_extract_septum(mesh, surf_folder, input_tags_setup)

    logger.info("Mapping surfaces")
    mapping_surfaces(mesh, surf_folder, input_tags_setup)

    logger.info("Removing the septum")
    remove_sept(mesh, surf_folder)

    logger.info("Preparing vtx files for UVCs")
    prepare_vtx_for_uvc(surf_folder)

    # extract surfacs for LA
    logger.info("Extracting the LA base")
    meshtool_extract_la_base(mesh, surf_folder_la, input_tags_setup)

    logger.info("Extracting the LA epi and LA endo")
    meshtool_extract_la_surfaces(mesh, surf_folder_la, input_tags_setup)

    logger.info("Mapping LA surfaces")
    mapping_surfaces_la(mesh, surf_folder_la, input_tags_setup)

    # extract surfaces for RA
    logger.info("Extracting the RA base")
    meshtool_extract_ra_base(mesh,surf_folder_ra,input_tags_setup)
    
    logger.info("Extracting the RA epi and RA endo")
    meshtool_extract_ra_surfaces(mesh,surf_folder_ra,input_tags_setup)
    
    logger.info("Mapping surfaces RA")
    mapping_surfaces_ra(mesh,surf_folder_ra,input_tags_setup)
    
    # Extracting the BiV mesh
    logger.info("Extracting the biventricular mesh ")
    meshtool_extract_biv(mesh,surf_folder,input_tags_setup)

    logger.info("Mapping vtx files from four-chamber mesh to BiV mesh ")
    meshtool_map_vtx(surf_folder)

    logger.info("Renaming files ")
    renaming_myo_files(surf_folder)

    # Extracting the LA mesh 
    logger.info("Extracting the left atrial mesh")
    meshtool_extract_la_for_UVCs(mesh,surf_folder_la,input_tags_setup)

    logger.info("Mapping vtx files from four-chamber mesh to left atrial mesh ")
    meshtool_map_vtx_la(surf_folder_la)

    logger.info("Copying blank files for LA apex and septum ")
    my_cp(f"{apex_septum_setup}/la.lvapex.vtx", f"{surf_folder_la}/la/la.lvapex.vtx", debug) 
    my_cp(f"{apex_septum_setup}/la.rvsept_pt.vtx",  f"{surf_folder_la}/la/la.rvsept_pt.vtx", debug)

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
    
    from common_4ch.create_commands import correct_fibres 
    if debug : 
        logger.debug(f"Correcting fibres in {mesh_path}")

    correct_fibres(f"{base_dir}/{mesh_path}/BiV") 


def main(args):

    logger.info("Starting 4ch docker container...")

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
        


   
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser(prog="docker run --rm --volume=/path/to/data:/data cemrg/uac:TAG",
                                           description="4ch docker container entrypoint",
                                           usage="%(prog)s [surfs|correctfibres|scalarmap|latfield|stim|labels|vis|getparfile] [options]",
                                           epilog="$ docker run cemrg/4ch:TAG MODE help \n# for specific help about the operation mode")
    input_parser.add_argument("operation",
                              metavar="mode_of_operation",
                              choices=["surfs", "correctfibres", "scalarmap", "latfield",
                                       "labels", "stim", "vis", "getparfile"],
                              type=str, help="Modes of operation [surfs|correctfibres|scalarmap|latfield|stim|labels|vis|getparfile]")

    input_parser.add_argument("help", nargs='?', type=bool, default=False, help="Help page specific to each mode")

    input_parser.add_argument("--input-tags-setup", metavar="input_tags_setup", nargs='?', type=str)
    input_parser.add_argument("--apex-septum-setup", metavar="apex_septum_setup", nargs='?', type=str)
    input_parser.add_argument("--mesh-path", metavar="meshPath", nargs='?', type=str) 
    input_parser.add_argument("--meshname", metavar="meshName", nargs='?', type=str, default="myocardium")
    input_parser.add_argument("--fch-path", metavar="fchPath", nargs='?', type=str)

    input_parser.add_argument("--vfib-a-endo", metavar="VFIBRES_alpha_endo", nargs='?', type=float, default=60)
    input_parser.add_argument("--vfib-a-epi", metavar="VFIBRES_alpha_epi", nargs='?', type=float, default=-60)
    input_parser.add_argument("--vfib-b-endo", metavar="VFIBRES_beta_endo", nargs='?', type=float, default=-65)
    input_parser.add_argument("--vfib-b-epi", metavar="VFIBRES_beta_epi", nargs='?', type=float, default=25)

    input_parser.add_argument("--output-dir", metavar="OutputName", nargs='?', type=str, default="")
    input_parser.add_argument("--debug", action='store_true', help="Only show command to run")

    input_parser.add_argument("--dev-base-dir", "-bdir", metavar="dir", nargs='?', default='/data', type=str, help="(only DEVs) Data path")
    input_parser.add_argument("--dev-code-dir", "-code",  metavar="dir", nargs='?', default='/code', type=str, help="(only DEVs) Code path")
    input_parser.add_argument("--dev-run-local", "-local",  action='store_true', help="(only DEVs) Run locally toggle")

    args = input_parser.parse_args()
    
    main(args)