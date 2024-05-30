import argparse
import os

DIR = os.path.dirname(os.path.realpath(__file__))
os.sys.path.append(os.path.join(DIR, ".."))

from common_4ch.file_utils import pjoin
from common_4ch.config import configure_logging
milog = configure_logging(log_name="cemrg/4ch")

def _surfs(args, debug=False, help=False): 
    """
    Mode of operation: surfs
    Extract surfaces from the mesh using meshtool

    Parameters:
    --par-folder: relative path to the folder containing the par files 
    --input-tags-setup: relative path to the json file containing the input tags
    --apex-septum-setup: relative path to the folder containing the apex and septum tags
    --meshname: name of the mesh file (default: meshing/myocardium_OUT/myocardium)

    Usage:
    surfs --par-folder parfiles --input-tags-setup file.json --apex-septum-setup apex_septum_setup
    """

    if help :
        print(_surfs.__doc__)
        return

    if debug :  milog.debug(f"Importing extract_surfs function")
    from common_4ch.process_handler import extract_surfs

    directory = args.base_dir
    input_tags_setup = pjoin(directory, args.par_folder, args.input_tags_setup)
    apex_septum_setup = pjoin(directory, args.par_folder, args.apex_septum_setup)

    extract_surfs(directory, input_tags_setup, apex_septum_setup, args.meshname, debug)

def _correct_fibres(args, debug=False, help=False) :
    """
    Mode of operation: correctfibres
    Correct the fibres

    Parameters:
    --mesh-path: subfolder to mesh

    Usage:
    correctfibres --mesh-path sub/dir/to/meshname ## do not use extension!

    """
    if help :
        print(_correct_fibres.__doc__)
        return
    
    if debug :  milog.debug(f"Importing correct_fibres function")
    from common_4ch.process_handler import correct_fibres 

    directory = args.base_dir
    mesh_path = args.mesh_path

    correct_fibres(f"{directory}/{mesh_path}") 

def _surf2vol(args, debug=False, help=False) :
    """
    Mode of operation: surf2vol
    Convert surfaces to volume, including fibres orientations using an endo-to-epi laplace field.
    

    Parameters:
    --atrium: name of the atrium {'la', 'ra'}
    --fibres-endo: name of the endocardial fibres file (without extension)
    --fibres-epi: name of the epicardial fibres file (without extension)
    --mesh-path: subfolder to UAC mesh folders (e.g 'atrial_fibres/UAC' )

    Usage:
    surf2vol --atrium la --fibres-endo Fibre_endo_l --fibres-epi Fibre_epi_l --mesh-path atrial_fibres/UAC  

    Important Assumptions: 
        + UAC files are stored inside the LA_endo and RA_endo folders.    
        + mesh is stored in f"{base_dir}/{mesh_path}/{atrium}/{meshname}", e.g /data/atrial_fibres/UAC/la/la

    """
    if help :
        print(_surf2vol.__doc__)
        return
    
    atrium = args.atrium

    uac_folder = f"{args.base_dir}/{args.mesh_path}"
    uac_output_folder = f"{uac_folder}/{atrium.upper()}_endo"
    
    mesh_path_no_ext = f"{uac_folder}/{atrium}/{args.meshname}"
    uac_mesh_path_no_ext = f"{uac_output_folder}/{args.fibres_endo}"
    endo_fibres_path = f"{uac_mesh_path_no_ext}.lon"
    epi_fibres_path =  f"{uac_output_folder}/{args.fibres_epi}.lon"
    endo_epi_laplace = f"{uac_output_folder}/{atrium}/endo_epi/phie.dat"
    outmeshname = f"{uac_folder}/{atrium}/{atrium}_fibres_l"

    if debug: milog.debug(f"Importing surf2vol function")
    from common_4ch.process_handler import surf_to_volume

    surf_to_volume(mesh_path_no_ext, uac_mesh_path_no_ext, endo_fibres_path, epi_fibres_path, endo_epi_laplace, outmeshname, debug)

def _laplace_prep(args, debug=False, help=False) :
    """
    Mode of operation: laplace_prep 
    Prepare the laplace field for the endo-to-epi fibres orientation

    Parameters:
    --atrium: name of the atrium {'la', 'ra'}
    --surf-endo: name of the endocardial surface file (without extension)
    --surf-epi: name of the epicardial surface file (without extension)
    --mesh-path: subfolder to UAC mesh folders (e.g 'atrial_fibres/UAC' )

    Usage:
    laplace_prep --atrium la --surf-endo LA_endo --surf-epi LA_epi --mesh-path atrial_fibres/UAC 

    Important Assumptions: 
        + UAC files are stored inside the LA_endo and RA_endo folders.    
        + surf files are stored in f"{base_dir}/{mesh_path}/{atrium}/{meshname}", e.g /data/atrial_fibres/UAC/la/la

    """

    if help :
        print(_laplace_prep.__doc__)
        return
    
    atrium = args.atrium

    uac_folder = f"{args.base_dir}/{args.mesh_path}"
    endo_surf_path = f"{uac_folder}/{atrium}/{args.surf_endo}"
    epi_surf_path = f"{uac_folder}/{atrium}/{args.surf_epi}"

    if debug: milog.debug(f"Importing laplace_prep function")
    from common_4ch.process_handler import laplace_preparation

    laplace_preparation(endo_surf_path, epi_surf_path, debug)

def _tags(args, debug=False, help=False) :
    """
    Mode of operation: tags
    Create the tags 

    Parameters:
    --data-subdir: name of the subfolder containing the BiV mesh (e.g surfaces_uvc)
    --mesh-path: subfolder to UAC mesh folders (e.g 'atrial_fibres' )
    --meshname: name of the mesh (without extension, e.g myocardium_fibres_l)
    --par-folder: name of the subfolder containing the paraview files (e.g parfiles)
    --input-tags-setup: name of the input tags setup JSON file (e.g tags_presim.json)
    --input-bb-settings: name of the input bb settings JSON file (e.g bachmann_bundle_fec_settings.json)

    Usage:
    tags --data-subdir surfaces_uvc --mesh-path atrial_fibres --meshname myocardium_fibres_l --par-folder parfiles --input-tags-setup tags_presim.json --input-bb-settings bachmann_bundle_fec_settings.json
    
    Important Assumptions: 
        + We assume there's a BiV, la, and ra folders inside the directory/data-subdir.
    """

    if help :
        print(_tags.__doc__)
        return
    
    presim_folder = f"{args.base_dir}/pre_simulation"
    biv_folder = f"{args.base_dir}/{args.data_subdir}/BiV"
    la_folder = f"{args.base_dir}/{args.data_subdir}/la"
    ra_folder = f"{args.base_dir}/{args.data_subdir}/ra"

    mesh_path_no_ext = f"{args.base_dir}/{args.mesh_path}/{args.meshname}"
    input_tags = pjoin(args.par_folder, args.input_tags_setup)
    input_bb = pjoin(args.par_folder, args.input_bb_settings)

    if debug :  milog.debug(f"Importing create_tags function")
    from common_4ch.process_handler import create_tags

    create_tags(presim_folder, biv_folder, la_folder, ra_folder, mesh_path_no_ext, input_tags, input_bb, debug)

def _surfs_presim(args, debug=False, help=False) :
    """
    Mode of operation: surfs_presim
    Create the surfaces for the pre-simulation

    Parameters:
    --data-subdir: name of the subfolder containing the BiV mesh (e.g surfaces_uvc)
    --par-folder: name of the subfolder containing the paraview files (e.g parfiles)
    --input-tags-setup: name of the input tags setup JSON file (e.g tags_presim.json)
    --map-settings: name of the map settings JSON file (e.g map_settings.json)
    --fch-apex: name of the file containing the apex coordinates (e.g fch_apex.txt)
    --fch-sa: name of the file containing the sa coordinates (e.g fch_sa.txt)

    Usage:
    surfs_presim --data-subdir surfaces_uvc --par-folder parfiles --input-tags-setup tags_presim.json --map-settings map_settings.json --fch-apex fch_apex.txt --fch-sa fch_sa.txt

    Important Assumptions:
        + We assume there are other subfolders with name {subfolder}_LA and {subfolder}_RA inside the directory. 
        
    """
    if help:
        print(_surfs_presim.__doc__)
        return
    
    la_folder = f"{args.base_dir}/{args.data_subdir}_LA"
    ra_folder = f"{args.base_dir}/{args.data_subdir}_RA"

    input_tags_setup = pjoin(args.par_folder, args.input_tags_setup)
    map_settings = pjoin(args.par_folder, args.map_settings)
    fch_apex = pjoin(args.par_folder, args.fch_apex)
    fch_sa = pjoin(args.par_folder, args.fch_sa)

    if debug :  milog.debug(f"Importing surf_presim function")
    from common_4ch.process_handler import surf_presim

    surf_presim(args.base_dir, la_folder, ra_folder, input_tags_setup, map_settings, fch_apex, fch_sa, args.dev_code_dir, debug)

def _fec(args, debug=False, help=False) :
    """
    Mode of operation: fec
    Split the fast endocardial conduction zone (FEC))
    This allows material/fibre properties to be assigned to the LV and RV independently

    Parameters:
    --meshname: name of the mesh (without extension, e.g myocardium)
    --mesh-path: subfolder to UAC mesh folders (e.g 'meshing/myocardium_OUT' )
    --par-folder: name of the subfolder containing the paraview files (e.g parfiles)
    --input-tags-setup: name of the input tags setup JSON file (e.g tags_presim.json)
    --lvrv-tags: name of the map settings JSON file (e.g map_settings.json)

    Usage:
    fec --meshname myocardium --mesh-path meshing/myocardium_OUT --par-folder parfiles --input-tags-setup tags_presim.json --lvrv-tags lvrv_tags.json 
        
    """
    if help:
        print(_fec.__doc__)
        return
    
    mesh_path_no_ext = f"{args.base_dir}/{args.mesh_path}/{args.meshname}"
    input_tags_setup = pjoin(args.par_folder, args.input_tags_setup)
    lvrv_tags = pjoin(args.par_folder, args.lvrv_tags)

    if debug :  milog.debug(f"Importing split_fec function")
    from common_4ch.process_handler import split_fec

    split_fec(mesh_path_no_ext, input_tags_setup, lvrv_tags, debug)

def _landmarks(args, debug=False, help=False) :
    """
    Mode of operation: landmarks
    Create the landmarks

    Parameters:
    --meshname: name of the mesh (without extension, e.g myocardium)
    --surface: name of the surface {'endo', 'epi'}
    --par-folder: name of the subfolder containing the paraview files (e.g parfiles)
    --input-tags-setup: name of the input tags setup JSON file (e.g tags_presim.json)
    --raa-apex-file: relative path to base_dir of the file containing the apex coordinates (e.g raa_apex.txt)
    --output-folder: relative path to base_dir of the output folder

    Usage:
    landmarks --meshname myocardium --mesh-path meshing/myocardium_OUT --par-folder parfiles --input-tags-setup tags_presim.json --lvrv-tags lvrv_tags.json 
        
    """
    if help:
        print(_landmarks.__doc__)
        return
    
    if debug :  milog.debug(f"Importing main_mesh function")
    from common_4ch.process_handler import main_mesh 

    input_tags = pjoin(args.base_dir, args.par_folder, args.input_tags_setup)
    raa_apex_file = pjoin(args.base_dir, args.raa_apex_file)
    if args.output_folder == "":
        output_folder = pjoin(args.base_dir, "atrial_fibres/UAC")
    
    output_folder = pjoin(args.base_dir, args.output_folder)

    # main_mesh(base_dir, meshname, surface, input_tags_setup, raa_apex_file, output_folder, debug=False)
    main_mesh(args.base_dir, args.meshname, args.surface, input_tags, raa_apex_file, output_folder, debug=False)

def main(args):

    milog.info("Starting 4ch docker container...")

    mode=args.operation
    myhelp=args.help
    mydebug=args.debug

    not_supported_yet = []
    if mode in not_supported_yet:
        milog.warning(f"\n\nMode [{mode}] is not supported yet")
        return 

    if mode == "surfs":
        _surfs(args, mydebug, myhelp)
    
    elif mode == "correctfibres": 
        _correct_fibres(args, mydebug, myhelp)

    elif mode == "surf2vol":
        _surf2vol(args, mydebug, myhelp)
    
    elif mode == "laplace_prep":
        _laplace_prep(args, mydebug, myhelp)
    
    elif mode == "tags": 
        _tags(args, mydebug, myhelp)
    
    elif mode == "presim": 
        _surfs_presim(args, mydebug, myhelp)

    elif mode == "fec":
        _fec(args, mydebug, myhelp)

    elif mode == "landmarks":
        _landmarks(args, mydebug, myhelp)
   
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser(prog="docker run --rm --volume=/path/to/data:/data cemrg/4ch:TAG",
                                           description="4ch docker container entrypoint",
                                           usage="%(prog)s [surfs|correctfibres|surf2vol|laplace_prep|landmarks|tags|presim|fec] [options]",
                                           epilog="$ docker run cemrg/4ch:TAG MODE help \n# for specific help about the operation mode")
    input_parser.add_argument("operation",
                              metavar="mode_of_operation",
                              choices=["surfs", "correctfibres", "surf2vol", "laplace_prep",
                                       "presim", "tags", "fec", "landmarks"],
                              type=str, help="Modes of operation [surfs|correctfibres|surf2vol|laplace_prep|tags|landmarks|presim|fec]")

    input_parser.add_argument("help", nargs='?', type=bool, default=False, help="Help page specific to each mode")

    common_group = input_parser.add_argument_group("Common options", description="Options available two or more modes")
    common_group.add_argument("--meshname", metavar="meshname", nargs='?', type=str, default="meshing/myocardium_OUT/myocardium")
    common_group.add_argument("--mesh-path", metavar="meshpath", nargs='?', type=str, help="Relative folder to base dir where mesh is stored") 
    common_group.add_argument("--data-subdir", metavar="datasubfolder", nargs='?', type=str, help="Subfolder containing data files (modes: tags, presim)")
    common_group.add_argument("--output", metavar="filename", nargs='?', type=str, default="")

    common_group.add_argument("--atrium", metavar="option", choices=['la', 'ra'], nargs='?', type=str, help="Atrium option (modes: surf2vol, laplace_prep)")
    common_group.add_argument("--file-endo", metavar="filename", nargs='?', type=str, help="File specific to endo (.lon, .surf, ...); modes: surf2vol, laplace_prep")
    common_group.add_argument("--file-epi", metavar="filename", nargs='?', type=str, help="File specific to endo (.lon, .surf, ...); modes: surf2vol, laplace_prep")

    parameter_files_group = input_parser.add_argument_group("Parameter files", description="Parameter files for the different modes")
    parameter_files_group.add_argument("--par-folder", metavar="parameter_files", nargs='?', type=str, help="Subfolder containing parameter files")
    parameter_files_group.add_argument("--apex-septum-setup", metavar="parameter_files", nargs='?', type=str)    
    parameter_files_group.add_argument("--input-tags-setup", metavar="parameter_files", nargs='?', type=str)
    parameter_files_group.add_argument("--bb-settings", metavar="parameter_files", nargs='?', type=str, help="Bachmann bundle settings file")
    parameter_files_group.add_argument("--map-settings", metavar="parameter_files", nargs='?', type=str, help="Map settings file")
    parameter_files_group.add_argument("--fch-apex", metavar="parameter_files", nargs='?', type=str, help="Apex file")
    parameter_files_group.add_argument("--fch-sa", metavar="parameter_files", nargs='?', type=str, help="Septal annulus file")
    parameter_files_group.add_argument("--lvrv-tags", metavar="parameter_files", nargs='?', type=str, help="JSON file with input tags settings for split FEC")

    laplace_prep_group = input_parser.add_argument_group("laplace_prep options", description="Options for laplace_prep mode")
    laplace_prep_group.add_argument("--surf-endo", metavar="filename", nargs='?', type=str, help="Endocardial surface file")
    laplace_prep_group.add_argument("--surf-epi", metavar="filename", nargs='?', type=str, help="Epicardial surface file")

    landmarks_group = input_parser.add_argument_group("landmarks options", description="Options for landmarks mode")
    landmarks_group.add_argument("--surface", choices=['endo', 'epi'], nargs='?', type=str, help="Surface option")
    landmarks_group.add_argument("--raa-apex-file", metavar="filename", nargs='?', type=str, help="RAA apex file")

    dev_group = input_parser.add_argument_group("Developer options", description="Options only available to developers")
    dev_group.add_argument("--debug", action='store_true', help="Only show command to run")
    dev_group.add_argument("--base-dir", "-bdir", metavar="dev", nargs='?', default='/data', type=str, help="(only DEVs) Data path")
    dev_group.add_argument("--code-dir", "-code",  metavar="dev", nargs='?', default='/code', type=str, help="(only DEVs) Code path")
    dev_group.add_argument("--run-local", "-local", action='store_true', help="(only DEVs) Run locally toggle")

    args = input_parser.parse_args()
    
    main(args)
