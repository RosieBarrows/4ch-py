import argparse
import logging

from common_4ch.config import configure_logging
milog = configure_logging(log_name="cemrg/4ch")

def _surfs(directory, input_tags_setup, apex_septum_setup, meshname="meshing/myocardium_OUT/myocardium", debug=False, help=False): 
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

    extract_surfs(directory, input_tags_setup, apex_septum_setup, meshname, debug=False)

def _correct_fibres(directory, mesh_path, debug=False, help=False) :
    """
    Mode of operation: correctfibres
    Correct the fibres

    Parameters:
    --mesh-path: subfolder to mesh

    Usage:
    correctfibres --mesh-path sub/dir/to/folder

    """
    if help :
        print(_correct_fibres.__doc__)
        return
    
    if debug :  milog.debug(f"Importing correct_fibres function")
    from common_4ch.process_handler import correct_fibres 

    correct_fibres(f"{directory}/{mesh_path}/BiV") 

def _surf2vol(directory, atrium, fibres_endo, fibres_epi, mesh_path, debug=False, help=False) :
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

    uac_folder = f"{directory}/{mesh_path}"
    uac_output_folder = f"{uac_folder}/{atrium.upper()}_endo"
    
    mesh_path_no_ext = f"{uac_folder}/{atrium}/{args.meshname}"
    uac_mesh_path_no_ext = f"{uac_output_folder}/{fibres_endo}"
    endo_fibres_path = f"{uac_mesh_path_no_ext}.lon"
    epi_fibres_path =  f"{uac_output_folder}/{fibres_epi}.lon"
    endo_epi_laplace = f"{uac_output_folder}/{atrium}/endo_epi/phie.dat"
    outmeshname = f"{uac_folder}/{atrium}/{output}"

    if debug: milog.debug(f"Importing surf2vol function")
    from common_4ch.process_handler import surf_to_volume

    surf_to_volume(mesh_path_no_ext, uac_mesh_path_no_ext, endo_fibres_path, epi_fibres_path, endo_epi_laplace, outmeshname, debug)

def _laplace_prep(directory, atrium, surf_endo, surf_epi, mesh_path, debug=False, help=False) :
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

    uac_folder = f"{directory}/{mesh_path}"
    endo_surf_path = f"{uac_folder}/{atrium}/{surf_endo}"
    epi_surf_path = f"{uac_folder}/{atrium}/{surf_epi}"

    if debug: milog.debug(f"Importing laplace_prep function")
    from common_4ch.process_handler import laplace_preparation

    laplace_preparation(endo_surf_path, epi_surf_path, debug)

def _tags(directory, subfolder, mesh_path, meshname, input_tags_setup, input_bb_settings, debug=False, help=False) :
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
    
    presim_folder = f"{directory}/pre_simulation"
    biv_folder = f"{directory}/{subfolder}/BiV"
    la_folder = f"{directory}/{subfolder}/la"
    ra_folder = f"{directory}/{subfolder}/ra"

    mesh_path_no_ext = f"{directory}/{mesh_path}/{meshname}"


    if debug :  milog.debug(f"Importing create_tags function")
    from common_4ch.process_handler import create_tags

    create_tags(presim_folder, biv_folder, la_folder, ra_folder, mesh_path_no_ext, input_tags_setup, input_bb_settings, debug)

def _surfs_presim(directory, subfolder, input_tags_setup, map_settings, fch_apex, fch_sa, code_d="/code", debug=False, help=False) :
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
    
    la_folder = f"{directory}/{subfolder}_LA"
    ra_folder = f"{directory}/{subfolder}_RA"

    if debug :  milog.debug(f"Importing surf_presim function")
    from common_4ch.process_handler import surf_presim

    surf_presim(directory, la_folder, ra_folder, input_tags_setup, map_settings, fch_apex, fch_sa, code_d, debug)

def _fec(directory, mesh_path, meshname, input_tags_setup, lvrv_tags, debug=False, help=False) :
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
    
    mesh_path_no_ext = f"{directory}/{mesh_path}/{meshname}"

    if debug :  milog.debug(f"Importing split_fec function")
    from common_4ch.process_handler import split_fec

    split_fec(mesh_path_no_ext, input_tags_setup, lvrv_tags, debug)


def main(args):

    milog.info("Starting 4ch docker container...")

    mode=args.operation
    myhelp=args.help
    output=args.output
    debug=args.debug
    base_dir=args.dev_base_dir # default=/data
    codes_d=args.dev_code_dir
    local=args.dev_run_local

    if mode == "surfs":
        par_folder=args.par_folder
        
        input_tags=pjoin(base_dir, par_folder, args.input_tags_setup)
        apex_septum=pjoin(base_dir, par_folder, args.apex_septum_setup)
        _surfs(base_dir, input_tags, apex_septum, args.meshname, debug=False, help=myhelp)
    
    elif mode == "correctfibres": 
        mesh_path=args.mesh_path
        _correct_fibres(base_dir, mesh_path, debug=False, help=myhelp)

    elif mode == "surf2vol":
        
        atrium=args.atrium
        mesh_path=args.mesh_path
        fibres_endo=args.file_endo
        fibres_epi=args.file_epi

        arg_list = [atrium, mesh_path, fibres_endo, fibres_epi]
        if None in arg_list:
            milog.error("Please provide the endo and epi fibre files")
            myhelp = True 
        else :
            if fibres_endo.endswith('.lon'):
                fibres_endo = fibres_endo[:-4]
            if fibres_epi.endswith('.lon'):
                fibres_epi = fibres_epi[:-4]

        _surf2vol(base_dir, atrium, fibres_endo, fibres_epi, mesh_path, debug=False, help=myhelp)
    
    elif mode == "laplace_prep":

        atrium=args.atrium
        mesh_path=args.mesh_path
        surf_endo=args.file_endo
        surf_epi=args.file_epi

        arg_list = [atrium, mesh_path, surf_endo, surf_epi]
        if None in arg_list:
            milog.error("Please provide the endo and epi surface files")
            myhelp = True
        else :
            if not surf_endo.endswith('.surf'): surf_endo += '.surf'
            if not surf_epi.endswith('.surf'): surf_epi += '.surf'

        _laplace_prep(base_dir, atrium, surf_endo, surf_epi, mesh_path, debug=False, help=myhelp)
    
    elif mode == "tags": 
        meshname=args.meshname # myocardium_fibres_l 
        mesh_path=args.mesh_path # atrial_fibres 
        biv_subfolder=args.data_subdir # surfaces_uvc
        par_folder=args.par_folder # parfiles
        input_tags=f"{par_folder}/{args.input_tags_setup}"
        input_bb_settings=f"{par_folder}/{args.input_bb_settings}"

        _tags(base_dir, biv_subfolder, mesh_path, meshname, input_tags, input_bb_settings, debug=False, help=myhelp)
    
    elif mode == "presim": 
        #_surfs_presim(directory, subfolder, input_tags_setup, map_settings, fch_apex, fch_sa, code_d="/code", debug=False, help=False)
        subfolder = args.data_subdir
        input_tags_setup = f"{args.par_folder}/{args.input_tags_setup}"
        map_settings = f"{args.par_folder}/{args.map_settings}"
        fch_apex = f"{args.par_folder}/{args.fch_apex}"
        fch_sa = f"{args.par_folder}/{args.fch_sa}"

        _surfs_presim(base_dir, subfolder, input_tags_setup, map_settings, fch_apex, fch_sa, codes_d, debug=False, help=myhelp)

    elif mode == "fec":
        subfolder = args.mesh_path
        meshname = args.meshname
        input_tags_setup = f"{args.par_folder}/{args.input_tags_setup}"
        lvrv_tags = f"{args.par_folder}/{args.lvrv_tags}"

        _fec(base_dir, subfolder, meshname, input_tags_setup, lvrv_tags, debug=False, help=myhelp)



   
if __name__ == '__main__':
    input_parser = argparse.ArgumentParser(prog="docker run --rm --volume=/path/to/data:/data cemrg/uac:TAG",
                                           description="4ch docker container entrypoint",
                                           usage="%(prog)s [surfs|correctfibres|surf2vol|laplace_prep|tags|presim|fec] [options]",
                                           epilog="$ docker run cemrg/4ch:TAG MODE help \n# for specific help about the operation mode")
    input_parser.add_argument("operation",
                              metavar="mode_of_operation",
                              choices=["surfs", "correctfibres", "surf2vol", "laplace_prep",
                                       "presim", "tags", "fec"],
                              type=str, help="Modes of operation [surfs|correctfibres|surf2vol|laplace_prep|tags|presim|fec]")

    input_parser.add_argument("help", nargs='?', type=bool, default=False, help="Help page specific to each mode")

    input_parser.add_argument("--par-folder", metavar="parameter_files", nargs='?', type=str, help="Subfolder containing parameter files")
    input_parser.add_argument("--input-tags-setup", metavar="parameter_files", nargs='?', type=str)
    input_parser.add_argument("--apex-septum-setup", metavar="parameter_files", nargs='?', type=str)
    input_parser.add_argument("--bb-settings", metavar="parameter_files", nargs='?', type=str, help="Bachmann bundle settings file")
    input_parser.add_argument("--map-settings", metavar="parameter_files", nargs='?', type=str, help="Map settings file")
    input_parser.add_argument("--fch-apex", metavar="parameter_files", nargs='?', type=str, help="Apex file")
    input_parser.add_argument("--fch-sa", metavar="parameter_files", nargs='?', type=str, help="Septal annulus file")
    input_parser.add_argument("--lvrv-tags", metavar="parameter_files", nargs='?', type=str, help="JSON file with input tags settings for split FEC")

    input_parser.add_argument("--meshname", metavar="meshname", nargs='?', type=str, default="meshing/myocardium_OUT/myocardium")
    input_parser.add_argument("--mesh-path", metavar="meshpath", nargs='?', type=str) 
    input_parser.add_argument("--data-subdir", metavar="datasubfolder", nargs='?', type=str, help="Subfolder containing data files (see help)")

    input_parser.add_argument("--atrium", metavar="option", choices=['la', 'ra'], nargs='?', type=str)
    input_parser.add_argument("--file-endo", metavar="filename", nargs='?', type=str, help="File specific to endo (.lon, .surf, ...)")
    input_parser.add_argument("--file-epi", metavar="filename", nargs='?', type=str, help="File specific to endo (.lon, .surf, ...)")

    input_parser.add_argument("--output", metavar="filename", nargs='?', type=str, default="")

    input_parser.add_argument("--debug", action='store_true', help="Only show command to run")
    input_parser.add_argument("--dev-base-dir", "-bdir", metavar="dev", nargs='?', default='/data', type=str, help="(only DEVs) Data path")
    input_parser.add_argument("--dev-code-dir", "-code",  metavar="dev", nargs='?', default='/code', type=str, help="(only DEVs) Code path")
    input_parser.add_argument("--dev-run-local", "-local", action='store_true', help="(only DEVs) Run locally toggle")

    args = input_parser.parse_args()
    
    main(args)
