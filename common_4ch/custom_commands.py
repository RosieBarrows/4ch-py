import logging
import numpy as np
import copy
import os
from tqdm import tqdm

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *
import common_4ch.file_utils as fu

FORMAT = '[%(funcName)s:%(levelname)s] %(message)s'
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(FORMAT)
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def my_cp(src, dst, debug=False):
    if debug:
        logger.debug(f"cp {src} {dst}")
    os.system(f"cp {src} {dst}")


def my_mkdir(folder_path, create=True, debug=False):
    if debug:
        logger.debug(f"mkdir -p {folder_path}")
    res = False
    if os.path.exists(folder_path):
        res = True
    elif create:
        os.makedirs(folder_path)
        res = True
    return res

def my_rm(file_path, debug=False):
    import glob 
    if debug:
        logger.debug(f"rm {file_path}")

    flie_list = glob.glob(file_path)
    for f in flie_list:
        try :
            os.remove(f)
        except OSError:
            logger.error(f"Error while deleting file {f}")

def correct_fibres(meshname):
    """Corrects fibre orientation in the mesh"""
    logger.info('Reading mesh...')
    lon = np.loadtxt(f"{meshname}.lon",dtype=float,skiprows=1)
    # pts = np.loadtxt(f"{meshname}.pts",dtype=float,skiprows=1)
    # elem = np.loadtxt(f"{meshname}.elem",dtype=int,skiprows=1,usecols=[1,2,3,4])

    logger.info('Reading element centres...')
    elemC = np.loadtxt(f"{meshname}+_elem_centres.pts",dtype=float,skiprows=1)

    to_correct = np.where(np.abs(lon[:,2])<1e-6)[0]
    logger.info(f"Found: {str(to_correct.shape[0])} elements to correct")

    lon_corrected = copy.deepcopy(lon)
    good_elements = np.setdiff1d(np.arange(elemC.shape[0]), to_correct, assume_unique=True)

    for idx in tqdm(to_correct, desc="Correcting fibre orientation..."):
        d = np.linalg.norm(elemC[good_elements,:] - elemC[idx,:],axis=1)
        closest = good_elements[np.where(d == np.min(d))[0]]
        lon_corrected[idx,:] = lon[closest,:]
    to_correct = np.where(np.abs(lon_corrected[:,2])<1e-6)[0]

    logger.info(f"Left to correct {str(to_correct.shape[0])}")
    logger.info('Correcting fibre orientation...')

    np.savetxt(f"{meshname}_corrected.lon",lon_corrected,fmt="%g",header='2',comments='')


def extract_surfs(base_dir, input_tags_setup, apex_septum_setup, meshname="myocardium", debug): 
    """Extract surfaces from the mesh using meshtool libraries"""

    mesh = f"{base_dir}/meshing/myocardium_OUT/{meshname}"
    surf_folder = f"{base_dir}/surfaces_uvc/"
    surf_folder_la = f"{base_dir}/surfaces_uvc_LA/"
    surf_folder_ra = f"{base_dir}/surfaces_uvc_RA/"

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
    meshtool_extract_ra_base(mesh, surf_folder_ra, input_tags_setup)

    logger.info("Extracting the RA epi and RA endo")
    meshtool_extract_ra_surfaces(mesh, surf_folder_ra, input_tags_setup)

    logger.info("Mapping surfaces RA")
    mapping_surfaces_ra(mesh, surf_folder_ra, input_tags_setup)

    # Extracting the BiV mesh
    logger.info("Extracting the biventricular mesh ")
    meshtool_extract_biv(mesh, surf_folder, input_tags_setup)

    logger.info("Mapping vtx files from four-chamber mesh to BiV mesh ")
    meshtool_map_vtx(surf_folder)

    logger.info("Renaming files ")
    renaming_myo_files(surf_folder)

    # Extracting the LA mesh
    logger.info("Extracting the left atrial mesh")
    meshtool_extract_la_for_UVCs(mesh, surf_folder_la, input_tags_setup)

    logger.info("Mapping vtx files from four-chamber mesh to left atrial mesh ")
    meshtool_map_vtx_la(surf_folder_la)

    logger.info("Copying blank files for LA apex and septum ")
    my_cp(f"{apex_septum_setup}/la.lvapex.vtx",f"{surf_folder_la}/la/la.lvapex.vtx", debug)
    my_cp(f"{apex_septum_setup}/la.rvsept_pt.vtx", f"{surf_folder_la}/la/la.rvsept_pt.vtx", debug)

def laplace_preparation(endo_surf_path, epi_surf_path, debug) : 
    """Prepares the mesh for the Laplace equation (outside of docker container)""" 
    logger.info("Converting surfaces to vtx files")

    def to_vtx(surf_path, debug) :
        if debug : logger.debug(f"Converting {surf_path} to vtx")
        surf = read_elem(surf_path, el_type='Tr', tags=False)
        vtx = surf2vtx(surf)

        logger.info(f"Writing {surf_path}.vtx")
        fu.write_vtx(f"{surf_path}.vtx", vtx, init_row=2) 
    
    to_vtx(endo_surf_path, debug)
    to_vtx(epi_surf_path, debug)

def surf_to_volume(mesh_path_no_ext, uac_mesh_path_no_ext, endo_fibres_path, epi_fibres_path, endo_epi_laplace, output, debug) :
    """Converts the surface mesh to a volume mesh using the UAC library"""
    meshname_3d = mesh_path_no_ext
    meshname_2d = uac_mesh_path_no_ext

    logger.info(f"Converting surface to volume mesh {meshname_3d}")

    logger.info("Mapping endo-epi Laplace solution from nodes to elements")
    laplace_endo2elem(meshname_3d, endo_epi_laplace)

    logger.info("Computing element centres on both meshes (3D/2D)")
    compute_elemCenters(meshname_3d, el_type="Tt")
    compute_elemCenters(meshname_2d, el_type="Tr")

    logger.info("Mapping fibres")
    map_fibres_3d(endo_epi_laplace[:-4]+"_el.dat",
                  meshname_3d+"_elemCenters.pts",
                  meshname_2d+"_elemCenters.pts",
                  endo_fibres_path,
                  epi_fibres_path,
                  output+".lon",
                  map_tr_tt_file=f"{meshname_3d}_endo_to_3d.map")

    logger.info("Finding sheet direction")
    find_transmural_direction(meshname_3d,
                              meshname_2d,
                              meshname_3d+"_elemCenters.pts",
                              meshname_2d+"_elemCenters.pts",
                              meshname_3d+"_transmural.lon")


    find_rotation_axes(output+".lon",
                       meshname_3d+"_transmural.lon",
                       meshname_3d+"_rotation_axes.lon")

    make_sheet_orthogonal(output+".lon",
                          meshname_3d+"_transmural.lon",
                          meshname_3d+"_rotation_axes.lon",
                          output+"_sheet.lon")

    logger.info("Converting to VTK for visualisation")
    my_cp(f"{meshname_3d}.elem", f"{output}.elem", debug)
    my_cp(f"{meshname_3d}.pts", f"{output}.pts", debug)

    cmd = f"meshtool convert -imsh={output}_sheet -omsh={output}_sheet.vtk"
    if debug: logger.info(f"COMMAND: {cmd}")
    os.system(cmd)

def create_tags(presim_folder, biv_folder, la_folder, ra_folder, mesh_path_no_ext, input_tags, bb_settings, debug) : 
    """Defines the tags for the mesh"""
    # presim_folder=f"{directory}/pre_simulation"
    my_mkdir(presim_folder, debug)

    input_tags_json = load_json(input_tags)
    bb_settings_json = load_json(bb_settings)

    msh="myocardium"

    logger.info("Setting up folder structure")
    my_cp(f"{mesh_path_no_ext}.elem", f"{presim_folder}/{msh}.elem", debug)
    my_cp(f"{mesh_path_no_ext}.pts", f"{presim_folder}/{msh}.pts", debug)
    my_cp(f"{mesh_path_no_ext}.lon", f"{presim_folder}/{msh}.lon", debug)

    msh_av = f"{msh}_av"
    msh_fec = f"{msh}_av_fec"
    msh_bb = f"{msh}_av_fec_bb"

    logger.info("Defining the AV separating plane")
    define_AV_separation(f"{presim_folder}/{msh}.elem", input_tags_json, input_tags_json["AV_plane"], f"{presim_folder}/{msh_av}.elem")

    logger.info("Defining the FEC layer")
    biv_mesh = f"{biv_folder}/BiV"
    zbiv_file = f"{biv_folder}/uvc/BiV.uvc_z.dat" 
    rhobiv_file = f"{biv_folder}/uvc/BiV.uvc_rho.dat"

    define_FEC(f"{presim_folder}/{msh_av}.elem", 
               biv_mesh, zbiv_file, rhobiv_file, 
               f"{presim_folder}/{msh_fec}.elem", 
               input_tags_json["FEC"], 
               include_septum=f"{biv_folder}/BiV.rvsept.surf", 
               FEC_height=bb_settings_json["FEC_height"])
    
    logger.info("Defining the BB area")
    la_mesh=f"{la_folder}/la"
    ra_mesh=f"{ra_folder}/ra"
    zla_file=f"{la_folder}/uvc/la.uvc_z.dat"
    zra_file=f"{ra_folder}/uvc/ra.uvc_z.dat"
    phila_file=f"{la_folder}/uvc/la.uvc_phi.dat"
    phira_file=f"{ra_folder}/uvc/ra.uvc_phi.dat"

    define_BB(f"{presim_folder}/{msh_fec}.elem", 
                la_mesh, ra_mesh, 
                zla_file, zra_file, phila_file, phira_file,
                bb_settings_json, 
                input_tags_json, 
                f"{presim_folder}/{msh_bb}.elem")
    
    my_cp(f"{presim_folder}/{msh}.pts", f"{presim_folder}/{msh_bb}.pts", debug)
    my_cp(f"{presim_folder}/{msh}.lon", f"{presim_folder}/{msh_bb}.lon", debug)

    cmd = f"meshtool convert -imsh={presim_folder}/{msh_bb} -omsh={presim_folder}/{msh_bb}.vtk"
    if debug: logger.info(f"COMMAND: {cmd}")
    os.system(cmd)

def surf_presim(directory, la_folder, ra_folder, input_tags_settings, map_settings, fch_apex, fch_sa, code_d, debug) :
    """Surfaces for pre-simulation"""
    input_tags_json = load_json(input_tags_settings)
    presim_folder=f"{directory}/pre_simulation"
    msh = "myocardium_AV_FEC_BB"

    msh_path = f"{presim_folder}/{msh}"

    logger.info("Extracting the surface for the pericardium boundary condition")
    meshtool_extract_peri(msh_path, presim_folder, input_tags_json)
    connected_component_to_surface(f"{presim_folder}/peri_surface_CC.part1", f"{presim_folder}/peri_surface.surf", f"{presim_folder}/epicardium_for_sim")
    surf2vtk(msh_path, f"{presim_folder}/epicardium_for_sim.surf", f"{presim_folder}/epicardium_for_sim.vtk") 

    my_rm(f"{presim_folder}/*CC*", debug)

    logger.info("Extracting the epi and LV/RV/LA/RA endo surfaces")
    my_mkdir(f"{presim_folder}/surfaces_simulation", debug) 
    meshtool_extract_epi_endo_surfs(msh_path,presim_folder,input_tags_json)

    my_rm(f"{presim_folder}/surfaces_simulation/surface_heart_CC.*", debug)

    logger.info("Extracting the surfaces of the rings")
    my_mkdir(f"{presim_folder}/surfaces_simulation/surfaces_rings", debug) 
    meshtool_extract_rings(msh_path, presim_folder, input_tags_json)

    logger.info("Setting up the pericardium scale...")
    set_pericardium(msh_path, presim_folder, directory)

    logger.info("Setting up the pericarsium for the atria ...")
    la_mesh=f"{la_folder}/la"
    uvcs=f"{la_mesh}/uvc/" 

    cmd=f"python3 -u {code_d}/motion_atria_BCs.py --mesh {la_mesh} --uvcs {uvcs} --chamber la --map_settings {map_settings}"
    if debug: logger.info(f"COMMAND: {cmd}")
    os.system(cmd)

    ra_mesh=f"{ra_folder}/ra"
    uvcs=f"{ra_mesh}/uvc/"

    cmd=f"python3 -u {code_d}/motion_atria_BCs.py --mesh {ra_mesh} --uvcs {uvcs} --chamber ra --map_settings {map_settings}"
    if debug: logger.info(f"COMMAND: {cmd}")
    os.system(cmd)

    logger.info("Combining the pericardium scaling maps for the ventricles and atria")
    combine_elem_dats(directory, presim_folder)

    logger.info("Setting up a folder with the simulation-ready mesh") 
    setup_sim(directory, presim_folder, fch_apex, fch_sa)

def split_fec(directory, input_tags_setup, lvrv_tags, mesh_path_no_ext, debug=False) : 
    """
    Split the fast endocardial conduction zone (FEC))
    This allows material/fibre properties to be assigned to the LV and RV independently
    """
    original_tags_json = load_json(input_tags_setup)
    new_tags_json = load_json(lvrv_tags)

    sim_dir = f"{directory}/sims_folder"
    mshname = "myocardium_AV_FEC_BB"
    separate_FEC_lvrv(mesh_path_no_ext, f"{sim_dir}/{mshname}.elem", 
                      f"{sim_dir}/LV_endo.surf", f"{sim_dir}/RV_endo.surf",
                      f"{sim_dir}/{mshname}_lvrv.elem", original_tags_json, new_tags_json)
    
    my_cp(f"{sim_dir}/{mshname}.pts", f"{sim_dir}/{mshname}_lvrv.pts", debug)
    my_cp(f"{sim_dir}/{mshname}.lon", f"{sim_dir}/{mshname}_lvrv.lon", debug)

    cmd = f"meshtool convert -imsh={sim_dir}/{mshname}_lvrv -omsh={sim_dir}/{mshname}_lvrv.vtk"
    if debug: logger.info(f"COMMAND: {cmd}")
    os.system(cmd)

    
