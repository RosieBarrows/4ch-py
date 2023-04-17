import logging
import numpy as np
import copy
import os
from tqdm import tqdm

from common_4ch.json_utils import *
from common_4ch.meshtools_utils import *

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
    os.system(cmd, debug)
