import os 
import json 

import common_4ch.file_utils as fu
import common_4ch.meshtools_utils as mu
import common_4ch.mesh_utils as mmu
from common_4ch.process_handler import correct_fibres

from common_4ch.config import configure_logging
milog = configure_logging(log_name=__name__) 

SDIR = {
    'mesh' : 'meshing/myocardium_OUT/myocardium',
    'uvc' : 'surfaces_uvc',
    'uvc_la' : 'surfaces_uvc_LA',
    'uvc_ra' : 'surfaces_uvc_RA', 
    'uac' : 'atrial_fibres/UAC'
}
FIBRE_ANGLES = {
    'alpha_endo' : '60',
    'alpha_epi' : '-60',
    'beta_endo' : '-65',
    'beta_epi' : '25'
}
NP=20

class CarpWrapper: 
    def __init__(self, heart_folder) -> None:
        fch_meshname = 'myocardium_bayer'
        _heart_folder = heart_folder
        _angles = FIBRE_ANGLES
        _debug = False
    
    # helper functions 
    def HEART(self, folder):
        return os.path.join(self._heart_folder, folder)
    
    def SURF(self, folder, which):
        if which not in SDIR.keys():
            msg = f"Invalid subfolder identifier: {which}"
            milog.error(msg)
            raise ValueError(msg)
        
        return os.path.join(self.HEART(SDIR[which]), folder)
    
    def SURF_FOLDER(self, folder) : 
        return self.SURF(folder, 'uvc')
    
    def SURF_FOLDER_LA(self, folder) :
        return self.SURF(folder, 'uvc_la')
    
    def SURF_FOLDER_RA(self, folder) :
        return self.SURF(folder, 'uvc_ra')
    
    def run_cmd(self, cmd) :
        milog.info(f"Running command:\n\t{cmd}\n")
        os.system(cmd)

    def set_fibres_angles(self, angles:dict):
        self._angles = angles
    
    def set_fibres_angles(self, alpha_endo, alpha_epi, beta_endo, beta_epi):
        self._angles = {
            'alpha_endo' : alpha_endo,
            'alpha_epi' : alpha_epi,
            'beta_endo' : beta_endo,
            'beta_epi' : beta_epi
        }

    def set_fch_meshname(self, name):
        self.fch_meshname = name

    def get_meshname_fibres(self, name = 'myocardium_bayer', ext=''):
        mshname = f'{name}_{self._angles["alpha_endo"]}_{self._angles["alpha_epi"]}'
        if ext != '':
            ext = f'.{ext}' if '.' not in ext else ext
            mshname += ext
        return mshname
        
    # wraps 
    def w_mguvc(self, model:str, input:str, output:str, tags:str, odir:str, np=NP, laplace_solution=True, id=''):
        id_str = f' --ID={id}' if id else ''

        cmd = f'mguvc {id_str} --model-name {model} --input-model {input} --output-model {output} --np {np} --tags-file {tags} --output-dir {odir}'
        cmd += ' --laplace-solution' if laplace_solution else ''

        self.run_cmd(cmd)

    def w_GlVTKConvert(self, meshname:str, nodedata_list:list, output:str, trim_names = True):
        nodedata = ' -n '.join([n for n in nodedata_list])

        cmd = f'GlVTKConvert -m {meshname} -n {nodedata} -o {output}'
        cmd += ' --trim-names' if trim_names else ''

        self.run_cmd(cmd)
    
    def w_GlRuleFibres(self, meshname:str, uvc_apba:str, uvc_epi:str, uvc_lv:str, uvc_rv:str, output:str, angles_dic=FIBRE_ANGLES, type='biv'):
        angles_str = ' '.join([f'-{k} {v}' for k,v in angles_dic.items()]) # -alpha_endo 60 -alpha_epi -60 -beta_endo -65 -beta_epi 25
        cmd = f'GlRuleFibres -m {meshname} -t {type} -a {uvc_apba} -e {uvc_epi} -l {uvc_lv} -r {uvc_rv} {angles_str} -o {output}'

        self.run_cmd(cmd)

    def w_GlElemCenters(self, meshname:str, output:str):
        cmd = f'GlElemCenters -m {meshname} -o {output}'

        self.run_cmd(cmd)

    def w_carp_pt(self, parfile:str, simID:str, meshname:str, stim0:str, stim1:str):
        cmd = f'carp.pt +F {parfile} -simID {simID} -meshname {meshname} -stimulus[0].vtx_file {stim0} -stimulus[1].vtx_file {stim1}'
        self.run_cmd(cmd)

    def w_igbextract(self, igb:str, out:str):
        cmd = f'igbextract {igb} -o ascii -f 0 -F 0 -O {out}'
        self.run_cmd(cmd)

    # main functions
    def main_uvc_process(self, input_tags: str, etags:str, apex_septum:str, msh_subpath =SDIR['mesh']): 
        meshname = self.HEART(msh_subpath)

        milog.info(" ## Calculating UVCs for the BiV mesh ##")
        biv_etags = self.SURF_FOLDER('BiV/etags.sh')
        model_name = self.SURF_FOLDER('BiV/BiV')
        odir = self.SURF_FOLDER('BiV/uvc/')
        
        fu.mycp(os.path.join(etags, 'etags.sh'), biv_etags)
        self.w_mguvc(model_name, 'biv', 'biv', NP, biv_etags, odir, laplace_solution=True)

        ndata_list = ['BiV/uvc/BiV.uvc_phi.dat', 'BiV/uvc/BiV.uvc_z.dat', 'BiV/uvc/BiV.uvc_ven.dat', 'BiV/uvc/BiV.uvc_rho.dat']
        ndata_list = [self.SURF_FOLDER(n) for n in ndata_list]
        output_name = self.SURF_FOLDER('BiV/uvc/uvc')
        self.w_GlVTKConvert(model_name, ndata_list, output_name, trim_names=True)

        ndata_list = ['BiV/uvc/BiV.sol_apba_lap.dat', 'BiV/uvc/BiV.sol_rvendo_lap.dat', 'BiV/uvc/BiV.sol_endoepi_lap.dat', 'BiV/uvc/BiV.sol_lvendo_lap.dat']
        ndata_list = [self.SURF_FOLDER(n) for n in ndata_list]
        output_name = self.SURF_FOLDER('BiV/uvc/laplace')
        self.w_GlVTKConvert(model_name, ndata_list, output_name, trim_names=True)

        milog.info(" ## Calculating UVCs for the la mesh ##")
        la_etags = self.SURF_FOLDER_LA('la/etags.sh')
        model_name = self.SURF_FOLDER_LA('la/la')
        odir = self.SURF_FOLDER_LA('la/uvc/')

        fu.mycp(os.path.join(etags, 'etags_la.sh'), la_etags)
        self.w_mguvc(model_name, 'lv', 'lv', NP, la_etags, odir, laplace_solution=True, id='$MESH/UVC_ek')

        ndata_list = ['la/uvc/la.uvc_phi.dat', 'la/uvc/la.uvc_z.dat', 'la/uvc/la.uvc_ven.dat', 'la/uvc/la.uvc_rho.dat']
        ndata_list = [self.SURF_FOLDER_LA(n) for n in ndata_list]
        output_name = self.SURF_FOLDER_LA('la/uvc/uvc')
        self.w_GlVTKConvert(model_name, ndata_list, output_name, trim_names=True)

        milog.info(" ## Calculating UVCs for the ra mesh ##")
        ra_etags = self.SURF_FOLDER_RA('ra/etags.sh')
        model_name = self.SURF_FOLDER_RA('ra/ra')
        odir = self.SURF_FOLDER_RA('ra/uvc/')

        fu.mycp(os.path.join(etags, 'etags_ra.sh'), ra_etags)
        self.w_mguvc(model_name, 'lv', 'lv', NP, ra_etags, output_name, laplace_solution=True, id='$MESH/UVC_ek')

        ndata_list = ['ra/uvc/ra.uvc_phi.dat', 'ra/uvc/ra.uvc_z.dat', 'ra/uvc/ra.uvc_ven.dat', 'ra/uvc/ra.uvc_rho.dat']
        ndata_list = [self.SURF_FOLDER_RA(n) for n in ndata_list]
        output_name = self.SURF_FOLDER_RA('ra/uvc/uvc')

        self.w_GlVTKConvert(model_name, ndata_list, output_name, trim_names=True)

    def main_fibres_process(self, msh_path, fch_path, fch_name, angles_dict=FIBRE_ANGLES) :
        milog.info("## Calculating fibres on the BiV mesh ##" ) 
        MSHPATH = lambda x : os.path.join(msh_path, x)
        self.set_fibres_angles(angles_dict)
        
        meshname = MSHPATH('BiV')
        uvc_apba = MSHPATH('uvc/BiV.sol_apba_lap.dat')
        uvc_endoepi = MSHPATH('uvc/BiV.sol_endoepi_lap.dat')
        uvc_lv = MSHPATH('uvc/BiV.sol_lvendo_lap.dat')
        uvc_rv = MSHPATH('uvc/BiV.sol_rvendo_lap.dat')
        output = MSHPATH(self.get_meshname_fibres('fibres_bayer'))

        self.w_GlRuleFibres(meshname, uvc_apba, uvc_endoepi, uvc_lv, uvc_rv, output, angles_dict, type='biv')

        milog.info("## Substituted biv.lon with new fibres ##")
        fu.mycp(output, MSHPATH('BiV.lon'))

        milog.info("## Generating elements centres ##")
        output = MSHPATH('BiV_elem_centres.pts')
        self.w_GlElemCenters(meshname, output)

        milog.info("## Correcting fibre orientation ##")
        correct_fibres(meshname)

        milog.info("## Substituted biv.lon with new fibres ##")
        fu.mycp(MSHPATH('BiV_corrected.lon'), MSHPATH('BiV.lon'))

        milog.info("## Converting for visualisation with Paraview ##")
        mu.convert_to_vtk(meshname, f'{MSHPATH("BiV_fibres")}')

        milog.info('## Copying myocardium mesh into surfaces folder ##') 
        original_mesh_files = self.HEART('meshing/myocardium_OUT/myocardium.*') 
        fu.mycp(original_mesh_files, fch_path) 

        milog.info('## Inserting in submesh ##')
        gen_fibres_mesh = os.path.join(fch_path, fch_name)
        self.set_fch_meshname(f'{fch_name}_bayer')
        inserted_mesh = os.path.join(fch_path, self.get_meshname_fibres(self.fch_meshname))
        mu.generate_fibres(gen_fibres_mesh, 2, gen_fibres_mesh)
        mu.insert_submesh(gen_fibres_mesh, meshname, inserted_mesh, 'carp_txt')
        mu.convert(inserted_mesh, 'carp_txt', inserted_mesh, 'vtk_bin')

        milog.info('## Done ##')

    def main_mesh_process(self, meshname, surface, input_tags_file, raa_apex, outdir=SDIR['uac']) :
        meshname = self.HEART(meshname)
        outdir = self.HEART(outdir)
        OUTD = lambda x : os.path.join(outdir, x)

        BI_FOLDER = OUTD('biatrial')
        LA_FOLDER = OUTD('la')
        RA_FOLDER = OUTD('ra')

        fu.mymkdir(BI_FOLDER)
        fu.mymkdir(LA_FOLDER)
        fu.mymkdir(RA_FOLDER)

        input_tags = fu.load_json(input_tags_file)

        mu.meshtool_extract_biatrial(meshname, outdir, input_tags)
        mu.meshtool_extract_LA(outdir, input_tags)
        mu.meshtool_extract_RA(outdir, input_tags)
        mu.meshtool_extract_surfaces(meshname, outdir, input_tags, export_sup_inf=True, rm_vv_from_aa=True, surface=surface)

        mu.export_vtk_meshes_caroline(outdir, raa_apex_file=None, surface=surface)

        fu.mymkdir(OUTD('Landmarks'))
        fu.mymkdir(OUTD('Landmarks/LA'))
        fu.mymkdir(OUTD('Landmarks/RA'))

        cp_files = [
            (OUTD("/LA_epi/prodLaLandmarks.txt"), OUTD("/Landmarks/LA/prodRaLandmarks.txt")),
            (OUTD("/LA_epi/prodLaRegion.txt"), OUTD("/Landmarks/LA/prodRaRegion.txt")),
            (OUTD("/RA_epi/prodRaLandmarks.txt"), OUTD("/Landmarks/RA/")),
            (OUTD("/RA_epi/prodRaRegion.txt"), OUTD("/Landmarks/RA/"))
        ]

        for src, dst in cp_files:
            fu.mycp(src, dst)

        fu.myrm(os.path.join(outdir, 'tmp'))

        mu.recompute_raa_base(outdir, raa_apex, OUTD('RA_epi/prodRaRegion.txt'), scale=0.001, surface=surface)
        mu.scale_landmarks(OUTD('RA_epi/prodRaLandmarks.txt'), scale=0.001)
        
        mu.scale_landmarks(OUTD('LA_epi/prodLaLandmarks.txt'), scale=0.001)
        mu.scale_landmarks(OUTD('LA_epi/prodLaRegion.txt'), scale=0.001)

        cp_files = [
            (OUTD("/LA_epi/prodLaLandmarks.txt"), OUTD("/LA_endo/prodLaLandmarks.txt")),
            (OUTD("/LA_epi/prodLaRegion.txt"), OUTD("/LA_endo/prodLaRegion.txt")),	
            (OUTD("/LA_epi/prodLaLandmarks.txt"), OUTD("/Landmarks/LA/prodRaLandmarks.txt")),
            (OUTD("/LA_epi/prodLaRegion.txt"), OUTD("/Landmarks/LA/prodRaRegion.txt")),	
            (OUTD("/RA_epi/prodRaLandmarks.txt"), OUTD("/RA_endo/prodRaLandmarks.txt")),
            (OUTD("/RA_epi/prodRaRegion.txt"), OUTD("/RA_endo/prodRaRegion.txt")),
            (OUTD("/RA_epi/prodRaLandmarks.txt"), OUTD("/Landmarks/RA/")),
            (OUTD("/RA_epi/prodRaRegion.txt"), OUTD("/Landmarks/RA/"))
        ]

        for src, dst in cp_files:
            fu.mycp(src, dst)

        fu.myrm(os.path.join(outdir, 'tmp'))

        convert_list = [
            (OUTD('LA_endo/LA_endo.vtk'), OUTD('LA_endo/LA_only')),
            (OUTD('LA_epi/LA_epi.vtk'), OUTD('LA_epi/LA_only')),
            (OUTD('RA_endo/RA_endo.vtk'), OUTD('RA_endo/RA_only')),
            (OUTD('RA_epi/RA_epi.vtk'), OUTD('RA_epi/RA_only'))
        ]

        for msh_vtk, msh_carp in convert_list:
            mu.convert(msh_vtk, '', msh_carp, 'carp_txt')

        milog.info('## Done ##')

    ATRIUM_CHOICES = ['la', 'ra']
    def main_laplace_process(self, atrium, parfiles_path) :
        def helper_surf_to_vtx(input_surf, output_vtx):
            _surf = fu.read_elem(input_surf, el_type='Tr', tags=False)
            _vtx = fu.surf2vtx(_surf)
            fu.write_vtx(output_vtx, _vtx, init_row=2)

        atrium = atrium.lower()
        if atrium not in self.ATRIUM_CHOICES:
            msg = f"Invalid atrium choice: {atrium}"
            milog.error(msg)
            raise ValueError(msg)
        
        fch_mesh = os.path.join(SDIR['uvc'], 'myocardium_bayer_60_-60') 
        UACFOLDER = self.HEART(SDIR['uac'])
        FCH = self.HEART(fch_mesh)
        
        meshname = os.path.join(UACFOLDER, f'{atrium}/{atrium}')
        endo = os.path.join(UACFOLDER, f'{atrium}/{atrium}_endo.surf')
        epi = os.path.join(UACFOLDER, f'{atrium}/{atrium}_epi.surf')
        outdir = os.path.join(UACFOLDER, f'{atrium}/endo_epi')

        milog.info(f'## Runing endo-epi Laplace on {meshname}')
        milog.info('## Converting surfaces to vtx files ##')

        helper_surf_to_vtx(endo, f'{endo}.vtx')
        helper_surf_to_vtx(epi, f'{epi}.vtx')

        milog.info('## Running Laplace solution endo to epi ##')
        self.w_carp_pt(parfiles_path, simID=outdir, meshname=meshname, stim0=endo, stim1=epi)

        milog.info('## Converting igb to .dat file ##')
        self.w_igbextract(f'{outdir}/phie.igb', f'{outdir}/phie.dat')

        milog.info('## Converting to VTK to visualize ##')
        self.w_GlVTKConvert(meshname=meshname, nodedata_list=[f'{outdir}/phie.dat'], output=f'{outdir}/laplace_endo_epi', trim_names=True)

        milog.info('## Done ##')

    def main_surf_to_volume_process(self, atrium) :
        atrium = atrium.lower()
        if atrium not in self.ATRIUM_CHOICES:
            msg = f"Invalid atrium choice: {atrium}"
            milog.error(msg)
            raise ValueError(msg)
        
        uac_folder = self.HEART(SDIR['uac'])
        
        meshname = os.path.join(uac_folder, f'{atrium}/{atrium}')
        meshname_uac = os.path.join(uac_folder, f'{atrium.upper()}_endo/Fibre_endo_l')
        endo_fibres = os.path.join(uac_folder, f'{atrium.upper()}_endo/Fibre_endo_l.lon')
        epi_fibres = os.path.join(uac_folder, f'{atrium.upper()}_endo/Fibre_epi_l.lon')
        endo_epi_laplace = os.path.join(uac_folder, atrium, 'endo_epi/phie.dat') 
        outmeshname = os.path.join(uac_folder, atrium, f'{atrium}_fibres_l')

        milog.info('## Mapping fibres from 2D to 3D ##')
        mmu.laplace_endo2elem(meshname, endo_epi_laplace)

        milog.info('## Computing element centers on both meshes 3D and 2D meshes ##')
        mmu.compute_elemCenters(meshname, el_type='Tt')
        mmu.compute_elemCenters(meshname_uac, el_type='Tr')

        milog.info('## Mapping fibres ##')
        endo_epi_laplace_el = f'{endo_epi_laplace[:-4]}_el.dat'
        meshname_elem_centers = f'{meshname}_elemCenters.pts'
        meshname_uac_elem_centers = meshname_uac + '_elemCenters.pts'
        outmsh_ext = f'{outmeshname}.lon'
        mmu.map_fibres_3d(endo_epi_laplace_el, meshname_elem_centers, meshname_uac_elem_centers, endo_fibres, epi_fibres, outmsh_ext, f'{meshname}_endo_to_3d.map')

        milog.info('## Finding sheet direction ##')
        msh_transmural = f'{meshname}_transmural.lon'
        msh_rot_axes = f'{meshname}_rotation_axes.lon'
        msh_sheet = f'{outmeshname}_sheet'
        mmu.find_transmural_direction(meshname, meshname_uac, meshname_elem_centers, meshname_uac_elem_centers, msh_transmural)
        mmu.find_rotation_axes(outmsh_ext, msh_transmural, msh_rot_axes)
        mmu.make_sheet_orthogonal(outmsh_ext, msh_transmural, msh_rot_axes, f'{msh_sheet}.lon')

        milog.info('## Converting to VTK for visualization ##')
        fu.mycp(f'{meshname}.elem', f'{msh_sheet}.elem')
        fu.mycp(f'{meshname}.pts', f'{msh_sheet}.pts')
        mu.convert_to_vtk(f'{msh_sheet}', f'{msh_sheet}')

        milog.info('## Done ##')

    def main_mapping_to_biatrial(self, atrium) : 
        atrium = atrium.lower()
        if atrium not in self.ATRIUM_CHOICES:
            msg = f"Invalid atrium choice: {atrium}"
            milog.error(msg)
            raise ValueError(msg)
        
        uac_folder = self.HEART(SDIR['uac'])
        uacp = lambda x : os.path.join(uac_folder, x)

        cp_files = [
            (uacp(f'{atrium}/{atrium}.nod'), uacp(f'{atrium}/{atrium}_fibres_l_sheet.nod')), 
            (uacp(f'{atrium}/{atrium}.eidx'), uacp(f'{atrium}/{atrium}_fibres_l_sheet.eidx')), 
            
        ]

        ## start here
        