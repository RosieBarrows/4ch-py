import numpy as np
import os
import trimesh
import sys
import copy
import meshio
from tqdm import tqdm
import time
import json
from scipy import spatial

from common_4ch.motion_volume import *
from common_4ch.file_utils import *
from common_4ch.linalg_utils import *

def compute_elemCenters(meshname,
						el_type):

	pts = read_pts(meshname+".pts") 
	elem = read_elem(meshname+".elem",el_type=el_type,tags=False) 
	n_tet = elem.shape[0] 

	cog_el = np.zeros((n_tet,3),dtype=float) 
	for i,el in enumerate(elem):
		pts_el = pts[el,:]
		cog_el[i,:] = np.mean(pts_el,axis=0)
	
	write_pts_caroline(cog_el,meshname+"_elemCenters.pts")

def laplace_endo2elem(meshname,
					  phi_nodes,
					  output_phi=None):

	if output_phi is None:
		output_phi = phi_nodes[:-4]+"_el.dat"

	elem = read_elem(meshname+".elem",el_type='Tt',tags=False)
	n_tet = elem.shape[0]	

	phi = np.loadtxt(phi_nodes,dtype=float)

	phi_el = np.empty((n_tet,),dtype=float)
	for i,tt in enumerate(elem):
		flag=0
		for k,node in enumerate(tt):
			if phi[tt[k]] == 1:
				phi_el[i] = 1
				flag=1
			if phi[tt[k]] == 0:
				phi_el[i] = 0
				flag=1
		if flag==0:
			phi_el[i]=np.mean(phi[tt])	

	np.savetxt(output_phi,phi_el,fmt="%g")

def tag_node2elem(meshname,
				  tag_file,
				  output_file):

	print("Mapping "+tag_file+" from nodes to elements...")

	elem = read_elem(meshname+".elem",el_type='Tr',tags=False)
	tags = np.loadtxt(tag_file,dtype=int)

	tags_el = np.zeros((elem.shape[0],),dtype=int)

	for i,tr in enumerate(elem):

		all_in = np.sum((tags[tr]!=0))
		
		if all_in==3:
			tags_el[i] = 1

	np.savetxt(output_file,tags_el,fmt="%d")

def map_fibres_3d(phi_el_file,
				  elemCenters_3d,
				  elemCenters_2d,
				  f_endo_file,
				  f_epi_file,
				  output_fibres,
				  map_tr_tt_file=None):
		
	phi = np.loadtxt(phi_el_file,dtype=float)
	cog_tt = read_pts(elemCenters_3d)
	cog_tr = read_pts(elemCenters_2d)

	f_endo = np.loadtxt(f_endo_file,dtype=float) 
	f_epi = np.loadtxt(f_epi_file,dtype=float) 

	save_mapping = False
	if map_tr_tt_file is not None:
		if os.path.exists(map_tr_tt_file):
			mapping_tr_tt = np.loadtxt(map_tr_tt_file,dtype=int)
		else:
			mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)
			save_mapping = True
	else:
		mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)

	if phi.shape[0] != cog_tt.shape[0]:
		raise Exception("The laplace solution is not defined on the elements.")

	f_tt = np.zeros((phi.shape[0],3),dtype=float)
	for i in range(phi.shape[0]):

		if (map_tr_tt_file is None) or save_mapping:
			distance = np.linalg.norm(cog_tr-cog_tt[i,:],axis=1)
			idx_closest_tr = np.argmin(distance)
			mapping_tr_tt[i]= idx_closest_tr
		else:
			idx_closest_tr = mapping_tr_tt[i]

		if phi[i] > 0.5:
			f_tt[i,:] = f_epi[idx_closest_tr,:]
		else:
			f_tt[i,:] = f_endo[idx_closest_tr,:]

	if save_mapping:
		np.savetxt(map_tr_tt_file,mapping_tr_tt,fmt="%d")

	write_lon(f_tt,output_fibres)	

def map_tags_3d(phi_el_file,
				elemCenters_3d,
				elemCenters_2d,
				tags_2d_file,
				output_tags,
				endo_tags=[],
				map_tr_tt_file=None):
		
	phi = np.loadtxt(phi_el_file,dtype=float)
	cog_tt = read_pts(elemCenters_3d)
	cog_tr = read_pts(elemCenters_2d)

	tags_2d = np.loadtxt(tags_2d_file,dtype=int)

	save_mapping = False
	if map_tr_tt_file is not None:
		if os.path.exists(map_tr_tt_file):
			mapping_tr_tt = np.loadtxt(map_tr_tt_file,dtype=int)
		else:
			mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)
			save_mapping = True
	else:
		mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)

	if phi.shape[0] != cog_tt.shape[0]:
		raise Exception("The laplace solution is not defined on the elements.")

	tags_tt = np.zeros((phi.shape[0],),dtype=float)
	for i in range(phi.shape[0]):

		if (map_tr_tt_file is None) or save_mapping:
			distance = np.linalg.norm(cog_tr-cog_tt[i,:],axis=1)
			idx_closest_tr = np.argmin(distance)
			mapping_tr_tt[i]= idx_closest_tr
		else:
			idx_closest_tr = mapping_tr_tt[i]

		if tags_2d[idx_closest_tr] in endo_tags:
			if phi[i]==0:
				tags_tt[i] = tags_2d[idx_closest_tr]
		else:
			tags_tt[i] = tags_2d[idx_closest_tr]

	if save_mapping:
		np.savetxt(map_tr_tt_file,mapping_tr_tt,fmt="%d")

	np.savetxt(output_tags,tags_tt)	

def find_transmural_direction_Gl(meshname,
							  phi_el_file,
							  output_file):

	cmd = "GlGradient -m "+meshname
	cmd += " -S elem_ctr -t elem_ctr "
	cmd += "-d "+phi_el_file+" -o "+output_file

	os.system(cmd)


def find_transmural_direction(meshname_3d,
							  meshname_2d,
							  elemCenters_3d,
				  			  elemCenters_2d,
				  			  output_file,
							  map_tr_tt_file=None):
	
	elem_3d = read_elem(meshname_3d+".elem",el_type="Tt",tags=False)

	cog_tt = read_pts(elemCenters_3d)
	cog_tr = read_pts(elemCenters_2d)

	normals_surf = compute_normal_from_surface(meshname_2d)
	normals_tets = np.zeros((elem_3d.shape[0],3),dtype=float)

	save_mapping = False
	if map_tr_tt_file is not None:
		if os.path.exists(map_tr_tt_file):
			mapping_tr_tt = np.loadtxt(map_tr_tt_file,dtype=int)
		else:
			mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)
			save_mapping = True
	else:
		mapping_tr_tt = np.zeros((cog_tt.shape[0],),dtype=int)

	for i in range(elem_3d.shape[0]):

		if (map_tr_tt_file is None) or save_mapping:
			distance = np.linalg.norm(cog_tr-cog_tt[i,:],axis=1)
			idx_closest_tr = np.argmin(distance)
			mapping_tr_tt[i]= idx_closest_tr
		else:
			idx_closest_tr = mapping_tr_tt[i]

		normals_tets[i,:] = normals_surf[idx_closest_tr,:]

	if save_mapping:
		np.savetxt(map_tr_tt_file,mapping_tr_tt,fmt="%d")

	write_lon(normals_tets,output_file)

def tag_connected_component_2d(meshname,
							   tag_file,
							   tag_outfile,
							   vtx_threshold=10):

	print("Removing components with less than "+str(vtx_threshold)+" point from "+tag_file+"...")

	pts = read_pts(meshname+".pts")
	elem = read_elem(meshname+".elem",el_type="Tr",tags=False)
	tags = np.loadtxt(tag_file,dtype=int)

	elem_tag = elem[np.where(tags!=0)[0],:]
	tag_vtx = surf2vtx(elem_tag)
	pts_tag = pts[tag_vtx,:]
	elem_tag_reindexed = reindex_surf(tag_vtx,elem_tag)

	mesh_tag = trimesh.Trimesh(vertices=pts_tag,
							   faces=elem_tag_reindexed,
							   process=False)

	cc_labels_tag = trimesh.graph.connected_component_labels(mesh_tag.edges)
	n_cc = np.max(cc_labels_tag)+1

	if n_cc > 1:
		cc_size = np.array([(cc_labels_tag == i).sum() for i in range(n_cc)])

		remove_cc_idx = np.where(cc_size<=vtx_threshold)[0]

		for cc in remove_cc_idx:
			vtx_cc = tag_vtx[np.where(cc_labels_tag==cc)[0]]
			for i,tr in enumerate(elem):
				inters = np.intersect1d(vtx_cc,tr)
				if inters.size:
					tags[i] = 0

	np.savetxt(tag_outfile,tags,fmt="%d")

def combine_tags(meshname,
				 tagfile_list,
				 tags,
				 default_tag,
				 output_file):

	if len(tagfile_list)!=len(tags):
		raise Exception("Tags and tag files do not match in size.")

	elem = read_elem(meshname+".elem",el_type="Tr",tags=False)

	tags_combined = np.zeros((elem.shape[0],),dtype=int)+default_tag

	for i,tag_file in enumerate(tagfile_list):
		tmp = np.loadtxt(tag_file,dtype=int)
		tags_combined[np.where(tmp==1)[0]] = tags[i]

	np.savetxt(output_file,tags_combined,fmt="%d")

def surf2vtk(msh_file,surf_file,output_name):
    print('Converting '+surf_file+' to vtk...' )

    surf = np.loadtxt(surf_file,dtype=int,skiprows=1,usecols=[1,2,3])
    vtx = surf2vtx(surf)

    # msh = meshio.read(msh_file)
    pts = np.loadtxt(msh_file+".pts",dtype=float,skiprows=1)

    pts_surf = pts[vtx,:]

    surf_reindexed = copy.deepcopy(surf)



    for i in range(surf.shape[0]):
            surf_reindexed[i,0] = np.where(vtx==surf[i,0])[0]
            surf_reindexed[i,1] = np.where(vtx==surf[i,1])[0]
            surf_reindexed[i,2] = np.where(vtx==surf[i,2])[0]

    cells = {"triangle": surf_reindexed}
    surf_vtk_msh = meshio.Mesh(
                            pts_surf,
                            cells
                            )

    meshio.write(output_name, surf_vtk_msh,file_format="vtk")

def surf2vtx(surf):
	return np.unique(surf.flatten())


def write_surf(surf,surf_file):
        f = open(surf_file,"w")
        f.write(str(int(surf.shape[0]))+"\n")
        for t in surf:
                f.write("Tr "+str(t[0])+" "+str(t[1])+" "+str(t[2])+"\n")
        f.close()

def write_vtx(vtx,vtx_file):
        f = open(vtx_file,"w")
        f.write(str(int(vtx.shape[0]))+"\n")
        f.write("intra\n")
        for v in vtx:
                f.write(str(v)+"\n")
        f.close()

def connected_component_to_surface(eidx_file,original_surface,output_surface):
	eidx = np.fromfile(eidx_file+".eidx", dtype=int, count=-1)
	nod = np.fromfile(eidx_file+".nod", dtype=int, count=-1)
	surf = np.loadtxt(original_surface,dtype=int,skiprows=1,usecols=[1,2,3])
	vtx = surf2vtx(surf)

	subsurf = surf[eidx,:]
	subvtx = vtx[nod]

	write_surf(subsurf,output_surface+".surf")
	write_vtx(subvtx,output_surface+".vtx")

def remove_sept(mesh,surf_folder):
	print("Converting rvendo.surf and rvsept.surf to rvendo.vtx and rvsept.vtx")
	rvendo_surf = np.loadtxt(surf_folder+"/tmp/myocardium.rvendo.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rvsept_surf = np.loadtxt(surf_folder+"/tmp/myocardium.rvsept.surf",dtype=int,skiprows=1,usecols=[1,2,3])

	rvendo_vtx = surf2vtx(rvendo_surf)
	rvsept_vtx = surf2vtx(rvsept_surf)

	print("Finding the free wall")
	freewall_vtx = np.setdiff1d(rvendo_vtx,rvsept_vtx)

	freewall_tr = []
	sept_tr = []

	for tr in tqdm (rvendo_surf, desc="Looping over triangles..."):
		if len(np.intersect1d(tr,rvsept_vtx))<3:
			freewall_tr.append(tr)
        # else:
        	# sept_tr.append(tr)

	print("Finished looping over triangles")

	freewall_tr = np.array(freewall_tr)

	write_surf(freewall_tr,surf_folder+"/tmp/myocardium.rvendo_nosept.surf")
	freewall_vtx = surf2vtx(freewall_tr)
	write_vtx(freewall_vtx,surf_folder+"/tmp/myocardium.rvendo_nosept.surf"+".vtx")

	surf2vtk(mesh,surf_folder+"/tmp/myocardium.rvendo_nosept.surf",surf_folder+"/tmp/myocardium.rvendo_nosept.surf"+".vtk")

def prepare_vtx_for_uvc(surf_folder):
	epi_surf = np.loadtxt(surf_folder+"/tmp/myocardium.epi.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	lvendo_surf = np.loadtxt(surf_folder+"/tmp/myocardium.lvendo.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rvendo_surf = np.loadtxt(surf_folder+"/tmp/myocardium.rvendo.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rvendo_nosept_surf = np.loadtxt(surf_folder+"/tmp/myocardium.rvendo_nosept.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rvsept_surf = np.loadtxt(surf_folder+"/tmp/myocardium.rvsept.surf",dtype=int,skiprows=1,usecols=[1,2,3])

	epi_surf_vtx = surf2vtx(epi_surf)
	lvendo_surf_vtx = surf2vtx(lvendo_surf)
	rvendo_surf_vtx = surf2vtx(rvendo_surf)
	rvendo_nosept_surf_vtx = surf2vtx(rvendo_nosept_surf)
	rvsept_surf_vtx = surf2vtx(rvsept_surf)

	write_vtx(epi_surf,surf_folder+"/tmp/myocardium.epi.surf.vtx")
	write_vtx(lvendo_surf,surf_folder+"/tmp/myocardium.lvendo.surf.vtx")
	write_vtx(rvendo_surf,surf_folder+"/tmp/myocardium.rvendo.surf.vtx")
	write_vtx(rvendo_nosept_surf,surf_folder+"/tmp/myocardium.rvendo_nosept.surf.vtx")
	write_vtx(rvsept_surf,surf_folder+"/tmp/myocardium.rvsept.surf.vtx")

def extract_tags(input_tags,
				 labels):
	
	tags_list = []
	for l in labels:
		if type(input_tags[l])==int:
			tags_list += [input_tags[l]]
		elif type(input_tags[l])==np.ndarray:
			tags_list += list(input_tags[l])
		else:
			raise Exception("Type not recognised.")

	return tags_list

def define_AV_separation(elem_file,
			  			 tags,
			  			 tag_AVplane,
			  			 elem_output):

	print('Reading mesh elem file...')
	elem = read_tets(elem_file)
	print('Done.')

	print('Intersecting ventricles and atria...')

	tags_list_vent_fec = extract_tags(tags,["LV","RV","FEC"])
	tags_list_atr_bb = extract_tags(tags,["LA","RA","BB"])

	VEIDX = np.where(np.isin(elem[:,-1],tags_list_vent_fec)==1)[0]
	VVTX = np.unique(elem[VEIDX,0:4].flatten())
	AEIDX = np.where(np.isin(elem[:,-1],tags_list_atr_bb)==1)[0]
	AVTX = np.unique(elem[AEIDX,0:4].flatten())
	AVPLANEVTX = np.intersect1d(VVTX,AVTX)
	print('Done.')

	v1 = np.isin(elem[:,0],AVPLANEVTX)
	v2 = np.isin(elem[:,1],AVPLANEVTX)
	v3 = np.isin(elem[:,2],AVPLANEVTX)
	v4 = np.isin(elem[:,3],AVPLANEVTX)

	ind_v1 =  np.where(v1==1)[0]
	ind_v2 =  np.where(v2==1)[0]
	ind_v3 =  np.where(v3==1)[0]
	ind_v4 =  np.where(v4==1)[0]

	tags_list_vent = extract_tags(tags,["LV","RV"])

	AVPLANE_EIDX = np.unique(np.concatenate((ind_v1,ind_v2,ind_v3,ind_v4)))
	AVPLANE_EIDX = np.intersect1d(AVPLANE_EIDX,np.where(np.isin(elem[:,-1],tags_list_vent)==1)[0])

	elem[AVPLANE_EIDX,-1] = tag_AVplane
	print('Done.')

	print('Saving output...')
	write_tets(elem_output,elem)
	print('Done.')

def define_FEC(elem_file,
			   biv_meshname,
			   Zbiv_file,
			   RHObiv_file,
			   elem_output,
			   tag_fec,
			   include_septum=None,
			   FEC_height=0.7):
	if (FEC_height>1.0) or (FEC_height<0.0):
		raise Exception("Max height for the FEC should be between 0 and 1.")

	print('Reading mesh elem file...')
	elem = read_tets(elem_file)
	print('Done.')

	print('Reading coordinates files...')
	Zbiv = np.loadtxt(Zbiv_file)
	RHObiv = np.loadtxt(RHObiv_file)
	print('Done.')

	print('Intersecting endo surfaces and bottom '+str(FEC_height*100)+'%...')
	endo_VTXbiv = np.where(RHObiv==0)[0]

	if include_septum is not None:
		print("include_septum is not None, so I am trying to read RV_endo.surf mapped onto the biv mesh.")
		rvendosurf=read_surf(include_septum)
		rvendovtx=surf2vtx(rvendosurf)
		septum = np.intersect1d(rvendovtx,np.where(RHObiv==1)[0])
		endo_VTXbiv = np.concatenate((endo_VTXbiv,rvendovtx))

	bottom70_VTXbiv = np.where(Zbiv<=FEC_height)[0]
	endo_bottom_VTXbiv = np.intersect1d(endo_VTXbiv,bottom70_VTXbiv)
	print('Done.')
	
	print('Mapping back to four-chamber...')	
	VTXbiv = np.fromfile(biv_meshname+".nod", dtype=int, count=-1)

	if (Zbiv.shape[0]!=VTXbiv.shape[0]) or (RHObiv.shape[0]!=VTXbiv.shape[0]):
		raise Exception("The biv mesh does not match the UVCs you gave.")

	endo_bottom_VTX = VTXbiv[endo_bottom_VTXbiv]
	print('Done.')

	print('Mapping from vtx to elements...')
	v1 = np.isin(elem[:,0],endo_bottom_VTX)
	v2 = np.isin(elem[:,1],endo_bottom_VTX)
	v3 = np.isin(elem[:,2],endo_bottom_VTX)
	v4 = np.isin(elem[:,3],endo_bottom_VTX)

	ind_v1 =  np.where(v1==1)[0]
	ind_v2 =  np.where(v2==1)[0]
	ind_v3 =  np.where(v3==1)[0]
	ind_v4 =  np.where(v4==1)[0]

	FEC_EIDX = np.unique(np.concatenate((ind_v1,ind_v2,ind_v3,ind_v4)))
	elem[FEC_EIDX,-1] = tag_fec
	print('Done.')

	print('Saving output...')
	write_tets(elem_output,elem)
	print('Done.')

def define_BB(elem_file,
			  la_meshname,
			  ra_meshname,
			  Zla_file,
			  Zra_file,
			  PHIla_file,
			  PHIra_file,
			  settings,
			  tags,
			  elem_output):

	print('Reading mesh elem file...')
	elem = read_tets(elem_file)
	print('Done.')

	print('Reading coordinates files...')
	Zla = np.loadtxt(Zla_file)
	Zra = np.loadtxt(Zra_file)
	PHIla = np.loadtxt(PHIla_file)
	PHIra = np.loadtxt(PHIra_file)
	print('Done.')

	print('Finding Bachmann Bundle...')
	# LA
	Zla_VTX = np.intersect1d(np.where(Zla<=settings["LA"]["z_max"])[0],
							  np.where(Zla>=settings["LA"]["z_min"])[0])
	PHIla_VTX = np.intersect1d(np.where(PHIla<=settings["LA"]["phi_max"])[0],
							  np.where(PHIla>=settings["LA"]["phi_min"])[0])
	ZPHIla_VTX = np.intersect1d(Zla_VTX,PHIla_VTX)
	# RA
	Zra_VTX = np.intersect1d(np.where(Zra<=settings["RA"]["z_max"])[0],
							  np.where(Zra>=settings["RA"]["z_min"])[0])
	PHIra_VTX = np.intersect1d(np.where(PHIra<=settings["RA"]["phi_max"])[0],
							  np.where(PHIra>=settings["RA"]["phi_min"])[0])
	ZPHIra_VTX = np.intersect1d(Zra_VTX,PHIra_VTX)
	print('Done.')

	print('Mapping back to four-chamber...')	
	VTXla = np.fromfile(la_meshname+".nod", dtype=int, count=-1)
	VTXra = np.fromfile(ra_meshname+".nod", dtype=int, count=-1)
	BB_VTX = np.concatenate((VTXla[ZPHIla_VTX],VTXra[ZPHIra_VTX]))
	print('Done.')

	print('Mapping from vtx to elements...')
	v1 = np.isin(elem[:,0],BB_VTX)
	v2 = np.isin(elem[:,1],BB_VTX)
	v3 = np.isin(elem[:,2],BB_VTX)
	v4 = np.isin(elem[:,3],BB_VTX)

	ind_v1 =  np.where(v1==1)[0]
	ind_v2 =  np.where(v2==1)[0]
	ind_v3 =  np.where(v3==1)[0]
	ind_v4 =  np.where(v4==1)[0]

	tags_list_atr = extract_tags(tags,["LA","RA"])

	BB_EIDX = np.unique(np.concatenate((ind_v1,ind_v2,ind_v3,ind_v4)))
	BB_EIDX = np.intersect1d(BB_EIDX,np.where(np.isin(elem[:,-1],tags_list_atr)==1)[0])
	elem[BB_EIDX,-1] = tags["BB"]
	print('Done.')

	print('Saving output...')
	write_tets(elem_output,elem)
	print('Done.')

def separate_FEC_lvrv(elem_file,
					  fec_elem_file,
					  LV_endo_surf,
					  RV_endo_surf,
					  elem_output,
					  original_tags_settings,
					  new_tags_settings):
	
	print('Reading original tags...')
	elem = read_tets(elem_file)
	original_tags = elem[:,-1]
	print('Done')

	print('Reading fec tags...')
	elem_fec = read_tets(fec_elem_file)
	fec_tags = elem_fec[:,-1]
	print('Done')

	print('Reading surfaces...')
	lv_endo = read_surf(LV_endo_surf)
	lv_endo_vtx = surf2vtx(lv_endo)
	rv_endo = read_surf(RV_endo_surf)
	rv_endo_vtx = surf2vtx(rv_endo)
	print('Done')

	lv_eidx = np.where(original_tags == new_tags_settings["LV"])[0]		
	rv_eidx = np.where(original_tags == new_tags_settings["RV"])[0]
	fec_eidx = np.where(fec_tags == original_tags_settings["FEC"])[0]

	print('Extracted LV, FEC and RV')

	new_tags = copy.deepcopy(fec_tags)
	for eidx in fec_eidx:
		intersection_rv = np.intersect1d(elem[eidx,0:-1],rv_endo_vtx)
		if original_tags[eidx]==1 and len(intersection_rv)==0:
			new_tags[eidx] = new_tags_settings["FEC_LV"]
		elif original_tags[eidx]==1 and len(intersection_rv)>0:
			new_tags[eidx] = new_tags_settings["FEC_SV"]
		else:
			new_tags[eidx] = new_tags_settings["FEC_RV"]
			
	elem_fec[:,-1] = new_tags
	write_tets(elem_output,elem_fec)

def separate_FEC_lvrv_sets(elem_file, fec_elem_file, LV_endo_surf, RV_endo_surf, elem_output, original_tags_settings, new_tags_settings):

    print('Reading original tags...')
    elem = read_tets(elem_file)
    original_tags = elem[:, -1]
    print('Done')

    print('Reading fec tags...')
    elem_fec = read_tets(fec_elem_file)
    fec_tags = elem_fec[:, -1]
    print('Done')

    print('Reading surfaces...')
    lv_endo = read_surf(LV_endo_surf)
    lv_endo_vtx = set(surf2vtx(lv_endo))
    rv_endo = read_surf(RV_endo_surf)
    rv_endo_vtx = set(surf2vtx(rv_endo))
    print('Done')

    lv_eidx = np.where(original_tags == new_tags_settings["LV"])[0]
    rv_eidx = np.where(original_tags == new_tags_settings["RV"])[0]
    fec_eidx = np.where(fec_tags == original_tags_settings["FEC"])[0]

    print('Extracted LV, FEC, and RV')

    new_tags = copy.deepcopy(fec_tags)

    for eidx in fec_eidx:
        intersection_rv = set(elem[eidx, 0:-1]).intersection(rv_endo_vtx)

        if original_tags[eidx] == 1 and not intersection_rv:
            new_tags[eidx] = new_tags_settings["FEC_LV"]
        elif original_tags[eidx] == 1 and intersection_rv:
            new_tags[eidx] = new_tags_settings["FEC_SV"]
        else:
            new_tags[eidx] = new_tags_settings["FEC_RV"]

    elem_fec[:, -1] = new_tags
    write_tets(elem_output, elem_fec)

def find_cog_surf(pts):
	cog = np.mean(pts,axis=0)

	return cog

def find_cog_vol(mesh_pts,mesh_elem,tag):
	pts_idx = []
	for i,e in enumerate(mesh_elem):
		if e[4] == int(tag):	
			pts_idx.append(int(e[1]))
			pts_idx.append(int(e[2]))
			pts_idx.append(int(e[3]))

	pts_idx = np.unique(pts_idx)
	pts = np.zeros((pts_idx.shape[0],3))

	for i,p in enumerate(pts_idx):
		pts[i] = mesh_pts[p]

	cog = np.mean(pts,axis=0)

	return cog

def calculate_dist(cog_surf,cog_vol):
	dist = np.zeros((1,3))

	for i,c in enumerate(cog_vol):
		dist[:,i] = c - cog_surf[i]

	dist = np.linalg.norm(dist)
	return dist

def query_outwards_surf(surf,pts,cog):
	is_outward = np.zeros((surf.shape[0],),dtype=int)
	for i,t in enumerate(surf):
		v0 = pts[t[1],:] - pts[t[0],:]
		v0 = v0/np.linalg.norm(v0)

		v1 = pts[t[2],:] - pts[t[0],:]
		v1 = v1/np.linalg.norm(v1)

		n = np.cross(v0,v1)
		n = n/np.linalg.norm(n)

		dot_prod = np.dot(cog-pts[t[0],:],n)

		if dot_prod<0:
			is_outward[i] = 1

	if np.sum(is_outward)/surf.shape[0]>0.7:
		return True
	else:
		return False

def calculate_dist(cog_surf,cog_vol):
	dist = np.zeros((1,3))

	for i,c in enumerate(cog_vol):
		dist[:,i] = c - cog_surf[i]

	dist = np.linalg.norm(dist)
	return dist

def query_outwards_surf(surf,pts,cog):
	is_outward = np.zeros((surf.shape[0],),dtype=int)
	for i,t in enumerate(surf):
		v0 = pts[t[1],:] - pts[t[0],:]
		v0 = v0/np.linalg.norm(v0)

		v1 = pts[t[2],:] - pts[t[0],:]
		v1 = v1/np.linalg.norm(v1)

		n = np.cross(v0,v1)
		n = n/np.linalg.norm(n)

		dot_prod = np.dot(cog-pts[t[0],:],n)

		if dot_prod<0:
			is_outward[i] = 1

	if np.sum(is_outward)/surf.shape[0]>0.7:
		return True
	else:
		return False
		
def read_pnts(filename):
    return np.loadtxt(filename, dtype=float, skiprows=1)

def write_pnts(filename, pts):
    assert len(pts.shape) == 2 and pts.shape[1] == 3
    with open(filename, 'w') as fp:
        fp.write('{}\n'.format(pts.shape[0]))
        for pnt in pts:
            fp.write('{0[0]}\t{0[1]}\t{0[2]}\n'.format(pnt))


def read_surface(filename):
    return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3))

def write_surface(filename, surfs):
    assert len(surfs.shape) == 2 and surfs.shape[1] == 3
    with open(filename, 'w') as fp:
        fp.write('{}\n'.format(surfs.shape[0]))
        for tri in surfs:
            fp.write('Tr {0[0]}\t{0[1]}\t{0[2]}\n'.format(tri))

def read_neubc(filename):
   return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(0,1,2,3,4,5)) 


def read_elems(filename):
  return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3,4,5))


def vector_cprod(vec1, vec2):
  return np.array([vec1[1]*vec2[2]-vec1[2]*vec2[1],
                   vec1[2]*vec2[0]-vec1[0]*vec2[2],
                   vec1[0]*vec2[1]-vec1[1]*vec2[0]])

def read_dat(filename):
    return np.loadtxt(filename, dtype=float, skiprows=0)

def read_vtx(filename):
    return np.loadtxt(filename, dtype=int, skiprows=2)

def vector_sprod(vec1, vec2):
  return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]


def create_csys(vec):
	vec0 = None
	vec1 = None

	if (vec[0] < 0.5) and (vec[1] < 0.5):
		tmp = math.sqrt(vec[1]*vec[1]+vec[2]*vec[2])
		vec1 = np.array([0.0, -vec[2]/tmp, vec[1]/tmp])
		vec0 = vector_cprod(vec, vec1)
	else:
		tmp = math.sqrt(vec[0]*vec[0]+vec[1]*vec[1])
		vec1 = np.array([vec[1]/tmp, -vec[0]/tmp, 0.0])
		vec0 = vector_cprod(vec, vec1)

	return [vec0, vec1, vec]

def set_pericardium(mesh,presimFolder,heartFolder):
	outfile="elem_dat_UVC_ek.dat"
	basename = presimFolder+"/myocardium_AV_FEC_BB"

	pnts = read_pnts(basename+'.pts')    
	elem = read_elems(basename+'.elem')
	print("Read mesh")

	datdir = heartFolder+"/surfaces_uvc/BiV/uvc"
	datname = os.path.join(datdir, "BiV.uvc_z")
	UVCpts = read_dat(datname+".dat")
	print('Found z coords')


	vtxdir = heartFolder+'/surfaces_uvc/BiV'
	vtxname = os.path.join(vtxdir, 'BiV')
	vtxsubmsh = np.fromfile(vtxname+".nod", dtype=int, count=-1)
	print('Found BiV mesh')

	pcatags = [1, 2]
	UVCptsmsh = np.zeros(len(pnts[:,0]))
	for i, ind in enumerate(vtxsubmsh):
		UVCptsmsh[ind] = UVCpts[i]

	UVCelem = []
	for i, elm in enumerate(elem):
		if ( elem[i,4] in pcatags ):
			UVCelem.append((UVCptsmsh[elm[0]]+UVCptsmsh[elm[1]]+UVCptsmsh[elm[2]]+UVCptsmsh[elm[3]])*0.25)
		else:
			UVCelem.append(0.0)

	np.savetxt(os.path.join(heartFolder+"/pre_simulation/", "UVC_elem.dat"), UVCelem, fmt='%.8f')

	print("Computing the data on the elements")
	p1 = 1.5266
	p2 = -0.37
	p3 = 0.4964
	p4 = 0
	th = 0.82

	elemdat = []
	for i, l in enumerate(UVCelem):
		if ( elem[i,4] in pcatags ):
			if (UVCelem[i] >= th):
				elemdat.append(0.0) 
			else: 
				x = UVCelem[i]
				x_m = th-x
				elemdat.append(p1*x_m**3 + p2*x_m**2 + p3*x_m + p4)
		else:
			elemdat.append(0.0)    

	np.savetxt(presimFolder+"/"+outfile, elemdat, fmt='%.8f') 
	print("Saved pericardium scale file")

	os.system("GlVTKConvert -m "+basename+" -e "+presimFolder+"/"+outfile+" -e "+presimFolder+"/UVC_elem.dat -F bin -o "+basename+"_elem_dat_UVC --trim-names")

def combine_elem_file(input_elemFiles, output_file, mode='max'):
	
	print('Reading files...')
	elemData = np.loadtxt(input_elemFiles[0],dtype=float)
	for i in range(1,len(input_elemFiles)):
		tmp = np.loadtxt(input_elemFiles[i],dtype=float)
		elemData = np.column_stack((elemData,tmp))
	print('Done.')

	if mode=='max':
		elemData_new = np.max(elemData,axis=1)
	else:
		raise Exception('Mode not recognised')

	print('Writing output...')
	np.savetxt(output_file,elemData_new)

def combine_rot_coords(presimFolder):
	combine_elem_file([presimFolder+"/elem_dat_UVC_ek.dat", 
				  presimFolder+"/elem_dat_UVC_ek_inc_la.dat",
				  presimFolder+"/elem_dat_UVC_ek_inc_ra.dat"], 
				  presimFolder+"/elem_dat_UVC_ek_combined.dat")