import os
import sys
import numpy as np
import copy
import meshio

from common_4ch.file_utils import *
from common_4ch.distance_utils import *
from common_4ch.mesh_utils import *

def mycp(src,dst):
	os.system(f"cp {src} {dst}")

def mymv(src,dst):
	os.system(f"mv {src} {dst}")

def big_msg(msg):
	length = len(msg)
	print("-"*length)
	print(msg)
	print("-"*length)

def extract_tags(tags_setup ,
				 labels):	
	tags_list = []
	for l in labels:
		if type(tags_setup[l])==int:
			tags_list += [tags_setup[l]]
		elif type(tags_setup[l])==np.ndarray:
			tags_list += list(tags_setup[l])
		else:
			raise Exception("Type not recognised.")

	return tags_list

def get_tags_from_setup(tags_setup , labels):
		""" Extract tags from setup dictionary """
		tags_list = extract_tags(tags_setup,labels)
		tags_list = [str(t) for t in tags_list]
		tags_list_string = ",".join(tags_list)

		return tags_list_string

def meshtool_extract_biatrial(meshname,
							  output_folder,
							  tags_setup,
							  overwrite=False):
	
	do = True
	if os.path.exists(output_folder+'/biatrial/biatrial.pts') and not overwrite:
		do = False
		print('Mesh already found. Not overwriting.')

	if do:
		tags_list_string = get_tags_from_setup(tags_setup,["LA","RA"])

		big_msg("Extracting biatrial mesh...")
		cmd = "meshtool extract mesh -msh="+meshname+" -submsh="+output_folder+"/biatrial/biatrial -tags="+tags_list_string
		print(cmd)
		os.system(cmd)

def meshtool_extract_LA(output_folder,
						tags_setup,
						overwrite=False):

	do = True
	if os.path.exists(output_folder+'/la/la.pts') and not overwrite:
		do = False
		print('Mesh already found. Not overwriting.')

	if do:
		tags_list_string = get_tags_from_setup(tags_setup,["LA"])
		
		big_msg("Extracting LA mesh from biatrial...")		
		cmd = "meshtool extract mesh -msh="+output_folder+"/biatrial/biatrial -submsh="+output_folder+"/la/la -tags="+tags_list_string
		os.system(cmd)

def meshtool_extract_RA(output_folder,
						tags_setup,
						overwrite=False):

	do = True
	if os.path.exists(output_folder+'/ra/ra.pts') and not overwrite:
		do = False
		print('Mesh already found. Not overwriting.')

	if do:
		tags_list_string = get_tags_from_setup(tags_setup,["RA"])
		
		big_msg("Extracting RA mesh from biatrial...")
		cmd = f"meshtool extract mesh -msh={output_folder}/biatrial/biatrial -submsh={output_folder}/ra/ra -tags={tags_list_string}"
		os.system(cmd)

def meshtool_extract_surfaces(meshname,
		 					  output_folder,
		 					  tags_setup,
		 					  export_sup_inf=False,
		 					  rm_vv_from_aa=False,
		 					  surface="epi"):
	
	def extract_surface(surf, tag_init, tag_end, conn=":"):
		""" Creates meshtool commands from function inputs """
		return f"meshtool extract surface -msh={meshname} -surf={output_folder}/tmp/{surf} -op={tag_init}{conn}{tag_end} -ofmt=vtk"

	tag_la_string = get_tags_from_setup(tags_setup,["LA"])
	tag_ra_string = get_tags_from_setup(tags_setup,["RA"])

	aux_la_rm_labels = ["mitral","RSPV","RIPV","LSPV","LIPV","LAA","PV_planes"]
	aux_ra_rm_labels = ["SVC","IVC","tricuspid","VC_planes"]
	if rm_vv_from_aa:
		aux_la_rm_labels.append("LV")
		aux_ra_rm_labels.append("RV")
	
	tag_la_rm_string = get_tags_from_setup(tags_setup,aux_la_rm_labels)
	tag_ra_rm_string = get_tags_from_setup(tags_setup,aux_ra_rm_labels)

	tag_mv_string = get_tags_from_setup(tags_setup,["mitral"])
	tag_tv_string = get_tags_from_setup(tags_setup,["tricuspid"])
	tag_rpv_string = get_tags_from_setup(tags_setup,["RSPV","RIPV"])
	tag_svc_string = get_tags_from_setup(tags_setup,["SVC"])
	tag_lpv_string = get_tags_from_setup(tags_setup,["LSPV","LIPV"])
	tag_ivc_string = get_tags_from_setup(tags_setup,["IVC"])
	tag_pv_planes_string = get_tags_from_setup(tags_setup,["PV_planes"])
	tag_vc_planes_string = get_tags_from_setup(tags_setup,["VC_planes"])
	tag_lv_string = get_tags_from_setup(tags_setup,["LV"])
	tag_rv_string = get_tags_from_setup(tags_setup,["RV"])

	os.system(f"mkdir -p {output_folder}/tmp")

	surf_params = [
		("mitral", "Extracting mitral ring...", tag_mv_string),
		("rpv", "Extracting right PV ring...", tag_rpv_string),
		("lpv", "Extracting left PV ring...", tag_lpv_string),
	]

	for surf_type, msg, my_tag_end in surf_params:
		big_msg(msg)
		os.system(extract_surface(surf_type, tag_la_string, my_tag_end))


	if export_sup_inf:

		tag_rspv_string = get_tags_from_setup(tags_setup,["RSPV"])
		tag_ripv_string = get_tags_from_setup(tags_setup,["RIPV"])
		tag_lspv_string = get_tags_from_setup(tags_setup,["LSPV"])
		tag_lipv_string = get_tags_from_setup(tags_setup,["LIPV"])
		tag_laa_string = get_tags_from_setup(tags_setup,["LAA"])

		surf_params = [
			("laa", "Extracting LAA ring...", tag_laa_string),
			("ripv", "Extracting right inferior PV ring...", tag_ripv_string),
			("lipv", "Extracting left inferior PV ring...", tag_lipv_string),
			("rspv", "Extracting right superior PV ring...", tag_rspv_string),
			("lspv", "Extracting left superior PV ring...", tag_lspv_string)
		]

		for surf, msg, my_tag_end in surf_params:
			big_msg(msg)
			os.system(extract_surface(surf, tag_la_string, my_tag_end))

	big_msg("Extracting PV planes...")
	os.system(extract_surface("pv_planes", tag_la_string, tag_pv_planes_string))

	if export_sup_inf and surface=="endo":

		if "RSPV_vp" not in tags_setup:
			raise Exception("If you want to compute the landmarks on the endocardium, you need to provide separate tags for the valve planes")
		
		tag_rspv_vp_string = get_tags_from_setup(tags_setup,["RSPV_vp"])
		tag_ripv_vp_string = get_tags_from_setup(tags_setup,["RIPV_vp"])
		tag_lspv_vp_string = get_tags_from_setup(tags_setup,["LSPV_vp"])
		tag_lipv_vp_string = get_tags_from_setup(tags_setup,["LIPV_vp"])
		tag_laa_vp_string = get_tags_from_setup(tags_setup,["LAA_vp"])

		surf_params = [
			("laa_vp", "Extracting LAA valve plane...", tag_laa_vp_string),
			("ripv_vp", "Extracting right inferior PV valve plane...", tag_ripv_vp_string),
			("lipv_vp", "Extracting left inferior PV valve plane...", tag_lipv_vp_string),
			("rspv_vp", "Extracting right superior PV valve plane...", tag_rspv_vp_string),
			("lspv_vp", "Extracting left superior PV valve plane...", tag_lspv_vp_string)
		]

		for surf, msg, my_tag_end in surf_params:
			big_msg(msg)
			os.system(extract_surface(surf, tag_la_string, my_tag_end))
	
	big_msg("Extracting LA surface...")
	os.system(extract_surface("la", tag_la_string, tag_la_rm_string, conn="-"))

	surf_params = [
		("la_lv", "Extracting LA and LV intersection...", tag_la_string, tag_lv_string),
		("tricuspid", "Extracting tricuspid ring...", tag_ra_string, tag_tv_string),
		("svc", "Extracting SVC ring...", tag_ra_string, tag_svc_string),
		("ivc", "Extracting IVC ring...", tag_ra_string, tag_ivc_string)
	]

	for surf, msg, my_tag_init, my_tag_end in surf_params:
		big_msg(msg)
		os.system(extract_surface(surf, my_tag_init, my_tag_end))


	if surface=="endo":

		if "SVC_vp" not in tags_setup:
			raise Exception("If you want to compute the landmarks on the endocardium, you need to provide separate tags for the valve planes")
		
		tag_svc_vp_string = get_tags_from_setup(tags_setup,["SVC_vp"])
		tag_ivc_vp_string = get_tags_from_setup(tags_setup,["IVC_vp"])

		surf_params = [
			("svc_vp", "Extracting SVC valve plane...", tag_svc_vp_string),
			("ivc_vp", "Extracting IVC valve plane...", tag_ivc_vp_string)
		]

		for surf, msg, my_tag_end in surf_params:
			big_msg(msg)
			os.system(extract_surface(surf, tag_ra_string, my_tag_end))

	surf_params = [
		("pv_planes", "Extracting PV planes...", tag_la_string, tag_pv_planes_string, ":"),
		("vc_planes", "Extracting SVC IVC planes...", tag_ra_string, tag_vc_planes_string, ":"),
		("ra_rv", "Extracting RA and RV intersection...", tag_ra_string, tag_rv_string, ":"),
		("ra", "Extracting RA surface...", tag_ra_string, tag_ra_rm_string, "-")	
	]

	for surf, msg, my_tag_init, my_tag_end, my_conn in surf_params:
		big_msg(msg)
		os.system(extract_surface(surf, my_tag_init, my_tag_end, conn=my_conn))

	big_msg("Extracting LA surface connected components...")
	os.system(f"meshtool extract unreachable -msh={output_folder}/tmp/la.surfmesh.vtk -submsh={output_folder}/tmp/la_cc -ofmt=carp_txt")

	big_msg("Extracting RA surface connected components...")
	os.system(f"meshtool extract unreachable -msh={output_folder}/tmp/ra.surfmesh.vtk -submsh={output_folder}/tmp/ra_cc -ofmt=carp_txt")

	# Finish extraction of labels, moving on to the next step
	def find_and_append(list_of_files, prefix) :
		elements = []
		ix = 0
		is_file = True 
		while is_file :
			elemname = f"{prefix}.part{ix}"
			if f"{elemname}.elem" in list_of_files :
				elements.append(elemname)
			else :
				is_file = False
			ix += 1
		return elements
	
	def sort_and_clean(cc_list) : 
		if len(cc_list) > 2 : 
			cc_size = np.zeros((len(cc_list),),dtype=int)
			for i,cc in enumerate(cc_list):
				surf = read_elem(f"{output_folder}/tmp/{cc}.elem",el_type="Tr",tags=False)
				cc_size[i] = surf.shape[0]

			cc_list_old = copy.deepcopy(cc_list)
			sorted_size = np.argsort(cc_size)
			cc_list[0] = cc_list_old[sorted_size[-1]]
			cc_list[1] = cc_list_old[sorted_size[-2]]

			for i in range(len(cc_list)-2):
				os.system(f"rm {output_folder}/tmp/{cc_list_old[sorted_size[i]]}.*")
		
		return cc_list
		
	def id_and_save_endo_epi(cc_list, atrium : str) : 
		pts0 = read_pts(f"{output_folder}/tmp/{cc_list[0]}.pts")
		surf0 = read_elem(f"{output_folder}/tmp/{cc_list[0]}.elem",el_type="Tr",tags=False)

		pts1 = read_pts(f"{output_folder}/tmp/{cc_list[1]}.pts")
		surf1 = read_elem(f"{output_folder}/tmp/{cc_list[1]}.elem",el_type="Tr",tags=False)

		cog0 = np.mean(pts0,axis=0)
		is_outward = np.zeros((surf0.shape[0],),dtype=int)
		for i,t in enumerate(surf0):
			v0 = pts0[t[1],:] - pts0[t[0],:]
			v0 = v0/np.linalg.norm(v0)

			v1 = pts0[t[2],:] - pts0[t[0],:]
			v1 = v1/np.linalg.norm(v1)

			n = np.cross(v0,v1)
			n = n/np.linalg.norm(n)

			dot_prod = np.dot(cog0-pts0[t[0],:],n)

			if dot_prod>0:
				is_outward[i] = 1
		
		if np.sum(is_outward)/surf0.shape[0]>0.7:
			print(f'{cc_list[0]} is the epicardium')
			print(f'{cc_list[1]} is the endocardium')
			endo = 0
			epi = 1
		else:
			print(f'{cc_list[1]} is the epicardium')
			print(f'{cc_list[0]} is the endocardium')	
			endo = 1
			epi = 0
		
		print(f'Renaming {atrium.upper()} connected components...')
		formats = ["nod","eidx","elem","lon","pts"]
		for f in formats:
			mymv(f"{output_folder}/tmp/{cc_list[endo]}.{f}", f"{output_folder}/tmp/{atrium.lower()}_endo.{f}")
			mymv(f"{output_folder}/tmp/{cc_list[epi]}.{f}", f"{output_folder}/tmp/{atrium.lower()}_epi.{f}")
	
	tmp_files = os.listdir(output_folder+"/tmp")
	la_cc = find_and_append(tmp_files, "la_cc")
	ra_cc = find_and_append(tmp_files, "ra_cc")

	print("Checking connected component size and keeping only the two biggest...")
	la_cc = sort_and_clean(la_cc)
	ra_cc = sort_and_clean(ra_cc)

	id_and_save_endo_epi(la_cc, "la")
	id_and_save_endo_epi(ra_cc, "ra")

def export_LA_vtk_msh(output_folder,
					  tags_setup):

	tets = read_elem(output_folder+"/la/la.elem",el_type="Tt",tags=False)
	pts = read_pts(output_folder+"/la/la.pts")
	
	surface_list = [output_folder+"/tmp/la.surf",
					output_folder+"/tmp/mitral.surf",
					output_folder+"/tmp/rpv.surf",
					output_folder+"/tmp/lpv.surf",
					output_folder+"/tmp/pv_planes.surf",
					output_folder+"/tmp/la_lv.surf"]

	surface_list_string = ','.join(surface_list)

	surface_bia_list = [output_folder+"/biatrial/la.surf",
						output_folder+"/biatrial/mitral.surf",
						output_folder+"/biatrial/rpv.surf",
						output_folder+"/biatrial/lpv.surf",
						output_folder+"/biatrial/pv_planes.surf",
						output_folder+"/biatrial/la_lv.surf"]
	surface_bia_list_string = ','.join(surface_bia_list)

	print('Mapping surfaces onto biatrial mesh...')
	cmd = "meshtool map -submsh="+output_folder+"/biatrial/biatrial -files="+surface_list_string+" -outdir="+output_folder+"/biatrial/ -mode=m2s"
	os.system(cmd)

	print('Mapping surfaces onto LA mesh...')
	cmd = "meshtool map -submsh="+output_folder+"/la/la -files="+surface_bia_list_string+" -outdir="+output_folder+"/la/ -mode=m2s"
	os.system(cmd)

	la_tr = read_elem(output_folder+"/la/la.surf",el_type="Tr",tags=False)

	la_endo_eidx = read_nod_eidx(output_folder+"/tmp/la_endo.eidx")
	la_epi_eidx = read_nod_eidx(output_folder+"/tmp/la_epi.eidx")

	la_endo_tr = la_tr[la_endo_eidx,:]
	la_epi_tr = la_tr[la_epi_eidx,:]

	write_surf_caroline(output_folder+"/la/la_epi.surf",la_epi_tr)
	write_surf_caroline(output_folder+"/la/la_endo.surf",la_endo_tr)

	mitral_tr = read_elem(output_folder+"/la/mitral.surf",el_type="Tr",tags=False)
	rpv_tr = read_elem(output_folder+"/la/rpv.surf",el_type="Tr",tags=False)
	lpv_tr = read_elem(output_folder+"/la/lpv.surf",el_type="Tr",tags=False)
	pv_planes_tr = read_elem(output_folder+"/la/pv_planes.surf",el_type="Tr",tags=False)
	la_lv_tr = read_elem(output_folder+"/la/la_lv.surf",el_type="Tr",tags=False)

	la_endo_tr = np.concatenate((la_endo_tr,pv_planes_tr),axis=0)
	la_epi_tr = np.concatenate((la_epi_tr,la_lv_tr),axis=0)

	la_surface_tr = np.concatenate((la_endo_tr,la_epi_tr,mitral_tr,rpv_tr,lpv_tr),axis=0)

	tets_tags = np.zeros((tets.shape[0],),dtype=int)+tags_setup["LA"]["body"]
	surf_tags = np.zeros((la_surface_tr.shape[0],),dtype=int)

	surf_tags[:la_endo_tr.shape[0]] = tags_setup["LA"]["endo"]
	surf_tags[la_endo_tr.shape[0]:la_endo_tr.shape[0]+la_epi_tr.shape[0]] = tags_setup["LA"]["epi"]
	surf_tags[la_endo_tr.shape[0]+la_epi_tr.shape[0]:la_endo_tr.shape[0]+la_epi_tr.shape[0]+mitral_tr.shape[0]] = tags_setup["LA"]["mitral"]
	surf_tags[la_endo_tr.shape[0]+la_epi_tr.shape[0]+mitral_tr.shape[0]:la_endo_tr.shape[0]+la_epi_tr.shape[0]+mitral_tr.shape[0]+rpv_tr.shape[0]] = tags_setup["LA"]["RPV"]
	surf_tags[la_endo_tr.shape[0]+la_epi_tr.shape[0]+mitral_tr.shape[0]+rpv_tr.shape[0]:] = tags_setup["LA"]["LPV"]

	msh = meshio.Mesh(pts,
    				  [('tetra',tets),('triangle',la_surface_tr)],
    				  cell_data={'Ids': [tets_tags,surf_tags]})	

	msh.write(output_folder+"/la/la_tag.vtu")

def export_vtk_meshes_caroline(output_folder,
							   raa_apex_file=None,
							   scale_factor=1.0,
							   surface="epi"):

	tets = read_elem(output_folder+"/la/la.elem",el_type="Tt",tags=False)
	pts = read_pts(output_folder+"/la/la.pts")

	surface_names = ["la", "mitral", "ripv", "rspv", "lipv", "lspv", "laa", "pv_planes", "la_lv"]
	if surface=="endo":
		surface_names += ["ripv_vp", "rspv_vp", "lipv_vp", "lspv_vp", "laa_vp"]
	
	surface_list = [os.path.join(output_folder, f"/tmp/{name}.surf") for name in surface_names]
	surface_list_string = ','.join(surface_list)

	surface_bia_list = [os.path.join(output_folder, f"/biatrial/{name}.surf") for name in surface_names]
	surface_bia_list_string = ','.join(surface_bia_list)

	print('Mapping surfaces onto biatrial mesh...')
	cmd = f"meshtool map -submsh={output_folder}/biatrial/biatrial -files={surface_list_string} -outdir={output_folder}/biatrial/ -mode=m2s"
	os.system(cmd)

	print('Mapping surfaces onto LA mesh...')
	cmd = f"meshtool map -submsh={output_folder}/la/la -files={surface_bia_list_string} -outdir={output_folder}/la/ -mode=m2s"
	os.system(cmd)

	la_tr = read_elem(output_folder+"/la/la.surf",el_type="Tr",tags=False)

	la_endo_eidx = read_nod_eidx(output_folder+"/tmp/la_endo.eidx")
	la_epi_eidx = read_nod_eidx(output_folder+"/tmp/la_epi.eidx")

	la_endo_tr = la_tr[la_endo_eidx,:]
	la_epi_tr = la_tr[la_epi_eidx,:]

	write_surf_caroline(output_folder+"/la/la_epi.surf",la_epi_tr)
	write_surf_caroline(output_folder+"/la/la_endo.surf",la_endo_tr)

	tets = read_elem(output_folder+"/ra/ra.elem",el_type="Tt",tags=False)
	pts = read_pts(output_folder+"/ra/ra.pts")
	
	surface_names = ["ra", "tricuspid", "svc", "ivc", "vc_planes", "ra_rv"]
	if surface=="endo":
		surface_names += ["svc_vp", "ivc_vp"]

	surface_list = [os.path.join(output_folder, f"/tmp/{name}.surf") for name in surface_names]
	surface_list_string = ','.join(surface_list)

	surface_bia_list = [os.path.join(output_folder, f"/biatrial/{name}.surf") for name in surface_names]
	surface_bia_list_string = ','.join(surface_bia_list)

	print('Mapping surfaces onto biatrial mesh...')
	cmd = f"meshtool map -submsh={output_folder}/biatrial/biatrial -files={surface_list_string} -outdir={output_folder}/biatrial/ -mode=m2s"
	os.system(cmd)

	print('Mapping surfaces onto RA mesh...')
	cmd = f"meshtool map -submsh={output_folder}/ra/ra -files={surface_bia_list_string} -outdir={output_folder}/ra/ -mode=m2s"
	os.system(cmd)

	ra_tr = read_elem(output_folder+"/ra/ra.surf",el_type="Tr",tags=False)

	ra_endo_eidx = read_nod_eidx(output_folder+"/tmp/ra_endo.eidx")
	ra_epi_eidx = read_nod_eidx(output_folder+"/tmp/ra_epi.eidx")

	ra_endo_tr = ra_tr[ra_endo_eidx,:]
	ra_epi_tr = ra_tr[ra_epi_eidx,:]

	ra_rv_tr = read_elem(output_folder+"/ra/ra_rv.surf",el_type="Tr",tags=False)
	ra_epi_tr = np.concatenate((ra_epi_tr,ra_rv_tr),axis=0)

	write_surf_caroline(output_folder+"/ra/ra_epi.surf",ra_epi_tr)
	write_surf_caroline(output_folder+"/ra/ra_endo.surf",ra_endo_tr)

	big_msg('Finding LA and RA landmarks...')

	find_landmarks(output_folder,
				   surface=surface,
				   scale_factor=scale_factor,
				   raa_apex_file=raa_apex_file)

	big_msg('Organising folders...')

	os.system("mkdir -p "+output_folder+"/LA_endo/")
	os.system("mkdir -p "+output_folder+"/LA_epi/")
	os.system("mkdir -p "+output_folder+"/RA_endo/")
	os.system("mkdir -p "+output_folder+"/RA_epi/")

	landmarks_dic = {
		'LA_endo' : [f"la/{fn}" for fn in ['prodLaLandmarks.txt','prodLaRegion.txt']], 
		'LA_epi' : [f"la/{fn}" for fn in ['prodLaLandmarks.txt','prodLaRegion.txt']],
		'RA_endo' : [f"ra/{fn}" for fn in ['prodRaLandmarks.txt','prodRaRegion.txt']],
		'RA_epi' : [f"ra/{fn}" for fn in ['prodRaLandmarks.txt','prodRaRegion.txt']]
	}

	for key in landmarks_dic.keys(): 
		for fn in landmarks_dic[key] : 
			mycp(os.path.join(output_folder,fn), os.path.join(output_folder,key+"/"))

	for a in ["la","ra"]:
		for l in ["endo","epi"]:
			os.system(f"meshtool convert -imsh={output_folder}/tmp/{a}_{l} -ofmt=vtk_polydata -omsh={output_folder}/{a.upper()}_{l}/{a.upper()}_{l}")

	# os.system("meshtool convert -imsh="+output_folder+"/tmp/la_endo -ofmt=vtk_polydata -omsh="+output_folder+"/LA_endo/LA_endo")
	# os.system("meshtool convert -imsh="+output_folder+"/tmp/ra_endo -ofmt=vtk_polydata -omsh="+output_folder+"/RA_endo/RA_endo")
	# mycp(os.path.join(output_folder,"la","prodLaLandmarks.txt"), os.path.join(output_folder,"LA_endo/"))
	# mycp(os.path.join(output_folder,"la","prodLaRegion.txt"), os.path.join(output_folder,"LA_endo/"))
	# mycp(os.path.join(output_folder,"ra","prodRaLandmarks.txt"), os.path.join(output_folder,"RA_endo/"))
	# mycp(os.path.join(output_folder,"ra","prodRaRegion.txt"), os.path.join(output_folder,"RA_endo/"))

	# os.system("meshtool convert -imsh="+output_folder+"/tmp/la_epi -ofmt=vtk_polydata -omsh="+output_folder+"/LA_epi/LA_epi")
	# os.system("meshtool convert -imsh="+output_folder+"/tmp/ra_epi -ofmt=vtk_polydata -omsh="+output_folder+"/RA_epi/RA_epi")
	# mycp(os.path.join(output_folder,"la","prodLaLandmarks.txt"), os.path.join(output_folder,"LA_epi/"))
	# mycp(os.path.join(output_folder,"la","prodLaRegion.txt"), os.path.join(output_folder,"LA_epi/"))
	# mycp(os.path.join(output_folder,"ra","prodRaLandmarks.txt"), os.path.join(output_folder,"RA_epi/"))
	# mycp(os.path.join(output_folder,"ra","prodRaRegion.txt"), os.path.join(output_folder,"RA_epi/"))

def recompute_raa_base(output_folder,
					   raa_apex_file,
					   landmarks_file,
					   scale=1.0,
					   surface="epi"):
	
	if os.path.exists(raa_apex_file):
			landmarks_raa_apex = np.loadtxt(raa_apex_file)
			landmark_raa_base = find_raa_base(output_folder,landmarks_raa_apex,surface=surface)
			landmarks_raa_apex = np.reshape(landmarks_raa_apex,(1,3))
			region_landmarks_raa = np.concatenate((landmarks_raa_apex,landmark_raa_base),axis=0)
	else:
		raise Exception("Cannot find apex file.")

	region_landmarks = np.loadtxt(landmarks_file,dtype=float,delimiter=',')
	region_landmarks[-2:,:] = region_landmarks_raa
	np.savetxt(landmarks_file,region_landmarks*scale,delimiter=',')

def scale_landmarks(landmarks_file,
					scale=1.0):

	landmarks = np.loadtxt(landmarks_file,dtype=float,delimiter=',')
	np.savetxt(landmarks_file,landmarks*scale,delimiter=',')

def export_RA_vtk_msh(output_folder,
					  tags_setup,
					  r_geodesic=1000.0):

	tets = read_elem(output_folder+"/ra/ra.elem",el_type="Tt",tags=False)
	pts = read_pts(output_folder+"/ra/ra.pts")
	
	surface_list = [output_folder+"/tmp/ra.surf",
					output_folder+"/tmp/tricuspid.surf",
					output_folder+"/tmp/svc.surf",
					output_folder+"/tmp/ivc.surf",
					output_folder+"/tmp/vc_planes.surf",
					output_folder+"/tmp/ra_rv.surf"]
	surface_list_string = ','.join(surface_list)

	surface_bia_list = [output_folder+"/biatrial/ra.surf",
						output_folder+"/biatrial/tricuspid.surf",
						output_folder+"/biatrial/svc.surf",
						output_folder+"/biatrial/ivc.surf",
						output_folder+"/biatrial/vc_planes.surf",
						output_folder+"/biatrial/ra_rv.surf"]
	surface_bia_list_string = ','.join(surface_bia_list)

	print('Mapping surfaces onto biatrial mesh...')
	cmd = "meshtool map -submsh="+output_folder+"/biatrial/biatrial -files="+surface_list_string+" -outdir="+output_folder+"/biatrial/ -mode=m2s"
	os.system(cmd)

	print('Mapping surfaces onto RA mesh...')
	cmd = "meshtool map -submsh="+output_folder+"/ra/ra -files="+surface_bia_list_string+" -outdir="+output_folder+"/ra/ -mode=m2s"
	os.system(cmd)

	ra_tr = read_elem(output_folder+"/ra/ra.surf",el_type="Tr",tags=False)

	ra_endo_eidx = read_nod_eidx(output_folder+"/tmp/ra_endo.eidx")
	ra_epi_eidx = read_nod_eidx(output_folder+"/tmp/ra_epi.eidx")

	ra_endo_tr = ra_tr[ra_endo_eidx,:]
	ra_epi_tr = ra_tr[ra_epi_eidx,:]

	write_surf_caroline(output_folder+"/ra/ra_epi.surf",ra_epi_tr)
	write_surf_caroline(output_folder+"/ra/ra_endo.surf",ra_endo_tr)

	tricuspid_tr = read_elem(output_folder+"/ra/tricuspid.surf",el_type="Tr",tags=False)
	svc_tr = read_elem(output_folder+"/ra/svc.surf",el_type="Tr",tags=False)
	ivc_tr = read_elem(output_folder+"/ra/ivc.surf",el_type="Tr",tags=False)
	vc_planes_tr = read_elem(output_folder+"/ra/vc_planes.surf",el_type="Tr",tags=False)
	ra_rv_tr = read_elem(output_folder+"/ra/ra_rv.surf",el_type="Tr",tags=False)

	ra_endo_tr = np.concatenate((ra_endo_tr,vc_planes_tr),axis=0)
	ra_epi_tr = np.concatenate((ra_epi_tr,ra_rv_tr),axis=0)

	ra_surface_tr = np.concatenate((ra_endo_tr,ra_epi_tr,tricuspid_tr,svc_tr,ivc_tr),axis=0)

	big_msg("Extracting points for SVC IVC geodesic...")
	idx_geodesic,anterior_posterior_tag = find_SVC_IVC_geodesic(output_folder,r_geodesic=r_geodesic)
	tags_ra_epi = np.zeros((ra_epi_tr.shape[0],),dtype=int)+tags_setup["RA"]["epi"]
	tags_ra_epi[idx_geodesic] = tags_setup["RA"]["roof_line"]

	tricuspid_ant_post_tag = np.zeros((tricuspid_tr.shape[0],),dtype=int)+tags_setup["RA"]["tricuspid_anterior"]
	tricuspid_ant_post_tag[np.where(anterior_posterior_tag==1)[0]] = tags_setup["RA"]["tricuspid_posterior"]

	tets_tags = np.zeros((tets.shape[0],),dtype=int)+tags_setup["RA"]["body"]
	surf_tags = np.zeros((ra_surface_tr.shape[0],),dtype=int)

	surf_tags[:ra_endo_tr.shape[0]] = tags_setup["RA"]["endo"]
	surf_tags[ra_endo_tr.shape[0]:ra_endo_tr.shape[0]+ra_epi_tr.shape[0]] = tags_ra_epi
	surf_tags[ra_endo_tr.shape[0]+ra_epi_tr.shape[0]:ra_endo_tr.shape[0]+ra_epi_tr.shape[0]+tricuspid_tr.shape[0]] = tricuspid_ant_post_tag
	surf_tags[ra_endo_tr.shape[0]+ra_epi_tr.shape[0]+tricuspid_tr.shape[0]:ra_endo_tr.shape[0]+ra_epi_tr.shape[0]+tricuspid_tr.shape[0]+svc_tr.shape[0]] = tags_setup["RA"]["SVC"]
	surf_tags[ra_endo_tr.shape[0]+ra_epi_tr.shape[0]+tricuspid_tr.shape[0]+svc_tr.shape[0]:] = tags_setup["RA"]["IVC"]

	msh = meshio.Mesh(pts,
    				  [('tetra',tets),('triangle',ra_surface_tr)],
    				  cell_data={'Ids': [tets_tags,surf_tags]})	

	msh.write(output_folder+"/ra/ra_tag.vtu")

def meshtool_extract_base(mesh,surf_folder,input_tags):

	tags_list_vent_string = get_tags_from_setup(input_tags, ["LV","RV"])
	tags_list_VPs_string = get_tags_from_setup(input_tags, ["MV","TV","AV","PV"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/myocardium.base -ofmt=vtk -op="+tags_list_vent_string+":"+tags_list_VPs_string)

def meshtool_extract_surfaces_lv_rv_epi(mesh,surf_folder,input_tags):

	tags_list_vent = extract_tags(input_tags,["LV","RV"])
	tags_list_vent_string = get_tags_from_setup(input_tags, ["LV","RV"])
	tags_list_VPs_string = get_tags_from_setup(input_tags, ["MV","TV","AV","PV"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/epi_endo -ofmt=vtk -op="+tags_list_vent_string+"-"+tags_list_VPs_string)
	os.system("meshtool extract unreachable -msh="+surf_folder+"/tmp/epi_endo.surfmesh -ifmt=vtk -ofmt=vtk -ofmt=carp_txt -submsh="+surf_folder+"/tmp/epi_endo_CC")

	tmp_files = os.listdir(surf_folder+"/tmp")
	epi_endo_CC = []
	i = 0
	isfile=True
	while isfile:
		if "epi_endo_CC.part"+str(i)+".elem" in tmp_files:
			epi_endo_CC.append("epi_endo_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the 3 biggest...")
	if len(epi_endo_CC)>3:
		CC_size = np.zeros((len(epi_endo_CC),),dtype=int)
		for i,CC in enumerate(epi_endo_CC):
			surf = read_elem(surf_folder+"/tmp/"+CC+".elem",el_type="Tr",tags=False)
			CC_size[i] = surf.shape[0]

		epi_endo_CC_old = copy.deepcopy(epi_endo_CC)
		sorted_size = np.argsort(CC_size)
		epi_endo_CC[0] = epi_endo_CC_old[sorted_size[-1]]
		epi_endo_CC[1] = epi_endo_CC_old[sorted_size[-2]]
		epi_endo_CC[2] = epi_endo_CC_old[sorted_size[-3]]

		for i in range(len(epi_endo_CC)-3):
			os.system("rm "+surf_folder+"/tmp/"+epi_endo_CC_old[sorted_size[i]]+".*")

	# Find CoGs of surfaces
	pts0 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[0]+".pts")
	surf0 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[0]+".elem",el_type="Tr",tags=False)

	pts1 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[1]+".pts")
	surf1 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[1]+".elem",el_type="Tr",tags=False)

	pts2 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[2]+".pts")
	surf2 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[2]+".elem",el_type="Tr",tags=False)

	cog0 = np.mean(pts0,axis=0)
	cog1 = np.mean(pts1,axis=0)
	cog2 = np.mean(pts2,axis=0)

	# Find CoG for LV blood pool
	mesh_pts = read_pts(mesh+".pts")
	mesh_elem = read_elem(mesh+".elem",el_type="Tt",tags=True)

	lv_pts_idx = []
	for i,e in enumerate(mesh_elem):
		if e[4] == int(tags_list_vent[0]):	
			lv_pts_idx.append(int(e[1]))
			lv_pts_idx.append(int(e[2]))
			lv_pts_idx.append(int(e[3]))

	lv_pts_idx = np.unique(lv_pts_idx)
	lv_pts = np.zeros((lv_pts_idx.shape[0],3))

	for i,p in enumerate(lv_pts_idx):
		lv_pts[i] = mesh_pts[p]

	cog_lv = np.mean(lv_pts,axis=0)


	# Finding distance between surface CoGs and LV_BP CoG
	dist0 = np.zeros((1,3))
	dist1 = np.zeros((1,3))
	dist2 = np.zeros((1,3))

	for i,c in enumerate(cog_lv):
		dist0[:,i] = c - cog0[i]
		dist1[:,i] = c - cog1[i]
		dist2[:,i] = c - cog2[i]

	dist0 = np.linalg.norm(dist0)
	dist1 = np.linalg.norm(dist1)
	dist2 = np.linalg.norm(dist2)

	######## Checking orientation of surface normals on surf0 ########
	epi = 'not_yet_found'
	is_outward = np.zeros((surf0.shape[0],),dtype=int)
	for i,t in enumerate(surf0):
		v0 = pts0[t[1],:] - pts0[t[0],:]
		v0 = v0/np.linalg.norm(v0)

		v1 = pts0[t[2],:] - pts0[t[0],:]
		v1 = v1/np.linalg.norm(v1)

		n = np.cross(v0,v1)
		n = n/np.linalg.norm(n)

		dot_prod = np.dot(cog0-pts0[t[0],:],n)

		if dot_prod<0:
			is_outward[i] = 1

	if np.sum(is_outward)/surf0.shape[0]>0.7:
		print(epi_endo_CC[0]+' is the epicardium')
		epi = 0 
		
		if dist1 < dist2:
			print(epi_endo_CC[1]+' is the LV endocardium')
			print(epi_endo_CC[2]+' is the RV endocardium')
			lv_endo = 1
			rv_endo = 2
		else:
			print(epi_endo_CC[1]+' is the RV endocardium')
			print(epi_endo_CC[2]+' is the LV endocardium')
			rv_endo = 1
			lv_endo = 2

	if epi == 'not_yet_found':
		######## Checking orientation of surface normals on surf1 ########
		is_outward = np.zeros((surf1.shape[0],),dtype=int)
		for i,t in enumerate(surf1):
			v0 = pts1[t[1],:] - pts1[t[0],:]
			v0 = v0/np.linalg.norm(v0)

			v1 = pts1[t[2],:] - pts1[t[0],:]
			v1 = v1/np.linalg.norm(v1)

			n = np.cross(v0,v1)
			n = n/np.linalg.norm(n)

			dot_prod = np.dot(cog1-pts1[t[0],:],n)

			if dot_prod<0:
				is_outward[i] = 1

		if np.sum(is_outward)/surf1.shape[0]>0.7:
			print(epi_endo_CC[1]+' is the epicardium')
			epi = 1 
	
			if dist0 < dist2:
				print(epi_endo_CC[0]+' is the LV endocardium')
				print(epi_endo_CC[2]+' is the RV endocardium')
				lv_endo = 0
				rv_endo = 2
			else:
				print(epi_endo_CC[0]+' is the RV endocardium')
				print(epi_endo_CC[2]+' is the LV endocardium')
				rv_endo = 0
				lv_endo = 2

	if epi == 'not_yet_found':
		######## Checking orientation of surface normals on surf2 ########
		is_outward = np.zeros((surf2.shape[0],),dtype=int)
		for i,t in enumerate(surf2):
			v0 = pts2[t[1],:] - pts2[t[0],:]
			v0 = v0/np.linalg.norm(v0)

			v1 = pts2[t[2],:] - pts2[t[0],:]
			v1 = v1/np.linalg.norm(v1)

			n = np.cross(v0,v1)
			n = n/np.linalg.norm(n)

			dot_prod = np.dot(cog2-pts2[t[0],:],n)

			if dot_prod<0:
				is_outward[i] = 1

		if np.sum(is_outward)/surf2.shape[0]>0.7:
			print(epi_endo_CC[2]+' is the epicardium')
			epi = 2 
	
			if dist0 < dist1:
				print(epi_endo_CC[0]+' is the LV endocardium')
				print(epi_endo_CC[1]+' is the RV endocardium')
				lv_endo = 0
				rv_endo = 1
			else:
				print(epi_endo_CC[0]+' is the RV endocardium')
				print(epi_endo_CC[1]+' is the LV endocardium')
				rv_endo = 0
				lv_endo = 1


	if epi == 'not_yet_found':
		raise Exception("Surfaces could not be identified. Program terminated.")

	print('Renaming connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[epi]+"."+f+" "+surf_folder+"/tmp/myocardium.epi."+f)
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[lv_endo]+"."+f+" "+surf_folder+"/tmp/myocardium.lvendo."+f)
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[rv_endo]+"."+f+" "+surf_folder+"/tmp/myocardium.rvendo."+f)

def meshtool_extract_septum(mesh,surf_folder,input_tags):

	tags_list_lv_string = get_tags_from_setup(input_tags, ["LV"])
	tags_list_remove_string = get_tags_from_setup(input_tags, ["RV","RA","PArt"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/myocardium.rvsept -ofmt=vtk -op="+tags_list_lv_string+"-"+tags_list_remove_string)
	os.system("meshtool extract unreachable -msh="+surf_folder+"/tmp/myocardium.rvsept.surfmesh -ifmt=vtk -ofmt=vtk -ofmt=carp_txt -submsh="+surf_folder+"/tmp/myocardium.rvsept_CC")

	tmp_files = os.listdir(surf_folder+"/tmp")
	rvsept_CC = []
	i = 0
	isfile=True
	while isfile:
		if "myocardium.rvsept_CC.part"+str(i)+".elem" in tmp_files:
			rvsept_CC.append("myocardium.rvsept_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the 2 biggest...")
	CC_size = np.zeros((len(rvsept_CC),),dtype=int)
	for i,CC in enumerate(rvsept_CC):
		surf = read_elem(surf_folder+"/tmp/"+CC+".elem",el_type="Tr",tags=False)
		CC_size[i] = surf.shape[0]
	
	rvsept_CC_old = copy.deepcopy(rvsept_CC)
	sorted_size = np.argsort(CC_size)
	rvsept_CC[0] = rvsept_CC_old[sorted_size[-1]]
	rvsept_CC[1] = rvsept_CC_old[sorted_size[-2]]

	if len(rvsept_CC)>2:
		for i in range(len(rvsept_CC)-2):
			pass
			# os.system("rm "+surf_folder+"/tmp/"+rvsept_CC_old[sorted_size[i]]+".*")

	print('Renaming connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+surf_folder+"/tmp/"+rvsept_CC[0]+"."+f+" "+surf_folder+"/tmp/lvepi."+f)
		os.system("mv "+surf_folder+"/tmp/"+rvsept_CC[1]+"."+f+" "+surf_folder+"/tmp/myocardium.rvsept."+f)

def meshtool_extract_la_base(mesh,surf_folder,input_tags):

	tags_list_la_string = get_tags_from_setup(input_tags, ["LA"])
	tags_list_lv_string = get_tags_from_setup(input_tags, ["LV"])
	tags_list_mv_string = get_tags_from_setup(input_tags, ["MV"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/la.base -ofmt=vtk -op="+tags_list_la_string+":"+tags_list_mv_string+","+tags_list_lv_string)

def meshtool_extract_la_surfaces(mesh,surf_folder,input_tags):

	tags_list_la_string = get_tags_from_setup(input_tags, ["LA"])
	tags_list_lv_string = get_tags_from_setup(input_tags, ["LV"])
	tags_list_VPs_string = get_tags_from_setup(input_tags, ["MV","TV","AV","PV","LSPV","LIPV","RSPV","RIPV","LAA","SVC","IVC"])
	tags_list_rings_string = get_tags_from_setup(input_tags, ["LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring","LAA_ring","SVC_ring","IVC_ring"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/epi_endo -ofmt=vtk -op="+tags_list_la_string+"-"+tags_list_lv_string+","+tags_list_VPs_string+","+tags_list_rings_string)
	os.system("meshtool extract unreachable -msh="+surf_folder+"/tmp/epi_endo.surfmesh -ifmt=vtk -ofmt=vtk -ofmt=carp_txt -submsh="+surf_folder+"/tmp/epi_endo_CC")

	tmp_files = os.listdir(surf_folder+"/tmp")
	epi_endo_CC = []
	i = 0
	isfile=True
	while isfile:
		if "epi_endo_CC.part"+str(i)+".elem" in tmp_files:
			epi_endo_CC.append("epi_endo_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the 2 biggest...")
	if len(epi_endo_CC)>2:
		CC_size = np.zeros((len(epi_endo_CC),),dtype=int)
		for i,CC in enumerate(epi_endo_CC):
			surf = read_elem(surf_folder+"/tmp/"+CC+".elem",el_type="Tr",tags=False)
			CC_size[i] = surf.shape[0]

		epi_endo_CC_old = copy.deepcopy(epi_endo_CC)
		sorted_size = np.argsort(CC_size)
		epi_endo_CC[0] = epi_endo_CC_old[sorted_size[-1]]
		epi_endo_CC[1] = epi_endo_CC_old[sorted_size[-2]]

		for i in range(len(epi_endo_CC)-2):
			os.system("rm "+surf_folder+"/tmp/"+epi_endo_CC_old[sorted_size[i]]+".*")

	# Find CoGs of surfaces
	pts0 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[0]+".pts")
	surf0 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[0]+".elem",el_type="Tr",tags=False)

	pts1 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[1]+".pts")
	surf1 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[1]+".elem",el_type="Tr",tags=False)

	cog0 = np.mean(pts0,axis=0)
	cog1 = np.mean(pts1,axis=0)

	######## Checking orientation of surface normals on surf0 ########
	is_outward = np.zeros((surf0.shape[0],),dtype=int)
	for i,t in enumerate(surf0):
		v0 = pts0[t[1],:] - pts0[t[0],:]
		v0 = v0/np.linalg.norm(v0)

		v1 = pts0[t[2],:] - pts0[t[0],:]
		v1 = v1/np.linalg.norm(v1)

		n = np.cross(v0,v1)
		n = n/np.linalg.norm(n)

		dot_prod = np.dot(cog0-pts0[t[0],:],n)

		if dot_prod<0:
			is_outward[i] = 1

	if np.sum(is_outward)/surf0.shape[0]>0.7:
		print(epi_endo_CC[0]+' is the epicardium')
		print(epi_endo_CC[1]+' is the endocardium')
		epi=0
		endo=1
	else:
		print(epi_endo_CC[0]+' is the endocardium')
		print(epi_endo_CC[1]+' is the epicardium')
		endo=0
		epi=1

	print('Renaming connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[epi]+"."+f+" "+surf_folder+"/tmp/la.epi."+f)
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[endo]+"."+f+" "+surf_folder+"/tmp/la.lvendo."+f)

def meshtool_extract_ra_base(mesh,surf_folder,input_tags):
	tags_list_ra_string = get_tags_from_setup(input_tags, ["RA"])
	tags_list_rv_string = get_tags_from_setup(input_tags, ["RV"])
	tags_list_tv_string = get_tags_from_setup(input_tags, ["TV"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/ra.base -ofmt=vtk -op="+tags_list_ra_string+":"+tags_list_tv_string+","+tags_list_tv_string)

def meshtool_extract_ra_surfaces(mesh,surf_folder,input_tags):
	tags_list_ra_string = get_tags_from_setup(input_tags, ["RA"])
	tags_list_rv_string = get_tags_from_setup(input_tags, ["RV"])
	tags_list_VPs_string = get_tags_from_setup(input_tags, ["MV","TV","AV","PV","LSPV","LIPV","RSPV","RIPV","LAA","SVC","IVC"])
	tags_list_rings_string = get_tags_from_setup(input_tags, ["LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring","LAA_ring","SVC_ring","IVC_ring"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/epi_endo -ofmt=vtk -op="+tags_list_ra_string+"-"+tags_list_rv_string+","+tags_list_VPs_string+","+tags_list_rings_string)
	os.system("meshtool extract unreachable -msh="+surf_folder+"/tmp/epi_endo.surfmesh -ifmt=vtk -ofmt=vtk -ofmt=carp_txt -submsh="+surf_folder+"/tmp/epi_endo_CC")

	tmp_files = os.listdir(surf_folder+"/tmp")
	epi_endo_CC = []
	i = 0
	isfile=True
	while isfile:
		if "epi_endo_CC.part"+str(i)+".elem" in tmp_files:
			epi_endo_CC.append("epi_endo_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the 2 biggest...")
	if len(epi_endo_CC)>2:
		CC_size = np.zeros((len(epi_endo_CC),),dtype=int)
		for i,CC in enumerate(epi_endo_CC):
			surf = read_elem(surf_folder+"/tmp/"+CC+".elem",el_type="Tr",tags=False)
			CC_size[i] = surf.shape[0]

		epi_endo_CC_old = copy.deepcopy(epi_endo_CC)
		sorted_size = np.argsort(CC_size)
		epi_endo_CC[0] = epi_endo_CC_old[sorted_size[-1]]
		epi_endo_CC[1] = epi_endo_CC_old[sorted_size[-2]]

		for i in range(len(epi_endo_CC)-2):
			os.system("rm "+surf_folder+"/tmp/"+epi_endo_CC_old[sorted_size[i]]+".*")

	# Find CoGs of surfaces
	pts0 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[0]+".pts")
	surf0 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[0]+".elem",el_type="Tr",tags=False)

	pts1 = read_pts(surf_folder+"/tmp/"+epi_endo_CC[1]+".pts")
	surf1 = read_elem(surf_folder+"/tmp/"+epi_endo_CC[1]+".elem",el_type="Tr",tags=False)

	cog0 = np.mean(pts0,axis=0)
	cog1 = np.mean(pts1,axis=0)

	######## Checking orientation of surface normals on surf0 ########
	is_outward = np.zeros((surf0.shape[0],),dtype=int)
	for i,t in enumerate(surf0):
		v0 = pts0[t[1],:] - pts0[t[0],:]
		v0 = v0/np.linalg.norm(v0)

		v1 = pts0[t[2],:] - pts0[t[0],:]
		v1 = v1/np.linalg.norm(v1)

		n = np.cross(v0,v1)
		n = n/np.linalg.norm(n)

		dot_prod = np.dot(cog0-pts0[t[0],:],n)

		if dot_prod<0:
			is_outward[i] = 1

	if np.sum(is_outward)/surf0.shape[0]>0.7:
		print(epi_endo_CC[0]+' is the epicardium')
		print(epi_endo_CC[1]+' is the endocardium')
		epi=0
		endo=1
	else:
		print(epi_endo_CC[0]+' is the endocardium')
		print(epi_endo_CC[1]+' is the epicardium')
		endo=0
		epi=1

	print('Renaming connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[epi]+"."+f+" "+surf_folder+"/tmp/ra.epi."+f)
		os.system("mv "+surf_folder+"/tmp/"+epi_endo_CC[endo]+"."+f+" "+surf_folder+"/tmp/ra.lvendo."+f)

def mapping_surfaces(mesh,surf_folder,input_tags):
	connected_component_to_surface(surf_folder+'/tmp/myocardium.epi',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/myocardium.epi')
	connected_component_to_surface(surf_folder+'/tmp/myocardium.lvendo',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/myocardium.lvendo')
	connected_component_to_surface(surf_folder+'/tmp/myocardium.rvendo',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/myocardium.rvendo')

	connected_component_to_surface(surf_folder+'/tmp/lvepi',surf_folder+'/tmp/myocardium.rvsept.surf',surf_folder+'/tmp/lvepi')
	connected_component_to_surface(surf_folder+'/tmp/myocardium.rvsept',surf_folder+'/tmp/myocardium.rvsept.surf',surf_folder+'/tmp/myocardium.rvsept')

	surf2vtk(mesh,surf_folder+'/tmp/myocardium.epi'+'.surf',surf_folder+'/tmp/myocardium.epi'+'.vtk')

def mapping_surfaces_la(mesh,surf_folder,input_tags):
	connected_component_to_surface(surf_folder+'/tmp/la.epi',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/la.epi')
	connected_component_to_surface(surf_folder+'/tmp/la.lvendo',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/la.lvendo')

	surf2vtk(mesh,surf_folder+'/tmp/la.epi'+'.surf',surf_folder+'/tmp/la.epi'+'.vtk')
	surf2vtk(mesh,surf_folder+'/tmp/la.lvendo'+'.surf',surf_folder+'/tmp/la.lvendo'+'.vtk')
	


def mapping_surfaces_ra(mesh,surf_folder,input_tags):
	connected_component_to_surface(surf_folder+'/tmp/ra.epi',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/ra.epi')
	connected_component_to_surface(surf_folder+'/tmp/ra.lvendo',surf_folder+'/tmp/epi_endo.surf',surf_folder+'/tmp/ra.lvendo')

	surf2vtk(mesh,surf_folder+'/tmp/ra.epi'+'.surf',surf_folder+'/tmp/ra.epi'+'.vtk')
	surf2vtk(mesh,surf_folder+'/tmp/ra.lvendo'+'.surf',surf_folder+'/tmp/ra.lvendo'+'.vtk')
	
def picking_apex(segmentation,mesh,surf_folder,seg_tags):
	tags_list_lv = extract_tags(seg_tags,["LV_BP"])
	tags_list_la = extract_tags(seg_tags,["LA_BP"])
	lv_cavity = str(tags_list_lv[0])
	base = str(tags_list_la[0])

	pts = mesh+'.pts'
	vtx = surf_folder+'/tmp/myocardium.epi.vtx'

	os.system("pickapex "+segmentation+" "+lv_cavity+" "+base+" "+pts+" "+vtx+" > "+surf_folder+"/myocardium.apex.vtx")

def meshtool_extract_biv(mesh,surf_folder,input_tags):
	tags_list_lv_string = get_tags_from_setup(input_tags, ["LV"])
	tags_list_rv_string = get_tags_from_setup(input_tags, ["RV"])

	os.system("meshtool extract mesh -msh="+mesh+" -submsh="+surf_folder+"/BiV/BiV -tags="+tags_list_lv_string+","+tags_list_rv_string)

def meshtool_map_vtx(surf_folder):
	os.system("meshtool map -submsh="+surf_folder+"/BiV/BiV"
				  			+" -files="+surf_folder+"/myocardium.apex.vtx,"
				  					   +surf_folder+"/tmp/myocardium.base.surf.vtx,"
				  					   +surf_folder+"/tmp/myocardium.epi.surf.vtx,"
				  					   +surf_folder+"/tmp/myocardium.lvendo.surf.vtx,"
				  					   +surf_folder+"/tmp/myocardium.rvendo.surf.vtx,"
				  					   +surf_folder+"/tmp/myocardium.rvendo_nosept.surf.vtx,"
				  					   +surf_folder+"/tmp/myocardium.rvsept.surf.vtx"
				  			+" -outdir="+surf_folder+"/BiV")

	os.system("meshtool map -submsh="+surf_folder+"/BiV/BiV"
				  			+" -files="+surf_folder+"/tmp/myocardium.base.surf,"
				  					   +surf_folder+"/tmp/myocardium.epi.surf,"
				  					   +surf_folder+"/tmp/myocardium.lvendo.surf,"
				  					   +surf_folder+"/tmp/myocardium.rvendo.surf,"
				  					   +surf_folder+"/tmp/myocardium.rvendo_nosept.surf,"
				  					   +surf_folder+"/tmp/myocardium.rvsept.surf" 
				  			+" -outdir="+surf_folder+"/BiV")

def renaming_myo_files(surf_folder):
	os.system("mv "+surf_folder+"/BiV/myocardium.base.surf "+surf_folder+"/BiV/BiV.base.surf")
	os.system("mv "+surf_folder+"/BiV/myocardium.apex.vtx "+surf_folder+"/BiV/BiV.apex.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.base.surf.vtx "+surf_folder+"/BiV/BiV.base.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.epi.surf "+surf_folder+"/BiV/BiV.epi.surf")
	os.system("mv "+surf_folder+"/BiV/myocardium.epi.surf.vtx "+surf_folder+"/BiV/BiV.epi.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.lvendo.surf "+surf_folder+"/BiV/BiV.lvendo.surf")
	os.system("mv "+surf_folder+"/BiV/myocardium.lvendo.surf.vtx "+surf_folder+"/BiV/BiV.lvendo.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvendo.surf "+surf_folder+"/BiV/BiV.rvendo.surf")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvendo.surf.vtx "+surf_folder+"/BiV/BiV.rvendo.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvendo_nosept.surf.vtx "+surf_folder+"/BiV/BiV.rvendo_nosept.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvendo_nosept.surf "+surf_folder+"/BiV/BiV.rvendo_nosept.surf")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvsept.surf.vtx "+surf_folder+"/BiV/BiV.rvsept.surf.vtx")
	os.system("mv "+surf_folder+"/BiV/myocardium.rvsept.surf "+surf_folder+"/BiV/BiV.rvsept.surf")

def meshtool_extract_la_for_UVCs(mesh,surf_folder,input_tags):
	tags_list_la = extract_tags(input_tags,["LA"])
	tags_list_la_string = str(tags_list_la[0])
	os.system("meshtool extract mesh -msh="+mesh+" -submsh="+surf_folder+"/la/la -tags="+tags_list_la_string)

def meshtool_map_vtx_la(surf_folder):
	os.system("meshtool map -submsh="+surf_folder+"/la/la"
							+" -files="+surf_folder+"/tmp/la.apex.vtx,"
									   +surf_folder+"/tmp/la.base.surf.vtx,"
									   +surf_folder+"/tmp/la.epi.vtx,"
									   +surf_folder+"/tmp/la.lvendo.vtx"
							+" -outdir="+surf_folder+"/la")

	os.system("meshtool convert -imsh="+surf_folder+"/la/la -omsh="+surf_folder+"/la/la -ofmt=vtk_bin")

def meshtool_extract_ra_for_UVCs(mesh,surf_folder,input_tags):
	tags_list_ra = extract_tags(input_tags,["RA"])
	tags_list_ra_string = str(tags_list_ra[0])
	os.system("meshtool extract mesh -msh="+mesh+" -submsh="+surf_folder+"/ra/ra -tags="+tags_list_ra_string)

def meshtool_map_vtx_ra(surf_folder):
	os.system("meshtool map -submsh="+surf_folder+"/ra/ra" 
						  " -files="+surf_folder+"/tmp/ra.apex.vtx,"
						  			+surf_folder+"/tmp/ra.base.surf.vtx,"
						  			+surf_folder+"/tmp/ra.epi.vtx,"
						  			+surf_folder+"/tmp/ra.lvendo.vtx"
						  " -outdir="+surf_folder+"/ra")

	os.system("meshtool convert -imsh="+surf_folder+"/ra/ra -omsh="+surf_folder+"/ra/ra -ofmt=vtk_bin")


def meshtool_extract_peri(mesh,presimFolder,input_tags):

	tags_list_peri_string = get_tags_from_setup(input_tags, ["LV","RV","LA","RA","BB","AV"])
	tags_list_not_peri_string = get_tags_from_setup(input_tags, ["Ao","PArt",
											 	  "MV","TV","AV","PV",
											 	  "LSPV","LIPV","RSPV","RIPV",
											 	  "LAA","SVC","IVC",
											 	  "LAA_ring","SVC_ring","IVC_ring",
											 	  "LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/peri_surface -ofmt=vtk -op="+tags_list_peri_string+"-"+tags_list_not_peri_string)
	os.system("meshtool extract unreachable -msh="+presimFolder+"/peri_surface.surfmesh -ifmt=vtk -ofmt=vtk -submsh="+presimFolder+"/peri_surface_CC")

def meshtool_extract_epi_endo_surfs(mesh,presimFolder,input_tags):
	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"surfaces_simulation/surface_heart -ofmt=carp_txt")
	os.system("meshtool extract unreachable -msh="+presimFolder+"surfaces_simulation/surface_heart.surfmesh -submsh="+presimFolder+"surfaces_simulation/surface_heart_CC -ofmt=carp_txt")

	tmp_files = os.listdir(presimFolder+"/surfaces_simulation/")
	surf_heart_CC = []
	i = 0
	isfile=True
	while isfile:
		if "surface_heart_CC.part"+str(i)+".elem" in tmp_files:
			surf_heart_CC.append("surface_heart_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the 5 biggest...")
	if len(surf_heart_CC)>5:
		CC_size = np.zeros((len(surf_heart_CC),),dtype=int)
		for i,CC in enumerate(surf_heart_CC):
			surf = read_elem(presimFolder+"/surfaces_simulation/"+CC+".elem",el_type="Tr",tags=False)
			CC_size[i] = surf.shape[0]

		surf_heart_CC_old = copy.deepcopy(surf_heart_CC)
		sorted_size = np.argsort(CC_size)
		surf_heart_CC[0] = surf_heart_CC_old[sorted_size[-1]]
		surf_heart_CC[1] = surf_heart_CC_old[sorted_size[-2]]
		surf_heart_CC[2] = surf_heart_CC_old[sorted_size[-3]]
		surf_heart_CC[3] = surf_heart_CC_old[sorted_size[-4]]
		surf_heart_CC[4] = surf_heart_CC_old[sorted_size[-5]]

		for i in range(len(surf_heart_CC)-5):
			print("Removing extraneous surfaces...")
			os.system("rm "+presimFolder+"/surfaces_simulation/"+surf_heart_CC_old[sorted_size[i]]+".*")

	surf0 = read_elem(presimFolder+"/surfaces_simulation/"+surf_heart_CC[0]+".elem",el_type="Tr",tags=False)
	surf1 = read_elem(presimFolder+"/surfaces_simulation/"+surf_heart_CC[1]+".elem",el_type="Tr",tags=False)
	surf2 = read_elem(presimFolder+"/surfaces_simulation/"+surf_heart_CC[2]+".elem",el_type="Tr",tags=False)
	surf3 = read_elem(presimFolder+"/surfaces_simulation/"+surf_heart_CC[3]+".elem",el_type="Tr",tags=False)
	surf4 = read_elem(presimFolder+"/surfaces_simulation/"+surf_heart_CC[4]+".elem",el_type="Tr",tags=False)
	
	surfs = [surf0,surf1,surf2,surf3,surf4]

	pts0 = read_pts(presimFolder+"/surfaces_simulation/"+surf_heart_CC[0]+".pts")
	pts1 = read_pts(presimFolder+"/surfaces_simulation/"+surf_heart_CC[1]+".pts")
	pts2 = read_pts(presimFolder+"/surfaces_simulation/"+surf_heart_CC[2]+".pts")
	pts3 = read_pts(presimFolder+"/surfaces_simulation/"+surf_heart_CC[3]+".pts")
	pts4 = read_pts(presimFolder+"/surfaces_simulation/"+surf_heart_CC[4]+".pts")

	pts = [pts0,pts1,pts2,pts3,pts4]

	# Find CoGs of surfaces
	cog0 = find_cog_surf(pts[0])
	cog1 = find_cog_surf(pts[1])
	cog2 = find_cog_surf(pts[2])
	cog3 = find_cog_surf(pts[3])
	cog4 = find_cog_surf(pts[4])

	surf_cogs = [cog0,cog1,cog2,cog3,cog4]

	print("Finding the centre of gravity of each blood pool...")
	mesh_pts = read_pts(mesh+".pts")
	mesh_elem = read_elem(mesh+".elem",el_type="Tt",tags=True)

	lv_tag = extract_tags(input_tags,["LV"])
	cog_lv = find_cog_vol(mesh_pts,mesh_elem,lv_tag[0])
	print("LV centre of gravity found.")

	rv_tag = extract_tags(input_tags,["RV"])
	cog_rv = find_cog_vol(mesh_pts,mesh_elem,rv_tag[0])
	print("RV centre of gravity found.")

	la_tag = extract_tags(input_tags,["LA"])
	cog_la = find_cog_vol(mesh_pts,mesh_elem,la_tag[0])
	print("LA centre of gravity found.")

	ra_tag = extract_tags(input_tags,["RA"])
	cog_ra = find_cog_vol(mesh_pts,mesh_elem,ra_tag[0])
	print("RA centre of gravity found.")

	print("Searching for the epicardium by checking the direction of surface normals...")
	surf_orientated_out = False
	i = 0
	while surf_orientated_out == False:
		surf_orientated_out = query_outwards_surf(surfs[i],pts[i],surf_cogs[i])
		if surf_orientated_out == True:
			print("		"+surf_heart_CC[i]+' is the epicardium')
			ID_epi = i
		i += 1

	remaining_surf_list = list(range(5))
	remaining_surf_list.remove(ID_epi)

	surf_dist = []
	print("Finding distance between surf CoGs and LV CoG...")
	for i in remaining_surf_list:
		dist = calculate_dist(surf_cogs[i],cog_lv)
		surf_dist.append(dist)

	idx_ID_lv = np.argmin(surf_dist)
	ID_lv = remaining_surf_list[idx_ID_lv]
	print("		"+surf_heart_CC[ID_lv]+' is the LV endocardium')
	remaining_surf_list.remove(ID_lv)

	surf_dist = []
	print("Finding distance between surf CoGs and RV CoG...")
	for i in remaining_surf_list:
		dist = calculate_dist(surf_cogs[i],cog_rv)
		surf_dist.append(dist)

	idx_ID_rv = np.argmin(surf_dist)
	ID_rv = remaining_surf_list[idx_ID_rv]
	print("		"+surf_heart_CC[ID_rv]+' is the RV endocardium')
	remaining_surf_list.remove(ID_rv)

	surf_dist = []
	print("Finding distance between surf CoGs and LA CoG...")
	for i in remaining_surf_list:
		dist = calculate_dist(surf_cogs[i],cog_la)
		surf_dist.append(dist)

	idx_ID_la = np.argmin(surf_dist)
	ID_la = remaining_surf_list[idx_ID_la]
	print("		"+surf_heart_CC[ID_la]+' is the LA endocardium')
	remaining_surf_list.remove(ID_la)

	print("And therefore...")
	ID_ra = remaining_surf_list[0]
	print("		"+surf_heart_CC[ID_ra]+' is the RA endocardium')

	connected_component_to_surface(presimFolder+"/surfaces_simulation/"+surf_heart_CC[ID_epi],
								   presimFolder+"/surfaces_simulation/surface_heart.surf",
								   presimFolder+"/surfaces_simulation/epicardium")
	surf2vtk(mesh,
			 presimFolder+"/surfaces_simulation/epicardium.surf",
			 presimFolder+"/surfaces_simulation/epicardium.surf.vtk")

	connected_component_to_surface(presimFolder+"/surfaces_simulation/"+surf_heart_CC[ID_lv],
								   presimFolder+"/surfaces_simulation/surface_heart.surf",
								   presimFolder+"/surfaces_simulation/LV_endo")
	surf2vtk(mesh,
			 presimFolder+"/surfaces_simulation/LV_endo.surf",
			 presimFolder+"/surfaces_simulation/LV_endo.surf.vtk")

	connected_component_to_surface(presimFolder+"/surfaces_simulation/"+surf_heart_CC[ID_rv],
								   presimFolder+"/surfaces_simulation/surface_heart.surf",
								   presimFolder+"/surfaces_simulation/RV_endo")
	surf2vtk(mesh,
			 presimFolder+"/surfaces_simulation/RV_endo.surf",
			 presimFolder+"/surfaces_simulation/RV_endo.surf.vtk")

	connected_component_to_surface(presimFolder+"/surfaces_simulation/"+surf_heart_CC[ID_la],
								   presimFolder+"/surfaces_simulation/surface_heart.surf",
								   presimFolder+"/surfaces_simulation/LA_endo")
	surf2vtk(mesh,
			 presimFolder+"/surfaces_simulation/LA_endo.surf",
			 presimFolder+"/surfaces_simulation/LA_endo.surf.vtk")

	connected_component_to_surface(presimFolder+"/surfaces_simulation/"+surf_heart_CC[ID_ra],
								   presimFolder+"/surfaces_simulation/surface_heart.surf",
								   presimFolder+"/surfaces_simulation/RA_endo")
	surf2vtk(mesh,
			 presimFolder+"/surfaces_simulation/RA_endo.surf",
			 presimFolder+"/surfaces_simulation/RA_endo.surf.vtk")

def meshtool_extract_rings(mesh,presimFolder,input_tags):
	print("Extracting the RSPV ring and RIPV ring for use as boundary conditions...")

	tags_list_rpv_rings_string = get_tags_from_setup(input_tags, ["RSPV_ring","RIPV_ring"])
	tags_list_other_string = get_tags_from_setup(input_tags,["LV","RV","LA","RA",
											   "Ao","PArt",
											   "MV","TV","AV","PV",
											   "LSPV","LIPV","RSPV","RIPV",
											   "LAA","SVC","IVC",
											   "LAA_ring","SVC_ring","IVC_ring",
											   "LSPV_ring","LIPV_ring"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/surfaces_simulation/surfaces_rings/RPVs -ofmt=vtk -op="+tags_list_rpv_rings_string+"-"+tags_list_other_string)

	print("Extracting the SVC ring for use as a boundary condition...")
	tags_list_svc_ring_string = get_tags_from_setup(input_tags,["SVC_ring"])
	tags_list_other_string = get_tags_from_setup(input_tags,["LV","RV","LA","RA",
											   "Ao","PArt",
											   "MV","TV","AV","PV",
											   "LSPV","LIPV","RSPV","RIPV",
											   "LAA","SVC","IVC",
											   "LAA_ring","IVC_ring",
											   "LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring"])

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/surfaces_simulation/surfaces_rings/SVC -ofmt=vtk -op="+tags_list_svc_ring_string+"-"+tags_list_other_string)

	print("Converting necessary surfs to vtx files...")
	rpvs_surface = np.loadtxt(presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rpvs_vtx = surf2vtx(rpvs_surface)
	write_vtx(rpvs_vtx,presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf.vtx")

	svc_surface = np.loadtxt(presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	svc_vtx = surf2vtx(svc_surface)
	write_vtx(svc_vtx,presimFolder+"/surfaces_simulation/surfaces_rings/SVC.surf.vtx")

def combine_elem_dats(heartFolder,presimFolder):
	la_map_dat=os.path.join(heartFolder, "/surfaces_uvc_LA/la/uvc/map_rotational_z.dat")
	ra_map_dat=os.path.join(heartFolder, "/surfaces_uvc_RA/ra/uvc/map_rotational_z.dat")

	mycp(la_map_dat,os.path.join(presimFolder, "map_rotational_z_la.dat"))
	mycp(ra_map_dat,os.path.join(presimFolder, "map_rotational_z_ra.dat"))

	os.system("meshtool interpolate node2elem "
						+"-omsh="+heartFolder+"/surfaces_uvc_LA/la/la "
						+"-idat="+presimFolder+"/map_rotational_z_la.dat "
						+"-odat="+presimFolder+"/map_rotational_z_la_e.dat")

	os.system("meshtool interpolate node2elem "
						+"-omsh="+heartFolder+"/surfaces_uvc_RA/ra/ra "
						+"-idat="+presimFolder+"/map_rotational_z_ra.dat "
						+"-odat="+presimFolder+"/map_rotational_z_ra_e.dat")

	os.system("meshtool insert data "
						+"-msh="+presimFolder+"/myocardium_AV_FEC_BB "
						+"-submsh="+heartFolder+"/surfaces_uvc_LA/la/la "
						+"-submsh_data="+presimFolder+"/map_rotational_z_la_e.dat "
						+"-odat="+presimFolder+"/elem_dat_UVC_ek_inc_la.dat "
						+"-mode=1")

	os.system("meshtool insert data "
						+"-msh="+presimFolder+"/myocardium_AV_FEC_BB "
						+"-submsh="+heartFolder+"/surfaces_uvc_RA/ra/ra "
						+"-submsh_data="+presimFolder+"/map_rotational_z_ra_e.dat "
						+"-odat="+presimFolder+"/elem_dat_UVC_ek_inc_ra.dat "
						+"-mode=1")

	combine_rot_coords(presimFolder)

	os.system("GlVTKConvert "
				+"-m "+presimFolder+"/myocardium_AV_FEC_BB "
				+"-e "+presimFolder+"/elem_dat_UVC_ek_combined.dat "
				+"-F bin "
				+"-o "+presimFolder+"/myocardium_AV_FEC_BB_elem_dat_UVC_combined "
				+"--trim-names")
