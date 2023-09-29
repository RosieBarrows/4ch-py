import os
import sys
import numpy as np
import copy
import meshio

from common_4ch.file_utils import *
from common_4ch.distance_utils import *
from common_4ch.mesh_utils import *

def extract_tags(tags_setup,
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

def meshtool_extract_biatrial(meshname,
							  output_folder,
							  tags_setup,
							  overwrite=False):
	
	do = True
	if os.path.exists(output_folder+'/biatrial/biatrial.pts') and not overwrite:
		do = False
		print('Mesh already found. Not overwriting.')

	if do:	
		tags_list = extract_tags(tags_setup,["LA","RA"])
		tags_list = [str(t) for t in tags_list]
		tags_list_string = ",".join(tags_list)	

		print("----------------------------")
		print("Extracting biatrial mesh...")
		print("----------------------------")	

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
		tags_list = extract_tags(tags_setup,["LA"])
		tags_list = [str(t) for t in tags_list]
		tags_list_string = ",".join(tags_list)
		
		print("--------------------------------------")
		print("Extracting LA mesh from biatrial...")
		print("--------------------------------------")
		
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
		tags_list = extract_tags(tags_setup,["RA"])
		tags_list = [str(t) for t in tags_list]
		tags_list_string = ",".join(tags_list)
		
		print("--------------------------------------")
		print("Extracting RA mesh from biatrial...")
		print("--------------------------------------")
		
		cmd = "meshtool extract mesh -msh="+output_folder+"/biatrial/biatrial -submsh="+output_folder+"/ra/ra -tags="+tags_list_string
		os.system(cmd)

def meshtool_extract_surfaces(meshname,
		 					  output_folder,
		 					  tags_setup,
		 					  export_sup_inf=False,
		 					  rm_vv_from_aa=False,
		 					  surface="epi"):

	tag_la = extract_tags(tags_setup,["LA"])
	tag_la = [str(t) for t in tag_la]
	tag_la_string = ",".join(tag_la)

	tag_ra = extract_tags(tags_setup,["RA"])
	tag_ra = [str(t) for t in tag_ra]
	tag_ra_string = ",".join(tag_ra)

	if not rm_vv_from_aa:
		tag_la_rm = extract_tags(tags_setup,["mitral","RSPV","RIPV","LSPV","LIPV","LAA","PV_planes"])
	else:
		tag_la_rm = extract_tags(tags_setup,["mitral","RSPV","RIPV","LSPV","LIPV","LAA","PV_planes","LV"])
	tag_la_rm = [str(t) for t in tag_la_rm]
	tag_la_rm_string = ",".join(tag_la_rm)

	if not rm_vv_from_aa:
		tag_ra_rm = extract_tags(tags_setup,["SVC","IVC","tricuspid","VC_planes"])
	else:
		tag_ra_rm = extract_tags(tags_setup,["SVC","IVC","tricuspid","VC_planes","RV"])
	tag_ra_rm = [str(t) for t in tag_ra_rm]
	tag_ra_rm_string = ",".join(tag_ra_rm)

	tag_mv = extract_tags(tags_setup,["mitral"])
	tag_mv = [str(t) for t in tag_mv]
	tag_mv_string = ",".join(tag_mv)

	tag_tv = extract_tags(tags_setup,["tricuspid"])
	tag_tv = [str(t) for t in tag_tv]
	tag_tv_string = ",".join(tag_tv)

	tag_rpv = extract_tags(tags_setup,["RSPV","RIPV"])
	tag_rpv = [str(t) for t in tag_rpv]
	tag_rpv_string = ",".join(tag_rpv)

	tag_svc = extract_tags(tags_setup,["SVC"])
	tag_svc = [str(t) for t in tag_svc]
	tag_svc_string = ",".join(tag_svc)

	tag_lpv = extract_tags(tags_setup,["LSPV","LIPV"])
	tag_lpv = [str(t) for t in tag_lpv]
	tag_lpv_string = ",".join(tag_lpv)

	tag_ivc = extract_tags(tags_setup,["IVC"])
	tag_ivc = [str(t) for t in tag_ivc]
	tag_ivc_string = ",".join(tag_ivc)

	tag_pv_planes = extract_tags(tags_setup,["PV_planes"])
	tag_pv_planes = [str(t) for t in tag_pv_planes]
	tag_pv_planes_string = ",".join(tag_pv_planes)

	tag_vc_planes = extract_tags(tags_setup,["VC_planes"])
	tag_vc_planes = [str(t) for t in tag_vc_planes]
	tag_vc_planes_string = ",".join(tag_vc_planes)

	tag_lv = extract_tags(tags_setup,["LV"])
	tag_lv = [str(t) for t in tag_lv]
	tag_lv_string = ",".join(tag_lv)

	tag_rv = extract_tags(tags_setup,["RV"])
	tag_rv = [str(t) for t in tag_rv]
	tag_rv_string = ",".join(tag_rv)

	os.system("mkdir -p "+output_folder+"/tmp")

	print("--------------------------")
	print("Extracting mitral ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/mitral -op="+tag_la_string+":"+tag_mv_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting right PV ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/rpv -op="+tag_la_string+":"+tag_rpv_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting left PV ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/lpv -op="+tag_la_string+":"+tag_lpv_string+" -ofmt=vtk"
	os.system(cmd)

	if export_sup_inf:

		tag_rspv = extract_tags(tags_setup,["RSPV"])
		tag_rspv = [str(t) for t in tag_rspv]
		tag_rspv_string = ",".join(tag_rspv)	

		tag_ripv = extract_tags(tags_setup,["RIPV"])
		tag_ripv = [str(t) for t in tag_ripv]
		tag_ripv_string = ",".join(tag_ripv)	

		tag_lspv = extract_tags(tags_setup,["LSPV"])
		tag_lspv = [str(t) for t in tag_lspv]
		tag_lspv_string = ",".join(tag_lspv)	

		tag_lipv = extract_tags(tags_setup,["LIPV"])
		tag_lipv = [str(t) for t in tag_lipv]
		tag_lipv_string = ",".join(tag_lipv)

		tag_laa = extract_tags(tags_setup,["LAA"])
		tag_laa = [str(t) for t in tag_laa]
		tag_laa_string = ",".join(tag_laa)

		print("------------------------------------------")
		print("Extracting LAA ring...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/laa -op="+tag_la_string+":"+tag_laa_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting right inferior PV ring...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ripv -op="+tag_la_string+":"+tag_ripv_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting left inferior PV ring...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/lipv -op="+tag_la_string+":"+tag_lipv_string+" -ofmt=vtk"
		os.system(cmd)

		print("------------------------------------------")
		print("Extracting right superior PV ring...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/rspv -op="+tag_la_string+":"+tag_rspv_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting left superior PV ring...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/lspv -op="+tag_la_string+":"+tag_lspv_string+" -ofmt=vtk"
		os.system(cmd)


	print("--------------------------")
	print("Extracting PV planes...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/pv_planes -op="+tag_la_string+":"+tag_pv_planes_string+" -ofmt=vtk"
	os.system(cmd)

	if export_sup_inf and surface=="endo":

		if "RSPV_vp" not in tags_setup:
			raise Exception("If you want to compute the landmarks on the endocardium, you need to provide separate tags for the valve planes")

		tag_rspv_vp = extract_tags(tags_setup,["RSPV_vp"])
		tag_rspv_vp = [str(t) for t in tag_rspv_vp]
		tag_rspv_vp_string = ",".join(tag_rspv_vp)	

		tag_ripv_vp = extract_tags(tags_setup,["RIPV_vp"])
		tag_ripv_vp = [str(t) for t in tag_ripv_vp]
		tag_ripv_vp_string = ",".join(tag_ripv_vp)	

		tag_lspv_vp = extract_tags(tags_setup,["LSPV_vp"])
		tag_lspv_vp = [str(t) for t in tag_lspv_vp]
		tag_lspv_vp_string = ",".join(tag_lspv_vp)	

		tag_lipv_vp = extract_tags(tags_setup,["LIPV_vp"])
		tag_lipv_vp = [str(t) for t in tag_lipv_vp]
		tag_lipv_vp_string = ",".join(tag_lipv_vp)

		tag_laa_vp = extract_tags(tags_setup,["LAA_vp"])
		tag_laa_vp = [str(t) for t in tag_laa_vp]
		tag_laa_vp_string = ",".join(tag_laa_vp)

		print("------------------------------------------")
		print("Extracting LAA valve plane...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/laa_vp -op="+tag_la_string+":"+tag_laa_vp_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting right inferior PV valve plane...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ripv_vp -op="+tag_la_string+":"+tag_ripv_vp_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting left inferior PV valve plane...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/lipv_vp -op="+tag_la_string+":"+tag_lipv_vp_string+" -ofmt=vtk"
		os.system(cmd)

		print("------------------------------------------")
		print("Extracting right superior PV valve plane...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/rspv_vp -op="+tag_la_string+":"+tag_rspv_vp_string+" -ofmt=vtk"
		os.system(cmd)	

		print("------------------------------------------")
		print("Extracting left superior PV valve plane...")
		print("------------------------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/lspv_vp -op="+tag_la_string+":"+tag_lspv_vp_string+" -ofmt=vtk"
		os.system(cmd)

	print("--------------------------")
	print("Extracting LA surface...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/la -op="+tag_la_string+"-"+tag_la_rm_string+" -ofmt=vtk"
	os.system(cmd)

	print("------------------------------------")
	print("Extracting LA and LV intersection...")
	print("------------------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/la_lv -op="+tag_la_string+":"+tag_lv_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting tricuspid ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/tricuspid -op="+tag_ra_string+":"+tag_tv_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting SVC ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/svc -op="+tag_ra_string+":"+tag_svc_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting IVC ring...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ivc -op="+tag_ra_string+":"+tag_ivc_string+" -ofmt=vtk"
	os.system(cmd)

	if surface=="endo":

		if "SVC_vp" not in tags_setup:
			raise Exception("If you want to compute the landmarks on the endocardium, you need to provide separate tags for the valve planes")

		tag_svc_vp = extract_tags(tags_setup,["SVC_vp"])
		tag_svc_vp = [str(t) for t in tag_svc_vp]
		tag_svc_vp_string = ",".join(tag_svc_vp)	

		tag_ivc_vp = extract_tags(tags_setup,["IVC_vp"])
		tag_ivc_vp = [str(t) for t in tag_ivc_vp]
		tag_ivc_vp_string = ",".join(tag_ivc_vp)	

		print("--------------------------")
		print("Extracting SVC valve plane...")
		print("--------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/svc_vp -op="+tag_ra_string+":"+tag_svc_vp_string+" -ofmt=vtk"
		os.system(cmd)

		print("--------------------------")
		print("Extracting IVC valve plane...")
		print("--------------------------")
		cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ivc_vp -op="+tag_ra_string+":"+tag_ivc_vp_string+" -ofmt=vtk"
		os.system(cmd)

	print("--------------------------")
	print("Extracting PV planes...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/pv_planes -op="+tag_la_string+":"+tag_pv_planes_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting SVC IVC planes...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/vc_planes -op="+tag_ra_string+":"+tag_vc_planes_string+" -ofmt=vtk"
	os.system(cmd)

	print("------------------------------------")
	print("Extracting RA and RV intersection...")
	print("------------------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ra_rv -op="+tag_ra_string+":"+tag_rv_string+" -ofmt=vtk"
	os.system(cmd)

	print("--------------------------")
	print("Extracting RA surface...")
	print("--------------------------")
	cmd = "meshtool extract surface -msh="+meshname+" -surf="+output_folder+"/tmp/ra -op="+tag_ra_string+"-"+tag_ra_rm_string+" -ofmt=vtk"
	os.system(cmd)

	print("-----------------------------------------------")
	print("Extracting LA surface connected components...")
	print("-----------------------------------------------")
	cmd = "meshtool extract unreachable -msh="+output_folder+"/tmp/la.surfmesh.vtk -submsh="+output_folder+"/tmp/la_cc -ofmt=carp_txt"
	os.system(cmd)

	print("-----------------------------------------------")
	print("Extracting RA surface connected components...")
	print("-----------------------------------------------")
	cmd = "meshtool extract unreachable -msh="+output_folder+"/tmp/ra.surfmesh.vtk -submsh="+output_folder+"/tmp/ra_cc -ofmt=carp_txt"
	os.system(cmd)

	tmp_files = os.listdir(output_folder+"/tmp")
	la_cc = []
	i = 0
	isfile=True
	while isfile:
		if "la_cc.part"+str(i)+".elem" in tmp_files:
			la_cc.append("la_cc.part"+str(i))
		else: 
			isfile = False
		i += 1

	ra_cc = []
	i = 0
	isfile=True
	while isfile:
		if "ra_cc.part"+str(i)+".elem" in tmp_files:
			ra_cc.append("ra_cc.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only the two biggest...")
	if len(la_cc)>2:
		cc_size = np.zeros((len(la_cc),),dtype=int)
		for i,cc in enumerate(la_cc):
			surf = read_elem(output_folder+"/tmp/"+cc+".elem",el_type="Tr",tags=False)
			cc_size[i] = surf.shape[0]

		la_cc_old = copy.deepcopy(la_cc)
		sorted_size = np.argsort(cc_size)
		la_cc[0] = la_cc_old[sorted_size[-1]]
		la_cc[1] = la_cc_old[sorted_size[-2]]

		for i in range(len(la_cc)-2):
			os.system("rm "+output_folder+"/tmp/"+la_cc_old[sorted_size[i]]+".*")

	if len(ra_cc)>2:
		cc_size = np.zeros((len(ra_cc),),dtype=int)
		for i,cc in enumerate(ra_cc):
			surf = read_elem(output_folder+"/tmp/"+cc+".elem",el_type="Tr",tags=False)
			cc_size[i] = surf.shape[0]

		ra_cc_old = copy.deepcopy(ra_cc)
		sorted_size = np.argsort(cc_size)
		ra_cc[0] = ra_cc_old[sorted_size[-1]]
		ra_cc[1] = ra_cc_old[sorted_size[-2]]

		for i in range(len(ra_cc)-2):
			os.system("rm "+output_folder+"/tmp/"+ra_cc_old[sorted_size[i]]+".*")

	pts0 = read_pts(output_folder+"/tmp/"+la_cc[0]+".pts")
	surf0 = read_elem(output_folder+"/tmp/"+la_cc[0]+".elem",el_type="Tr",tags=False)

	pts1 = read_pts(output_folder+"/tmp/"+la_cc[1]+".pts")
	surf1 = read_elem(output_folder+"/tmp/"+la_cc[1]+".elem",el_type="Tr",tags=False)

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
		print(la_cc[0]+' is the epicardium')
		print(la_cc[1]+' is the endocardium')
		endo = 0
		epi = 1
	else:
		print(la_cc[1]+' is the epicardium')
		print(la_cc[0]+' is the endocardium')	
		endo = 1
		epi = 0

	print('Renaming LA connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+output_folder+"/tmp/"+la_cc[endo]+"."+f+" "+output_folder+"/tmp/la_endo."+f)
		os.system("mv "+output_folder+"/tmp/"+la_cc[epi]+"."+f+" "+output_folder+"/tmp/la_epi."+f)

	pts0 = read_pts(output_folder+"/tmp/"+ra_cc[0]+".pts")
	surf0 = read_elem(output_folder+"/tmp/"+ra_cc[0]+".elem",el_type="Tr",tags=False)

	pts1 = read_pts(output_folder+"/tmp/"+ra_cc[1]+".pts")
	surf1 = read_elem(output_folder+"/tmp/"+ra_cc[1]+".elem",el_type="Tr",tags=False)

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
		print(ra_cc[0]+' is the epicardium')
		print(ra_cc[1]+' is the endocardium')
		endo = 0
		epi = 1
	else:
		print(ra_cc[1]+' is the epicardium')
		print(ra_cc[0]+' is the endocardium')	
		endo = 1
		epi = 0

	print('Renaming RA connected components...')
	formats = ["nod","eidx","elem","lon","pts"]
	for f in formats:
		os.system("mv "+output_folder+"/tmp/"+ra_cc[endo]+"."+f+" "+output_folder+"/tmp/ra_endo."+f)
		os.system("mv "+output_folder+"/tmp/"+ra_cc[epi]+"."+f+" "+output_folder+"/tmp/ra_epi."+f)

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
	
	surface_list = [output_folder+"/tmp/la.surf",
					output_folder+"/tmp/mitral.surf",
					output_folder+"/tmp/ripv.surf",
					output_folder+"/tmp/rspv.surf",
					output_folder+"/tmp/lipv.surf",
					output_folder+"/tmp/lspv.surf",
					output_folder+"/tmp/laa.surf",
					output_folder+"/tmp/pv_planes.surf",
					output_folder+"/tmp/la_lv.surf"]
	if surface=="endo":
		surface_list += [output_folder+"/tmp/ripv_vp.surf",
						 output_folder+"/tmp/rspv_vp.surf",
						 output_folder+"/tmp/lipv_vp.surf",
						 output_folder+"/tmp/lspv_vp.surf",
						 output_folder+"/tmp/laa_vp.surf"]
	surface_list_string = ','.join(surface_list)

	surface_bia_list = [output_folder+"/biatrial/la.surf",
						output_folder+"/biatrial/mitral.surf",
						output_folder+"/biatrial/ripv.surf",
						output_folder+"/biatrial/rspv.surf",
						output_folder+"/biatrial/lipv.surf",
						output_folder+"/biatrial/lspv.surf",
						output_folder+"/biatrial/laa.surf",
						output_folder+"/biatrial/pv_planes.surf",
						output_folder+"/biatrial/la_lv.surf"]
	if surface=="endo":
		surface_bia_list += [output_folder+"/biatrial/ripv_vp.surf",
						 	 output_folder+"/biatrial/rspv_vp.surf",
						 	 output_folder+"/biatrial/lipv_vp.surf",
						 	 output_folder+"/biatrial/lspv_vp.surf",
						 	 output_folder+"/biatrial/laa_vp.surf"]
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

	tets = read_elem(output_folder+"/ra/ra.elem",el_type="Tt",tags=False)
	pts = read_pts(output_folder+"/ra/ra.pts")
	
	surface_list = [output_folder+"/tmp/ra.surf",
					output_folder+"/tmp/tricuspid.surf",
					output_folder+"/tmp/svc.surf",
					output_folder+"/tmp/ivc.surf",
					output_folder+"/tmp/vc_planes.surf",
					output_folder+"/tmp/ra_rv.surf"]

	if surface=="endo":
		surface_list += [output_folder+"/tmp/svc_vp.surf",
						 output_folder+"/tmp/ivc_vp.surf"]
	surface_list_string = ','.join(surface_list)

	surface_bia_list = [output_folder+"/biatrial/ra.surf",
						output_folder+"/biatrial/tricuspid.surf",
						output_folder+"/biatrial/svc.surf",
						output_folder+"/biatrial/ivc.surf",
						output_folder+"/biatrial/vc_planes.surf",
						output_folder+"/biatrial/ra_rv.surf"]
	if surface=="endo":
		surface_bia_list += [output_folder+"/biatrial/svc_vp.surf",
						 	 output_folder+"/biatrial/ivc_vp.surf"]
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

	ra_rv_tr = read_elem(output_folder+"/ra/ra_rv.surf",el_type="Tr",tags=False)
	ra_epi_tr = np.concatenate((ra_epi_tr,ra_rv_tr),axis=0)

	write_surf_caroline(output_folder+"/ra/ra_epi.surf",ra_epi_tr)
	write_surf_caroline(output_folder+"/ra/ra_endo.surf",ra_endo_tr)

	print('-----------------------------------------------------')
	print('Finding LA and RA landmarks...')
	print('-----------------------------------------------------')

	find_landmarks(output_folder,
				   surface=surface,
				   scale_factor=scale_factor,
				   raa_apex_file=raa_apex_file)

	print('-----------------------------------------------------')
	print('Organising folders...')
	print('-----------------------------------------------------')

	os.system("mkdir -p "+output_folder+"/LA_endo/")
	os.system("mkdir -p "+output_folder+"/LA_epi/")
	os.system("mkdir -p "+output_folder+"/RA_endo/")
	os.system("mkdir -p "+output_folder+"/RA_epi/")

	os.system("meshtool convert -imsh="+output_folder+"/tmp/la_endo -ofmt=vtk_polydata -omsh="+output_folder+"/LA_endo/LA_endo")
	os.system("meshtool convert -imsh="+output_folder+"/tmp/ra_endo -ofmt=vtk_polydata -omsh="+output_folder+"/RA_endo/RA_endo")
	os.system("cp "+output_folder+"/la/prodLaLandmarks.txt "+output_folder+"/LA_endo/")
	os.system("cp "+output_folder+"/la/prodLaRegion.txt "+output_folder+"/LA_endo/")
	os.system("cp "+output_folder+"/la/prodLaLandmarks.txt "+output_folder+"/LA_epi/")
	os.system("cp "+output_folder+"/la/prodLaRegion.txt "+output_folder+"/LA_epi/")

	os.system("meshtool convert -imsh="+output_folder+"/tmp/la_epi -ofmt=vtk_polydata -omsh="+output_folder+"/LA_epi/LA_epi")
	os.system("meshtool convert -imsh="+output_folder+"/tmp/ra_epi -ofmt=vtk_polydata -omsh="+output_folder+"/RA_epi/RA_epi")
	os.system("cp "+output_folder+"/ra/prodRaLandmarks.txt "+output_folder+"/RA_endo/")
	os.system("cp "+output_folder+"/ra/prodRaRegion.txt "+output_folder+"/RA_endo/")
	os.system("cp "+output_folder+"/ra/prodRaLandmarks.txt "+output_folder+"/RA_epi/")
	os.system("cp "+output_folder+"/ra/prodRaRegion.txt "+output_folder+"/RA_epi/")

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

	print("---------------------------------------------")
	print("Extracting points for SVC IVC geodesic...")
	print("---------------------------------------------")
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

	tags_list_vent = extract_tags(input_tags,["LV","RV"])
	tags_list_vent = [str(t) for t in tags_list_vent]
	tags_list_vent_string = ",".join(tags_list_vent)

	tags_list_VPs = extract_tags(input_tags,["MV","TV","AV","PV"])
	tags_list_VPs = [str(t) for t in tags_list_VPs]
	tags_list_VPs_string = ",".join(tags_list_VPs)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/myocardium.base -ofmt=vtk -op="+tags_list_vent_string+":"+tags_list_VPs_string)

def meshtool_extract_surfaces_lv_rv_epi(mesh,surf_folder,input_tags):

	tags_list_vent = extract_tags(input_tags,["LV","RV"])
	tags_list_vent = [str(t) for t in tags_list_vent]
	tags_list_vent_string = ",".join(tags_list_vent)

	tags_list_VPs = extract_tags(input_tags,["MV","TV","AV","PV","PArt"])
	tags_list_VPs = [str(t) for t in tags_list_VPs]
	tags_list_VPs_string = ",".join(tags_list_VPs)

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

	tags_list_lv = extract_tags(input_tags,["LV"])
	tags_list_lv = [str(t) for t in tags_list_lv]
	tags_list_lv_string = ",".join(tags_list_lv)

	tags_list_remove = extract_tags(input_tags,["RV","RA","PArt"])
	tags_list_remove = [str(t) for t in tags_list_remove]
	tags_list_remove_string = ",".join(tags_list_remove)

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
	tags_list_la = extract_tags(input_tags,["LA"])
	tags_list_la = [str(t) for t in tags_list_la]
	tags_list_la_string = ",".join(tags_list_la)

	tags_list_lv = extract_tags(input_tags,["LV"])
	tags_list_lv = [str(t) for t in tags_list_lv]
	tags_list_lv_string = ",".join(tags_list_lv)

	tags_list_mv = extract_tags(input_tags,["MV"])
	tags_list_mv = [str(t) for t in tags_list_mv]
	tags_list_mv_string = ",".join(tags_list_mv)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/la.base -ofmt=vtk -op="+tags_list_la_string+":"+tags_list_mv_string+","+tags_list_lv_string)

def meshtool_extract_la_surfaces(mesh,surf_folder,input_tags):
	tags_list_la = extract_tags(input_tags,["LA"])
	tags_list_la = [str(t) for t in tags_list_la]
	tags_list_la_string = ",".join(tags_list_la)

	tags_list_lv = extract_tags(input_tags,["LV"])
	tags_list_lv = [str(t) for t in tags_list_lv]
	tags_list_lv_string = ",".join(tags_list_lv)

	tags_list_VPs = extract_tags(input_tags,["MV","TV","AV","PV","LSPV","LIPV","RSPV","RIPV","LAA","SVC","IVC"])
	tags_list_VPs = [str(t) for t in tags_list_VPs]
	tags_list_VPs_string = ",".join(tags_list_VPs)

	tags_list_rings = extract_tags(input_tags,["LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring","LAA_ring","SVC_ring","IVC_ring"])
	tags_list_rings = [str(t) for t in tags_list_rings]
	tags_list_rings_string = ",".join(tags_list_rings)

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
	tags_list_ra = extract_tags(input_tags,["RA"])
	tags_list_ra = [str(t) for t in tags_list_ra]
	tags_list_ra_string = ",".join(tags_list_ra)

	tags_list_rv = extract_tags(input_tags,["RV"])
	tags_list_rv = [str(t) for t in tags_list_rv]
	tags_list_rv_string = ",".join(tags_list_rv)

	tags_list_tv = extract_tags(input_tags,["TV"])
	tags_list_tv = [str(t) for t in tags_list_tv]
	tags_list_tv_string = ",".join(tags_list_tv)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+surf_folder+"/tmp/ra.base -ofmt=vtk -op="+tags_list_ra_string+":"+tags_list_tv_string+","+tags_list_tv_string)

def meshtool_extract_ra_surfaces(mesh,surf_folder,input_tags):
	tags_list_ra = extract_tags(input_tags,["RA"])
	tags_list_ra = [str(t) for t in tags_list_ra]
	tags_list_ra_string = ",".join(tags_list_ra)

	tags_list_rv = extract_tags(input_tags,["RV"])
	tags_list_rv = [str(t) for t in tags_list_rv]
	tags_list_rv_string = ",".join(tags_list_rv)

	tags_list_VPs = extract_tags(input_tags,["MV","TV","AV","PV","LSPV","LIPV","RSPV","RIPV","LAA","SVC","IVC"])
	tags_list_VPs = [str(t) for t in tags_list_VPs]
	tags_list_VPs_string = ",".join(tags_list_VPs)

	tags_list_rings = extract_tags(input_tags,["LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring","LAA_ring","SVC_ring","IVC_ring"])
	tags_list_rings = [str(t) for t in tags_list_rings]
	tags_list_rings_string = ",".join(tags_list_rings)

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
	tags_list_lv = extract_tags(input_tags,["LV"])
	tags_list_lv = [str(t) for t in tags_list_lv]
	tags_list_lv_string = ",".join(tags_list_lv)

	tags_list_rv = extract_tags(input_tags,["RV"])
	tags_list_rv = [str(t) for t in tags_list_rv]
	tags_list_rv_string = ",".join(tags_list_rv)

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

	tags_list_peri = extract_tags(input_tags,["LV","RV","LA","RA","BB","AV_plane"])
	tags_list_peri = [str(t) for t in tags_list_peri]
	tags_list_peri_string = ",".join(tags_list_peri)

	tags_list_not_peri = extract_tags(input_tags,["Ao","PArt",
											 	  "MV","TV","AV","PV",
											 	  "LSPV","LIPV","RSPV","RIPV",
											 	  "LAA","SVC","IVC",
											 	  "LAA_ring","SVC_ring","IVC_ring",
											 	  "LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring"])
	tags_list_not_peri = [str(t) for t in tags_list_not_peri]
	tags_list_not_peri_string = ",".join(tags_list_not_peri)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/peri_surface -ofmt=carp_txt -op="+tags_list_peri_string+"-"+tags_list_not_peri_string)
	os.system("meshtool extract unreachable -msh="+presimFolder+"/peri_surface.surfmesh -ifmt=vtk -ofmt=carp_txt -submsh="+presimFolder+"/peri_surface_CC")

	tmp_files = os.listdir(presimFolder)
	peri_surface_CC = []
	i = 0
	isfile=True
	while isfile:
		if "peri_surface_CC.part"+str(i)+".elem" in tmp_files:
			peri_surface_CC.append("peri_surface_CC.part"+str(i))
		else: 
			isfile = False
		i += 1

	print("Checking connected component size and keeping only biggest...")
	if len(peri_surface_CC)>1:
		CC_size = np.zeros((len(peri_surface_CC),),dtype=int)
		for i,CC in enumerate(peri_surface_CC):
			surf = read_elem(presimFolder+CC+".elem",el_type="Tr",tags=False)
			CC_size[i] = surf.shape[0]

		peri_surface_CC_old = copy.deepcopy(peri_surface_CC)
		sorted_size = np.argsort(CC_size)
		peri_surface_CC[0] = peri_surface_CC_old[sorted_size[-1]]

		for i in range(len(peri_surface_CC)-1):
			print("Removing extraneous surfaces...")
			os.system("rm "+presimFolder+peri_surface_CC_old[sorted_size[i]]+".*")


		connected_component_to_surface(presimFolder+peri_surface_CC[0],
							   presimFolder+"/peri_surface.surf",
							   presimFolder+"/epicardium_for_sim")
	else:
		connected_component_to_surface(presimFolder+"peri_surface_CC.part0",
							   presimFolder+"/peri_surface.surf",
							   presimFolder+"/epicardium_for_sim")

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
	tags_list_rpv_rings = extract_tags(input_tags,["RSPV_ring","RIPV_ring"])
	tags_list_rpv_rings = [str(t) for t in tags_list_rpv_rings]
	tags_list_rpv_rings_string = ",".join(tags_list_rpv_rings)

	tags_list_other = extract_tags(input_tags,["LV","RV","LA","RA",
											   "Ao","PArt",
											   "MV","TV","AV","PV",
											   "LSPV","LIPV","RSPV","RIPV",
											   "LAA","SVC","IVC",
											   "LAA_ring","SVC_ring","IVC_ring",
											   "LSPV_ring","LIPV_ring"])
	tags_list_other = [str(t) for t in tags_list_other]
	tags_list_other_string = ",".join(tags_list_other)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/surfaces_simulation/surfaces_rings/RPVs -ofmt=vtk -op="+tags_list_rpv_rings_string+"-"+tags_list_other_string)

	print("Extracting the SVC ring for use as a boundary condition...")
	tags_list_svc_ring = extract_tags(input_tags,["SVC_ring"])
	tags_list_svc_ring = [str(t) for t in tags_list_svc_ring]
	tags_list_svc_ring_string = ",".join(tags_list_svc_ring)

	tags_list_other = extract_tags(input_tags,["LV","RV","LA","RA",
											   "Ao","PArt",
											   "MV","TV","AV","PV",
											   "LSPV","LIPV","RSPV","RIPV",
											   "LAA","SVC","IVC",
											   "LAA_ring","IVC_ring",
											   "LSPV_ring","LIPV_ring","RSPV_ring","RIPV_ring"])
	tags_list_other = [str(t) for t in tags_list_other]
	tags_list_other_string = ",".join(tags_list_other)

	os.system("meshtool extract surface -msh="+mesh+" -surf="+presimFolder+"/surfaces_simulation/surfaces_rings/SVC -ofmt=vtk -op="+tags_list_svc_ring_string+"-"+tags_list_other_string)

	print("Converting necessary surfs to vtx files...")
	rpvs_surface = np.loadtxt(presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	rpvs_vtx = surf2vtx(rpvs_surface)
	write_vtx(rpvs_vtx,presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf.vtx")

	svc_surface = np.loadtxt(presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf",dtype=int,skiprows=1,usecols=[1,2,3])
	svc_vtx = surf2vtx(svc_surface)
	write_vtx(svc_vtx,presimFolder+"/surfaces_simulation/surfaces_rings/SVC.surf.vtx")

def combine_elem_dats(heartFolder,presimFolder):
	la_map_dat=heartFolder+"/surfaces_uvc_LA/la/uvc/map_rotational_z.dat"
	ra_map_dat=heartFolder+"/surfaces_uvc_RA/ra/uvc/map_rotational_z.dat"

	os.system("cp "+la_map_dat+" "+presimFolder+"/map_rotational_z_la.dat")
	os.system("cp "+ra_map_dat+" "+presimFolder+"/map_rotational_z_ra.dat")

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
