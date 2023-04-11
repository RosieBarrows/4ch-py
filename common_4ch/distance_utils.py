import os
import sys
import numpy as np
import meshio
from math import cos,sin,sqrt
import math
import copy

import pyvista as pv

from common_4ch.file_utils import *

def calculate_rotation(ax, theta):

    R = np.zeros((3,3))  
    R[0,0] = ax[0]**2 + cos(theta) * (1 - ax[0]**2);
    R[0,1] = (1 - cos(theta)) * ax[0] * ax[1] - ax[2] * sin(theta);
    R[0,2] = (1 - cos(theta)) * ax[0] * ax[2] + ax[1] * sin(theta);
    R[1,0] = (1 - cos(theta)) * ax[0] * ax[1] + ax[2] * sin(theta);
    R[1,1] = ax[1]**2 + cos(theta) * (1 - ax[1]**2);
    R[1,2] = ( 1 - cos(theta)) * ax[1] * ax[2] - ax[0] * sin(theta);
    R[2,0] = ( 1 - cos(theta)) * ax[0] * ax[2] - ax[1] * sin(theta);
    R[2,1] = ( 1 - cos(theta)) * ax[1] * ax[2] + ax[0] * sin(theta);
    R[2,2] = ax[2]**2 + cos(theta) * (1 - ax[2]**2);z

    return R

def compute_path_length(path):

	pts = path.points

	d = 0
	for i in range(pts.shape[0]-1):
		d += np.linalg.norm(pts[i+1]-pts[i])

	return d

def find_midway_point(path):

	pts = path.points

	d = compute_path_length(path)
	d2 = 0.5*d

	d_tmp_array = [0.]
	d_tmp = 0.
	for i in range(pts.shape[0]-1):
		d_tmp += np.linalg.norm(pts[i+1]-pts[i])
		d_tmp_array.append(d_tmp) 
	d_tmp_array = np.array(d_tmp_array)
	vtx_midway = np.argmin(np.abs(d_tmp_array-d2))

	return pts[vtx_midway,:]

def find_closest_point(pv_mesh,
					   idx,
					   vtx):

	d = None
	path = None
	closest = None
	for i,v in enumerate(vtx):
		path_tmp = pv_mesh.geodesic(idx,v)	
		d_tmp = compute_path_length(path_tmp)	

		if (d is not None):
			if d_tmp < d:
				d = d_tmp
				path = path_tmp
				closest = v
			else:
				pass
		else:
			d = d_tmp
			path = path_tmp
			closest = v

	return closest,path

def find_landmarks(output_folder,
				   surface="epi",
				   scale_factor=1.0,
				   raa_apex_file=None):

	surfaces_la = read_surfaces(output_folder,
				  "la",
				  surface)

	la_roof_landmarks = find_roof_points_LSPV_RSPV(surfaces_la,
												   scale_factor=scale_factor)

	la_lspv_rspv_posterior_landmarks,la_region_landmarks_pv = find_LSPV_RSPV_posterior_points(surfaces_la,la_roof_landmarks)

	la_laa_septal_posterior_landmarks,laa_region_landmarks = find_LAA_septal_posterior_points(surfaces_la,
																	np.concatenate((la_roof_landmarks,la_lspv_rspv_posterior_landmarks),axis=0),
																	scale_factor=scale_factor)
	
	landmarks = np.concatenate((la_roof_landmarks,la_laa_septal_posterior_landmarks),axis=0)
	landmarks = np.concatenate((landmarks,la_lspv_rspv_posterior_landmarks),axis=0)

	landmarks_visualise = np.zeros((7,3),dtype=float)
	landmarks_visualise[0,:] = landmarks[0,:]
	landmarks_visualise[1:,:] = landmarks
	write_pts_caroline(landmarks_visualise,output_folder+'/la/landmarks.pts')
	np.savetxt(output_folder+"/la/prodLaLandmarks.txt",landmarks,delimiter=',')

	region_landmarks = np.zeros((6,3),dtype=float)
	region_landmarks[:4,:] = la_region_landmarks_pv
	region_landmarks[4:,:] = laa_region_landmarks
	np.savetxt(output_folder+"/la/prodLaRegion.txt",region_landmarks,delimiter=',')

	landmarks_visualise = np.zeros((7,3),dtype=float)
	landmarks_visualise[0,:] = region_landmarks[0,:]
	landmarks_visualise[1:,:] = region_landmarks
	write_pts_caroline(landmarks_visualise,output_folder+'/la/landmarks_regions.pts')

	# ------------------------------------------------------------------------- #

	surfaces_ra = read_surfaces(output_folder,
				  "ra",
				  surface)

	landmarks_septum = find_septal_point(surfaces_ra)

	landmarks_svc_ivc_posterior, landmarks_regions_roof = find_SVC_IVC_roof_points(surfaces_ra,landmarks_septum)

	landmarks = np.concatenate((landmarks_septum,landmarks_svc_ivc_posterior),axis=0)

	ivc_posterior_landmarks,region_landmarks_ivc = find_IVC_posterior_points(surfaces_ra,landmarks,
																			 scale_factor=scale_factor)
	
	if raa_apex_file is not None:
		if os.path.exists(raa_apex_file):
			landmarks_raa_apex = np.loadtxt(raa_apex_file)
			landmark_raa_base = find_raa_base(surfaces_ra,landmarks_raa_apex)
			region_landmarks_raa = np.concatenate((landmarks_raa_apex,landmark_raa_base),axis=0)
		else:
			raise Exception("Cannot find apex file.")
	else:
		region_landmarks_raa = find_raa_points(surfaces_ra,landmarks_septum[0,:])


	landmarks_rearranged = np.zeros((6,3),dtype=float)
	landmarks_rearranged[0,:] = landmarks_septum[1,:]
	landmarks_rearranged[1,:] = ivc_posterior_landmarks[0,:]
	landmarks_rearranged[2,:] = landmarks_septum[0,:]
	landmarks_rearranged[3,:] = ivc_posterior_landmarks[1,:]
	landmarks_rearranged[4,:] = landmarks_svc_ivc_posterior[0,:]
	landmarks_rearranged[5,:] = landmarks_svc_ivc_posterior[1,:]

	landmarks_visualise[0,:] = landmarks_rearranged[0,:]
	landmarks_visualise[1:,:] = landmarks_rearranged
	write_pts_caroline(landmarks_visualise,output_folder+'/ra/landmarks.pts')
	np.savetxt(output_folder+"/ra/prodRaLandmarks.txt",landmarks_rearranged,delimiter=',')

	region_landmarks = np.zeros((6,3),dtype=float)
	region_landmarks[:2,:] = region_landmarks_ivc
	region_landmarks[2:4,:] = landmarks_regions_roof
	region_landmarks[4:,:] = region_landmarks_raa
	np.savetxt(output_folder+"/ra/prodRaRegion.txt",region_landmarks,delimiter=',')

	landmarks_visualise = np.zeros((7,3),dtype=float)
	landmarks_visualise[0,:] = region_landmarks[0,:]
	landmarks_visualise[1:,:] = region_landmarks
	write_pts_caroline(landmarks_visualise,output_folder+'/ra/landmarks_regions.pts')

def read_surfaces(output_folder,
				  chamber,
				  surface):

	pts = read_pts(output_folder+"/"+chamber+"/"+chamber+".pts")
	surf = read_elem(output_folder+"/"+chamber+"/"+chamber+"_"+surface+".surf",
						 el_type="Tr",
						 tags=False)
	if chamber=="la":

		surfaces_names = ["lspv","rspv","lipv","ripv","laa"]
		if surface == "endo":
			surfaces_list = [s+"_vp" for s in surfaces_names]
		elif surface=="epi":
			surfaces_list = copy.deepcopy(surfaces_names)
		else:
			raise Exception("Surface not recognised. Choose between endo or epi")

		surfaces_list += ["mitral"]
		surfaces_names += ["mitral"]

	elif chamber=="ra":

		surfaces_names = ["svc","ivc"]
		if surface == "endo":
			surfaces_list = [s+"_vp" for s in surfaces_names]
		elif surface=="epi":
			surfaces_list = copy.deepcopy(surfaces_names)
		else:
			raise Exception("Surfacve not recognised. Choose between endo or epi")

		surfaces_list += ["tricuspid"]
		surfaces_names += ["tricuspid"]

	else:
		raise Exception("Chamber not recognised. Choose between la or ra")

	surf_vtx = surf2vtx(surf)
	surf_pts = pts[surf_vtx,:]
	surf_reindexed = reindex_surf(surf_vtx,surf)	

	surfaces_dct = {}
	surfaces_dct["mesh"] = {}
	surfaces_dct["mesh"]["pts"] = pts 

	surfaces_dct["surface"] = {}
	surfaces_dct["surface"]["pts"] = surf_pts 
	surfaces_dct["surface"]["tr"] = surf_reindexed 
	surfaces_dct["surface"]["vtx"] = surf_vtx 

	surfaces_vtx = []
	for i,s in enumerate(surfaces_list):
		
		surf_tmp = read_elem(output_folder+"/"+chamber+"/"+s+".surf",el_type="Tr",tags=False)
		tmp_vtx = surf2vtx(surf_tmp)

		surfaces_dct[surfaces_names[i]] = {}
		if s not in ["mitral","tricuspid"]:
			tmp_vtx_surf = np.intersect1d(surf_vtx,tmp_vtx)
			tmp_vtx_surf = reindex_vtx(tmp_vtx_surf,surf_vtx)	

			surfaces_dct[surfaces_names[i]]["vtx"] = tmp_vtx_surf 

		else:
			surfaces_dct[surfaces_names[i]]["vtx"] = tmp_vtx 
			surfaces_dct[surfaces_names[i]]["surf"] = surf_tmp 

	return surfaces_dct

def find_point_along_direction(points,
							   vtx,
							   n,
							   cog,
							   mode="furthest"):
	
	if mode=="furthest":
		dot_prod = 0.
		idx = None
		for v in vtx:
			dot_prod_tmp = np.dot(points[v,:]-cog,n)	

			if dot_prod_tmp>dot_prod:
				idx = v
				dot_prod = dot_prod_tmp

	elif mode=="closest":		
		dot_prod = 1e10
		for v in vtx:
			dot_prod_tmp = np.dot(points[v,:]-cog,n)	

			if dot_prod_tmp<dot_prod and dot_prod_tmp>=0:
				idx = v
				dot_prod = dot_prod_tmp
	return idx

def find_SVC_IVC_geodesic(output_folder,
						  r_geodesic=1000.0):	
	
	surfaces_dct = read_surfaces(output_folder,
		   				     "ra",
		   				     "epi")

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	tricuspid_vtx = surfaces_dct["tricuspid"]["vtx"]
	svc_vtx = surfaces_dct["svc"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	cog_tricuspid = np.mean(pts[tricuspid_vtx,:],axis=0)
	cog_svc = np.mean(surf_pts[svc_vtx,:],axis=0)
	cog_ivc = np.mean(surf_pts[ivc_vtx,:],axis=0)
	v_svc_ivc = cog_svc-cog_ivc
	v_svc_ivc = v_svc_ivc/np.linalg.norm(v_svc_ivc)
	v_tricuspid_ivc = cog_tricuspid-cog_ivc

	cog_tricuspid_pj = cog_ivc+v_svc_ivc*np.dot(v_svc_ivc,v_tricuspid_ivc)

	long_axis = cog_tricuspid_pj-cog_tricuspid
	long_axis = long_axis/np.linalg.norm(long_axis)

	proj_ra_epi = np.dot(surf_pts-cog_tricuspid,long_axis)
	vtx_roof = np.where(proj_ra_epi==np.max(proj_ra_epi))[0][0]

	pv_mesh = pv.PolyData(surf_pts,
				   		  np.hstack([ np.full(surf.shape[0], 3)[:,None] , surf ]))

	closest_svc,path_svc  = find_closest_point(pv_mesh,	
											   vtx_roof,
											   svc_vtx)

	closest_ivc,path_ivc  = find_closest_point(pv_mesh,	
											   vtx_roof,
											   ivc_vtx)

	path = pv_mesh.geodesic(closest_svc,closest_ivc)
	path_points = path.points

	lines = np.zeros((path_points.shape[0]-1,2),dtype=int)
	lines[:,0] = np.arange(path_points.shape[0]-1)
	lines[:,1] = np.arange(path_points.shape[0]-1)+1
	msh_geod = meshio.Mesh(path_points,
    				  [('line',lines)])	
	msh_geod.write(output_folder+"/ra/svc_ivc_geodesic.vtu")

	anterior_posterior_tag = split_tricuspid_ring(surfaces_dct["mesh"]["pts"],
						 surfaces_dct["tricuspid"]["surf"],
						 surf_pts[closest_svc,:],
						 surf_pts[closest_ivc,:],
						 cog_tricuspid)

	idx_geodesic = []
	for p in path_points:
		d = np.linalg.norm(surf_pts-p,axis=1)
		idx_geodesic += list(np.where(d<=r_geodesic)[0])
	idx_geodesic = np.unique(np.array(idx_geodesic))

	eidx_geodesic = []
	for i,t in enumerate(surf):
		if len(np.intersect1d(t,idx_geodesic))>0:
			eidx_geodesic.append(i)
	eidx_geodesic = np.array(eidx_geodesic)

	return eidx_geodesic,anterior_posterior_tag

def split_tricuspid_ring(pts,
						 tricuspid_surf,
						 svc_pts,
						 ivc_pts,
						 cog_tricuspid):

	v1 = svc_pts-cog_tricuspid
	v1 = v1/np.linalg.norm(v1)

	v2 = ivc_pts-cog_tricuspid
	v2 = v2/np.linalg.norm(v1)

	n_plane = np.cross(v1,v2)
	n_plane = n_plane/np.linalg.norm(n_plane)

	ant_post_tag = np.zeros((tricuspid_surf.shape[0],),dtype=int)
	for i,t in enumerate(tricuspid_surf):
		dot_prod = np.dot(pts[t[0],:]-cog_tricuspid,n_plane)

		if dot_prod>0:
			ant_post_tag[i] = 1	

	return ant_post_tag

def find_roof_points_LSPV_RSPV(surfaces_dct,
							   scale_factor=1.0):

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	mitral_vtx = surfaces_dct["mitral"]["vtx"]
	lspv_vtx = surfaces_dct["lspv"]["vtx"]
	rspv_vtx = surfaces_dct["rspv"]["vtx"]

	cog_mitral = np.mean(pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(surf_pts[lspv_vtx,:],axis=0)
	cog_rspv = np.mean(surf_pts[rspv_vtx,:],axis=0)

	v_lspv_rspv = cog_lspv-cog_rspv
	v_lspv_rspv = v_lspv_rspv/np.linalg.norm(v_lspv_rspv)
	v_mitral_rpsv = cog_mitral-cog_rspv

	cog_mitral_pj = cog_rspv+v_lspv_rspv*np.dot(v_lspv_rspv,v_mitral_rpsv)

	long_axis = cog_mitral_pj-cog_mitral
	long_axis = long_axis/np.linalg.norm(long_axis)

	# making sure the whole roof is included by moving this point down
	cog_mitral_pj = cog_mitral_pj-5000.0*scale_factor*long_axis

	proj_la_epi = np.dot(surf_pts-cog_mitral,long_axis)
	vtx_roof = np.where(proj_la_epi==np.max(proj_la_epi))[0][0]

	pv_mesh = pv.PolyData(surf_pts,
				   		  np.hstack([ np.full(surf.shape[0], 3)[:,None] , surf ]))

	closest_lspv,path_lspv  = find_closest_point(pv_mesh,	
											   vtx_roof,
											   lspv_vtx)

	closest_rspv,path_rspv  = find_closest_point(pv_mesh,	
											   vtx_roof,
											   rspv_vtx)

	pts_lspv_roof = surf_pts[closest_lspv,:]
	pts_rspv_roof = surf_pts[closest_rspv,:]

	landmarks = np.row_stack((pts_lspv_roof,pts_rspv_roof))

	return landmarks

def find_LSPV_RSPV_posterior_points(surfaces_dct,
								    lspv_rspv_roof_points):
	
	pts_lspv_roof = lspv_rspv_roof_points[0,:]
	pts_rspv_roof = lspv_rspv_roof_points[1,:]

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	mitral_vtx = surfaces_dct["mitral"]["vtx"]
	lspv_vtx = surfaces_dct["lspv"]["vtx"]
	rspv_vtx = surfaces_dct["rspv"]["vtx"]
	lipv_vtx = surfaces_dct["lipv"]["vtx"]
	ripv_vtx = surfaces_dct["ripv"]["vtx"]

	cog_mitral = np.mean(pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(surf_pts[lspv_vtx,:],axis=0)
	cog_rspv = np.mean(surf_pts[rspv_vtx,:],axis=0)

	v1 = cog_lspv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	v2 = cog_rspv-cog_mitral
	v2 = v2/np.linalg.norm(v1)

	n_plane = np.cross(v1,v2)
	n_plane = n_plane/np.linalg.norm(n_plane)

	print('Finding posterior LSPV rim...')
	idx_posterior_lpsv = find_point_along_direction(surf_pts,lspv_vtx,-n_plane,cog_mitral)

	print('Finding posterior RSPV rim...')
	idx_posterior_rpsv = find_point_along_direction(surf_pts,rspv_vtx,-n_plane,cog_mitral)

	pts_lspv_posterior = surf_pts[idx_posterior_lpsv,:]
	pts_rspv_posterior = surf_pts[idx_posterior_rpsv,:]

	landmarks = np.concatenate((pts_lspv_posterior,pts_rspv_posterior),axis=0)
	landmarks = np.reshape(landmarks,(2,3))

	# -------------------------------------------------------------------------------------- #
	region_landmarks = np.zeros((4,3),dtype=float)

	v_rspv_across = cog_rspv-pts_rspv_posterior	
	v_rspv_across = v_rspv_across/np.linalg.norm(v_rspv_across)

	print('Finding anterior RSPV rim...')
	idx_anterior_rspv = find_point_along_direction(surf_pts,rspv_vtx,v_rspv_across,pts_rspv_posterior)
	region_landmarks[0,:] = surf_pts[idx_anterior_rspv,:]

	print('Finding posterior RIPV rim...')
	idx_posterior_risv = find_point_along_direction(surf_pts,ripv_vtx,-n_plane,cog_mitral)
	region_landmarks[1,:] = surf_pts[idx_posterior_risv,:]

	print('Finding posterior LIPV rim...')
	idx_posterior_lipv = find_point_along_direction(surf_pts,lipv_vtx,-n_plane,cog_mitral)
	region_landmarks[2,:] = surf_pts[idx_posterior_lipv,:]

	v_lspv_across = cog_lspv-pts_lspv_posterior	
	v_lspv_across = v_lspv_across/np.linalg.norm(v_lspv_across)

	print('Finding anterior LSPV rim...')
	idx_anterior_lspv = find_point_along_direction(surf_pts,lspv_vtx,v_lspv_across,pts_lspv_posterior)
	region_landmarks[3,:] = surf_pts[idx_anterior_lspv,:]

	return landmarks,region_landmarks

def find_LAA_septal_posterior_points(surfaces_dct,
								     landmarks,
							   		 scale_factor=1.0):

	pts_lspv_roof = landmarks[0,:]
	pts_rspv_roof = landmarks[1,:]
	pts_lspv_posterior = landmarks[2,:]
	pts_rspv_posterior = landmarks[3,:]

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	mitral_vtx = surfaces_dct["mitral"]["vtx"]
	lspv_vtx = surfaces_dct["lspv"]["vtx"]
	rspv_vtx = surfaces_dct["rspv"]["vtx"]
	lipv_vtx = surfaces_dct["lipv"]["vtx"]
	ripv_vtx = surfaces_dct["ripv"]["vtx"]
	laa_vtx = surfaces_dct["laa"]["vtx"]

	cog_mitral = np.mean(pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(surf_pts[lspv_vtx,:],axis=0)
	cog_rspv = np.mean(surf_pts[rspv_vtx,:],axis=0)
	cog_ripv = np.mean(surf_pts[ripv_vtx,:],axis=0)
	cog_laa = np.mean(surf_pts[laa_vtx,:],axis=0)

	v1 = cog_lspv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	v2 = cog_rspv-cog_mitral
	v2 = v2/np.linalg.norm(v1)

	n_plane = np.cross(v1,v2)
	n_plane = n_plane/np.linalg.norm(n_plane)

	long_axis = cog_rspv-cog_mitral
	long_axis = long_axis/np.linalg.norm(long_axis)

	print('Finding most posterior LAA point...')
	idx_posterior_laa = find_point_along_direction(surf_pts,laa_vtx,-n_plane,cog_laa,mode="furthest")
	pts_laa_posterior = surf_pts[idx_posterior_laa,:]

	n_distance_posterior = np.dot(pts_lspv_posterior-pts_lspv_roof,n_plane)
	pts_laa_posterior_pj = pts_laa_posterior+n_plane*n_distance_posterior
	cog_mitral_pj = cog_mitral+n_plane*n_distance_posterior

	distance_laa = np.linalg.norm(surf_pts-pts_laa_posterior_pj,axis=1)
	landmark_laa_posterior = surf_pts[np.where(distance_laa==np.min(distance_laa))[0][0]]

	print("Finding septal posterior point...")
	v1 = cog_ripv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	n_septum = np.cross(v1,v2)
	n_septum = n_septum/np.linalg.norm(n_septum)

	vtx_septum_area = []
	for i in range(surf_pts.shape[0]):
		dot_prod_vertical = np.dot(surf_pts[i,:]-landmark_laa_posterior,long_axis)
		if dot_prod_vertical>5000.0*scale_factor:
			vtx_septum_area.append(i)
	vtx_septum_area = np.array(vtx_septum_area)

	idx_posterior_septum = find_point_along_direction(surf_pts,vtx_septum_area,-n_septum,landmark_laa_posterior)

	landmark_sept_posterior = surf_pts[idx_posterior_septum,:]

	landmarks = np.concatenate((landmark_laa_posterior,landmark_sept_posterior),axis=0)
	landmarks = np.reshape(landmarks,(2,3))

	# --------------------------------------------------------------------------------- #
	region_landmarks = np.zeros((2,3),dtype=float)
	region_landmarks[0,:] = pts_laa_posterior

	dot_prod = 0.
	print('Finding most anterior LAA point...')
	idx_anterior_laa = find_point_along_direction(surf_pts,laa_vtx,n_plane,cog_mitral)
	region_landmarks[1,:] = surf_pts[idx_anterior_laa,:]

	return landmarks,region_landmarks

def find_septal_point(surfaces_dct):

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	tricuspid_vtx = surfaces_dct["tricuspid"]["vtx"]
	svc_vtx = surfaces_dct["svc"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	cog_tricuspid = np.mean(pts[tricuspid_vtx,:],axis=0)
	cog_svc = np.mean(surf_pts[svc_vtx,:],axis=0)
	cog_ivc = np.mean(surf_pts[ivc_vtx,:],axis=0)
	v_svc_ivc = cog_svc-cog_ivc
	v_svc_ivc = v_svc_ivc/np.linalg.norm(v_svc_ivc)
	v_tricuspid_ivc = cog_tricuspid-cog_ivc
	v_tricuspid_ivc = v_tricuspid_ivc/np.linalg.norm(v_tricuspid_ivc)

	cog_svc_tricuspid = 0.5*(cog_tricuspid+cog_svc)

	n_septum = np.cross(v_tricuspid_ivc,v_svc_ivc)
	n_septum = n_septum/np.linalg.norm(n_septum)

	vtx_septum_area = []
	for i in range(surf_pts.shape[0]):
		dot_prod_sept = np.dot(surf_pts[i,:]-cog_svc_tricuspid,n_septum)
		if dot_prod_sept<0.:
			vtx_septum_area.append(i)
	vtx_septum_area = np.array(vtx_septum_area)	

	distances = np.linalg.norm(surf_pts[vtx_septum_area,:]-cog_svc_tricuspid,axis=1)
	idx_septum = vtx_septum_area[np.where(distances==np.min(distances))[0][0]]
	landmark_septum = surf_pts[idx_septum,:]

	pv_mesh = pv.PolyData(surf_pts,
				   		  np.hstack([ np.full(surf.shape[0], 3)[:,None] , surf ]))
	closest_svc,path = find_closest_point(pv_mesh,idx_septum,svc_vtx)

	landmark_svc_septum = surf_pts[closest_svc,:]

	landmarks = np.row_stack((landmark_septum,landmark_svc_septum))
	landmarks = np.reshape(landmarks,(2,3))

	return landmarks

def find_SVC_IVC_roof_points(surfaces_dct,
							 landmarks_septum):

	landmark_septum = landmarks_septum[0,:]
	landmark_svc_septum = landmarks_septum[1,:]

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	tricuspid_vtx = surfaces_dct["tricuspid"]["vtx"]
	svc_vtx = surfaces_dct["svc"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	cog_tricuspid = np.mean(pts[tricuspid_vtx,:],axis=0)
	cog_svc = np.mean(surf_pts[svc_vtx,:],axis=0)
	cog_ivc = np.mean(surf_pts[ivc_vtx,:],axis=0)

	v_svc_ivc = cog_svc-cog_ivc
	v_svc_ivc = v_svc_ivc/np.linalg.norm(v_svc_ivc)
	v_tricuspid_ivc = cog_tricuspid-cog_ivc
	d_tricuspid_ivc = np.linalg.norm(v_tricuspid_ivc)
	v_tricuspid_ivc = v_tricuspid_ivc/d_tricuspid_ivc

	n_septum = np.cross(v_tricuspid_ivc,v_svc_ivc)
	n_septum = n_septum/np.linalg.norm(n_septum)

	print('Finding posterior SVC rim...')
	idx_posterior_svc = find_point_along_direction(surf_pts,svc_vtx,n_septum,cog_tricuspid)
	landmark_svc_posterior = surf_pts[idx_posterior_svc,:]

	cog_tricuspid_pj = cog_ivc+v_svc_ivc*np.dot(v_svc_ivc,d_tricuspid_ivc*v_tricuspid_ivc)

	long_axis = cog_tricuspid_pj-cog_tricuspid
	long_axis = long_axis/np.linalg.norm(long_axis)

	proj_ra_epi = np.dot(surf_pts-cog_tricuspid,long_axis)
	vtx_roof = np.where(proj_ra_epi==np.max(proj_ra_epi))[0][0]

	pv_mesh = pv.PolyData(surf_pts,
				   		  np.hstack([ np.full(surf.shape[0], 3)[:,None] , surf ]))

	closest_ivc,path_ivc  = find_closest_point(pv_mesh,	
											   vtx_roof,
											   ivc_vtx)
	landmark_ivc_posterior = surf_pts[closest_ivc,:]

	path_svc_ivc = pv_mesh.geodesic(idx_posterior_svc,closest_ivc)

	landmarks = np.row_stack((landmark_svc_posterior,landmark_ivc_posterior))
	landmarks = np.reshape(landmarks,(2,3))

	# ----------------------------------------------------------------------- #
	landmarks_regions = np.zeros((2,3),dtype=float)

	midway_point = find_midway_point(path_svc_ivc)

	landmarks_regions[0,:] = midway_point
	landmarks_regions[1,:] = landmark_svc_posterior

	return landmarks,landmarks_regions

def find_IVC_posterior_points(surfaces_dct,
							  landmarks,
							  scale_factor=1.0):

	landmark_septum = landmarks[0,:]
	landmark_svc_septum = landmarks[1,:]
	landmark_svc_posterior = landmarks[2,:]
	landmark_ivc_top = landmarks[3,:]

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	surf_vtx = surfaces_dct["surface"]["vtx"]
	tricuspid_vtx = surfaces_dct["tricuspid"]["vtx"]
	svc_vtx = surfaces_dct["svc"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	tricuspid_vtx_epi = np.intersect1d(surf_vtx,tricuspid_vtx)
	tricuspid_vtx_epi = reindex_vtx(tricuspid_vtx_epi,surf_vtx)	

	cog_tricuspid = np.mean(surf_pts[tricuspid_vtx_epi,:],axis=0)
	cog_svc = np.mean(surf_pts[svc_vtx,:],axis=0)
	cog_ivc = np.mean(surf_pts[ivc_vtx,:],axis=0)

	v_svc_ivc = cog_svc-cog_ivc
	v_svc_ivc = v_svc_ivc/np.linalg.norm(v_svc_ivc)
	v_tricuspid_ivc = cog_tricuspid-cog_ivc
	d_tricuspid_ivc = np.linalg.norm(v_tricuspid_ivc)
	v_tricuspid_ivc = v_tricuspid_ivc/d_tricuspid_ivc

	v_ivc_across = cog_ivc-landmark_ivc_top	
	v_ivc_across = v_ivc_across/np.linalg.norm(v_ivc_across)

	idx_ivc_bottom = find_point_along_direction(surf_pts,ivc_vtx,v_ivc_across,landmark_ivc_top)
	landmarks_ivc_bottom = surf_pts[idx_ivc_bottom,:]

	distances = np.linalg.norm(surf_pts[tricuspid_vtx_epi	,:]-landmarks_ivc_bottom,axis=1)
	idx_tricuspid = tricuspid_vtx_epi[np.where(distances==np.min(distances))[0][0]]
	shortest = np.min(distances)

	midway = 0.5*(surf_pts[idx_tricuspid,:]+landmarks_ivc_bottom)
	distances = np.linalg.norm(surf_pts-midway,axis=1)
	idx_tricuspid_posterior = np.where(distances==np.min(distances))[0][0]
	landmark_tricuspid_posterior = surf_pts[idx_tricuspid_posterior,:]

	landmarks = np.row_stack((landmarks_ivc_bottom,landmark_tricuspid_posterior))
	landmarks = np.reshape(landmarks,(2,3))

	# ----------------------------------------------------------------------------------- #
	region_landmarks = np.zeros((2,3),dtype=float)
	region_landmarks[0,:] = landmarks_ivc_bottom

	v_tricuspid_svc = cog_svc-cog_tricuspid
	v_tricuspid_svc = v_tricuspid_svc/np.linalg.norm(v_tricuspid_svc)

	v_tricuspid_septum = landmark_septum-cog_tricuspid
	v_tricuspid_septum = v_tricuspid_septum/np.linalg.norm(v_tricuspid_septum)

	n_raa = np.cross(v_tricuspid_septum,v_tricuspid_svc)
	n_raa = n_raa/np.linalg.norm(n_raa)

	septal_point_CS = landmark_septum-n_raa*10000.0*scale_factor
	distances = np.linalg.norm(surf_pts-septal_point_CS,axis=1)
	region_landmarks[1,:] = surf_pts[np.where(distances==np.min(distances))[0][0],:]

	return landmarks,region_landmarks

def find_raa_base(output_folder,
				  raa_apex,
				  raa_length=10000.0,
				  surface="epi"):
	
	surfaces_dct = read_surfaces(output_folder,"ra",surface)

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	surf_vtx = surfaces_dct["surface"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	cog_ra = np.mean(surf_pts,axis=0)

	v_ivc_raa = raa_apex-cog_ra
	v_ivc_raa = v_ivc_raa/np.linalg.norm(v_ivc_raa)

	raa_tip_pj = raa_apex-v_ivc_raa*raa_length
	distances = np.linalg.norm(surf_pts-raa_tip_pj,axis=1)
	landmark_raa_base = surf_pts[np.where(distances==np.min(distances))[0][0],:]

	landmarks = np.reshape(landmark_raa_base,(1,3))

	return landmarks

def find_raa_points(surfaces_dct,
					septal_point,
					raa_length=10000.0):

	pts = surfaces_dct["mesh"]["pts"]
	surf_pts = surfaces_dct["surface"]["pts"]
	surf = surfaces_dct["surface"]["tr"]
	surf_vtx = surfaces_dct["surface"]["vtx"]
	tricuspid_vtx = surfaces_dct["tricuspid"]["vtx"]
	svc_vtx = surfaces_dct["svc"]["vtx"]
	ivc_vtx = surfaces_dct["ivc"]["vtx"]

	cog_ra = np.mean(surf_pts,axis=0)

	tricuspid_vtx_epi = np.intersect1d(surf_vtx,tricuspid_vtx)
	tricuspid_vtx_epi = reindex_vtx(tricuspid_vtx_epi,surf_vtx)	

	cog_tricuspid = np.mean(pts[tricuspid_vtx,:],axis=0)
	cog_svc = np.mean(surf_pts[svc_vtx,:],axis=0)
	cog_ivc = np.mean(surf_pts[ivc_vtx,:],axis=0)	

	v_tricuspid_svc = cog_svc-cog_tricuspid
	v_tricuspid_svc = v_tricuspid_svc/np.linalg.norm(v_tricuspid_svc)	

	v_tricuspid_septum = septal_point-cog_tricuspid
	v_tricuspid_septum = v_tricuspid_septum/np.linalg.norm(v_tricuspid_septum)	

	n_raa = np.cross(v_tricuspid_septum,v_tricuspid_svc)
	n_raa = n_raa/np.linalg.norm(n_raa)

	print('Finding RAA tip...')
	vtx_raa_tip = find_point_along_direction(surf_pts,np.arange(surf_pts.shape[0]),n_raa,cog_ivc)

	landmark_raa_tip = surf_pts[vtx_raa_tip,:]

	v_ivc_raa = landmark_raa_tip-cog_ivc
	v_ivc_raa = v_ivc_raa/np.linalg.norm(v_ivc_raa)

	raa_tip_pj = landmark_raa_tip-v_ivc_raa*raa_length
	distances = np.linalg.norm(surf_pts-raa_tip_pj,axis=1)
	landmark_raa_base = surf_pts[np.where(distances==np.min(distances))[0][0],:]

	landmarks = np.row_stack((landmark_raa_tip,landmark_raa_base))
	landmarks = np.reshape(landmarks,(2,3))

	return landmarks