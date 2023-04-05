import os
import sys
import numpy as np
import copy
import meshio

from common_4ch.file_utils import *
from common_4ch.distance_utils import *

import os
import sys
import numpy as np
import meshio
from math import cos,sin,sqrt
import math

from py_atrial_fibres.file_utils import *

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
    R[2,2] = ax[2]**2 + cos(theta) * (1 - ax[2]**2);

    return R

def surf_find_roof_points_LSPV_RSPV(la_pts,
								    la_surf,
								    lspv_vtx,
								    rspv_vtx,
								    mitral_vtx):

	cog_mitral = np.mean(la_pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(la_pts[lspv_vtx,:],axis=0)
	cog_rspv = np.mean(la_pts[rspv_vtx,:],axis=0)
	v_lspv_rspv = cog_lspv-cog_rspv
	v_lspv_rspv = v_lspv_rspv/np.linalg.norm(v_lspv_rspv)
	v_mitral_rpsv = cog_mitral-cog_rspv

	cog_mitral_pj_ = cog_rspv+v_lspv_rspv*np.dot(v_lspv_rspv,v_mitral_rpsv)
	distance_roof = np.linalg.norm(la_pts-cog_mitral_pj_,axis=1)
	cog_mitral_pj = la_pts[np.where(distance_roof==np.min(distance_roof))[0][0],:]

	distance_lspv = np.linalg.norm(la_pts[lspv_vtx,:]-cog_mitral_pj,axis=1)
	lspv_vtx_roof = lspv_vtx[np.where(distance_lspv==np.min(distance_lspv))[0][0]]

	distance_rspv = np.linalg.norm(la_pts[rspv_vtx,:]-cog_mitral_pj,axis=1)
	rspv_vtx = rspv_vtx[np.where(distance_rspv==np.min(distance_rspv))[0][0]]

	landmark_lspv_roof = la_pts[lspv_vtx_roof,:]
	landmark_rspv_roof = la_pts[rspv_vtx,:]

	return landmark_lspv_roof,landmark_rspv_roof

def surf_find_LSPV_RSPV_posterior_points(la_pts,
								    la_surf,
								    lspv_vtx,
								    lipv_vtx,
								    rspv_vtx,
								    ripv_vtx,
								    mitral_vtx):
	
	cog_mitral = np.mean(la_pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(la_pts[lspv_vtx,:],axis=0)
	cog_rspv = np.mean(la_pts[rspv_vtx,:],axis=0)

	v1 = cog_lspv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	v2 = cog_rspv-cog_mitral
	v2 = v2/np.linalg.norm(v1)

	n_plane = np.cross(v1,v2)
	n_plane = n_plane/np.linalg.norm(n_plane)

	dot_prod = 0.
	print('Finding posterior LSPV rim...')
	for v in lspv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)

		if dot_prod_tmp<dot_prod:
			vtx_posterior_lpsv = v
			dot_prod = dot_prod_tmp

	dot_prod = 0.
	print('Finding posterior RSPV rim...')
	for v in rspv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)
		
		if dot_prod_tmp<dot_prod:
			vtx_posterior_rpsv = v
			dot_prod = dot_prod_tmp

	landmark_lspv_posterior = la_pts[vtx_posterior_lpsv,:]
	landmark_rspv_posterior = la_pts[vtx_posterior_rpsv,:]

	# -------------------------------------------------------------------------------------- #
	pv_region_landmarks = np.zeros((4,3),dtype=float)

	v_rspv_across = cog_rspv-landmark_rspv_posterior	
	v_rspv_across = v_rspv_across/np.linalg.norm(v_rspv_across)

	print('Finding anterior RSPV rim...')
	dot_prod = 0.
	for v in rspv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-landmark_rspv_posterior,v_rspv_across)

		if dot_prod_tmp>dot_prod:
			dot_prod = dot_prod_tmp
			vtx_rspv_anterior = v
	pv_region_landmarks[0,:] = la_pts[vtx_rspv_anterior,:]

	dot_prod = 0.
	print('Finding posterior RIPV rim...')
	for v in ripv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)
		
		if dot_prod_tmp<dot_prod:
			vtx_posterior_ripv = v
			dot_prod = dot_prod_tmp
	pv_region_landmarks[1,:] = la_pts[vtx_posterior_ripv,:]

	dot_prod = 0.
	print('Finding posterior LIPV rim...')
	for v in lipv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)
		
		if dot_prod_tmp<dot_prod:
			vtx_posterior_lipv = v
			dot_prod = dot_prod_tmp
	pv_region_landmarks[2,:] = la_pts[vtx_posterior_lipv,:]

	v_lspv_across = cog_lspv-landmark_lspv_posterior	
	v_lspv_across = v_lspv_across/np.linalg.norm(v_lspv_across)

	print('Finding anterior LSPV rim...')
	dot_prod = 0.
	for v in lspv_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-landmark_lspv_posterior,v_lspv_across)

		if dot_prod_tmp>dot_prod:
			dot_prod = dot_prod_tmp
			vtx_lspv_anterior = v
	pv_region_landmarks[3,:] = la_pts[vtx_lspv_anterior,:]

	return landmark_lspv_posterior,landmark_rspv_posterior,pv_region_landmarks

def surf_find_LAA_septal_posterior_points(la_pts,
	  								      la_surf,
	  								      lspv_vtx,
	  								      lipv_vtx,
	  								      rspv_vtx,
	  								      ripv_vtx,
	  								      laa_vtx,
	  								      mitral_vtx,
	  								      landmark_lspv_posterior,
	  								      landmark_lspv_roof,
	  								      scale_factor=0.001):

	cog_mitral = np.mean(la_pts[mitral_vtx,:],axis=0)
	cog_lspv = np.mean(la_pts[lspv_vtx,:],axis=0)
	cog_ripv = np.mean(la_pts[ripv_vtx,:],axis=0)
	cog_laa = np.mean(la_pts[laa_vtx,:],axis=0)
	cog_rspv = np.mean(la_pts[rspv_vtx,:],axis=0)

	v1 = cog_lspv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	v2 = cog_rspv-cog_mitral
	v2 = v2/np.linalg.norm(v1)

	n_plane = np.cross(v1,v2)
	n_plane = n_plane/np.linalg.norm(n_plane)

	long_axis = cog_rspv-cog_mitral
	long_axis = long_axis/np.linalg.norm(long_axis)

	dot_prod = 1e10
	print('Finding most posterior LAA point...')
	for v in laa_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)
		
		if dot_prod_tmp<dot_prod:
			vtx_posterior_laa = v
			dot_prod = dot_prod_tmp

	landmark_laa_posterior_tmp = la_pts[vtx_posterior_laa,:]

	distance_ant_post = np.abs(np.dot(la_pts-cog_mitral,n_plane))
	distance_range = np.max(distance_ant_post)-np.min(distance_ant_post)
	distance_th = distance_range*0.01 # 1% of distance range
	vtx_ant_post_interface = np.where(distance_ant_post<distance_th)[0]

	distance_long_axis = np.linalg.norm(la_pts[laa_vtx,:]-cog_mitral,axis=1)
	vtx_laa_lowest = laa_vtx[np.where(distance_long_axis==np.min(distance_long_axis))[0][0]]
	landmark_laa_lowest = la_pts[vtx_laa_lowest,:]

	distance_laa_posterior = np.linalg.norm(la_pts[vtx_ant_post_interface,:]-landmark_laa_lowest,axis=1)
	vtx_posterior_free_wall = vtx_ant_post_interface[np.where(distance_laa_posterior==np.min(distance_laa_posterior))[0][0]]
	landmark_laa_posterior = la_pts[vtx_posterior_free_wall,:]

	# n_distance_posterior = np.dot(landmark_lspv_posterior-landmark_lspv_roof,n_plane)
	# landmark_laa_posterior_pj = landmark_laa_posterior+n_plane*n_distance_posterior
	# cog_mitral_pj = cog_mitral+n_plane*n_distance_posterior

	# distance_laa = np.linalg.norm(la_pts-landmark_laa_posterior_pj,axis=1)
	# landmark_laa_posterior = la_pts[np.where(distance_laa==np.min(distance_laa))[0][0]]

	print("Finding septal posterior point...")
	v1 = cog_ripv-cog_mitral
	v1 = v1/np.linalg.norm(v1)

	n_septum = np.cross(v1,v2)
	n_septum = n_septum/np.linalg.norm(n_septum)

	dot_prod = 1e10
	for i in range(la_pts.shape[0]):
		dot_prod_tmp = np.dot(la_pts[i,:]-landmark_laa_posterior_tmp,n_septum)
		dot_prod_vertical = np.dot(la_pts[i,:]-landmark_laa_posterior_tmp,long_axis)

		if dot_prod_tmp<dot_prod and dot_prod_vertical>5000.0*scale_factor:
			vtx_posterior_septum = i
			dot_prod = dot_prod_tmp

	landmark_sept_posterior = la_pts[vtx_posterior_septum,:]

	# --------------------------------------------------------------------------------- #

	dot_prod = 0.
	print('Finding most anterior LAA point...')
	for v in laa_vtx:
		dot_prod_tmp = np.dot(la_pts[v,:]-cog_mitral,n_plane)
		
		if dot_prod_tmp>dot_prod:
			vtx_anterior_laa = v
			dot_prod = dot_prod_tmp
	region_landmark_laa_anterior = la_pts[vtx_anterior_laa,:]

	return landmark_laa_posterior,landmark_sept_posterior,region_landmark_laa_anterior

def extract_boundaries(foldername,
					   clipper_name,
					   input_tags):

	la_pts,la_faces,la_tags = read_surf_polydata(foldername+"/clean-Labelled-refined.vtk")
	mv_center,mv_radius = read_clipper(foldername+"/"+clipper_name)

	mv_radius_new = mv_radius+2.0 

	la_eidx = np.where(la_tags==input_tags["LA"])[0]
	la_vtx = np.unique(la_faces[la_eidx,:].flatten())

	lipv_eidx = np.where(la_tags==input_tags["LIPV"])[0]
	lipv_vtx = np.intersect1d(np.unique(la_faces[lipv_eidx,:].flatten()),la_vtx)

	lspv_eidx = np.where(la_tags==input_tags["LSPV"])[0]
	lspv_vtx = np.intersect1d(np.unique(la_faces[lspv_eidx,:].flatten()),la_vtx)

	ripv_eidx = np.where(la_tags==input_tags["RIPV"])[0]
	ripv_vtx = np.intersect1d(np.unique(la_faces[ripv_eidx,:].flatten()),la_vtx)

	rspv_eidx = np.where(la_tags==input_tags["RSPV"])[0]
	rspv_vtx = np.intersect1d(np.unique(la_faces[rspv_eidx,:].flatten()),la_vtx)

	laa_eidx = np.where(la_tags==input_tags["LAA"])[0]
	laa_vtx = np.intersect1d(np.unique(la_faces[laa_eidx,:].flatten()),la_vtx)

	print('Increasing radius to find mitral valve ring...')

	distance = np.linalg.norm(la_pts[la_vtx,:]-mv_center,axis=1)
	mitral_vtx = la_vtx[np.where(distance<=mv_radius_new)[0]]

	landmark_lspv_roof,landmark_rspv_roof = surf_find_roof_points_LSPV_RSPV(la_pts,
																		    la_faces,
																		    lspv_vtx,
																		    rspv_vtx,
																		    mitral_vtx)

	landmark_lspv_posterior,landmark_rspv_posterior,pv_region_landmarks = surf_find_LSPV_RSPV_posterior_points(la_pts,
																											    la_faces,
																											    lspv_vtx,
																											    lipv_vtx,
																											    rspv_vtx,
																											    ripv_vtx,
																											    mitral_vtx)

	landmark_laa_posterior,landmark_sept_posterior,region_landmark_laa_anterior = surf_find_LAA_septal_posterior_points(la_pts,
																			 	  								        la_faces,
																			 	  								        lspv_vtx,
																			 	  								        lipv_vtx,
																			 	  								        rspv_vtx,
																			 	  								        ripv_vtx,
																			 	  								        laa_vtx,
																			 	  								        mitral_vtx,
																			 	  								        landmark_lspv_posterior,
																			 	  								        landmark_lspv_roof)

	landmarks = np.zeros((4,3),dtype=float)
	landmarks[0,:] = landmark_laa_posterior
	landmarks[1,:] = landmark_sept_posterior
	landmarks[2,:] = landmark_lspv_posterior
	landmarks[3,:] = landmark_rspv_posterior


	np.savetxt(foldername+"/prodLaLandmarks.txt",landmarks,delimiter=',')

	