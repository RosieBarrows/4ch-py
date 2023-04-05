import numpy as np
import os
import copy

from common.file_utils import *

def rotation_matrix(axis,theta):
    R = np.zeros((3,3), dtype=float)

    R[0][0] = axis[0]**2 + np.cos(theta)*(1-axis[0]**2)
    R[0][1] = (1 - np.cos(theta)) * axis[0] * axis[1] - axis[2] * np.sin(theta)
    R[0][2] = (1 - np.cos(theta)) * axis[0] * axis[2] + axis[1] * np.sin(theta)
    R[1][0] = (1 - np.cos(theta)) * axis[0] * axis[1] + axis[2] * np.sin(theta)
    R[1][1] = axis[1]**2 + np.cos(theta)*(1-axis[1]**2)
    R[1][2] = (1 - np.cos(theta)) * axis[1] * axis[2] - axis[0] * np.sin(theta)
    R[2][0] = (1 - np.cos(theta)) * axis[0] * axis[2] - axis[1] * np.sin(theta)
    R[2][1] = (1 - np.cos(theta)) * axis[1] * axis[2] + axis[0] * np.sin(theta)
    R[2][2] = axis[2]**2 + np.cos(theta)*(1-axis[2]**2)

    return R

def normalise_vectors(v):

	v_norm = copy.deepcopy(v)

	for i,vv in enumerate(v):
		v_norm[i,:] = vv/np.linalg.norm(vv)

	return v_norm

def find_rotation_axes(v1_file,
					   v2_file,
					   output_file):
	
	v1 = read_lon(v1_file)
	v1_norm = normalise_vectors(v1)

	v2 = read_lon(v2_file)
	v2_norm = normalise_vectors(v2)

	rotation_axes = copy.deepcopy(v1)
	for i in range(v1.shape[0]):
		rotation_axes[i,:] = np.cross(v1_norm[i,:],v2_norm[i,:])
	ax_norm = normalise_vectors(rotation_axes)

	write_lon(ax_norm,output_file)

def make_sheet_orthogonal(fibres_file,
						  transmural_d_file,
						  rotation_axes_file,
						  output_file):

	fibres = read_lon(fibres_file)
	transmural = read_lon(transmural_d_file)
	transmural = normalise_vectors(transmural)
	rot_ax = read_lon(rotation_axes_file)

	sheet = np.zeros((fibres.shape[0],3), dtype=float)
	for i,n in enumerate(rot_ax):

	    alpha = np.arccos(np.dot(fibres[i,:],transmural[i,:]))
	    theta = np.pi/2 - alpha
	    sheet[i,:] = np.matmul(rotation_matrix(rot_ax[i,:], theta),transmural[i,:])	

	sheet = normalise_vectors(sheet)
	fibres_sheet = np.column_stack((fibres,sheet))

	write_lon(fibres_sheet,output_file)

def compute_normal_from_surface(meshname):

	pts = read_pts(meshname+".pts")
	elem = read_elem(meshname+".elem",el_type="Tr",tags=False)

	normals = np.zeros((elem.shape[0],3),dtype=float)
	for i,tr in enumerate(elem):
		v10 = pts[tr[1],:]-pts[tr[0],:]
		v20 = pts[tr[2],:]-pts[tr[0],:]

		v10 = v10/np.linalg.norm(v10)
		v20 = v20/np.linalg.norm(v20)

		n = np.cross(v10,v20)
		normals[i,:] = -n/np.linalg.norm(n)

	np.savetxt(meshname+"_normals.txt",normals,fmt="%g")

	return normals





