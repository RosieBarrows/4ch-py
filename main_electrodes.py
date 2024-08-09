import argparse
import os
import json
import numpy as np
from functools import reduce

from SIMULATION_library.electrode_utils import *
from common_4ch.file_utils import *



def create_SAN(heartFolder):
	"""
	Function to extract the sinoatrial node vtx from the atrial fibres code.
	"""

	SAN_tags_pts_ra_endo = np.loadtxt(os.path.join(heartFolder, "atrial_fibres", "UAC", "RA_endo", "MappedScalar_SAN.dat"), dtype=float)
	ra_endo_pts = read_pts(os.path.join(heartFolder, "atrial_fibres", "UAC", "RA_endo", "RA_only.pts"))
	fourch_pts = read_pts(os.path.join(heartFolder, "pre_simulation", "myocardium_AV_FEC_BB.pts"))

	fourch_pts = np.array(fourch_pts)  # Convert to NumPy array for efficient calculations

	def euclidean_distances(pts1, pts2):
		return np.linalg.norm(pts1[:, np.newaxis] - pts2, axis=2)

	SAN_vtx_ra_endo = np.where(SAN_tags_pts_ra_endo > 0)

	SAN_pts = ra_endo_pts[SAN_vtx_ra_endo]

	# Calculate Euclidean distances for all SAN_pts against fourch_pts
	distances = euclidean_distances(SAN_pts, fourch_pts)

	# Find the indices with the minimum distance
	SAN_vtx_fourch = np.argmin(distances, axis=1)

	write_vtx(filename=os.path.join(heartFolder, "pre_simulation", "SAN.vtx"),
				vtx=SAN_vtx_fourch)


def find_electrode_UVC_cylinder(uvc,
					   fascicles_settings):

	v0 = fascicles_settings["v"]
	rho0 = fascicles_settings["rho"]
	radius_rho = fascicles_settings["radius_rho"]

	rho = uvc[:,1] 
	v = uvc[:,3] 

	surface_v = reduce(np.intersect1d,
		    			(np.where(v==v0)[0],
       					 np.where(rho<=rho0)[0],
						 np.where(rho>=rho0-radius_rho)[0]
						 )
					   )
	
	distance = np.linalg.norm(uvc[surface_v,:][:,[0,2]]-np.array([fascicles_settings["z"],fascicles_settings["phi"]]),axis=1)
	
	electrode_idx = surface_v[np.where(distance<=fascicles_settings["radius_phi"])[0]]

	electrode_idx = np.unique(electrode_idx)

	return electrode_idx

def define_fascicles_UVC(bivnod_file,
							  bivuvc_file,
							  fascicles_settings,
							  output_folder):


	biv_nod = read_nod_eidx(bivnod_file)

	uvc = np.loadtxt(bivuvc_file,dtype=float)
	z = uvc[:,0]

	if z.shape[0]!=biv_nod.shape[0]:
		raise Exception(".nod file and UVCs do not match in size.")

	fascicles_list = list(fascicles_settings.keys())

	print('=============================================')
	print('Finding fascicles...')
	print('=============================================')

	for f in fascicles_list:
				print(f)
				electrode_idx = find_electrode_UVC_cylinder(uvc                = uvc,
												fascicles_settings = fascicles_settings[f])

				write_vtx(output_folder+"/"+f+".vtx", biv_nod[electrode_idx],init_row=2)

def main(args):

	# We need the combined UVCs without the first line

	bivuvc_file           = os.path.join(args.heartFolder,"surfaces_uvc","BiV","uvc","BiV.uvc")
	bivuvc_file_no_header = os.path.join(args.heartFolder,"surfaces_uvc","BiV","uvc","BiV_no_header.uvc")

	cmd = "tail -n +2 " + bivuvc_file + " > " + bivuvc_file_no_header
	os.system(cmd)

	# We read the fascicle settings
	f = open(args.fascicles_settings,"r")
	fascicles_settings = json.load(f)
	f.close()

	# ### Ventricular electrodes

	define_fascicles_UVC(bivnod_file        = os.path.join(args.heartFolder,"surfaces_uvc","BiV","BiV.nod"),
					 bivuvc_file        = bivuvc_file_no_header,
					 fascicles_settings = fascicles_settings,
					 output_folder      = os.path.join(args.heartFolder,"pre_simulation"))
	
	combine_electrode(vtx_file_list = [os.path.join(args.heartFolder,"pre_simulation",electrode) for electrode in ["LVsept.vtx","LVpost.vtx","LVant.vtx"]],
					  vtx_file_out  = os.path.join(args.heartFolder,"pre_simulation","fascicles_lv.vtx"))
	combine_electrode(vtx_file_list = [os.path.join(args.heartFolder,"pre_simulation",electrode) for electrode in ["RVsept.vtx","RVmod.vtx"]],
					  vtx_file_out  = os.path.join(args.heartFolder,"pre_simulation","fascicles_rv.vtx"))
	
	# Atrial electrodes

	create_SAN(heartFolder = args.heartFolder)

	# For quality check
	vtx2pts(vtx_file = os.path.join(args.heartFolder,"pre_simulation","fascicles_lv.vtx"),
			pts_file = os.path.join(args.heartFolder,"pre_simulation","myocardium_AV_FEC_BB.pts"))
	vtx2pts(vtx_file = os.path.join(args.heartFolder,"pre_simulation","fascicles_rv.vtx"),
			pts_file = os.path.join(args.heartFolder,"pre_simulation","myocardium_AV_FEC_BB.pts"))
	vtx2pts(vtx_file = os.path.join(args.heartFolder,"pre_simulation","SAN.vtx"),
			pts_file = os.path.join(args.heartFolder,"pre_simulation","myocardium_AV_FEC_BB.pts"))
					  
					 
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

	parser.add_argument('--heartFolder', type=str, default=None,
						help='Provide path to the heart folder')
	parser.add_argument('--fascicles_settings', type=str, default=None,
						help='Provide path to the fascicle settings file')
	
	args = parser.parse_args()

	main(args)