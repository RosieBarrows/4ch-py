import argparse
import os
import json
import numpy as np
from functools import reduce

# from SIMULATION_library.electrode_utils import *
from common_4ch.file_utils import *

def float_to_tuple(f, n=0):
	mult = 10**n
	return tuple(mult*np.floor(x/mult) for x in f)  
def create_SAN(heartFolder):
	"""
	Function to extract the sinoatrial node vtx from the atrial fibres code.
	"""

	SAN_tags_pts_ra_endo = np.loadtxt(os.path.join(heartFolder,"atrial_fibres","UAC","RA_endo","MappedScalar_SAN.dat"),dtype=float)
	ra_endo_pts          = read_pts(os.path.join(heartFolder,"atrial_fibres","UAC","RA_endo","RA_only.pts"))
	fourch_pts           = read_pts(os.path.join(heartFolder,"pre_simulation","myocardium_AV_FEC_BB.pts"))

	finished = False
	N=0
	count = 0
	while not finished:
		fourch_dict = {float_to_tuple(elem, N): index for index, elem in enumerate(fourch_pts)}
		SAN_vtx_ra_endo = np.where(SAN_tags_pts_ra_endo > 0)

		SAN_pts = ra_endo_pts[SAN_vtx_ra_endo,:]
		SAN_vtx_fourch = np.array([fourch_dict.get(float_to_tuple(san_elem, N)) for san_elem in SAN_pts[0]])

		## get number of None values
		n_none = np.sum(SAN_vtx_fourch == None)
		if n_none == 0 or count > 3:
			finished = True
		else:
			print('[Warning] Could not find SAN points. Increasing the number of decimals...')
			count+=1
			N+=1
		
		if count > 3:
			raise Exception("Could not find SAN points")

	output_filename = os.path.join(heartFolder,"pre_simulation","SAN.vtx")
	write_vtx(filename = output_filename, vtx = SAN_vtx_fourch)

	if os.path.exists(output_filename):
		print(f'SAN vtx file created: {output_filename}')


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

def combine_electrode(vtx_file_list, vtx_file_out):
	
	vtx_out = []
	for v in vtx_file_list:
		vtx_tmp = read_vtx(v)
		vtx_out += list(vtx_tmp)

	vtx_out = np.array(vtx_out)
	vtx_out = np.unique(vtx_out)

	write_vtx(vtx_file_out,vtx_out,init_row=2)

def write_pts(pts,filename):

	print('Writing '+filename+'...')

	assert pts.shape[1] == 3
	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(pts.shape[0]))
		for pnt in pts:
			fp.write('{0[0]} {0[1]} {0[2]}\n'.format(pnt))
	fp.close()

def vtx2pts(vtx_file, pts_file):

	vtx = read_vtx(vtx_file)
	pts = read_pts(pts_file)

	pts_vtx = pts[vtx,:]

	write_pts(pts_vtx,vtx_file+".pts")


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