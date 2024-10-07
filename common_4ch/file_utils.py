import numpy as np
import copy
import json
import meshio
import vtk
import os
import sys
from vtk.util import numpy_support

from common_4ch.config import configure_logging
milog = configure_logging(log_name=__name__)

def chooseplatform(): 
	return sys.platform


class NumpyEncoder(json.JSONEncoder):
	def default(self, obj):
		if isinstance(obj, np.integer):
			return int(obj)
		elif isinstance(obj, np.floating):
			return float(obj)
		elif isinstance(obj, np.ndarray):
			return obj.tolist()
		return json.JSONEncoder.default(self, obj)

def mycp(src,dst, debug=False):
	if chooseplatform() == 'win32':
		cmd = f"copy {src} {dst}"
	else:
		cmd = f"cp {src} {dst}"
	if debug:
		milog.info(cmd)
	os.system(cmd)

def mymv(src,dst, debug=False):
	if chooseplatform() == 'win32':
		cmd = f"move {src} {dst}"
	else:
		cmd = f"mv {src} {dst}"
	if debug:
		milog.info(cmd)
	os.system(cmd)

def mymkdir(path, full_path=False, debug=False):
	# check if path exists
	if os.path.exists(path):
		milog.info(f"Path [{path}] already exists")
		return
	
	flag = "-p" if full_path else ""
	cmd = f"mkdir {flag} {path}"
	if debug:
		milog.info(cmd)
	os.system(cmd)

def myrm(path, debug=False):
	if chooseplatform() == 'win32':
		cmd = f"rmdir /s /q {path}"
	else:
		cmd = f"rm -rf {path}"
	if debug:
		milog.info(cmd)
	os.system(cmd)

def big_msg(msg):
	length = len(msg)
	milog.info("-"*length)
	milog.info(msg)
	milog.info("-"*length)

def pjoin(*paths) :
	return os.path.join(*paths)
	
def read_pts(filename):
	milog.info(f'Reading {filename}...')
	return np.loadtxt(filename, dtype=float, skiprows=1)

def read_elem(filename,el_type='Tt',tags=True):
	milog.info(f'Reading {filename}...')
	cols_notags_dic = {'Tt':(1,2,3,4),'Tr':(1,2,3),'Ln':(1,2)}
	try: 
		cols = cols_notags_dic[el_type]
		if tags:
			# add tags column (largest + 1)
			cols += (cols[-1]+1,)

		return np.loadtxt(filename, dtype=int, skiprows=1, usecols=cols)
	
	except KeyError:
		error_msg = f"element type not recognised. Accepted: {list(cols_notags_dic.keys())}"
		milog.error(error_msg)
		raise Exception(error_msg)

def read_lon(filename):
	milog.info(f'Reading {filename}...')

	return np.loadtxt(filename, dtype=float, skiprows=1)

def get_polydata(filename) -> vtk.vtkPolyData:
	milog.info(f'Reading {filename}...')

	polydata_reader = vtk.vtkPolyDataReader()

	polydata_reader.SetFileName(filename)
	polydata_reader.Update()

	return polydata_reader.GetOutput()

def read_surf_polydata(filename):

	polydata = get_polydata(filename)

	points = numpy_support.vtk_to_numpy(polydata.GetPoints().GetData())

	faces = numpy_support.vtk_to_numpy(polydata.GetPolys().GetData())
	faces = faces.reshape((int(faces.shape[0]/4),4))
	faces = faces[:,1:]

	tags = numpy_support.vtk_to_numpy(polydata.GetCellData().GetArray('elemTag'))

	milog.info('Done.')

	return points,faces,tags

def read_clipper(filename):

	polydata = get_polydata(filename)

	points = numpy_support.vtk_to_numpy(polydata.GetPoints().GetData())
	center = np.mean(points,axis=0)
	radius = np.linalg.norm(points[0,:]-center)

	milog.info('Done.')

	return center,radius

def read_nod_eidx(filename):
	idx = np.fromfile(filename, dtype=int, count=-1)

	return idx

def surf2vtx(surf):

	vtx = np.unique(surf.flatten())

	return vtx

def reindex_vtx(vtx,vtx_surf):
	
	vtx = np.intersect1d(vtx,vtx_surf)
	vtx_reindexed = copy.deepcopy(vtx)
	for i,v in enumerate(vtx):
		vtx_reindexed[i] = np.where(vtx_surf==v)[0]

	return vtx_reindexed

def reindex_surf(vtx_surf,
				 surf):
	
	surf_reindexed = copy.deepcopy(surf)
	for i,t in enumerate(surf):
		surf_reindexed[i,0] = np.where(vtx_surf==t[0])[0][0]
		surf_reindexed[i,1] = np.where(vtx_surf==t[1])[0][0]
		surf_reindexed[i,2] = np.where(vtx_surf==t[2])[0][0]

	return surf_reindexed

def write_surf_caroline(filename,surf):
	milog.info(f'Writing {filename}...')
	milog.info("Check you are using the correct write_surf function")

	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(surf.shape[0]))
		for t in surf:
			fp.write('Tr {} {} {}\n'.format(int(t[0]),int(t[1]),int(t[2])))

def write_pts_caroline(pts,filename):

	milog.info(f'Writing {filename}...')

	assert pts.shape[1] == 3
	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(pts.shape[0]))
		for pnt in pts:
			fp.write('{0[0]} {0[1]} {0[2]}\n'.format(pnt))

def write_elem_caroline(elem,
			   tags,
			   filename,
			   el_type='Tt'):
	milog.info(f'Writing {filename}...')
	elem_shape_dic = {'Tt':4,'Tr':3,'Ln':2}
	try:
		elem_shape = elem_shape_dic[el_type]
		assert elem.shape[1] == elem_shape
	except KeyError:
		error_msg = f"element type not recognised. Accepted: {list(elem_shape_dic.keys())}"
		milog.error(error_msg)
		raise Exception(error_msg)
	
	assert elem.shape[0] == tags.shape[0]

	with open(filename, 'w') as fe:
		fe.write('{}\n'.format(elem.shape[0]))
		for i,el in enumerate(elem):
			fe.write(el_type)
			for e in el:
				fe.write(' '+str(e))
			fe.write(' '+str(tags[i]))
			fe.write('\n')

def write_tets_ln(filename,
				  elem,
				  ln):

	milog.info(f'Writing {filename}...')

	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(elem.shape[0] + ln.shape[0]))
		for el in elem:
			fp.write('Tt {} {} {} {} {}\n'.format(int(el[0]),int(el[1]),int(el[2]),int(el[3]),int(el[4])))
		for cn in ln:
			fp.write('Ln {} {} {}\n'.format(int(cn[0]),int(cn[1]),int(cn[2])))

def write_lon(lon,filename):
	milog.info(f'Writing {filename}...')

	assert lon.shape[1] % 3 == 0
	with open(filename, 'w') as fl:
		fl.write('{}\n'.format(int(lon.shape[1]/3)))
		for ll in lon:
			for i,l in enumerate(ll):
				if i==len(ll)-1:
					fl.write(str(l)+'\n')
				else:
					fl.write(str(l)+' ')

def write_vtx(filename, vtx,init_row=2):
	milog.info(f'Writing {filename}...')

	print(f'\n\n\n vtx: {vtx.shape}')

	with open(filename, 'w') as fd:
		if init_row==2:
			fd.write('{}\n'.format(vtx.shape[0]))
			fd.write('intra\n')
			for v in vtx:
				fd.write('{}\n'.format(int(v)))

def save_json(dct, filename):
	with open(filename, "w") as f:
		json.dump(dct, f, cls=NumpyEncoder, indent=4)
	return

def numpy_hook(dct):
	for key, value in dct.items():
		if isinstance(value, list):
			value = np.array(value)
			dct[key] = value
	return dct

def load_json(filename):
	milog.info(f'Reading {filename}...')

	dct = {}
	with open(filename, "r") as f:
		dct = json.load(f, object_hook=numpy_hook)
	return dct

def read_vtx(filename,init_row=2):
	milog.info(f'Reading {filename}...')

	return np.loadtxt(filename, dtype=int, skiprows=init_row)


def setup_sim(heartFolder,presimFolder,electrodes_paths_array):
	sims_folder = pjoin(heartFolder,"sims_folder")
	mymkdir(sims_folder) 

	surf_sim_folder = pjoin(presimFolder,"surfaces_simulation")
	surf_rings_folder = pjoin(surf_sim_folder,"surfaces_rings")
	
	mycp(pjoin(presimFolder,"elem_dat_UVC_ek_combined.dat"), pjoin(sims_folder, "pericardium_scale.dat"))
	fnames = [
		( presimFolder, "epicardium_for_sim.surf"),
		( surf_sim_folder, "LA_endo.surf"),
		( surf_sim_folder, "LV_endo.surf"),
		( surf_sim_folder, "RA_endo.surf"),
		( surf_sim_folder, "RV_endo.surf"),
		( surf_rings_folder, "RPVs.surf"),
		( surf_rings_folder, "SVC.surf"),
	]	

	for folder, fname in fnames:
		mycp(pjoin(folder, fname), pjoin(sims_folder, fname))

	for electrode in electrodes_paths_array:
		mycp(electrode, pjoin(sims_folder, os.path.basename(electrode)))

