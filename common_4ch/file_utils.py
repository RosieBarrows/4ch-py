import numpy as np
import copy
import json
import meshio
import vtk
import os
from vtk.util import numpy_support

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
	cmd = f"cp {src} {dst}"
	if debug:
		print(cmd)
	os.system(cmd)

def mymv(src,dst, debug=False):
	cmd = f"mv {src} {dst}"
	if debug:
		print(cmd)
	os.system(cmd)

def mymkdir(path, full_path=False, debug=False):
	# check if path exists
	if os.path.exists(path):
		print(f"Path [{path}] already exists")
		return
	
	flag = "-p" if full_path else ""
	cmd = f"mkdir {flag} {path}"
	if debug:
		print(cmd)
	os.system(cmd)

def myrm(path, debug=False):
	cmd = f"rm -rf {path}"
	if debug:
		print(cmd)
	os.system(cmd)

def big_msg(msg):
	length = len(msg)
	print("-"*length)
	print(msg)
	print("-"*length)

def pjoin(*paths) :
	return os.path.join(*paths)
	
def read_pts(filename):
	print('Reading '+filename+'...')
	return np.loadtxt(filename, dtype=float, skiprows=1)

def read_elem(filename,el_type='Tt',tags=True):
	print('Reading '+filename+'...')

	if el_type=='Tt':
		if tags:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3,4,5))
		else:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3,4))
	elif el_type=='Tr':
		if tags:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3,4))
		else:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3))
	elif el_type=='Ln':
		if tags:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3))
		else:
			return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2))
	else:
		raise Exception('element type not recognised. Accepted: Tt, Tr, Ln')

def read_lon(filename):
	print('Reading '+filename+'...')

	return np.loadtxt(filename, dtype=float, skiprows=1)

def read_surf_polydata(filename):

	print('Reading '+filename+'...')

	polydata_reader = vtk.vtkPolyDataReader()

	polydata_reader.SetFileName(filename)
	polydata_reader.Update()

	polydata = polydata_reader.GetOutput()

	points = numpy_support.vtk_to_numpy(polydata.GetPoints().GetData())

	faces = numpy_support.vtk_to_numpy(polydata.GetPolys().GetData())
	faces = faces.reshape((int(faces.shape[0]/4),4))
	faces = faces[:,1:]

	tags = numpy_support.vtk_to_numpy(polydata.GetCellData().GetArray('elemTag'))

	print('Done.')

	return points,faces,tags

def read_clipper(filename):

	print('Reading center and radius for clipper '+filename+'...')

	polydata_reader = vtk.vtkPolyDataReader()

	polydata_reader.SetFileName(filename)
	polydata_reader.Update()

	polydata = polydata_reader.GetOutput()

	points = numpy_support.vtk_to_numpy(polydata.GetPoints().GetData())
	center = np.mean(points,axis=0)
	radius = np.linalg.norm(points[0,:]-center)

	print('Done.')

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
	print('Writing '+filename+'...')
	print("Check you are using the correct write_surf function")

	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(surf.shape[0]))
		for t in surf:
			fp.write('Tr {} {} {}\n'.format(int(t[0]),int(t[1]),int(t[2])))

def write_pts_caroline(pts,filename):

	print('Writing '+filename+'...')

	assert pts.shape[1] == 3
	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(pts.shape[0]))
		for pnt in pts:
			fp.write('{0[0]} {0[1]} {0[2]}\n'.format(pnt))
	fp.close()

def write_elem_caroline(elem,
			   tags,
			   filename,
			   el_type='Tt'):
	print('Writing '+filename+'...')

	if el_type=='Tt':
		assert elem.shape[1] == 4
	elif el_type=='Tr':
		assert elem.shape[1] == 3
	elif el_type=='Ln':
		assert elem.shape[1] == 2
	else:
		raise Exception('element type not recognised. Accepted: Tt, Tr, Ln')

	assert elem.shape[0] == tags.shape[0]

	with open(filename, 'w') as fe:
		fe.write('{}\n'.format(elem.shape[0]))
		for i,el in enumerate(elem):
			fe.write(el_type)
			for e in el:
				fe.write(' '+str(e))
			fe.write(' '+str(tags[i]))
			fe.write('\n')
	fe.close()

def write_tets_ln(filename,
				  elem,
				  ln):

	print('Writing '+filename+'...')

	with open(filename, 'w') as fp:
		fp.write('{}\n'.format(elem.shape[0] + ln.shape[0]))
		for el in elem:
			fp.write('Tt {} {} {} {} {}\n'.format(int(el[0]),int(el[1]),int(el[2]),int(el[3]),int(el[4])))
		for cn in ln:
			fp.write('Ln {} {} {}\n'.format(int(cn[0]),int(cn[1]),int(cn[2])))

def write_lon(lon,filename):
	print('Writing '+filename+'...')

	assert lon.shape[1] % 3 == 0
	with open(filename, 'w') as fl:
		fl.write('{}\n'.format(int(lon.shape[1]/3)))
		for ll in lon:
			for i,l in enumerate(ll):
				if i==len(ll)-1:
					fl.write(str(l)+'\n')
				else:
					fl.write(str(l)+' ')
	fl.close()

def write_vtx(filename, vtx,init_row=2):
	print('Writing '+filename+'...')

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
	print('Reading '+filename+'...')

	dct = {}
	with open(filename, "r") as f:
		dct = json.load(f, object_hook=numpy_hook)
	return dct

def read_vtx(filename,init_row=2):
	print('Reading '+filename+'...')

	return np.loadtxt(filename, dtype=int, skiprows=init_row)


def setup_sim(heartFolder,presimFolder,electrodes_paths_array):
	os.system("mkdir "+heartFolder+"/sims_folder")

	os.system("cp "+presimFolder+"/elem_dat_UVC_ek_combined.dat "+heartFolder+"/sims_folder/pericardium_scale.dat")
	os.system("cp "+presimFolder+"/epicardium_for_sim.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/LA_endo.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/LV_endo.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/RA_endo.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/RV_endo.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/surfaces_rings/RPVs.surf.vtx "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/surfaces_rings/SVC.surf "+heartFolder+"/sims_folder")
	os.system("cp "+presimFolder+"/surfaces_simulation/surfaces_rings/SVC.surf.vtx "+heartFolder+"/sims_folder")

	for electrode in electrodes_paths_array:
		os.system("cp "+electrode + " " +heartFolder+"/sims_folder")

