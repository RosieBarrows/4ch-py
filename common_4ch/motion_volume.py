import numpy as np
import copy
import meshio
import os

um3_to_mL = 1e-12

def read_pts(ptsname):
	return np.loadtxt(ptsname,dtype=float,skiprows=1)

def read_tets(elemname):
	return np.loadtxt(elemname,dtype=int,usecols=[1,2,3,4,5],skiprows=1)

def write_tets(filename,elem):
    with open(filename, 'w') as fp:
        fp.write('{}\n'.format(elem.shape[0]))
        for el in elem:
            fp.write('Tt {} {} {} {} {}\n'.format(int(el[0]),int(el[1]),int(el[2]),int(el[3]),int(el[4])))
	
def write_surf(surf,surf_file):
	f = open(surf_file,"w")
	f.write(str(int(surf.shape[0]))+"\n")
	for t in surf:
		f.write("Tr "+str(t[0])+" "+str(t[1])+" "+str(t[2])+"\n")
	f.close()

def write_pts(pts,pts_file):
	f = open(pts_file,"w")
	f.write(str(int(pts.shape[0]))+"\n")
	for p in pts:
		f.write(str(p[0])+" "+str(p[1])+" "+str(p[2])+"\n")
	f.close()

def neubc2surf(neubc_file,surf_file):
	surf = np.loadtxt(neubc_file,dtype=int,usecols=[0,1,2],skiprows=1)
	write_surf(surf,surf_file)

def track_mesh(dofFile,refMesh,nFrames,outputFolder):

	frames = np.arange(0,100,int(100/nFrames),dtype=int)
	for i,n in enumerate(frames):
		cmd="mirtk-transform-points "+refMesh+" "+outputFolder+"/transformed-"+str(i)+".vtk -dofin "+dofFile+" -St "+str(n)+" -binary"
		os.system(cmd)

def transform_mesh_nii(mshName,mshoutName,transform_vector):
	msh = meshio.read(mshName)
	msh.points = transform_vector*msh.points

	msh = meshio.write(mshoutName,msh,file_format="vtk42")

def meshtool_vtkASCII(filename):
	cmd="meshtool convert -imsh="+filename+" -omsh="+filename+" -ifmt=vtk_bin -ofmt=vtk"
	os.system(cmd)

def read_surf(surfname):
	return np.loadtxt(surfname,dtype=int,usecols=[1,2,3],skiprows=1)

def read_neubc(surfname):
	return [np.loadtxt(surfname,dtype=int,usecols=[0,1,2],skiprows=1),np.loadtxt(surfname,dtype=int,usecols=[3],skiprows=1)]

def surf2vtx(surf):
	return np.unique(surf.flatten())

def surf2vtk_tags(msh,surf,tags_array,output_name):

	vtx = surf2vtx(surf)
	pts_surf = msh.points[vtx,:]

	surf_reindexed = copy.deepcopy(surf)
	for i in range(surf.shape[0]):
		surf_reindexed[i,0] = np.where(vtx==surf[i,0])[0]
		surf_reindexed[i,1] = np.where(vtx==surf[i,1])[0]
		surf_reindexed[i,2] = np.where(vtx==surf[i,2])[0]

	cells = {"triangle": surf_reindexed}
	cell_data = {"elemTag":[tags_array]}
	surf_vtk_msh = meshio.Mesh(
    				pts_surf,
    				cells,
    				cell_data=cell_data
				)

	meshio.write(output_name, surf_vtk_msh,file_format="vtu")

def compute_enclosed_volume(msh_name,surf_name):

	print("Reading mesh...")
	pts = read_pts(msh_name+".pts")
	print("Done.")

	surf = read_surf(surf_name)

	v02 = np.array([pts[surf[:,0],0]-pts[surf[:,2],0],
					pts[surf[:,0],1]-pts[surf[:,2],1],
					pts[surf[:,0],2]-pts[surf[:,2],2]],dtype=float)
	v01 = np.array([pts[surf[:,0],0]-pts[surf[:,1],0],
					pts[surf[:,0],1]-pts[surf[:,1],1],
					pts[surf[:,0],2]-pts[surf[:,1],2]],dtype=float)
	cr = np.cross(v02,v01,axisa=0,axisb=0)
	area = 0.5*np.sqrt(np.power(cr[:,0],2)+np.power(cr[:,1],2)+np.power(cr[:,2],2))
	totalArea = np.sum(area)

	crNorm = np.sqrt(np.power(cr[:,0],2)+np.power(cr[:,1],2)+np.power(cr[:,2],2))
	zMean = (pts[surf[:,0],2]+pts[surf[:,1],2]+pts[surf[:,2],2])/3.;
	nz = -cr[:,2]/crNorm
	volume = area*zMean*nz

	return sum(volume)*um3_to_mL

def compute_volume_endo(endo_vtkBasename,nFrames,tags,output_file,transform=0.001,header="time LV RV LA RA"):

	volume_transient = np.zeros((nFrames,len(tags)+1),dtype=float)
	for n in range(nFrames):
		print("Reading "+endo_vtkBasename+str(n)+".vtu")

		endoMsh = meshio.read(endo_vtkBasename+str(n)+".vtu",file_format="vtu")
		cells = endoMsh.get_cells_type("triangle")
		pts = endoMsh.points
		tags_array = endoMsh.cell_data["elemTag"][0]

		for i,t in enumerate(tags):
			cells_endo = cells[np.where(tags_array==t)[0],:]
			v02 = np.array([pts[cells_endo[:,0],0]-pts[cells_endo[:,2],0],
							pts[cells_endo[:,0],1]-pts[cells_endo[:,2],1],
							pts[cells_endo[:,0],2]-pts[cells_endo[:,2],2]],dtype=float)
			v01 = np.array([pts[cells_endo[:,0],0]-pts[cells_endo[:,1],0],
							pts[cells_endo[:,0],1]-pts[cells_endo[:,1],1],
							pts[cells_endo[:,0],2]-pts[cells_endo[:,1],2]],dtype=float)
			cr = np.cross(v02,v01,axisa=0,axisb=0)
			area = 0.5*np.sqrt(np.power(cr[:,0],2)+np.power(cr[:,1],2)+np.power(cr[:,2],2))
			totalArea = np.sum(area)

			crNorm = np.sqrt(np.power(cr[:,0],2)+np.power(cr[:,1],2)+np.power(cr[:,2],2))
			zMean = (pts[cells_endo[:,0],2]+pts[cells_endo[:,1],2]+pts[cells_endo[:,2],2])/3.;
			nz = -cr[:,2]/crNorm
			volume = area*zMean*nz
			volume_transient[n,i+1] = sum(volume)
	volume_transient[:,0] = range(nFrames)
	np.savetxt(output_file,volume_transient*transform,fmt="%.1f",header=header)


def tracked_to_volume(vtkBasename,nFrames,surfaces,tags,outBasename):

	surf_dict = {}
	nT = 0
	for s in surfaces:
		surf_dict[s] = read_surf(s)
		nT += surf_dict[s].shape[0]

	surf_array = np.zeros((nT,3),dtype=int)
	tags_array = np.zeros((nT,),dtype=int)
	shift = 0
	for i,s in enumerate(surfaces):
		surf_array[shift:shift+surf_dict[s].shape[0],:] = surf_dict[s]
		tags_array[shift:shift+surf_dict[s].shape[0]] = tags[i]
		shift += surf_dict[s].shape[0]
	surf_vtx = surf2vtx(surf_array)

	msh0 = meshio.read(vtkBasename+"0.vtk")
	surf2vtk_tags(msh0,surf_array,tags_array,outBasename+"0.vtu")
	surf0 = meshio.read(outBasename+"0.vtu")			
	for n in range(nFrames):
		print('Extracting surface for frame '+str(n)+'...')
		filename = vtkBasename+str(n)
		# meshtool_vtkASCII(filename)
		msh = meshio.read(filename+".vtk")

		pts_surf = msh.points[surf_vtx,:]
		surf0.points = pts_surf
		meshio.write(outBasename+str(n)+".vtu", surf0,file_format="vtu")


		
