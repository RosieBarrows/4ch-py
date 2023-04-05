import numpy as np
from common.file_utils import *
import os
import sys
import copy
import math
import scipy.spatial
import h5py

def LIFEX_to_CARP(msh_name,
				  carp_meshname):
	
	with h5py.File(msh_name+".h5", "r") as f:

		pts_lifex = f['nodes'][()]*1000.0
		f0 = f['f0'][()]
		s0 = f['s0'][()]

	f.close()

	print('Correcting NaN values...')
	idx_notnan = np.where(~np.isnan(f0[:,0]))[0]
	idx_nan = np.where(np.isnan(f0[:,0]))[0]
	kdt_pts_lifex = scipy.spatial.cKDTree(pts_lifex[idx_notnan,:])
	nearest_dist_pts, nearest_ind_pts = kdt_pts_lifex.query(pts_lifex[idx_nan,:], k=1)
	for i,idx in enumerate(idx_nan):
		f0[idx,:] = f0[idx_notnan[nearest_ind_pts[i]],:]
		s0[idx,:] = s0[idx_notnan[nearest_ind_pts[i]],:]

	print('Mapping fibres and sheet directions from nodes to elements...')

	cmd="GlElemCenters -m "+carp_meshname+" -f 0 -p 1 -o "+carp_meshname+"_elCentres.vpts"
	os.system(cmd)

	pts = read_pts(carp_meshname+".pts")
	elem_centres = read_pts(carp_meshname+"_elCentres.vpts")

	print('Mapping nodes order from lifex to carp...')
	kdt_pts_lifex = scipy.spatial.cKDTree(pts_lifex)
	nearest_dist_pts, nearest_ind_pts = kdt_pts_lifex.query(pts, k=1)
	f0_remapped = np.zeros(f0.shape,dtype=float)
	s0_remapped = np.zeros(s0.shape,dtype=float)
	pts_test = np.zeros(pts.shape,dtype=float)
	for j in range(f0.shape[0]):
		pts_test[j,:] = pts_lifex[nearest_ind_pts[j],:]
		f0_remapped[j,:] = f0[nearest_ind_pts[j],:]
		s0_remapped[j,:] = s0[nearest_ind_pts[j],:]

	fibres = np.concatenate((f0_remapped,s0_remapped),axis=1)

	kdt_pts = scipy.spatial.cKDTree(pts)
	nearest_dist, nearest_ind = kdt_pts.query(elem_centres, k=1)

	fibres_elem = np.zeros((elem_centres.shape[0],6),dtype=float)
	for i in range(fibres_elem.shape[0]):
		fibres_elem[i,:] = fibres[nearest_ind[i],:] 

	write_lon(fibres_elem,carp_meshname+"_fibres.lon")
	os.system("cp "+carp_meshname+".pts "+carp_meshname+"_fibres.pts")
	os.system("cp "+carp_meshname+".elem "+carp_meshname+"_fibres.elem")

def define_stimulus(msh_name,
					stimulus_centre,
					stim_name,
					stim_radius=1500.0):

	stimulus_centre = np.loadtxt(stimulus_centre,dtype=int)
	pts = read_pts(msh_name+".pts")
	kdt_pts = scipy.spatial.cKDTree(pts)	

	vtx_stim = kdt_pts.query_ball_point(pts[stimulus_centre], stim_radius)

	write_vtx(stim_name,np.array(vtx_stim))

def combine_la_ra_biatrial(la_mshname,
						   ra_mshname,
						   biatrial_mshname,
						   output_mshname):

	print("Combining fibres la and ra into biatrial mesh...")

	cmd="meshtool insert meshdata -msh="+biatrial_mshname+" -imsh="+la_mshname+" -op=1 -outmsh="+output_mshname
	os.system(cmd)

	cmd="meshtool insert meshdata -msh="+output_mshname+" -imsh="+ra_mshname+" -op=1 -outmsh="+output_mshname
	os.system(cmd)

def write_simulation_setup(labels,
						   tags,
						   CV,
						   anisotropy,
						   output_file):

	if (len(labels)!=len(tags)) or len(labels)!=len(CV) or len(labels)!=len(anisotropy):
		raise Exception("The length of the supplied inputs does not match.")

	dct = {}
	for i,l in enumerate(labels):
		dct[l] = {}
		dct[l]["tags"] = tags[i]
		dct[l]["CV"] = CV[i]
		dct[l]["anisotropy"] = anisotropy[i]

	save_json(dct,output_file)

def get_nod_file(msh_name,
				 submsh_name,
				 tags_list,
				 tag_file):
	print('Using meshtool to get .nod file...')

	tags_list_str = []
	for t in tags_list:
		tags_list_str.append(str(t))
	tag_str = ','.join(tags_list_str)

	cmd = "meshtool extract mesh -msh="+msh_name+" "+"-submsh="+submsh_name+" -tags="+tag_str+" -tag_file="+tag_file+" -norm"
	print(cmd)
	os.system(cmd)

def check_init_stim_folder(init_folder,
	    			  	    msh_name,
	    			        nod_file,
	    			        results_folder):
	
	files = os.listdir(init_folder)
	init_files = []
	results_files = []
	filenames = []
	for f in files:
		if f[-5:] == '.init':
			print(f)
			filenames.append(f[:-5])
			init_files.append(init_folder+'/'+f[:-5])
			results_files.append(results_folder+'/'+f[:-5])
	print('Found '+str(len(init_files))+' init files to check.')

	for i,f in enumerate(init_files):
		check_init_stim(f+'.init',
				  	    msh_name,
				        nod_file,
				        results_files[i]+'.dat')


def check_init_stim(init_file,
			  	    msh_name,
			        nod_file,
			        act_file):

	print('Checking file '+init_file+'...')

	count = 0

	file1 = open(init_file, 'r')
	lines = file1.readlines()
	count = 0

	vtx = []
	for line in lines:
		if count==2:
			line_txt =line.strip()
			line_txt_list = line_txt.split(' ')
			N_vtx = int(line_txt_list[0])
			print('Checking '+line_txt_list[0]+' node(s)')

		if count>2 and count<=N_vtx+2:
			line_txt =line.strip()
			line_txt_list = line_txt.split(' ')
			vtx.append(int(line_txt_list[0]))

		count += 1

	nod = read_nod(nod_file)
	stim = np.intersect1d(nod,vtx)

	if len(stim)==0:
		print('No stimulus in domain. Writing fake activation file...')
		pts = read_pts(msh_name+'.pts')
		N_pts = pts.shape[0]
		act = -np.ones((N_pts,),dtype=float)
		np.savetxt(act_file,act,fmt="%.1f")
	else:
		print('Found stimulus in domain.')
	

def CARPtoInit(setup_file,
			   stim_time,
			   stimulus,
			   init_file,
			   vPS=4.0,
			   aDelay=10.0,
			   rDelay=3.0):

	setup = load_json(setup_file)
	CV = []
	tags = []
	anisotropy = []
	for l in setup:
		tags_tmp = setup[l]["tags"]
		for t in tags_tmp:
			tags.append(t)
			CV.append(setup[l]["CV"])
			anisotropy.append(setup[l]["anisotropy"])

	vtx = []
	stim_times = []
	for i,vtxFile in enumerate(stimulus):
		temp = read_vtx(vtxFile)
		if len(temp.shape)==0:
			vtx.append(temp)
			stim_times.append(stim_time[i])
		else:
			for v in temp:
				vtx.append(v)
				stim_times.append(stim_time[i])

	# ----------------------------------------------------------------------------------- #
	f = open(init_file,'w')

	f.write('vf:1.000000 vs:1.000000 vn:1.000000 vPS:%f\n' % vPS)
	f.write('retro_delay:%f antero_delay:%f\n' % (rDelay,aDelay))

	f.write('%d %d\n' % (int(len(vtx)), int(len(tags))))

	for i,v in enumerate(vtx):
		f.write('%d %f\n' % (int(v),stim_times[i]))

	for i,t in enumerate(tags):
		f.write('%d %f %f %f\n' % (int(t),CV[i],CV[i]*anisotropy[i],CV[i]*anisotropy[i]))

	f.close()

def get_tag_list(setup_file):

	setup = load_json(setup_file)
	tags = []
	for l in setup:
		tags_tmp = setup[l]["tags"]
		for t in tags_tmp:
			tags.append(t)
	return tags

def run_ekbatch(mesh_name,
			    init_folder,
			    results_folder,
				tag_file=None,
				tags=None,
				overwrite=False):

	print("WARNING: make sure ekbatch is in your searchpath.")

	files = os.listdir(init_folder)
	init_files = []
	filenames = []
	for f in files:
		if f[-5:] == '.init':
			print(f)
			filenames.append(f[:-5])
			init_files.append(init_folder+'/'+f[:-5])
	print('Found '+str(len(init_files))+' init files.')

	if not overwrite:
		print('You chose not to overwrite... Removing simulations you have run already. Left:')
		init_files_not_exist = []
		for f in filenames:
			if not os.path.exists(results_folder+'/'+f+'.dat'):
				print(f)
				init_files_not_exist.append(init_folder+'/'+f)
		init_files = init_files_not_exist

	if len(init_files)>0:
		init_files_list = ','.join(init_files)	

		cmd = "ekbatch "+mesh_name+" "+init_files_list	

		if tags is not None:
			print('Using tags:')
			tags_list_str = []
			for t in tags:
				print(t)
				tags_list_str.append(str(t))
			tag_str = ','.join(tags_list_str)
			cmd += " "+tag_str	

		if tag_file is not None:
			print('Using tags file '+tag_file)
			cmd += " "+tag_file
		
		os.system(cmd)	

		os.system("mkdir "+results_folder)
		cmd = "mv "+init_folder+'/*.dat '+results_folder
		os.system(cmd)
	else:
		print('It looks like all simulations were run. If you want to overwrite, set overwrite=True.')