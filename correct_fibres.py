#!/usr/bin/env python2

import os
import sys
import numpy as np
import copy
# from tqdm import tqdm

mshName=sys.argv[1]

print('Reading mesh...')
lon = np.loadtxt(mshName+'.lon',dtype=float,skiprows=1)
pts = np.loadtxt(mshName+'.pts',dtype=float,skiprows=1)
elem = np.loadtxt(mshName+'.elem',dtype=int,skiprows=1,usecols=[1,2,3,4])
print('Done...')

print('Reading element centres...')
elemC = np.loadtxt(mshName+'_elem_centres.pts',dtype=float,skiprows=1)
print('Done...')

to_correct = np.where(np.abs(lon[:,2])==0)[0]
print('Found '+str(to_correct.shape[0])+' elements to correct')

lon_corrected = copy.deepcopy(lon)
good_elements = np.setdiff1d(np.arange(elemC.shape[0]), to_correct, assume_unique=True)

# for idx in tqdm(to_correct, desc="Correcting fibre orientation..."):
for idx in to_correct:
	d = np.linalg.norm(elemC[good_elements,:] - elemC[idx,:],axis=1)
	closest = good_elements[np.where(d == np.min(d))[0]]
	lon_corrected[idx,:] = lon[closest,:]
to_correct = np.where(np.abs(lon_corrected[:,2])==0)[0]

print('Left to correct '+str(to_correct.shape[0]))

np.savetxt(mshName+'_corrected.lon',lon_corrected,fmt="%g",header='2',comments='')
