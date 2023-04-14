import logging
import numpy as np
import copy
import os
from tqdm import tqdm

def my_log(funcname, msg) : 
    print(f"[{funcname}] : {msg}")

def correct_fibres(meshname):
    """Corrects fibre orientation in the mesh"""
    fun_name = "correct_fibres"
    my_log(fun_name, 'Reading mesh...')
    lon = np.loadtxt(f"{meshname}.lon",dtype=float,skiprows=1)
    # pts = np.loadtxt(f"{meshname}.pts",dtype=float,skiprows=1)
    # elem = np.loadtxt(f"{meshname}.elem",dtype=int,skiprows=1,usecols=[1,2,3,4])

    my_log(fun_name, 'Reading element centres...')
    elemC = np.loadtxt(f"{meshname}+_elem_centres.pts",dtype=float,skiprows=1)

    to_correct = np.where(np.abs(lon[:,2])<1e-6)[0]
    my_log(fun_name, f"Found: {str(to_correct.shape[0])} elements to correct")

    lon_corrected = copy.deepcopy(lon)
    good_elements = np.setdiff1d(np.arange(elemC.shape[0]), to_correct, assume_unique=True)

    for idx in tqdm(to_correct, desc="Correcting fibre orientation..."):
        d = np.linalg.norm(elemC[good_elements,:] - elemC[idx,:],axis=1)
        closest = good_elements[np.where(d == np.min(d))[0]]
        lon_corrected[idx,:] = lon[closest,:]
    to_correct = np.where(np.abs(lon_corrected[:,2])<1e-6)[0]

    my_log(fun_name, f"Left to correct {str(to_correct.shape[0])}")
    my_log(fun_name, 'Correcting fibre orientation...')

    np.savetxt(f"{meshname}_corrected.lon",lon_corrected,fmt="%g",header='2',comments='')

