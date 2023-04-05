import os
import sys
import json
import time
import argparse
# import pandas
import copy
import numpy as np
import meshio

def main(args):

  json_settings = args.map_settings
  f_input = open(json_settings,"r")
  settings = json.load(f_input)
  f_input.close()

  # read UVCs z and phi
  uvcs_z = np.loadtxt(args.uvcs+'/'+args.chamber+'.uvc_z.dat',dtype=float)
  uvcs_phi = np.loadtxt(args.uvcs+'/'+args.chamber+'.uvc_phi.dat',dtype=float)

  map_phi = np.zeros((uvcs_z.shape[0],))
  for i in range(uvcs_phi.shape[0]):
    if ((uvcs_phi[i] < settings[args.chamber]["begin_interp_aorta_side"]) or (uvcs_phi[i] > settings[args.chamber]["end_interp_not_aorta"])):
      map_phi[i] = 1.
    elif (uvcs_phi[i] > settings[args.chamber]["end_interp_aorta_side"]) and (uvcs_phi[i] < settings[args.chamber]["begin_interp_not_aorta"]):
      map_phi[i] = 0.
    elif (uvcs_phi[i] >= settings[args.chamber]["begin_interp_aorta_side"]) and (uvcs_phi[i] <= settings[args.chamber]["end_interp_aorta_side"]):
      map_phi[i] = ( uvcs_phi[i] - settings[args.chamber]["begin_interp_aorta_side"] )/( settings[args.chamber]["end_interp_aorta_side"] - settings[args.chamber]["begin_interp_aorta_side"] ) * ( 0.0-1.0 ) + 1.0
    elif (uvcs_phi[i] >= settings[args.chamber]["begin_interp_not_aorta"]) and (uvcs_phi[i] <= settings[args.chamber]["end_interp_not_aorta"]):
      map_phi[i] = ( uvcs_phi[i] - settings[args.chamber]["begin_interp_not_aorta"] )/( settings[args.chamber]["end_interp_not_aorta"] - settings[args.chamber]["begin_interp_not_aorta"] ) * ( 1.0-0.0 ) + 0.0

  map_z = np.zeros((uvcs_z.shape[0],))
  for i in range(uvcs_z.shape[0]):
    if (uvcs_z[i] < settings["Iz"]["Iz_1"]):
      map_z[i] = 1.
    elif (uvcs_z[i] > settings["Iz"]["Iz_2"]):
      map_z[i] = 0.
    elif (uvcs_z[i] >= settings["Iz"]["Iz_1"]) and (uvcs_z[i] <= settings["Iz"]["Iz_2"]):
      map_z[i] = ( uvcs_z[i] - settings["Iz"]["Iz_1"] )/( settings["Iz"]["Iz_2"] - settings["Iz"]["Iz_1"] ) * ( 0.0-1.0 ) + 1.0

  map_combined = np.multiply(map_phi,map_z)

  f=open(args.uvcs + '/map_rotational.dat','w')
  for m in map_phi:
    f.write(str(m)+'\n')
  f.close()

  f=open(args.uvcs + '/map_z.dat','w')
  for m in map_z:
    f.write(str(m)+'\n')
  f.close()

  f=open(args.uvcs + '/map_rotational_z.dat','w')
  for m in map_combined:
    f.write(str(m)+'\n')
  f.close()

  cmd = "GlVTKConvert -m "+args.mesh+"/"+args.chamber+" -n "+args.uvcs+"/map_rotational_z.dat -o "+args.uvcs+"/map_rotational"
  os.system(cmd)

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

  parser.add_argument('--chamber', type=str, default='la',
                      help='choose chamber to process (la or ra)')

  parser.add_argument('--mesh', type=str, default='',
                      help='give la or ra CARP mesh name')

  parser.add_argument('--uvcs', type=str, default=None,
                      help='path to the uvc files')

  parser.add_argument('--map_settings', type=str, default=None,
                      help='path to the atria_map_settings file')

  args = parser.parse_args()

  if not os.path.isdir(args.uvcs):
    raise RuntimeError('UVCs directory "{}" does not exist!'.format(args.uvcs))
  main(args)
