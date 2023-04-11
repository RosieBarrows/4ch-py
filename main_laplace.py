import os
import sys

from common_4ch.file_utils import *

import argparse
import warnings

def main(args):

	warnings.warn("The following executables should be in your PATH: carp.pt, igbextract, GlVTKConvert.")

	print('=======================================================')
	print('Running endo-epi Laplace on '+args.meshname+'...')
	print('=======================================================')

	meshname = args.meshname
	endo = args.endo	
	epi = args.epi	

	print('=========================================')
	print('Converting surfaces to vtx files...')
	print('=========================================')
	endo_surf = read_elem(endo,el_type='Tr',tags=False)
	endo_vtx = surf2vtx(endo_surf)
	write_vtx(endo+".vtx",endo_vtx,init_row=2)	

	epi_surf = read_elem(epi,el_type='Tr',tags=False)
	epi_vtx = surf2vtx(epi_surf)
	write_vtx(epi+".vtx",epi_vtx,init_row=2)	

	print('=========================================')
	print('Running Laplace solution endo to epi...')
	print('=========================================')

	cmd = "carp.pt +F "+args.parfile+" -simID "+args.outdir+" -meshname "+meshname+" -stimulus[0].vtx_file "+endo+" -stimulus[1].vtx_file "+epi
	os.system(cmd)	

	print('=========================================')
	print('Converting igb to .dat file...')
	print('=========================================')

	cmd = "igbextract "+args.outdir+"/phie.igb -o ascii -f 0 -F 0 -O "+args.outdir+"/phie.dat"
	os.system(cmd)	

	print('=========================================')
	print('Converting to VTK to visualize...')
	print('=========================================')

	cmd = "GlVTKConvert -m "+meshname+" -n "+args.outdir+"/phie.dat -o "+args.outdir+"/laplace_endo_epi --trim-names" 
	os.system(cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--meshname', type=str, default=None,
                        help='Provide path to the carp mesh')

    parser.add_argument('--endo', type=str, default=None,
                        help='Provide path to the surf file for the endocardium')

    parser.add_argument('--epi', type=str, default=None,
                        help='Provide path to the surf file for the epicardium')

    parser.add_argument('--parfile', type=str, default="./parfiles/carpf_laplace_endo_epi.par",
                        help='Provide the par file for the Laplace solution')

    parser.add_argument('--outdir', type=str, default=None,
                        help='Provide the output directory for the Laplace solution')

    args = parser.parse_args()

    main(args)