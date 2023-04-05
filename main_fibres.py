import argparse
import warnings
import os

def main(args):

	# os.system("clear")

	heartFolder = args.heartFolder
	mshPath = args.mshPath
	fchPath = args.fchPath
	fchName = args.fchName
	alpha_endo = args.alpha_endo
	alpha_epi = args.alpha_epi
	beta_endo = args.beta_endo
	beta_epi = args.beta_epi

	# ----------------------------------------------------------------------------------------------
	# Calculating fibres on the BiV mesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Calculating fibres on the BiV mesh ##")
	os.system("GlRuleFibers -m "+mshPath+"/BiV --type biv "
						  +"-a "+mshPath+"/uvc/BiV.sol_apba_lap.dat "
						  +"-e "+mshPath+"/uvc/BiV.sol_endoepi_lap.dat "
						  +"-l "+mshPath+"/uvc/BiV.sol_lvendo_lap.dat "
						  +"-r "+mshPath+"/uvc/BiV.sol_rvendo_lap.dat "
						  +"--alpha_endo "+alpha_endo+" "
						  +"--alpha_epi "+alpha_epi+" "
						  +"--beta_endo "+beta_endo+" "
						  +"--beta_epi "+beta_epi+" "
						  +"-o "+mshPath+"/fibres_bayer_"+alpha_endo+"_"+alpha_epi+".lon")

	# ----------------------------------------------------------------------------------------------
	# Substituted biv.lon with new fibres
	# ----------------------------------------------------------------------------------------------
	print(" ## Substituted biv.lon with new fibres ##")
	os.system("cp "+mshPath+"/fibres_bayer_"+alpha_endo+"_"+alpha_epi+".lon "+mshPath+"/BiV.lon")

	# ----------------------------------------------------------------------------------------------
	# Generating elements centres
	# ----------------------------------------------------------------------------------------------
	print(" ## Generating elements centres ##")
	os.system("GlElemCenters -m "+mshPath+"/BiV -o "+mshPath+"/BiV_elem_centres.pts")

	# ----------------------------------------------------------------------------------------------
	# Correcting fibre orientation
	# ----------------------------------------------------------------------------------------------
	print(" ## Correcting fibre orientation ##")
	os.system("./correct_fibres.py "+mshPath+"/BiV")

	# ----------------------------------------------------------------------------------------------
	# Substituted biv.lon with new fibres
	# ----------------------------------------------------------------------------------------------
	print(" ## Substituted biv.lon with new fibres ##")
	os.system("cp "+mshPath+"/BiV_corrected.lon "+mshPath+"/BiV.lon")

	# ----------------------------------------------------------------------------------------------
	# Converting for visualisation with Paraview
	# ----------------------------------------------------------------------------------------------
	print(" ## Converting for visualisation with Paraview ##")
	os.system("meshtool convert -imsh="+mshPath+"/BiV -ifmt=carp_txt -ofmt=vtk -omsh="+mshPath+"/BiV_fibres")
 
 	# ----------------------------------------------------------------------------------------------
	# Copying myocardium mesh into surfaces folder
	# ----------------------------------------------------------------------------------------------
	print(" ## Copying myocardium mesh into surfaces folder ##")
	os.system("cp "+heartFolder+"/meshing/myocardium_OUT/myocardium.* "+fchPath)

	# ----------------------------------------------------------------------------------------------
	# Inserting in submesh
	# ----------------------------------------------------------------------------------------------
	print(" ## Inserting in submesh ##")
	os.system("meshtool generate fibres -msh="+fchPath+"/"+fchName+" "
							 +"-op=2 "
							 +"-outmsh="+fchPath+"/"+fchName)

	os.system("meshtool insert submesh -submsh="+mshPath+"/BiV "
							 +"-msh="+fchPath+"/"+fchName+" "
							 +"-outmsh="+fchPath+"/"+fchName+"_bayer_"+alpha_endo+"_"+alpha_epi+" "
							 +"-ofmt=carp_txt")

	os.system("meshtool convert -imsh="+fchPath+"/"+fchName+"_bayer_"+alpha_endo+"_"+alpha_epi+" "
					  +"-omsh="+fchPath+"/"+fchName+"_bayer_"+alpha_endo+"_"+alpha_epi+" "
					  +"-ofmt=vtk_bin")
 


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser.add_argument('--heartFolder', type=str, default=None,
                        help='Provide path to the heart folder')

    parser.add_argument('--mshPath', type=str, default=None,
                        help='Provide path to the mesh to which fibres will be added')

    parser.add_argument('--fchPath', type=str, default=None,
                        help='?')

    parser.add_argument('--fchName', type=str, default=None,
                        help='Provide name of the four chamber mesh')

    parser.add_argument('--CARPFOLDER', type=str, default=None,
                        help='Provide path to CARP')

    parser.add_argument('--alpha_endo', type=str, default=None,
                        help='Provide angle of fibres on endo')

    parser.add_argument('--alpha_epi', type=str, default=None,
                        help='Provide angle of fibres on epi')

    parser.add_argument('--beta_endo', type=str, default=None,
                        help='Provide angle of sheet direction on endo')

    parser.add_argument('--beta_epi', type=str, default=None,
                        help='Provide angle of sheet direction on epi')

    args = parser.parse_args()

    main(args)
