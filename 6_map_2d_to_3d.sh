#!/bin/bash

heart_folder=$(cat /data/Dropbox/4ch-py/parfiles/heartFolder.txt)

FCH="${heart_folder}/surfaces_uvc/myocardium_bayer_60_-60"
UACFOLDER="${heart_folder}/atrial_fibres/UAC/"

# ---------------------------------------------------------------------------------

CMD="python3 main_laplace.py --meshname ${UACFOLDER}/la/la 
							--endo ${UACFOLDER}/la/la_endo.surf
							--epi ${UACFOLDER}/la/la_epi.surf
							--outdir ${UACFOLDER}/la/endo_epi/"
eval $CMD

CMD="python3 main_laplace.py --meshname ${UACFOLDER}/ra/ra 
							--endo ${UACFOLDER}/ra/ra_endo.surf
							--epi ${UACFOLDER}/ra/ra_epi.surf
							--outdir ${UACFOLDER}/ra/endo_epi/"
eval $CMD

# ---------------------------------------------------------------------------------

CMD="python3 main_surf_to_volume.py --meshname ${UACFOLDER}/la/la
								   --meshname_uac ${UACFOLDER}/LA_endo/Fibre_endo_l
								   --endo_fibres ${UACFOLDER}/LA_endo/Fibre_endo_l.lon
								   --epi_fibres ${UACFOLDER}/LA_endo/Fibre_epi_l.lon
								   --endo_epi_laplace ${UACFOLDER}/la/endo_epi/phie.dat
								   --outmeshname ${UACFOLDER}/la/la_fibres_l"
eval $CMD

# ---------------------------------------------------------------------------------

CMD="python3 main_surf_to_volume.py --meshname ${UACFOLDER}/ra/ra
								   --meshname_uac ${UACFOLDER}/RA_endo/Fibre_endo_l
								   --endo_fibres ${UACFOLDER}/RA_endo/Fibre_endo_l.lon
								   --epi_fibres ${UACFOLDER}/RA_endo/Fibre_epi_l.lon
								   --endo_epi_laplace ${UACFOLDER}/ra/endo_epi/phie.dat
								   --outmeshname ${UACFOLDER}/ra/ra_fibres_l"
eval $CMD

# ---------------------------------------------------------------------------------

# CMD="python main_surf_to_volume_tags.py --meshname ${UACFOLDER}/la/la_fibres_l_sheet
# 								   --meshname_uac ${UACFOLDER}/LA_endo/Fibre_endo_l
# 								   --tagfiles ${UACFOLDER}/LA_endo/MappedScalar_BB_LA.dat
# 								   --tag_labels BB
# 								   --default_tag_label LA
# 								   --endo_epi_laplace ${UACFOLDER}/la/endo_epi/phie.dat
# 								   --outmeshname ${UACFOLDER}/la/la_fibres_l_tags"
# eval $CMD

# CMD="python main_surf_to_volume_tags.py --meshname ${UACFOLDER}/ra/ra_fibres_l_sheet
# 								   --meshname_uac ${UACFOLDER}/RA_endo/Fibre_endo_l
# 								   --tagfiles ${UACFOLDER}/RA_endo/MappedScalar_BB.dat ${UACFOLDER}/RA_endo/MappedScalar_CT.dat ${UACFOLDER}/RA_endo/MappedScalar_PM.dat
# 								   --tag_labels BB CT PM
# 								   --tag_labels_endo CT PM
# 								   --default_tag_label RA
# 								   --endo_epi_laplace ${UACFOLDER}/ra/endo_epi/phie.dat
# 								   --outmeshname ${UACFOLDER}/ra/ra_fibres_l_tags"
# eval $CMD

# ---------------------------------------------------------------------------------

echo "Mapping to biatrial mesh..."

CMD="cp ${UACFOLDER}/la/la.nod ${UACFOLDER}/la/la_fibres_l_sheet.nod;cp ${UACFOLDER}/la/la.eidx ${UACFOLDER}/la/la_fibres_l_sheet.eidx"
eval $CMD

CMD="cp ${UACFOLDER}/ra/ra.nod ${UACFOLDER}/ra/ra_fibres_l_sheet.nod;cp ${UACFOLDER}/ra/ra.eidx ${UACFOLDER}/ra/ra_fibres_l_sheet.eidx"
eval $CMD

CMD="meshtool insert submesh -submsh=${UACFOLDER}/la/la_fibres_l_sheet -msh=${UACFOLDER}/biatrial/biatrial -ofmt=carp_txt -outmsh=${UACFOLDER}/biatrial/biatrial_la_fibres"
eval $CMD

CMD="meshtool convert -imsh=${UACFOLDER}/biatrial/biatrial_la_fibres -omsh=${UACFOLDER}/biatrial/biatrial_la_fibres.vtk"
eval $CMD

CMD="meshtool insert submesh -submsh=${UACFOLDER}/ra/ra_fibres_l_sheet -msh=${UACFOLDER}/biatrial/biatrial_la_fibres -ofmt=carp_txt -outmsh=${UACFOLDER}/biatrial/biatrial_fibres_l"
eval $CMD

CMD="meshtool convert -imsh=${UACFOLDER}/biatrial/biatrial_fibres_l -omsh=${UACFOLDER}/biatrial/biatrial_fibres_l.vtk"
eval $CMD

rm ${UACFOLDER}/biatrial/biatrial_la_fibres*

echo "Mapping to four-chamber mesh..."

CMD="cp ${UACFOLDER}/biatrial/biatrial.nod ${UACFOLDER}/biatrial/biatrial_fibres_l.nod;cp ${UACFOLDER}/biatrial/biatrial.eidx ${UACFOLDER}/biatrial/biatrial_fibres_l.eidx"
eval $CMD

CMD="meshtool insert submesh -submsh=${UACFOLDER}/biatrial/biatrial_fibres_l -msh=${FCH} -ofmt=carp_txt -outmsh=${heart_folder}/atrial_fibres/myocardium_fibres_l"
eval $CMD

CMD="meshtool convert -imsh=${heart_folder}/atrial_fibres/myocardium_fibres_l -omsh=${heart_folder}/atrial_fibres/myocardium_fibres_l.vtk"
eval $CMD