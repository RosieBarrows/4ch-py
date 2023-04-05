#!/bin/bash

# =============================================================================
# SMOOTHED HEART FOUR CHAMBERS
# (c) images          Steven Niederer
# (c) segmentation    SN
# (c) meshing         AJP
# (c) fibers & sheets JB, AJP
# =============================================================================


# ---SEGMENTATION LABELS-------------------------------------------------------
# #
T_UNUSED=200             # unused tag
# T_BATH=0                # default bath tag
#
# ---Left atrium---------------------------------------------------------------
T_LA=200                 # left atrial wall
# T_LARA=22               # combined left/right atrial wall
T_LABP=200               # left atrial blood pool

# T_LINFPULMVEIN=         # left inferior pulmonary vein
T_LINFPULMVEINCUT=200    # left inferior pulmonary vein (cut)
# T_LINFPULMVEINBP=34     # left inferior pulmonary vein blood pool
# T_LINFPULMVEINBPCUT=35  # left inferior pulmonary vein blood pool (cut)

# T_LSUPPULMVEIN=11         # left superior pulmonary vein
T_LSUPPULMVEINCUT=200    # left superior pulmonary vein (cut)
# T_LSUPPULMVEINBP=36     # left superior pulmonary vein blood pool
# T_LSUPPULMVEINBPCUT=37  # left superior pulmonary vein blood pool (cut)

# T_RINFPULMVEIN=8         # right inferior pulmonary vein
T_RINFPULMVEINCUT=200      # right inferior pulmonary vein (cut)
# T_RINFPULMVEINBP=38     # right inferior pulmonary vein blood pool
# T_RINFPULMVEINBPCUT=    # right inferior pulmonary vein blood pool (cut)

# T_RSUPPULMVEIN=         # right superior pulmonary vein
T_RSUPPULMVEINCUT=200      # right superior pulmonary vein (cut)
# T_RSUPPULMVEINBP=39     # right superior pulmonary vein blood pool
# T_RSUPPULMVEINBPCUT=    # right superior pulmonary vein blood pool (cut)
#
# ---Right atrium--------------------------------------------------------------
T_RA=200                 # right atrial wall
T_RABP=200               # right atrial blood pool
# T_CORONARYSINUS=        # coronary sinus (closing hole in RA appendage)
# T_SINUSNODE=            # sinus node
# T_AVNODE=               # AV node
# T_BACHMBUNDLES=         # Bachmann's bundles
# T_PECTINATES=           # pectinate muscles
# T_FO=                   # fossa ovalis
# T_FORIM=                # rim of fossa ovalis
# T_LOB=                  # line of block
#
# ---Left ventricle
T_LV=3                  # left ventricular tissue
# T_LVENDO=61             # left ventricular tissue (trabeculae)
T_LVBP=200              # left ventricular blood pool
T_AORTA=200              # aorta
T_AORTABP=200            # aortic blood pool
T_MITRALVV=200           # mitral valve
T_AORTICVV=200          # aortic valve
#
# ---Right ventricle-----------------------------------------------------------
T_RV=200                 # right ventricle
T_RVBP=200               # right ventricular blood pool
T_VCINF=200             # vena cava inferior
T_VCSUP=200              # vena cava superior
# T_VCINFBP=91            # vena cava inferior blood pool
# T_VCINFBPCUT=
# T_VCSUPBP=92            # vena cava superior blood pool
# T_VCBP=93               # venae cavae blood pool
T_PULMARTERY=200         # pulmonary artery
T_PULMARTERYBP=200       # pulmonary artery blood pool
T_TRICUSPVV=200          # tricuspic valve
T_PULMVV=200             # pulmonic valve
#
# ---Torso---------------------------------------------------------------------
# T_TORSO=401             # torso
# T_BOWELGAS=402          # bowel gas
# T_LUNGS=407
# #
# # ---Mechanics-----------------------------------------------------------------
# T_APICALCUSHION=12      # apical cushion      (CARDIOPROFF)
# T_APICALCUSHIONBASE=13  # apical cushion base (CARDIOPROOF)
# T_AORTA_TERMINI=11      #                     (CARDIOPROOF)
# T_BASALCUSHION=14       #                     (PUSHCART)
# T_BASALCUSHIONBASE=15   #                     (PUSHCART)
# T_LVLID=41              #                     (PUSHCART)
# T_RVLID=46              #                     (PUSHCART)
# #
# # ---Postprocessing------------------------------------------------------------
# T_OFFSETAPEX=100
# T_OFFSETMID=200
# T_OFFSETBASE=300
# T_LVENDO=25
# T_LVMID=50
# T_LVEPI=75
