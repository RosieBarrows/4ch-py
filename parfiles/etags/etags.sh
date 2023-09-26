#!/bin/bash

# =============================================================================
# SMOOTHED HEART FOUR CHAMBERS
# (c) images          Steven Niederer
# (c) segmentation    SN
# (c) meshing         AJP
# (c) fibers & sheets JB, AJP
# =============================================================================


## CHANGE ONLY THESE LABELS SO THAT THEY MATCH YOUR MESH LABELS
## ============================================================
T_LV=1                  # left ventricular tissue
T_RV=2                  # right ventricular tissue


## DO NOT CHANGE ANY OF THE LABELS BELOW
## ==============================

# ---SEGMENTATION LABELS-------------------------------------------------------
# #
T_UNUSED=200             # unused tag
#
# ---Left atrium---------------------------------------------------------------
T_LA=200                 # left atrial wall
T_LABP=200               # left atrial blood pool
T_LINFPULMVEINCUT=200    # left inferior pulmonary vein (cut)
T_LSUPPULMVEINCUT=200    # left superior pulmonary vein (cut)
T_RINFPULMVEINCUT=200      # right inferior pulmonary vein (cut)
T_RSUPPULMVEINCUT=200      # right superior pulmonary vein (cut)
#
# ---Right atrium--------------------------------------------------------------
T_RA=200                 # right atrial wall
T_RABP=200               # right atrial blood pool

# ---Left ventricle
T_LVBP=200              # left ventricular blood pool
T_AORTA=200              # aorta
T_AORTABP=200            # aortic blood pool
T_MITRALVV=200           # mitral valve
T_AORTICVV=200          # aortic valve
#
# ---Right ventricle-----------------------------------------------------------
T_RVBP=200               # right ventricular blood pool
T_VCINF=200             # vena cava inferior
T_VCSUP=200              # vena cava superior
T_PULMARTERY=200         # pulmonary artery
T_PULMARTERYBP=200       # pulmonary artery blood pool
T_TRICUSPVV=200          # tricuspic valve
T_PULMVV=200             # pulmonic valve

