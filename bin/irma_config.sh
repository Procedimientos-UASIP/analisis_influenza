# HEADER
PARAM_FILE_NAME="FLU-avian"
PARAM_FILE_AUTHOR="S. Shepard"
PARAM_FILE_VERSION="1.0"
PARAM_FILE_DATE="2016-01-28"

# CONSENSUS REFINEMENT & READ SELECTION
INS_T=0.25                              # threshold for insertion refinement
DEL_T=0.60                              # threshold for deletion refinement
SKIP_E=0                                # skip reference elongation
MIN_LEN=125
GRID_ON=0

#SINGLE_LOCAL_PROC=28
TMP="/home/admcenapa/DAVID_RENDON/TEMP/IRMA/"
REF_SET="/home/admcenapa/DAVID_RENDON/TEMP/CPA-03682-25R2_TEST/ENSAMBLE/IRMA/ref_mod.fasta"

# STAGES
MATCH_PROG="BLAT"
SORT_PROG="LABEL BLAT"
ALIGN_PROG="SAM BLAT"
ASSEM_PROG="SSW"
