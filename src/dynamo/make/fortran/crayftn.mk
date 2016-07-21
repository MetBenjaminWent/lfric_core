##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Compiler specific details: Cray Fortran and C/C++
#
$(info ** Chosen Cray Fortran compiler)

FFLAGS_COMPILER           = -h omp
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O3
FFLAGS_DEBUG              = -Gfast
FFLAGS_WARNINGS           = -m 0

DEPRULE_FLAGS = -moduleobjects

