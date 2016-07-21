##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Compiler specific details: Portland Fortran with GCC C/C++
#
$(info ** Chosen Portland Fortran compiler)

ifdef PGFORTRAN_VERSION
  ifeq ($(shell test $(PGFORTRAN_VERSION) -lt 015007000; echo $$?), 0)
    $(error PGFortran is too old. Must be at least 15.7)
  endif
endif

FFLAGS_COMPILER           = -mp
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O4
FFLAGS_DEBUG              = -g -traceback
FFLAGS_RUNTIME            = -Mbounds -Mchkptr -Mchkstk

LDFLAGS_COMPILER = -g