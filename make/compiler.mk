##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Identify the compiler and bring in the correct specifics.
# This file is intended to be "include"d from other make files.
#
$(info Make version $(MAKE_VERSION))
ifeq ($(filter else-if,$(value .FEATURES)),)
  $(error The build system requires else-if support from GMake)
endif

EMPTY :=
SPACE := $(EMPTY) # Comment to highlight space character.

# The default value of FC is almost always "f77" which is of no use to us.
# An empty FC is also of no use.
ifneq "$(or $(filter default, $(origin FC)), $(filter x, x$(FC)))" ""
  $(error The FC environment variable must be set to a Fortran compiler command)
endif

FORTRAN_COMPILER := $(shell basename $(FC))
# Some compiler names are sometimes extended so we special-case them.
ifeq ( $(findstring gfortran,$(FORTRAN_COMPILER)), gfortran )
  FORTRAN_COMPILER = gfortran
else ifeq ( $(findstring xlf,$(FORTRAN_COMPILER)), xlf )
  FORTRAN_COMPILER = xlf
endif

# Attempt to identify Cray systems...
#
ifdef PE_ENV

    CRAY_ENVIRONMENT = true

    ifeq '$(PE_ENV)' 'CRAY'
        FORTRAN_COMPILER = crayftn
    else ifeq '$(PE_ENV)' 'INTEL'
        FORTRAN_COMPILER = ifort
    else ifeq '$(PE_ENV)' 'GNU'
        FORTRAN_COMPILER = gfortran
    else ifeq '$(PE_ENV)' 'PGI'
        FORTRAN_COMPILER = pgfortran
    else
        $(error Unrecognised Cray programming environment)
    endif

endif
