##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

export IGNORE_DEPENDENCIES = netcdf MPI ESMF
export EXTERNAL_DYNAMIC_LIBRARIES = netcdf netcdff netcdf_c++ esmf
export EXTERNAL_STATIC_LIBRARIES = 

export ROOT = ../..

include $(ROOT)/include.mk

ifdef GFORTRAN_VERSION
  ifeq ($(shell test $(GFORTRAN_VERSION) -lt 040601; echo $$?), 0)
    $(error GFortran is too old. Must be at least 4.6.1)
  endif
endif

OBJ_DIR = $(BUILD_DIR)/dynamo

.PHONY: all
all: modules

.PHONY: modules
modules: applications
	$(MAKE) -f dynamo.mk modules OBJ_DIR=$(OBJ_DIR)

.PHONY: applications
applications:
	$(MAKE) -f dynamo.mk applications OBJ_DIR=$(OBJ_DIR)

.PHONY: clean
clean:
	$(MAKE) -f dynamo.mk clean OBJ_DIR=$(OBJ_DIR)
