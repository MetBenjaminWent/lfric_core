##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

$(info UM physics specific compile options)

science/%.o science/%.mod: export FFLAGS := $(FFLAGS) $(FFLAGS_UM_PHYSICS)
