##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

PGFORTRAN_VERSION := $(shell pgfortran --version \
                             | awk -F "[. -]" '/[0-9]+\.[0-9]+-[0-9]+/ { printf "%03i%03i%03i", $$(2), $$(3), $$(4) }' )

$(info ** Portland Fortran version $(PGFORTRAN_VERSION))

F_MOD_DESTINATION_ARG = -module$(SPACE)
