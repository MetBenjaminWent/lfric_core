##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

LYX_DOCUMENTS = formulation/dynamo_formulation.lyx glossary/Glossary.lyx
LYX_OUTPUTS = $(patsubst %.lyx,%.pdf,$(LYX_DOCUMENTS))

DOXYGEN ?= doxygen

TOP := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
$(info $(TOP))

.PHONY: most
most: api lyx

.PHONY: api
api:
	$(DOXYGEN) Doxyfile

.PHONY: lyx
lyx: $(LYX_OUTPUTS)

%.pdf: %.lyx
	env TEXINPUTS=$(TOP)/tex//: lyx -e pdf $<

.PHONY: gource
gource:
	$(MAKE) -C visualisation

.PHONY: clean
clean:
	-rm -rf docs/api/*
	-rm $(LYX_OUTPUTS)

