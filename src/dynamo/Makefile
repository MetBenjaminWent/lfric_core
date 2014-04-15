##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

ROOT = ../..
include $(ROOT)/Makefile.inc

export SRC        = $(wildcard *.F90)
export OBJ_DIR    = $(BUILD_DIR)/dynamo
export DEP_FILE   = dependencies.mk

.PHONY: all
all: $(OBJ_DIR)/$(DEP_FILE)
	$(MAKE) -f dynamo.mk all

$(OBJ_DIR)/$(DEP_FILE): $(SRC) $(OBJ_DIR)
	@$(TOOL_DIR)/DependencyAnalyser $(OBJ_DIR) $(OBJ_DIR) $(OBJ_DIR) .

$(OBJ_DIR):
	@echo "Creating $@"
	@mkdir -p $@

.PHONY: clean
clean:
	$(MAKE) -f dynamo.mk clean
