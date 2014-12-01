##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

include ../../include.mk

AR ?= ar

BIN_DIR  = $(ROOT)/bin
DATABASE = $(OBJ_DIR)/dependencies.db

include $(OBJ_DIR)/programs.mk $(OBJ_DIR)/dependencies.mk

ALL_SRC     = $(shell find . -name "*.[Ff]90")
TOUCH_FILES = $(patsubst ./%.f90,$(OBJ_DIR)/%.t,$(patsubst ./%.F90,$(OBJ_DIR)/%.t,$(ALL_SRC)))
PROGRAMS    = $(patsubst %.o,%,$(notdir $(PROG_OBJS)))
ALL_MODULES = $(filter-out $(PROG_OBJS),$(patsubst %.t,%.o,$(TOUCH_FILES)))

.SECONDEXPANSION:

.PHONY: applications
applications: $(patsubst %,$(BIN_DIR)/%,$(PROGRAMS))

.PHONY: modules
modules: $(OBJ_DIR)/modules.a

$(BIN_DIR)/%: $(OBJ_DIR)/%.x | $(BIN_DIR)
	@echo "Installing $@"
	$(Q)cp $(OBJ_DIR)/$(notdir $@).x $@

# Directories

$(BIN_DIR):
	@echo "Creating $@"
	$(Q)mkdir -p $@

SUBDIRS = $(shell find * -type d -prune)
OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%/,$(SUBDIRS))

$(OBJ_DIR) $(OBJ_SUBDIRS):
	@echo "Creating $@"
	$(Q)mkdir -p $@

# Build Rules

INCLUDE_ARGS = -I $(OBJ_DIR) $(patsubst %,-I $(OBJ_DIR)/%,$(SUBDIRS))

$(OBJ_DIR)/%.o $(OBJ_DIR)/%.mod: %.F90 | $(OBJ_DIR)
	@echo "Compile $<"
	$(Q)$(FCOM) $(CPPFLAGS) $(FFLAGS) \
	            $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	            $(INCLUDE_ARGS) -c -o $(basename $@).o $<

$(OBJ_DIR)/%.o $(OBJ_DIR)/%.mod: %.f90 | $(OBJ_DIR)
	@echo "Compile $<"
	$(Q)$(FCOM) $(CPPFLAGS) $(FFLAGS) \
	            $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	            $(INCLUDE_ARGS) -c -o $(basename $@).o $<

$(OBJ_DIR)/modules.a: $(ALL_MODULES)
	$(Q)$(AR) -r $@ $^

$(OBJ_DIR)/%.x: $$($$(shell echo $$* | tr a-z A-Z)_OBJS)
	@echo "Linking $@"
	$(Q)$(FCOM) $(FFLAGS) $(LDFLAGS) -o $@ $^

# Dependencies

$(OBJ_DIR)/programs.mk: | $(OBJ_DIR)
	$(TOOL_DIR)/ProgramObjects -database $(DATABASE) $@

$(OBJ_DIR)/dependencies.mk: $(TOUCH_FILES) | $(OBJ_DIR) $(OBJ_SUBDIRS)
	$(Q)$(TOOL_DIR)/DependencyRules -database $(DATABASE) $@

IGNORE_DEPENDENCIES := $(patsubst %,-ignore %,$(IGNORE_DEPENDENCIES))
export LDFLAGS += $(patsubst %, -l%,$(EXTERNAL_LIBRARIES))

$(OBJ_DIR)/%.t: %.F90 | $(OBJ_SUBDIRS)
	@echo Analysing $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_DEPENDENCIES) \
	                               $(DATABASE) $< && touch $@

$(OBJ_DIR)/%.t: %.f90 | $(OBJ_SUBDIRS)
	@echo Analysing $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_DEPENDENCIES) \
	                                   $(DATABASE) $< && touch $@

# Special Rules

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
	-rm -f $(BIN_DIR)/*
