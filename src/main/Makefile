ROOT = ../..
include $(ROOT)/Makefile.inc

export SRC        = $(wildcard *.f90)
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
