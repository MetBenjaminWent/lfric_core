include ../../Makefile.inc

BIN_DIR   = ../../bin

EXE  = dynamo
-include $(OBJ_DIR)/$(EXE).mk

.PHONY: all
all: $(BIN_DIR)/$(EXE)

$(BIN_DIR)/$(EXE): $(OBJ_DIR)/$(EXE) | $(BIN_DIR)
	@echo "Installing $@"
	@cp $< $(BIN_DIR)

# Directories
$(BIN_DIR):
	@echo "Creating $@"
	@mkdir -p $@

$(OBJ_DIR):
	@echo "creating $@"
	@mkdir -p $@

# Rules
$(OBJ_DIR)/%.mod: $(OBJ_DIR)/%.o
	@echo "Require $@"

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	@echo "Compile $<"
	@$(FC95) $(F95FLAGS) \
	         $(F95_MOD_DESTINATION_ARG) -I $(OBJ_DIR) \
	         -c -o $@ $<

$(OBJ_DIR)/$(EXE): $($(shell echo $(EXE) | tr a-z A-Z)_OBJS)
	@echo "Linking $@"
	@$(FC95) $(F95FLAGS) -o $@ $^

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
	-rm -f $(BIN_DIR)/dynamo

-include $(OBJ_DIR)/$(DEP_FILE)
