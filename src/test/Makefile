##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

ROOT = ../..
include $(ROOT)/Makefile.inc

OBJ_DIR = $(BUILD_DIR)/dynamo

TEST_SRC     = $(shell find . -name '*.pf')
TEST_OBJS    = $(patsubst %.pf,$(OBJ_DIR)/%.o,$(TEST_SRC))

SRC        = $(wildcard *.[Ff]90)
OBJS       = $(patsubst %.[Ff]90,$(OBJ_DIR)/%.o,$(SRC))

DRIVER_OBJ = $(OBJ_DIR)/driver.o

TARGET_OBJ = $(patsubst %_test.pf,$(OBJ_DIR)/%.o,$(TEST_SRC))

.PHONY: test
test: pfunit $(OBJ_DIR)/test
	$(OBJ_DIR)/test

include pfunit.mk

$(OBJ_DIR)/test: $(TEST_OBJS) $(OBJS) $(DRIVER_OBJ) $(TARGET_OBJ)
	$(FC) -L$(PFUNIT_INSTALL)/lib -o $@ $^ -lpfunit

$(DRIVER_OBJ): $(PFUNIT_INSTALL)/include/driver.F90 testSuites.inc
	$(FC) $(FFLAGS) -c -I$(PFUNIT_INSTALL)/mod -I. -o $@ $<

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.F90
	$(FC) $(FFLAGS) -c $(F_MOD_DESTINATION_ARG) -I$(PFUNIT_INSTALL)/mod \
	    -o $@ $<

.PRECIOUS: $(OBJ_DIR)/%.F90
$(OBJ_DIR)/%.F90: %.pf
	$(PFUNIT_INSTALL)/bin/pFUnitParser.py $< $@

.PHONY: clean
clean: cleanpfunit
	-rm -rf $(OBJ_DIR)
