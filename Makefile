##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

# A trio of targets are offered for various levels of debugging support and
# speed of execution.
#
.PHONY: full-debug
full-debug: OPTIMISATION?=NONE
full-debug: SYMBOLS?=YES
full-debug: CHECKS?=YES
full-debug: build

.PHONY: fast-debug
fast-debug: OPTIMISATION?=SAFE
fast-debug: SYMBOLS?=YES
fast-debug: build

.PHONY: production
production: OPTIMISATION?=RISKY
production: SYMBOLS?=YES
production: build

.PHONY: build
build:
	$(MAKE) -C src/dynamo OPTIMISATION=$(OPTIMISATION) \
	                      SYMBOLS=$(SYMBOLS) CHECKS=$(CHECKS)
	$(MAKE) -C src/test

# The test target allows tests to be run on their own. Since there is no way
# to know which of the above builds is desired this target does not list a
# build as a prerequisite. If you try to run it without building first it will
# fail.
#
.PHONY: test
test:
	$(MAKE) -C src/test

# Build the projects documentation. This includes both API and design documents.
.PHONY: doc docs
doc docs:
	$(MAKE) -C Docs

# Clean only the dynamo build. This leaves pFUnit alone.
#
.PHONY: clean
clean:
	$(MAKE) -C src/dynamo clean
	$(MAKE) -C src/test clean

# Clean both dynamo and pFUnit builds.
#
.PHONY: clean-all
clean-all:
	$(MAKE) -C src/dynamo clean ALL=1
	$(MAKE) -C src/test clean ALL=1
