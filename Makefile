##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

.PHONY: test
test: build
	$(MAKE) -C src/test

.PHONY: build
build:
	$(MAKE) -C src/dynamo

.PHONY: doc docs
doc docs:
	$(MAKE) -C Docs

.PHONY: clean
clean:
	$(MAKE) -C src/dynamo clean
	$(MAKE) -C src/test clean

.PHONY: clean-all
clean-all:
	$(MAKE) -C src/dynamo clean ALL=1
	$(MAKE) -C src/test clean ALL=1
