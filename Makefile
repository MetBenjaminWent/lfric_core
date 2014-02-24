.PHONY: test
test: all
	$(MAKE) -C src/test

.PHONY: all
all:
	$(MAKE) -C src/main

.PHONY: run
run: test
	$(MAKE) -C run

.PHONY: doc docs
doc docs:
	$(MAKE) -C Docs docs

.PHONY: clean
clean:
	$(MAKE) -C src/main clean
	$(MAKE) -C src/test clean
