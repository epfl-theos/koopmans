# List of available tasks
.PHONY: help install tests clean

install:
	@# install testcode
	git submodule init
	git submodule update --remote --merge

tests:
ifeq ($(wildcard ./testcode/bin/testcode.py),)
	@$(MAKE) install
endif
ifneq ($(wildcard ./test_06/converged.json),)
	@rm test_06/converged.json 2> /dev/null
endif
	testcode/bin/testcode.py $(TC) --verbose

KEEP_FILES := bench* *.json 
KEEP_TESTS := $(patsubst %,! -name %,${KEEP_FILES})

CLEAN_FILES := test.out.* test.err.*
CLEAN_TESTS := $(patsubst %,-o -name %,${CLEAN_FILES})

CLEAN_TEST_DIRS := -name ki -o -name pkipz -o -name kipz

CLEAN_TEST_SUBDIRS := -name init -o -name calc_alpha -o -name final -o \
                      -name TMP-CP -o -name ecutwfc*


clean:
	find test_?? -type f ${KEEP_TESTS} ${CLEAN_TESTS} | xargs -i rm -v {}
	find ${CLEAN_TEST_DIRS} | xargs -i rm -vr {}
	find ${CLEAN_TEST_SUBDIRS} | xargs -i rm -vrf {}

help:
	@echo ''
	@echo 'python_KI QC tests'
	@echo ''
	@echo 'To run:'
	@echo '> make tests'
	@echo 'To remove all test outputs in the test subdirectories:'
	@echo '> make clean'
	@echo 'To run only test_01:'
	@echo '> make tests TC="-c test_01"'
	@echo 'To install testcode:'
	@echo '> make install'
	@echo ''
	@echo 'For information on testcode see http://testcode.readthedocs.org'
	@echo ''

