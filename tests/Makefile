# List of available tasks
.PHONY: install tests clean help

install:
	@# install testcode
	git submodule init
	git submodule update

tests:
ifeq ($(wildcard ./testcode/bin/testcode.py),)
	@$(MAKE) install
endif
	testcode/bin/testcode.py $(TC) --verbose

KEEP_FILES := bench* *.json 
KEEP_TESTS := $(patsubst %,! -name %,${KEEP_FILES})
CLEAN_FILES := test.out.* test.err.*
CLEAN_TESTS := $(patsubst %,-o -name %,${CLEAN_FILES})

clean:
	find test_?? -type f ${KEEP_TESTS} ${CLEAN_TESTS} | xargs -i rm -v {}
	find -name 'init' -o -name 'calc_alpha' -o -name 'final' | xargs -i rm -vr {}

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

