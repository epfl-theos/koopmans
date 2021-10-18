# List of available tasks
.PHONY: help install submodules configure_4 configure_6 configure espresso_4 espresso_6 espresso workflow tests mock_tests clean clean_espresso clean_tests

MPIF90 = "mpif90"

help:
	@echo '  _                                                '
	@echo ' | | _____   ___  _ __  _ __ ___   __ _ _ __  ___  '
	@echo ' | |/ / _ \ / _ \| '"'"'_ \| '"'"'_ ` _ \ / _` | '"'"'_ \/ __| '
	@echo ' |   < (_) | (_) | |_) | | | | | | (_| | | | \__ \ '
	@echo ' |_|\_\___/ \___/| .__/|_| |_| |_|\__,_|_| |_|___/ '
	@echo '                 |_|                               '
	@echo ''
	@echo ' For performing Koopmans spectral functional calculations with Quantum ESPRESSO'
	@echo ' See README.rst for more information'
	@echo ''
	@echo ' To install'
	@echo ' > make install'
	@echo ' (For more details, see README.rst)'
	@echo ''
	@echo ' To run the test suite'
	@echo ' > make tests'
	@echo ''
	@echo ' To test without calling QE'
	@echo ' > make mock_tests'
	@echo ''
	@echo ' To clean the repository (removing compilations of QE and all test outputs)'
	@echo ' > make clean'
	@echo ''

install: submodules espresso workflow

submodules:
	git submodule init
	git submodule update

configure_4:
	cd quantum_espresso/cp_koopmans; ./configure MPIF90=$(MPIF90);

configure_6:
	cd quantum_espresso/qe_koopmans; ./configure MPIF90=$(MPIF90);

configure: configure_4 configure_6

espresso_4:
	@(cd quantum_espresso/cp_koopmans; $(MAKE) kcp)

espresso_6:
	@(cd quantum_espresso/qe_koopmans; $(MAKE) kc)

espresso: configure_4 espresso_4 configure_6 espresso_6

workflow:
	python3 -m pip install --upgrade pip
	python3 -m pip install -e . -e ase/

clean: clean_espresso clean_tests

clean_espresso:
	@(cd quantum_espresso/cp_koopmans; $(MAKE) veryclean)
	@(cd quantum_espresso/qe_koopmans; $(MAKE) veryclean)

tests:
	python3 -m pytest -m "standard" tests/

mock_tests:
	python3 -m pytest -m "mock" tests/

KEEP_TESTS := -name '*.py' -o -name '*.json'
TEST_DIRS := ki; pkipz; kipz; init; calc_alpha; final; TMP-CP; ecutwfc; charged; neutral
CLEAN_TEST_DIRS := $(pathsubst %,;, -o name %)

clean_tests:
	find tests/test_??/ -path */*/input_files -prune -false -o \
		-type f ! \( ${KEEP_TESTS} \) | xargs -i rm -v {}
	find tests/test_??/* -type d ${CLEAN_TEST_DIRS} | grep -v input_files | xargs -i rm -vrf {}

