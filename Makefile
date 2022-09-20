# List of available tasks
.PHONY: help install submodules configure_4 configure_7 configure espresso_4 espresso_7 espresso_utils espresso espresso_4_install espresso_7_install espresso_utils_install espresso_install workflow tests clean clean_espresso clean_tests benchmarks

MPIF90 = "mpif90"
PREFIX ?= /usr/local

default: submodules espresso workflow

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
	@echo ' To install everything'
	@echo ' > make'
	@echo ' > sudo make install'
	@echo ' (For more details, see README.rst)'
	@echo ''
	@echo ' To run the test suite'
	@echo ' > make tests'
	@echo ''
	@echo ' To clean the repository (removing compilations of QE and all test outputs)'
	@echo ' > make clean'
	@echo ''

submodules:
	git submodule init
	git submodule update

configure_4:
	cd quantum_espresso/kcp; ./configure MPIF90=$(MPIF90);

configure_7:
	cd quantum_espresso/q-e; ./configure MPIF90=$(MPIF90);

configure: configure_4 configure_7

espresso_4:
	@(cd quantum_espresso/kcp && $(MAKE) kcp)

espresso_4_install:
	@(cd quantum_espresso/kcp && $(MAKE) install PREFIX=$(PREFIX))

espresso_7:
	@(cd quantum_espresso/q-e && $(MAKE) pw kcw ph)

espresso_7_install:
	@(cd quantum_espresso/q-e && $(MAKE) install PREFIX=$(PREFIX))

espresso_utils:
	@(cd quantum_espresso/utils && $(MAKE) all)

espresso_utils_install:
	@(cd quantum_espresso/utils && $(MAKE) install PREFIX=$(PREFIX))

espresso: configure_4 espresso_4 configure_7 espresso_7 espresso_utils

install: espresso_4_install espresso_7_install espresso_utils_install

workflow:
	python3 -m pip install --upgrade pip
	python3 -m pip install -e .

clean: clean_espresso clean_tests

clean_espresso:
	@(cd quantum_espresso/kcp; $(MAKE) veryclean)
	@(cd quantum_espresso/q-e; $(MAKE) veryclean)
	@(cd quantum_espresso/utils; $(MAKE) clean)

tests:
	python3 -m pytest -m "not tutorials" tests/

clean_tests:
	rm -rf tests/tmp

benchmark:
	python3 -m pytest --generate_benchmark
