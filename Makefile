# List of available tasks
.PHONY: help install submodules configure_4 configure_7 configure espresso_4 espresso_7 espresso_utils espresso workflow tests clean clean_espresso clean_tests benchmarks

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
	@echo ' To clean the repository (removing compilations of QE and all test outputs)'
	@echo ' > make clean'
	@echo ''

install: submodules espresso workflow

submodules:
	git submodule init
	git submodule update

configure_4:
	cd quantum_espresso/kcp; ./configure MPIF90=$(MPIF90);

configure_7:
	cd quantum_espresso/q-e; ./configure MPIF90=$(MPIF90);

configure: configure_4 configure_7

espresso_4:
	@(cd quantum_espresso/kcp; $(MAKE) kcp)

espresso_7:
	@(cd quantum_espresso/q-e; $(MAKE) pw kcw)

espresso_utils:
	@(cd quantum_espresso/utils; $(MAKE) all)

espresso: configure_4 espresso_4 configure_7 espresso_7 espresso_utils

workflow:
	python3 -m pip install --upgrade pip
	python3 -m pip install -e . -e ase/

clean: clean_espresso clean_tests

clean_espresso:
	@(cd quantum_espresso/kcp; $(MAKE) veryclean)
	@(cd quantum_espresso/q-e; $(MAKE) veryclean)
	@(cd quantum_espresso/utils; $(MAKE) clean)

tests:
	python3 -m pytest -m "not tutorials" tests/

clean_tests:
	rm -r tests/tmp

tests:
	python3 -m pytest --generate_benchmark
