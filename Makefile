# List of available tasks
.PHONY: help install submodules configure espresso workflow tests mock_tests clean clean_espresso clean_tests

MPIF90 = "mpiifort"
FFT_LIBS="-lfftw3"
BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
LAPACK_LIBS=" "

help:
	@echo '  _                                                '
	@echo ' | | _____   ___  _ __  _ __ ___   __ _ _ __  ___  '
	@echo ' | |/ / _ \ / _ \| '"'"'_ \| '"'"'_ ` _ \ / _` | '"'"'_ \/ __| '
	@echo ' |   < (_) | (_) | |_) | | | | | | (_| | | | \__ \ '
	@echo ' |_|\_\___/ \___/| .__/|_| |_| |_|\__,_|_| |_|___/ '
	@echo '                 |_|                               '
	@echo ''
	@echo ' Python module for performing KI and KIPZ calculations'
	@echo ' For more information, see README.rst'
	@echo ''
	@echo ' To run the test suite'
	@echo ' > make tests'
	@echo ''
	@echo ' To test without calling QE:'
	@echo ' > make mock_tests'
	@echo ''
	@echo ' To remove all test outputs in the test subdirectories:'
	@echo ' > make clean'
	@echo ''

install: submodules configure espresso workflow

submodules:
	git submodule init
	git submodule update --remote --merge

configure:
	cd quantum_espresso/cp_koopmans; ./configure MPIF90=$(MPIF90) FFT_LIBS=$(FFT_LIBS) BLAS_LIBS=$(BLAS_LIBS) LAPACK_LIBS=$(LAPACK_LIBS);
	cd quantum_espresso/qe_koopmans; ./configure MPIF90=$(MPIF90);

espresso:
	@(cd quantum_espresso/cp_koopmans; $(MAKE) kcp)
	@(cd quantum_espresso/qe_koopmans; $(MAKE) kc)

workflow:
	python3 -m pip install --upgrade pip
	python3 -m pip install . ase/

clean: clean_espresso clean_tests

clean_espresso:
	@(cd quantum_espresso/cp_koopmans; $(MAKE) clean)
	@(cd quantum_espresso/qe_koopmans; $(MAKE) clean)

tests:
	python3 -m pytest -m "not mock" tests/

mock_tests:
	python3 -m pytest -m "mock" tests/

KEEP_TESTS := -name '*.py' -o -name '*.json'
TEST_DIRS := ki; pkipz; kipz; init; calc_alpha; final; TMP-CP; ecutwfc; charged; neutral
CLEAN_TEST_DIRS := $(pathsubst %,;, -o name %)

clean_tests:
	find tests/test_??/ -path */*/input_files -prune -false -o \
		-type f ! \( ${KEEP_TESTS} \) | xargs -i rm -v {}
	find tests/test_??/* -type d ${CLEAN_TEST_DIRS} | grep -v input_files | xargs -i rm -vrf {}

