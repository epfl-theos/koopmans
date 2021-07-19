# List of available tasks
.PHONY: help tests mock_tests

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

tests:
	python3 -m pytest -m "not mock" tests/

mock_tests:
	python3 -m pytest -m "mock" tests/

KEEP_TESTS := -name '*.py' -o -name '*.json'
TEST_DIRS := ki; pkipz; kipz; init; calc_alpha; final; TMP-CP; ecutwfc; charged; neutral
CLEAN_TEST_DIRS := $(pathsubst %,;, -o name %)

clean:
	find tests/test_??/ -path */*/input_files -prune -false -o \
		-type f ! \( ${KEEP_TESTS} \) | xargs -i rm -v {}
	find tests/test_??/* -type d ${CLEAN_TEST_DIRS} | grep -v input_files | xargs -i rm -vrf {}

