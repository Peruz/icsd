# Build, package, test, and clean
PROJECT=icsd
TESTDIR=tmp-test-dir-with-unique-name
PYTEST_ARGS=-v --color yes --cov-config ../.coveragerc --cov icsd --cov-append --cov-report term-missing --cov-report xml tests
STYLE_CHECK_FILES=setup.py $(PROJECT) examples

help:
	@echo "Commands:"
	@echo ""
	@echo "  install   install in editable mode"	
	@echo "  test      run the test suite (including doctests) and report coverage"
	@echo "  format    run isort and black to automatically format the code"
	@echo "  check     run code style and quality checks (black, isort and flake8)"
	@echo "  clean     clean up build and generated files"
	@echo "  doc       build sphinx doc"
	@echo "  publish   publish to pypi"
	@echo ""

install:
	pip install --no-deps -e .

test: test_coverage

test_coverage:
	# Run a tmp folder to make sure the tests are run on the installed version
	#mkdir -p $(TESTDIR)
	#cd $(TESTDIR); NUMBA_DISABLE_JIT=1 MPLBACKEND='agg' pytest $(PYTEST_ARGS) $(PROJECT)
	#cp $(TESTDIR)/.coverage* .
	#rm -rvf $(TESTDIR)
	#pytest -v --color yes --cov-config --cov icsd --cov-append --cov-report term-missing --cov-report xml tests
	pytest -v --color yes --cov-config .coveragerc --cov icsd --cov-append --cov-report term-missing --cov-report xml tests
	
format: isort black

check: isort-check black-check flake8

black:
	black $(STYLE_CHECK_FILES)

black-check:
	black --check $(STYLE_CHECK_FILES)

isort:
	isort $(STYLE_CHECK_FILES)

isort-check:
	isort --check $(STYLE_CHECK_FILES)

flake8:
	flake8 $(STYLE_CHECK_FILES)

clean:
	find . -name "*.pyc" -exec rm -v {} \;
	find . -name ".coverage.*" -exec rm -v {} \;
	rm -rvf build dist MANIFEST *.egg-info __pycache__ .coverage .cache .pytest_cache $(PROJECT)/_version.py
	rm -rvf $(TESTDIR) dask-worker-space 
	
doc: 
	sphinx-build ./docsrc/ _sphinx_build
	
publish:
	python3 setup.py sdist bdist_wheel
        twine upload dist/*
        grayskull pypi icsd

        
        
        
