# ICSD: Inversion of current source density
**for Geophysicists**

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![PyPI version](https://badge.fury.io/py/icsd.svg)](https://badge.fury.io/py/icsd)
[![Docs](https://github.com/Peruz/icsd/actions/workflows/documentation.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/documentation.yml)
[![Tests](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml)

<!--
![sphinx doc](https://github.com/Peruz/icsd/actions/workflows/sphinx_doc.yml/badge.svg)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/tesspy.svg)](https://anaconda.org/conda-forge/tesspy)
-->

`iCSD` is a python library for current source density inversion.


## Installation
You can install ``iCSD`` from PyPI using pip:
```
pip install iCSD
```

### Developement environment for icsd

Create a new environment, activate it and install the package:
```shell
conda env create -f environment.yml
conda activate icsd_env
pip install -e .
```

Update iCSD:
```shell
conda env update --file environment.yml --prune
```

### Run Tests

```shell
pytest -v --color yes --cov-config .coveragerc --cov icsd --cov-append --cov-report term-missing --cov-report xml tests
```


## Dependencies

`iCSD`'s dependencies are: `?`, `?`

## Documentation
The official documentation is hosted on **[ReadTheDocs](?)**.

## Examples

Some examples of laboratory acquisitions and iCSD (based on cotton experiments)


