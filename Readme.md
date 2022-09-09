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
You can install ``iCSD`` from PyPI using pip (**Not Recommended**):
```
pip install iCSD
```

and from conda (**Recommended**):
```
conda install -c conda-forge iCSD
```


### Creating a new environment for icsd

`iCSD` depends on `?`, which could make the installation sometimes tricky because of the conflicts with the current packages. Therefore, we recommend creating a new clean environment and installing the dependencies from the conda-forge channel.


Create a new environment:
```shell
conda create -n icsd_env -c conda-forge
```

Activate this environment:
```shell
conda activate icsd_env
```

Install iCSD from conda-forge channel:
```shell
conda install -c conda-forge iCSD
```



### Developement environment for icsd

Create a new environment:
```shell
conda env create -f environment.yml
```

Activate this environment:
```shell
conda activate icsd_env
```

Install iCSD from setup.cfg:
```shell
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

some examples of laboratory acquisitions and iCSD (based on cotton experiments)


