<h2 align="center">pyGeoCSD</h2>

<p align="center">
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![PyPI version](https://badge.fury.io/py/icsd.svg)](https://badge.fury.io/py/icsd)
[![Docs](https://github.com/Peruz/icsd/actions/workflows/documentation.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/documentation.yml)
[![Tests](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml)
</p>

<!--
![sphinx doc](https://github.com/Peruz/icsd/actions/workflows/sphinx_doc.yml/badge.svg)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/tesspy.svg)](https://anaconda.org/conda-forge/tesspy)
-->

## About 

`pyGeoCSD` is a python library for current source density inversion **for Geophysicists**. 



### Installation
You can install ``pyGeoCSD`` from PyPI using pip:
```
pip install pyGeoCSD
```

#### Developement environment for pyGeoCSD

Create a new environment, activate it and install the package:
```shell
conda env create -f environment.yml
conda activate pyGeoCSD_env
pip install -e .
```

Update pyGeoCSD:
```shell
conda env update --file environment.yml --prune
```

#### Run Tests

```shell
pytest -v --color yes --cov-config .coveragerc --cov pyGeoCSD --cov-append --cov-report term-missing --cov-report xml tests
```


#### Dependencies

`pyGeoCSD`'s dependencies are limited to numpy, scipy, matplotlib, pyvista and kneed.

## Documentation
The official documentation is hosted on **[ReadTheDocs](?)**.

### Examples

Some examples of laboratory acquisitions and pyGeoCSD (based on cotton experiments)

<!-- ![Drag Racing](Dragster.jpg) -->
