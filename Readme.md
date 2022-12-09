<h2 align="center">pyGeoCSD</h2>

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![PyPI version](https://badge.fury.io/py/icsd.svg)](https://badge.fury.io/py/icsd)
[![Docs](https://github.com/Peruz/icsd/actions/workflows/documentation.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/documentation.yml)
[![Tests](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml/badge.svg)](https://github.com/Peruz/icsd/actions/workflows/tests_package.yml)

<!--
![sphinx doc](https://github.com/Peruz/icsd/actions/workflows/sphinx_doc.yml/badge.svg)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/tesspy.svg)](https://anaconda.org/conda-forge/tesspy)
-->

## About 

`pyGeoCSD` is a python library for current source density inversion **for Geophysicists**. 

## Documentation
The official documentation is hosted on **[ReadTheDocs](?)**.

## Installation

`pyGeoCSD`'s dependencies are limited to numpy, scipy, matplotlib, pyvista and kneed.


You can install ``pyGeoCSD`` from PyPI using pip:
```
pip install pyGeoCSD
```

Developement environment for pyGeoCSD: create a new environment, activate it and install the package:
```shell
conda env create -f environment.yml
conda activate pyGeoCSD_env
pip install -e .
```

Update pyGeoCSD:
```shell
conda env update --file environment.yml --prune
```

#### Run tests

```shell
pytest -v --color yes --cov-config .coveragerc --cov pyGeoCSD --cov-append --cov-report term-missing --cov-report xml tests
```



## Getting involved

- [Contributing Guide](https://github.com/Peruz/icsd/blob/main/CONTRIBUTING.md)
to see how you can help and give feedback.

- This project is released with a [Code of Conduct](https://github.com/Peruz/icsd/blob/main//CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
[`LICENSE.txt`](https://github.com/Peruz/icsd/blob/main/LICENSE.md).


