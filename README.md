embvr
======
[![Travis-CI Build Status](https://travis-ci.org/DominikMueller64/embvr.svg?branch=master)](https://travis-ci.org/DominikMueller64/embvr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.377030.svg)](https://doi.org/10.5281/zenodo.377030)

embvr is a minimal [R](http://www.r-project.org) package for the estimation of
expected maximum haploid breeding values. The package contains a fast C++ routine and an
[Rcpp](http://www.rcpp.org/)-wrapper for exposing it to R.

[//]: # (TODO: Add reference to publication.)

---

### Installation

You can install embvr from its [GitHub repository](http://github.com/DominikMueller64/embvr).
You first need to install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install embvr with 

```r
devtools::install_github("DominikMueller64/embvr")
```

Windows users need to make sure that [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
is installed.

---

### Example

Please refer to the documentation of the function `embv` for an example of usage.


---

### Author

Dominik Mueller

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>
