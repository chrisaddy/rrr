# rrr
[![Build Status](https://travis-ci.org/chrisaddy/rrr.svg?branch=master)](https://travis-ci.org/chrisaddy/rrr)

R package for reduced-rank regression.

### Install latest release from CRAN

```{r}
install.packages("rrr")
```

### Install development version from github

```{r}
devtools::install_github("chrisaddy/rrr")
```

Please send any issues, with minimally-reproducible code to, and any feature requests or suggestions to the package [github issues portal](https://github.com/chrisaddy/rrr/issues).

The rrr package provides four categories of reduced-rank regression functions:
* reduced-rank regression
* principal component analysis
* canonical variate analysis
* multiclass linear discriminant analysis

The framework of this package was based on Alan Izenman's book *Modern Multivariate Statistical Techniques: Regression, Classification, and Manifold Learning* and the Fall 2016 course -- taught by him and based on the text -- at Temple University.

Special thanks to Dr. Izenman -- his insights were invaluable in keeping the focus of this package on reduced-rank regression.

The package author would also like to acknowledge the work of Charles Miller who authored the original reduced-rank regression code for the text.
