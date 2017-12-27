# natural

The `natural` package contains implementations of two methods of estimating error variance in a high-dimensional linear model.
The details of the method can be found in 
[Yu, Bien (2017) *Estimating the error variance in a high-dimensional linear model*](https://arxiv.org/abs/1712.02412).

To install `natural` from [github](http://github.com), type in R console
```R
devtools::install_github("hugogogo/natural", build_vignettes = TRUE)
```
Note that the installation above requires using R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
(which can be installed using `install.packages("devtools")`).

To see the vignette for examples of how to use the package, type in R console (after installation of the package),
```R
vignette("using_natural")
```
