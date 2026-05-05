
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- The goal of MRpms is to ... -->

## Installation

You can install the development version of MRpms from
[GitHub](https://github.com/) with:

    # install.packages("remotes")
    remotes::install_github("Brigantium/MRpms")

## About this package

The main goal of this package is to performs modal regression over data
sets with an, as to date, unidimensional covariable. That is: to detect
modes of the response variable $Y$’s distribution conditional to
$X = x$, where $X$ is the covaribale.

<!-- incluir aquí varias de las referencias o, incluso, nuestro propio paper -->

The basic function included in this package is `PMS`, which performs the
estimaton over a given matrix sample using the *partial mean shift*
algorithm:

``` r
library(MRpms)
data(twosines) # a data set with two sines shifted by 3 and a normal error of sd = 0.3 
modas <- PMS(twosines, h1 = 0.5, h2 = 0.5) 
plot(twosines) 
plot(modas, pch = 19, col = 2) # a previous plot is needed
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-2-1.png" alt="Results of *partial mean shift* over the twosines data set." width="100%" />
<p class="caption">

Results of *partial mean shift* over the twosines data set.
</p>

</div>

In addition, is possible to compute confidence regions too —blue are
pointwise intervals, grey is an uniform confidence region, both to
conditional modes—:

``` r

# invisible is just to avoid a big console output
invisible(ConfPMS(twosines, modas = modas, h1 = 0.5, h2 = 0.5, B = 10))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="Confidence regions for the twosines data set. In blue, the pointwise confidence intervals for modes in each $X$ point. In grey, an uniform confidence region for all the modal mainfolds." width="100%" />
<p class="caption">

Confidence regions for the twosines data set. In blue, the pointwise
confidence intervals for modes in each $X$ point. In grey, an uniform
confidence region for all the modal mainfolds.
</p>

</div>

And, last, with `PredPMS` an uniform prediction region is computed for
the response variable $Y$ —in grey—:

``` r

invisible(PredPMS(muestra = twosines, h1 = 0.5, h2 = 0.5))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="Uniform prediccion intervals to ensure, at $95\%$, the value of Y known $X = x$." width="100%" />
<p class="caption">

Uniform prediccion intervals to ensure, at $95\%$, the value of Y known
$X = x$.
</p>

</div>

Furthermore, we programed a bandwith selector by two different methods:
a “cross validation”-like proposed by Zhou, H. and Huang, X. (2019); and
other build-up the size of prediction intervals,Chen et tal (2016). To
reduce computability times, the function performs a Nedler-Mead
heuristic algorithm.

<!-- ``` r -->

<!-- # install.packages("pak") -->

<!-- pak::pak("Brigantium/MRpms") -->

<!-- ``` -->

<!-- ## Example -->

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- library(MRpms) -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
