
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- The goal of MRpms is to ... -->

## Installation

You can install the development version of MRpms from
[GitHub](https://github.com/) with:

    # install.packages("remotes")
    remotes::install_github("Brigantium/MRpms")

## About this package

The main goal of this package is to perform modal regression over data
sets with a, as of now, one-dimensional covariate. That is: to detect
local modes of the response variable $Y$’s distribution conditional on
$X = x$, where $X$ is the covariate.

<!-- incluir aquí varias de las referencias o, incluso, nuestro propio paper -->

The basic function included in this package is `PMS`, which performs the
estimation over a given sample using the *partial mean shift* algorithm:

``` r
library(MRpms)
data(twosines) # a data set with two sines at different heights and a normal error with sd = 0.3 
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

A better description of this function can be consulted in the manual:
`?PMS`.

In addition, it is possible to compute confidence sets for the local
modes too —blue lines are pointwise sets, grey area is a uniform
confidence set—:

``` r

# invisible is just to avoid a big console output
invisible(ConfPMS(twosines, modas = modas, h1 = 0.5, h2 = 0.5, B = 10))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="Confidence sets for conditional modes on the twosines data set. In blue, pointwise confidence sets for local modes conditional on a mesh of $X$ points. In grey, a uniform confidence set for all the modal mainfolds." width="100%" />
<p class="caption">

Confidence sets for conditional modes on the twosines data set. In blue,
pointwise confidence sets for local modes conditional on a mesh of $X$
points. In grey, a uniform confidence set for all the modal mainfolds.
</p>

</div>

And, last, with `PredPMS` a uniform prediction set is computed for the
response variable $Y$ —in grey—:

``` r

invisible(PredPMS(muestra = twosines, h1 = 0.5, h2 = 0.5))
```

<div class="figure" style="text-align: center">

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="$95\%$ uniform prediction set for $Y$." width="100%" />
<p class="caption">

$95\%$ uniform prediction set for $Y$.
</p>

</div>

Furthermore, we programmed a bandwith selector, `bwselector`, which
implements two different methods: one taking into account the number of
estimated modes and the distance from the estimated set of modes to the
sample points, proposed by Zhou, H. and Huang, X. (2019); and other
based on the size of uniform prediction sets, as seen in Chen et
al. (2016). To reduce computing times, the function performs a
Nedler-Mead heuristic search.

In case the user wants a better control over the estimation, like
choosing the maximum number of modes or the points where they will be
estimated, we also created the ‘mallador’ function. This function can be
used to build a mesh with either equidistant or user provided $X$
points. Every function admits a `malla` argument, but this must be
constructed by `mallador` function.

``` r

# equidistant
malla <- mallador(twosines, k = 10, len = 100) # k: maximum number of modes. 
                                               # len: number of covariate points where the estimation will be computed

# user provided
# malla <- mallador(twosines, x.malla = c(0,0.1,0.2,0.4,0.55, 0.6,...), k = 10)
```

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
