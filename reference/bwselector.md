# Bandwith selector

Performs an heuristic search based on the Nedler-Mead algorithm to find
appropiate values for the bandwith parameters `h1` and `h2`. (Warning:
this is a computationally requiring method which might take some time to
finish).

## Usage

``` r
bwselector(
  data,
  malla = mallador(data, x.malla = X),
  dim.y = 1,
  eps = 1e-08,
  conf.level = 0.95,
  method = "CV",
  eps2 = 0.001,
  maxiter = 100,
  maxiter.tol = 10,
  alpha = 1,
  gamma = 2,
  rho = 1/2,
  sigma = 1/2,
  initial.hs
)
```

## Arguments

- data:

  Matrix or data frame containing the sample. Currently, the response
  values are needed to be in the last column.

- malla:

  Matrix containing the initial points used to compute the modes using
  the Partial Mean-Shift algorithm. If not provided, `mallador`function
  is called. It must be a MRpms_malla class object.

- dim.y:

  Response dimension.

- eps:

  Convergence tolerance for Partial Mean-Shift algorithm.

- conf.level:

  Confidence level for prediction sets (if `method = P`).

- method:

  `CV` (by default) for Zhou and Huang selector or `P` for Chen et al.
  selector based on the size of prediction sets.

- eps2:

  Convergence tolerance for heuristic algorithm.

- maxiter:

  Maximum number of iterations for heuristic algorithm.

- maxiter.tol:

  Maximum number of consecutive iterations without improvement in
  heuristic algorithm before termination.

- alpha:

  Parameter of reflection in the Nelder-Mead algorithm.

- gamma:

  Parameter of expansion in the Nelder-Mead algorithm.

- rho:

  Parameter of contraction in the Nelder-Mead algorithm.

- sigma:

  Parameter of shrink in the Nelder-Mead algorithm.

- initial.h:

  Candidates to construct the initial simplex in the Nelder-Mead
  algorithm. Currently, there are needed only three diferent candidates
  saves by rows.

## Value

A vector with two elements, `h1` and `h2`, first one been for covariable
and second one for response.

## Details

This function allows a Nelder-Mead based search of smooth parameter.
There is two diferent measurments:

- If `method = CV` then the algorithm uses the cross-validation like
  formula suggested by Zhou and Huang as degree of fitting,

  \$\$f(\mathbf h) = \displaystyle \frac{1}{n}\sum\_{i=1}^n
  d^2\left(\widehat{M}\_{n,\mathbf{h},-i}\left(X_i\right),Y_i\right)\widehat{N}\_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right),\$\$

  where \\\widehat M\_{n,\mathbf{h},-i}(X_i)\\ is the estimated modes in
  \\X_i\\ without consider \\X_i\\ and
  \\\widehat{N}\_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right)\\ is
  the number of modes estimated in \\X_i\\.

- If `method = P` then uses the volume of prediction sets instead as
  measurement of fitting, \$\$f(\mathbf h) = \displaystyle
  Vol\left(\widehat{\mathcal{P}}\_{1-\alpha,\mathbf{h}}\right)=\widehat{\varepsilon}\_{1-\alpha,\mathbf{h}}\int\_{x\in
  D}\widehat{N}\_{\mathbf{h}}(x)\\ dx,\$\$ where \\\mathcal {\widehat
  P}\_{1-\alpha, \mathbf {h}}\\ is the estimation of the prediction set,
  i.e. a set for which \\\mathbb{P}(Y\in \mathcal P \geq 1-\alpha)\\;
  \\\hat \varepsilon\_{1-\alpha,\mathbf h}\\ is the estimation of the
  \\1-\alpha\\ quantile of the random variable \\d(Y,M(X))\\, computed
  as the sample quantile of the sample \\\\d(Y_i,\widehat{M}\_{n,\mathbf
  h}(X_i))\\\colon i \in \\1,\dots,n\\\\\\; and \\\widehat{N}\_{\mathbf
  h}(x)\\ is the estimated number of modes condicional to \\X = x\\.

### **Nelder-Mead algorithm:**

Selected the ranking method, the algorithm preforms a Nelder-Mead
heuristic strategy to minimize \\f\\. This proccess starts with 3
different candidates to \\\mathbf h\\: \\\mathbf h_1\\, \\\mathbf h_2\\
and \\\mathbf h_3\\. In each one computes \\f(\mathbf h_i)\\. Now, the
idea is the same in each iteration:

1.  Firstly, the \\f(\mathbf h_i)\\ values are ordered, without lost of
    generality: \\f(\mathbf h_1) \leq f(\mathbf h_2) \leq f(\mathbf
    h_3)\\.

2.  the centroid of \\\mathbf h_1\\ and \\\mathbf h_2\\ is computed
    (\\\mathbf h_0\\), i.e., the mean point in coordinate terms —plus
    certain tolerance to avoid zeros—.

    NOTE: *the original Nelder-Mead algorithm allows only unrestricted
    spaces, in other words the entirety of \\\mathbb R^n\\. To bypass
    this, we previously transform through a natural logarithm function
    the espace where the smooth parameters lives,
    \\(\mathbb{R}^{+})^n\\, to \\\mathbb R^n\\. So, we will work always
    with \\\log \mathbf h\\ except when evaluate \\f\\.*

3.  **Reflection**: we compute a new candidate: \\\mathbf h_r = \mathbf
    h_0 + \alpha(\mathbf h_0 - \mathbf h_3)\\, with \\\alpha \> 0\\, and
    compute\\f(\mathbf h_r)\\:

    - if \\f(\mathbf h_1) \< f(\mathbf h_r) \< f(\mathbf h_2)\\ (not
      enough improvement): replace \\\mathbf h_3\\ with \\\mathbf h_r\\
      and go to (1.).

    - if \\f(\mathbf h_r) \< f(\mathbf h_1)\\ (improvement): go to (4.).

    - if \\f(\mathbf h_2) \<f(\mathbf h_r)\\ (no improvement): go to
      (5.).

    NOTE: *we use \\\alpha = 1\\*.

4.  **Expansion**: we compute a new candidate: \\\mathbf h_e = \mathbf
    h_0 + \gamma(\mathbf h_r - \mathbf h_0)\\, with \\\gamma \>1\\, then

    - if \\f(\mathbf h_e) \< f(\mathbf h_r)\\ then replace \\\mathbf
      h_3\\ with \\\mathbf h_e\\ and go to (1.).

    - if \\f(\mathbf h_r) \< f(\mathbf h_e)\\ then replace \\\mathbf
      h_3\\ with \\\mathbf h_r\\ and go to (1.).

    NOTE: *we use \\\gamma = 2\\*.

5.  **Contraction**: two situations:

    - if \\f(\mathbf h_r) \< f(\mathbf h_3)\\ then compute a new
      candidate: \\\mathbf h_c = \mathbf h_0 + \rho (\mathbf h_r -
      \mathbf h_0)\\, with \\0\<\rho \leq 0.5\\.

    - if \\f(\mathbf h_3) \< f(\mathbf h_r)\\ then compute a new
      condadte: \\\mathbf h_c = \mathbf h_0 + \rho (\mathbf h_3 -
      \mathbf h_0)\\, with \\0 \< \rho \leq 0.5\\.

    In case of improvement, i.e. \\f(\mathbf h_c) \< \min\\f(\mathbf
    h_3),f(\mathbf h_r)\\\\, replace \\\mathbf h_3\\ with \\\mathbf
    h_c\\. In other case, go to (6.).

    NOTE: *we use \\\rho = 1/2\\*.

6.  **Shrink**: Fixed \\\mathbf h_1\\, replace the others with \\\mathbf
    h_i = \mathbf h_1 + \sigma (\mathbf h_i - \mathbf h_1)\\. Returns to
    (1.).

    NOTE: *we take \\\sigma = 1/2\\*.

**Exit criteria:** The heuristic strategy ends when anyone of this
statements is true:

- The algorithm performs `maxiter` iterations or there is no improvement
  in `maxiter.tol` consecutive iterations.

- The standard deviation of \\f(\mathbf h_1),\\ f(\mathbf h_2)\\ and
  \\f(\mathbf h_3)\\, reescalated by the value of \\f(\mathbf h_1)\\ is
  less than a certain tolerance `eps`.

Further details could be reviewed in Luersen, Le Riche and Guyon (2004).

## References

Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L.
(2016). *Nonparametric modal regression*. The Annals of Statistics,
**44**(2), 489–514.

Zhou, H. and Huang, X. (2019). *Bandwidth selection for nonparametric
modal regression*. Communication in Statistics - Simulation and
Computation, **48**(4), 968–984.

Luersen, M., Le Riche, R. and Guyon, F. (2004). *A constrained,
globalized, and bounded Nelder–Mead method for engineering
optimization.* Structural and Multidisciplinary Optimization, **27**(1).
43–54.

## Examples

``` r
h <- bwselector(twosines)
modas <- PMS(twosines, h1 = h[1], h2 = h[2])

plot(modas,twosines, pch = 19, col = "red")

```
