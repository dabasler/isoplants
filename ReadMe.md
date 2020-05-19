
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!--# isoplants <img src="man/figures/logo.png" align="right" height="200" width="200"/> -->
Modelling stable isotope composition of plant tissue
----------------------------------------------------

#### *Note*

This package is under development and may yet produce unstable results

Description
-----------

{isoplants} is a lightweight R package to model the isotopic ratios in in plant tissue and analyze the sensitivity to changes in environmental conditions. It currently focuses on stable Oxygen isotopes but will integrate Carbon and possibly Hydrogen isotopes in future releases. It uses the R package [tealeaves](https://CRAN.R-project.org/package=tealeaves) to integrate the calculation of leaf temperatures.

Get isoplants
-------------

From GitHub

``` r
install.packages("devtools")
devtools::install_github("dabasler/isoplants")
```

And load isoplants

``` r
library("isoplants")
```

Vignette
--------

The {isoplants} package allows to calculate isotopic composition of plant tissue in respone to environmental factors. This vignette shows a basic example to usse the {isoplants} package:

-   run a minimum worked example using default parameters
-   replace default parameters
-   include uncertainty estimates for specific parameters

Minimum worked example
----------------------

You can use the models with the default parameter settings using the `get_default_parameters()` function and `plant18O_model()`. Basic information about the inputparmeters can be displayed by calling `get_parameter_definition()`.

``` r

library(magrittr)
library(isoplants)

# Get all default parameters (the default parameters will use include parameters for the peclet model)
df_parameter<-get_default_parameters()
result<-plant18O_model(df_parameter)
result %>% knitr::kable()
```

|        ea|        ei|  ea\_ei|        vpd|        eq|        ek|   gs|          E|    D|         pn|   D18O\_e|    d18O\_e|  D18O\_lw|  d18O\_lw|   D18O\_c|   d18O\_c|  d18O\_pt|  d18O\_pt|
|---------:|---------:|-------:|----------:|---------:|---------:|----:|----------:|----:|----------:|---------:|----------:|---------:|---------:|---------:|---------:|---------:|---------:|
|  1.642629|  2.346613|     0.7|  0.7039839|  9.806829|  25.42857|  0.4|  0.0019851|    0|  0.5429544|  10.57811|  -9.633452|  8.162601|  -11.8374|  31.89756|  11.25961|  11.25961|  11.25961|

``` r


# The output can also be limited to a specific variable
plant18O_model(df_parameter,output = 'd18O_c')
#> [1] 11.25961
```

Parameter Checking
------------------

The package includes the `check_parameters()` function which provides (1) a list of selected model options, (2) basic checking of provided parameters and, (3) a list of any errors found.

``` r

check<-check_parameters(df_parameter) # Run a check on the specified parameters

check$model_options %>% knitr::kable()
```

<table>
<colgroup>
<col width="100%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">model_options</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Provided 1 set(s) for 12 parameters: Tair, RH, P, d18O_sw, Tleaf, rb, gs, Lm, px, pex, ecp, D18O_wv</td>
</tr>
<tr class="even">
<td align="left">Mixing: peclet mixing</td>
</tr>
<tr class="odd">
<td align="left">Other: D18O of water vapor over soil water assumed to be in equilibrium</td>
</tr>
</tbody>
</table>

``` r
check$model_parameters %>% knitr::kable()
```

|     | group       | name     |     lower|  upper|     default| unit            | description                                 | check  |  valid|  invalid|   na| range      | type       |
|-----|:------------|:---------|---------:|------:|-----------:|:----------------|:--------------------------------------------|:-------|------:|--------:|----:|:-----------|:-----------|
| 1   | environment | Tair     |  -30.0000|     50|   20.000000| \[deg C\]       | Air Temperature                             | passed |      1|        0|    0| 20.000000  | constant   |
| 2   | environment | RH       |    0.0000|    100|   70.000000| \[%\]           | relative humidity                           | passed |      1|        0|    0| 70.000000  | constant   |
| 4   | environment | P        |   40.0000|    110|  101.325000| \[kPa\]         | barometric pressure                         | passed |      1|        0|    0| 101.325000 | constant   |
| 6   | environment | d18O\_sw |  -30.0000|      0|  -20.000000| \[permil\]      | d18O soil water                             | passed |      1|        0|    0| -20.000000 | constant   |
| 7   | environment | D18O\_wv |  -15.0000|      0|   -9.806829| \[permil\]      | D18O of water vapor over soil water         | passed |      1|        0|    0| -9.806829  | calculated |
| 9   | leaf        | Tleaf    |  -30.0000|     50|   20.000000| \[deg C\]       | absolute Leaf temperature                   | passed |      1|        0|    0| 20.000000  | constant   |
| 10  | leaf        | rb       |    0.4000|      6|    1.000000| \[m2 s mol-1\]  | boundary resistance                         | passed |      1|        0|    0| 1.000000   | constant   |
| 11  | leaf        | gs       |    0.0000|      2|    0.400000| \[mol m-2 s-1\] | stomatal conductance                        | passed |      1|        0|    0| 0.400000   | constant   |
| 13  | leaf        | Lm       |    0.0001|      2|    0.030000| \[m\]           | Peclet-model scaled path length             | passed |      1|        0|    0| 0.030000   | constant   |
| 15  | leaf        | px       |    0.0000|      1|    0.400000| \[\]            | conversion parameter for cellulose          | passed |      1|        0|    0| 0.400000   | constant   |
| 16  | leaf        | pex      |    0.0000|      1|    1.000000| \[\]            | conversion parameter for cellulose          | passed |      1|        0|    0| 1.000000   | constant   |
| 17  | leaf        | ecp      |  -10.0000|     10|    0.000000| \[permil\]      | offset from cellulose to bulk leaf material | passed |      1|        0|    0| 0.000000   | constant   |

``` r
check$errors %>% knitr::kable() # should be empty
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
</td>
</tr>
</tbody>
</table>
Replace parameters
------------------

The input parameter data.frame object can have multiple rows. Here we calculate d18O\_lt for different values of RH.

``` r
# Run model on default parameters, manipulating one parameter (RH)
RH<-seq(10,90,10)
# PArameters can be set using the set_parameter function
df_parameter<-get_default_parameters()
df_parameter<-set_parameters(data.frame(RH),df_parameter)
# OR by assigning them directly to a parametertable with the appropriate number of rows
df_parameter<-get_default_parameters(n=length(RH), mode='peclet')
df_parameter$RH<-RH
#run model
result<-plant18O_model(df_parameter)
plot(df_parameter$RH,result$d18O_pt,xlab='RH',ylab='d18O_lt',pch=19,las=1)
```

<img src="man/figures/README-cpar-1.png" width="100%" />

Parameter uncertainty
---------------------

The included fuction `get_randomized_parameters()` samples n parametersets from parameter specific distributions in order to estimate the effect of uncertainty of input parameters. `set_parameters()` combines fixed parameters with a set of varying parameters.

``` r
# Run model on default parameters with uncertainty on multiple parameter

# Get default parameters
df_parameter<-get_default_parameters(mode='peclet')
# Defines ranges and distributions for parameters (see ?get_randomized_parameters for detailed information on the imput format and further options)
rnd_par <- data.frame(name= c('RH',   'd18O_sw', 'px'   ),
                      pdist=c('norm', 'norm'   , 'unif' ),
                      pdm=  c(70,      -10     , 0.4    ),
                      pdv=  c(20,       5      , 0.2    ),
                      stringsAsFactors = FALSE
                      )

n<-1000 # Number of samples

df_parameter <- set_parameters(get_randomized_parameters(n,rnd_par) , df_parameter)
result<-plant18O_model(df_parameter)
```

``` r
hist(result$d18O_pt,main = 'd18O_pt',xlab='d18O plant tissue')
boxplot(result$d18O_pt,main = 'd18O_pt')
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

Vignette Sensititivy
--------------------

The models in the {isoplants} package are include various environmental and leaf specific input parameters. The package includes functions to analyze the sensitivity to changes in the input variables This vignette shows some basic example to do some seinsitivity analysis with {isoplants}:

Direct sensitivity analysis
---------------------------

For example, we can look at how d18O of cellulose is affected by differnt pex values under a set of differnet RH conditions.

``` r
library(isoplants)

pex<-seq(0.1,0.9,0.01)
lf18Opar<-get_default_parameters(n=length(pex))

lf18Opar$pex<-pex
result<-plant18O_model(lf18Opar)
```

``` r

plot(lf18Opar$pex,result$d18O_c,ylim=c(0,25),xlim=c(0,1),type='l',las=1,xlab='pex',ylab='d18O cellulose')
rhs<-seq(10,90,10)
col<-rainbow(length(rhs))
for (i in 1:length(rhs)){
  lf18Opar$RH<-rhs[i]
  result<-plant18O_model(lf18Opar)
  lines(lf18Opar$pex,result$d18O_c,col=col[i])
}
legend('bottomleft',as.character(rhs),col=col,lty=1,title = 'RH (%)',bty='n')
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Sobol incices
-------------

Sobol incices are a Global form of sensitivity indices whereby the variance of the output of the model is decomposed into fractions which can be attributed to individual input parameters or from the interaction between parameters. Here we analyze the sensitivity of d18O cellulose to changes in the input parameters.

``` r
# Sensitivity analysis (with sobol indices)

library('sensitivity')
library('boot')
library('ggplot2')

# Define parameter ranges (min-max)
rnd_par <- data.frame(
  name  = c("Tair","DTleaf","RH",  "P", "d18O_sw","D18O_wv","rb", "Lm", "gsref", "px", "pex"),
  pdist = 'minmax',
  pdm   = c(     5,     -5,  30,   90.0,       -30,     -20, 0.4, 0.001,    0.01 ,  0.8,  0.4),
  pdv    = c(    30,      10,  95, 102.0 ,       0,       0,   5,    0.5,     0.6 ,  0.8,  0.4),
  stringsAsFactors = FALSE
)
#Initialize two random parametersets
n<-10000
X1 <- get_randomized_parameters(n,rnd_par)
X2 <- get_randomized_parameters(n,rnd_par)

x <- sobol2007(model =plant18O_model, X1 = X1, X2 = X2, nboot = 1000, output = "d18O_c")
```

``` r
ggplot(x)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Fixing some parameters

The `plant18O_model()` model allows to pass some parameters with the `addpar` argument. These will be internally joined with the main parameters, but allow for the exclusion of these parameters in certain kinds of analysis. Here we analyze the sensitivity of d18O leafwater to changes in the input parameters (excluding variation in soil water, and the parameters for the steps after the calculation of leaf water).

``` r
# Calculate Sobol indices with a set of fixed parameters (set to default values)
fixpar <- get_default_parameters()[,c('d18O_sw','Lm','px','pex','ecp')]
n<-10000
# Using the previously defined parameter ranges (rnd_par), excluding the fixed parameters
X1 <- get_randomized_parameters(n,rnd_par)[,!names(fixpar) %in% names (rnd_par)]
X2 <- get_randomized_parameters(n,rnd_par)[,!names(fixpar) %in% names (rnd_par)]
X1<-X1[,!names(X1) %in% names (fixpar)]
X2<-X2[,!names(X2) %in% names (fixpar)]

x <- sobol2007(model =plant18O_model, X1 = X1, X2 = X2, nboot = 1000, output = "d18O_lw",addpar=fixpar)
```

``` r
ggplot(x)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

### Parmeter gradients

In order to analyze the importance of multiple input parameters across a gradient of a certain parameter, we can use the `scan_sensitivity()` function

``` r
n<-10000
X1 <- get_randomized_parameters(n,rnd_par)
X2 <- get_randomized_parameters(n,rnd_par)
fixpar <- get_default_parameters()[,c('d18O_sw','P','px','pex')]
scanpar<-"Tair"
ts<-seq(0,30,2)
ss<-scan_sensitivity(X1,X2,fixpar,ts,scanpar,outpar="D18O_lw")
#> [1] 0
#> [1] 2
#> [1] 4
#> [1] 6
#> [1] 8
#> [1] 10
#> [1] 12
#> [1] 14
#> [1] 16
#> [1] 18
#> [1] 20
#> [1] 22
#> [1] 24
#> [1] 26
#> [1] 28
#> [1] 30
```

``` r
plot_sensitivity(ss$ds,scanpar,'normalized first-order indices','D18O_lw',parcol = 1)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
plot_sensitivity(ss$dt,scanpar,'normalized total indices','',parcol = 1,legend = 2)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Contributors
------------

-   [David Basler](https://github.com/dabasler)

<!--
## Comments and contributions

#I welcome comments, criticisms, and especially contributions!
#GitHub issues are the preferred way to report bugs, ask questions, or request new features.
#You can submit issues here:
#https://github.com/dabasler/isoplants/issues
-->
Meta
----

<!--
# Please [report any issues or bugs](https://github.com/dabasler/isoplants/issues).
-->
-   License: MIT

<!-- * Get citation information for `isoplants` in R doing `citation(package = 'isoplants')` 
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
-->
