
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EpiInvert

<!-- badges: start -->
<!-- badges: end -->

EpiInvert estimates time varying epidemic reproduction numbers and
restored incidence curves by inverting a renewal equation through a
variational model as described in [PNAS,
2021](https://www.pnas.org/doi/10.1073/pnas.2105112118) and [Biology,
2022](https://www.mdpi.com/2079-7737/11/4/540). EpiInvert also corrects
the administrative weekly bias in the daily registration of cases and
the bias introduced by the festive days. EpiInvert can manage daily
incidence data and weekly aggregated incidence data. This version of the
package also includes EpiInvertForecast, a learning method for the short
time forecast of the restored incidence curve.

## Vignettes

-   [EpiInvert](https://ctim.ulpgc.es/covid19/EpiInvertVignette.html) :
    A detailed description of EpiInvert R package functionalities.

-   [Rt Comparison](https://ctim.ulpgc.es/covid19/RtComparison.html) : A
    comparative analysis of the methods : EpiInvert,
    [EpiEstim](https://CRAN.R-project.org/package=EpiEstim),
    [Wallinga-Teunis](https://academic.oup.com/aje/article/160/6/509/79472)
    and [EpiNow2](https://CRAN.R-project.org/package=EpiNow2).

-   [EpiInvertForecast](https://ctim.ulpgc.es/covid19/EpiInvertForecast.html)
    : a learning method for the short time forecast of the restored
    incidence curve.

## EpiInvert Installation

You can install the development version of EpiInvert from
[GitHub](https://github.com/) with:

``` r
 install.packages("devtools")
 devtools::install_github("lalvarezmat/EpiInvert")
```

## Example

We attach some required packages

``` r
library(EpiInvert)
library(ggplot2)
library(dplyr)
library(grid)
```

Loading data on COVID-19 daily incidence up to 2022-05-05 for
[France](https://www.santepubliquefrance.fr/dossiers/coronavirus-covid-19/coronavirus-chiffres-cles-et-evolution-de-la-covid-19-en-france-et-dans-le-monde),
[Germany](https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4),
[the USA](https://ourworldindata.org/coronavirus-source-data) and [the
UK](https://ourworldindata.org/coronavirus-source-data):

``` r
data(incidence)
tail(incidence)
#>           date   FRA    DEU    USA    UK
#> 828 2022-04-30 49482  11718  23349     0
#> 829 2022-05-01 36726   4032  16153     0
#> 830 2022-05-02  8737 113522  81644    32
#> 831 2022-05-03 67017 106631  61743 35518
#> 832 2022-05-04 47925  96167 114308 16924
#> 833 2022-05-05 44225  85073  72158 12460
```

Loading some festive days for the same countries:

``` r
data(festives)
head(festives)
#>          USA        DEU        FRA         UK
#> 1 2020-01-01 2020-01-01 2020-01-01 2020-01-01
#> 2 2020-01-20 2020-04-10 2020-04-10 2020-04-10
#> 3 2020-02-17 2020-04-13 2020-04-13 2020-04-13
#> 4 2020-05-25 2020-05-01 2020-05-01 2020-05-08
#> 5 2020-06-21 2020-05-21 2020-05-08 2020-05-25
#> 6 2020-07-03 2020-06-01 2020-05-21 2020-06-21
```

Executing EpiInvert using Germany data:

``` r
res <- EpiInvert(incidence$DEU,"2022-05-05",festives$DEU)
```

Plotting the results:

``` r
EpiInvert_plot(res)
```

<img src="man/figures/README-fig1-1.png" width="100%" style="display: block; margin: auto;" />

For a detailed description of EpiInvert outcomes see the [EpiInvert
vignette](https://ctim.ulpgc.es/covid19/EpiInvertVignette.html).
