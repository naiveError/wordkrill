
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Wordkrill: Extending Wordfish into the Multidimensional Political Space

<!-- badges: start -->
<!-- badges: end -->

Wordkrill fits Wordkrill models (generalisation of Wordfish models to
multiple dimensions) to document-feature matrices. For additional
information, see Riesch (2025): ‘Wordkrill: Extending Wordfish into the
multidimensional political space’. <https://arxiv.org/abs/2506.20275>.

## Installation

You can install the development version of wordkrill from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("naiveError/wordkrill")
```

## Example 1

This simple example, which analyses speeches from the Irish Dáil in 2009
(included in the `quanteda.textmodels` package), demonstrates the use of
the package.

``` r
library(wordkrill)
# install.packages("quanteda.textmodels")
library(quanteda.textmodels)

toks_irish <- tokens(data_corpus_irishbudget2010, remove_punct = TRUE)
dfmat_irish <- dfm(toks_irish)
params_int_krill2irish <- initialize_urkrill_kD(dfmat_irish, K=2, ort=T)
tmod_wk2irish <- wordkrill_constr(dtm=dfmat_irish, dir=c(1,2), K=2,
                                  epsilon = 0.1, 
                                  control=list(tol=1e-08, sigma=9,
                                            startparams=params_int_krill2irish,
                                            conv.check='ll'),
                                  verbose=T, asympSE=T)
print(summary(object = tmod_wk2irish, top_n = 15))
```

## Example 2

Alternatively, the application can also be visualised using simulated
data:

``` r
library(wordkrill)

simulation <- sim_wordkrill(docs = 10, vocab =20, K=2, doclen =5000, 
                            dist = "normal")

epsilon_est <- compute_epsilon(simulation$Y, alpha = 0.95, norm = T, K=2)

test_model <- wordkrill_constr(dtm=as.dfm(simulation$Y), dir=c(10,1), K=2, 
                          epsilon = epsilon_est,
                 control=list(tol=1e-08, 
                              sigma=9, 
                              startparams=NULL, 
                              conv.check='ll'), 
                 asympSE=TRUE)
print(summaryl(object = test_model, top_n = 15))
```
