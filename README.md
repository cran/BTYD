# BTYD 2.4.2

A patch for the `pnbd.LL()` function in the original [BTYD](https://CRAN.R-project.org/package=BTYD) package, 
first proposed by [Theo Strinopoulos](https://github.com/theofilos). Now on [CRAN](https://CRAN.R-project.org/package=BTYD), 
so you can get it the usual way, with `install.packages()`.

## Justification

In its original version the Pareto/NBD (pnbd) part of the BTYD package failed for me, as it did for Theo, for reasons he explained 
[here](https://github.com/theofilos/BTYD). So, I implemented his fix, rebuilt from source, and then it worked. This is the short
version. The long version is that in the process of implementing the fix I made some changes to the choice of optimization routine 
(now using `optimx` as opposed to `optim`), Gaussian hypergeometric (now you have the option of using the `hypergeo` package) and 
I found some opportunities for refactoring functions defined in both the `pnbd` and the `bgnbd` (Beta-Geometric Negative Binomial) 
groups. I explained these changes are explained in separate documentation (fix_pnbd.html). Finally, I rebuilt BTYD
and checked that it would return the same numbers as BTYD 2.4 did when running the examples in the original BTYD vignette (see 
threeway_walkthrough.R).

## How to install from here

First `git clone`. 

If you use [`devtools`](https://devtools.r-lib.org), then at the R console just call `document(); build(); install(build_vignettes = TRUE)`. Done.

Otherwise, build the standard way in two steps:

1. At the command line, build and check the source tarball:

```
R CMD build BTYD
R CMD check BTYD_2.4.2.tar.gz
```

2. Then, at the R console, install it:

```
install.packages("BTYD2_2.4.2.tar.gz", repos = NULL, type = "source")
# works only in RStudio:
.rs.restartR()
```
