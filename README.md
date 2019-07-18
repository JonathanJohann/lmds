# lmds

The **lmds** package implements local-multidimensional scaling [1] in C++ with an R wrapper and provides the dimensionality reduction quality metrics mentioned in [2].

To install the package, run the following
```R
install.packages("devtools")
library("devtools")
devtools::install_github("JonathanJohann/lmds")

```

Example Usage:
```R
library(lmds)
library(MASS)

n = 500

# Run LMDS
X <- MASS::mvrnorm(n = n, mu = c(1,-1,1,-1), Sigma = diag(4))
Y <- lmds::lmds(X, d = 2)
print(Y$X)

# Compare LMDS output with the original X
metric <- lmds::spectral(X = X, Y = Y$X, n = n-1)
print(metric)
```

*Code was developed using R version 3.6.0*

# References
[1] Chen, Lisha, and Andreas Buja. "Local multidimensional scaling for nonlinear dimension reduction, graph drawing, and proximity analysis." Journal of the American Statistical Association 104.485 (2009): 209-219.

[2] Johannemann, Jonathan, and Robert Tibshirani. "Spectral Overlap and a Comparison of Parameter-Free, Dimensionality Reduction Quality Metrics." arXiv preprint arXiv:1907.01974, 2019
