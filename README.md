# lmds

The **lmds** package implements local-multidimensional scaling [1] in C++ with an R wrapper and provides the quality metrics introduced in [2].

To install the package, run the following
```R
install.packages(devtools)
library(devtools)
devtools::install_github("JonathanJohann/lmds")

```

Example Usage:
```R
library(lmds)
library(MASS)

n = 500

# Run LMDS
X <- MASS::mvrnorm(n = n, mu = c(1,-1,1), Sigma = diag(3))
Y <- lmds::lmds(X)
print(Y$X)

# Compare LMDS output with the original X
metric <- lmds::spectral(X = X, Y = Y$X, n = n-1)
print(metric)
```

# References
[1] Chen, Lisha, and Andreas Buja. "Local multidimensional scaling for nonlinear dimension reduction, graph drawing, and proximity analysis." Journal of the American Statistical Association 104.485 (2009): 209-219.
[2] *Coming soon*
