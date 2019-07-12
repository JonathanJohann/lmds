

library(lmds)
library(Rtsne)
library(stats)
library(MASS)
library(CVXR)

source("utils.R")
source("get_sims.R")

# Choose the nonlinear dimensionality reduction algorithm and dataset and run! That's it.
simulate("lmds","two lines")
