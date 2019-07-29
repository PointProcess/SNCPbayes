rm(list = ls())

library(SNCPbayes)

setwd("~/pkg/treedata/")

#Load the the data with filename "treesPP" in ShotSliceSampler
data(list = "treesPP", package = "SNCPbayes")

#Plot the locations of trees in the first dataframe in the list
plot(trees[[1]]$x,trees[[1]]$y)
