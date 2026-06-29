# ============================================================================
# dependencies.R — package dependencies
# Source this first; loads the packages used across the library.
# Part of the BayeSPeer function library; logic unchanged from the original implementation.
# ============================================================================

library(MCMCpack)
library(extraDistr)

library(TopKLists)
library("expm")  # For matrix calculations used in the geometric mean

library(tidyverse)
library(reshape2)
library(matrixStats) #colVar
