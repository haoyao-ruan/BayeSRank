# Reproducibility: simulation studies

This folder contains the scripts used to produce the **simulation results and
figures** reported in the manuscript *"Bias Meets Bayes: Bayesian Aggregation
for Peer and Self Ranking"* (BayeSPeer). It is organised as a numbered
pipeline.

> The real-world classroom analyses and MCMC convergence diagnostics from the
> paper are intentionally **not** included here, and no real student data is
> distributed with this repository (see the privacy note below).

## Pipeline

| Stage | Folder | What it does |
|-------|--------|--------------|
| Functions | `functions/` | Shared function library sourced by every script: the BayeSRank/BayeSPeer Gibbs sampler, the half-Cauchy prior variant, the two-step tuning wrapper, and the competing methods (BIRRA, BiGER, plus wrappers for RRA, MC1/2/3, CEMC, Stuart, summary statistics). |
| 01 | `01_simulate_data/` | Generate simulated peer/self ranking data sets across sample sizes, ranking quality (rho), bias mean (muBeta) and bias variance (sigmaBeta). Variants: default, repetitions, triangular, robustness, mixed-rho. |
| 02 | `02_run_methods/` | Run BayeSPeer and all competing methods over the simulated data (parallelised with `doParallel`/`foreach`), with and without self-ranks removed. |
| 03 | `03_process_results/` | Merge and reshape the raw result objects into tidy data frames / `.rds` files for plotting. |
| 04 | `04_figures/` | Produce the manuscript figures: bar plots, correlation bubble plots, top-k line plots, performance "envelope" plots, and robustness panels (gamma / t / triangular). |

## How the scripts find their inputs

These scripts were written to run from a **single flat working directory** and
therefore use bare relative paths, e.g.:

```r
source("00_functions.r")
load("RData/02_compare_1000_seed42_42.RData")
```

To reproduce, two things are needed that are **not** shipped in this repository:

1. **The function library on the path.** Either copy the files from
   `analysis/functions/` next to the script you are running, or change the
   `source(...)` lines to point at `analysis/functions/`. Note the core sampler
   file is `BayeSRank_core.R` here but several scripts `source("bayesrank_core.r")`
   / `source("BayeSRank_core.r")` — adjust the name/case to match your system.
2. **The intermediate `RData/` objects.** Stage 02 produces large `.RData`
   result files (multiple GB) that stages 03–04 load. These are excluded from
   version control (see `.gitignore`). Re-run stages 01–02 to regenerate them,
   or request them from the authors.

The scripts also assume the packages they `library()` at the top are installed,
including: `MCMCpack`, `extraDistr`, `matrixStats`, `data.table`,
`doParallel`, `foreach`, `progressr`, `RobustRankAggreg`, `TopKLists`,
`tidyverse`, `patchwork`, `scales`, `GGally`, `readxl`.

## Privacy note

The manuscript's real-world application uses identifiable classroom grade and
peer-ranking spreadsheets. Those data files, and the scripts that read them,
are **not** part of this repository. Only fully synthetic, generated data is
used in everything published here. See `../examples/example_synthetic.R` for a
runnable end-to-end demo on simulated data.
