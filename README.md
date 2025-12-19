# BayeSRank

**BayeSRank** is an R implementation of a bias-aware Bayesian rank aggregation framework designed for settings with peer- and self-evaluations.\
The method explicitly models self-evaluation bias and uncertainty through a latent utility formulation, producing aggregated rankings with posterior uncertainty quantification.

This repository contains the full implementation of the Gibbs sampler, data-generating mechanisms used in simulation studies, and evaluation metrics reported in the accompanying manuscript.

------------------------------------------------------------------------

## Repository Structure

-   `R/`
    -   `utils.R` — Rank-format conversions and initialization helpers
    -   `BayeSRank_core.R` — BayeSRank Gibbs sampler and two-step tuning
    -   `data_generation.R` — Simulation data generators
    -   `metrics.R` — Ranking accuracy metrics (Top-1, Top-3)
-   `BayeSRank.Rproj`
-   `README.md`

### Bayesian Rank Aggregation

-   **`BayeSRank()`**\
    Main Gibbs sampler returning posterior summaries, aggregated rankings, and accuracy metrics.

-   **`bayesrank_2step()`**\
    Two-stage procedure that first runs a pilot chain with diffuse hyperparameters to estimate scale parameters, then fits the final model using empirically tuned hyperpriors.

### Data Generation (Simulations)

-   **`gen_data_dir()`** — Normal latent effects
-   **`gen_data_t()`** — Heavy-tailed (t-distributed) latent effects
-   **`gen_data_gamma()`** — Skewed positive latent effects

Each generator simulates latent utilities, converts them to observed ranks, and regenerates data if rankings are degenerate.

### Rank Utilities and Diagnostics

-   Conversion between rank matrices, ordinal lists, and numeric rankings
-   Diagnostics comparing self-ranks to peer-based ranks
-   Initialization routines for latent utility matrices

### Evaluation Metrics

-   Top-1 and Top-3 accuracy

------------------------------------------------------------------------

## Example Usage

\`\`\`r 
source("R/BayeSRank_core.R") 
source("R/data_generation.R")
source("R/utils.R")
source("R/metrics.R")

set.seed(42)

dat \<- gen_data_dir( n = 15, mu_beta = 0.3, sigma_beta = 0.5, sigma_epsilon = 1 )

res \<- BayeSRank_2step( n = 5, r = dat\$rank, true_rank = dat\$true_rank)

res$post_mean_mus
res$corr_spearman
\`\`\`