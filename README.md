# BayeSRank

**BayeSRank** is an R implementation of **BayeSPeer**, a bias-aware Bayesian
rank-aggregation framework for settings with peer- and self-evaluations.
The method explicitly models self-evaluation bias and uncertainty through a
latent-utility formulation, producing aggregated rankings with posterior
uncertainty quantification.

It accompanies the manuscript *"Bias Meets Bayes: Bayesian Aggregation for Peer
and Self Ranking."* This repository provides the Gibbs sampler, the
data-generating mechanisms, the evaluation metrics, and the function library
(including the competing aggregators) used in the study.

---

## Repository structure

```
BayeSPeer/
├── R/                          # Core method API (clean, user-facing)
│   ├── utils.R                 #   rank-format conversions, init & tuning helpers
│   ├── data_generation.R       #   synthetic data generators
│   ├── metrics.R               #   Top-1 / Top-3 accuracy
│   └── BayeSRank_core.r        #   BayeSRank() sampler + bayesrank_2step() wrapper
├── analysis/
│   └── functions/              # Full function library used in the study
│       ├── 00_functions.R      #   combined sampler + helpers + competing-method wrappers
│       ├── 00_half_cauchy.R    #   half-Cauchy prior variant of the sampler
│       ├── BayeSRank_core.R    #   two-step BayeSRank wrapper
│       ├── BIRRA.R             #   competing method: BIRRA
│       └── biger.R             #   competing method: BiGER
├── examples/
│   └── example_synthetic.R     # runnable end-to-end demo on simulated data
├── BayeSPeer.Rproj
├── .gitignore
└── README.md
```

> **Note on file names.** The core file is `R/BayeSRank_core.r` (lower-case
> `.r`). On case-sensitive systems use that exact name in `source()`.

---

## Core API (`R/`)

### Bayesian rank aggregation
- **`BayeSRank()`** — main Gibbs sampler; returns posterior summaries of the
  latent performance scores, self-evaluation bias, and variance components.
- **`bayesrank_2step()`** — two-stage procedure that first runs a pilot chain
  with diffuse hyper-priors to estimate scale parameters, then fits the final
  model with empirically tuned hyper-priors.

### Data generation (simulations)
- **`gen_data_dir()`** — Normal latent effects
- **`gen_data_t()`** — heavy-tailed (t-distributed) latent effects, rescaled to unit variance
- **`gen_data_gamma()`** — skewed positive latent effects (unit variance)

Each generator simulates latent utilities, converts them to observed ranks, and
regenerates data if the resulting rankings are degenerate.

### Rank utilities & evaluation
- Conversion between rank matrices, ordinal lists, and numeric rankings
- Self-rank vs. peer-rank diagnostics; latent-utility initialization
- Top-1 and Top-3 accuracy metrics

---

## Quick start

```r
source("R/utils.R")
source("R/data_generation.R")
source("R/metrics.R")
source("R/BayeSRank_core.r")

set.seed(42)

# Simulate a small peer/self ranking data set
dat <- gen_data_dir(n = 8, mu_beta = 1.0, sigma_beta = sqrt(0.5), sigma_epsilon = 1.0)

# Fit BayeSRank with the two-step tuned procedure
res <- bayesrank_2step(rank_matrix = dat$rank, M_itrns = 2000, m_burn = 1000, seed = 42)

res$post_mean_mus                     # posterior mean latent scores
agg_rank <- rank(-res$post_mean_mus)  # aggregated ranking (rank 1 = best)

calculate_top1(agg_rank, dat$true_rank)
calculate_top3(agg_rank, dat$true_rank)
```

A fuller, commented version is in
[`examples/example_synthetic.R`](examples/example_synthetic.R).

---

## Function library (`analysis/functions/`)

`R/` holds the clean, user-facing API. `analysis/functions/` holds the complete
function library used in the study, including:

- `00_functions.R` — the combined implementation (sampler, data generators,
  metrics, and wrappers for the competing aggregators);
- `00_half_cauchy.R` — a half-Cauchy prior variant of the variance components;
- `BIRRA.R`, `biger.R` — the competing Bayesian aggregators (BIRRA, BiGER).

These reproduce the method and its competitors; the simulation, processing, and
figure scripts used to generate the manuscript's results are not distributed
here.

---

## Dependencies

Core API: `MCMCpack`, `extraDistr`, `matrixStats`.
Function library additionally uses: `expm`, `reshape2`, `tidyverse`,
`TopKLists`, `RobustRankAggreg` (for the competing aggregators).

```r
install.packages(c(
  "MCMCpack", "extraDistr", "matrixStats", "expm", "reshape2",
  "tidyverse", "TopKLists", "RobustRankAggreg"
))
```

---

## Data & privacy

No real data is distributed with this repository. The manuscript's real-world
classroom application relies on identifiable student grade and peer-ranking
records, which are **not** included here for privacy reasons. All code shipped
in this repository runs on fully synthetic, simulated data.
