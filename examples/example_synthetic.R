# ============================================================================
# BayeSRank / BayeSPeer - minimal synthetic example
# ----------------------------------------------------------------------------
# A fully self-contained demo that uses only simulated data (no real student
# records). Run this from the repository root.
# ============================================================================

# --- Load the core method API -----------------------------------------------
source("R/utils.R")            # rank-format conversions, init & tuning helpers
source("R/data_generation.R")  # synthetic data generators
source("R/metrics.R")          # Top-1 / Top-3 accuracy
source("R/BayeSRank_core.R")   # BayeSRank() sampler + bayesrank_2step() wrapper
# NOTE: the core file ships as BayeSRank_core.r (lower-case extension). On
# case-sensitive systems (Linux/GitHub Actions) use the exact name:
# source("R/BayeSRank_core.r")

set.seed(42)

# --- 1. Simulate a small peer/self ranking data set -------------------------
# n            : number of individuals (= number of rankers)
# mu_beta      : mean self-evaluation bias
# sigma_beta   : sd of self-evaluation bias across individuals
# sigma_epsilon: noise sd in the latent utilities
dat <- gen_data_dir(
  n             = 8,
  mu_beta       = 1.0,
  sigma_beta    = sqrt(0.5),
  sigma_epsilon = 1.0
)

# dat$rank      : n x n matrix of observed ranks (column j = ranker j's ranking)
# dat$true_rank : the ground-truth ranking used to generate the data
print(dat$rank)

# --- 2. Fit BayeSRank with the two-step (empirically tuned) procedure --------
res <- bayesrank_2step(
  rank_matrix = dat$rank,
  M_itrns     = 2000,   # total MCMC iterations
  m_burn      = 1000,   # burn-in
  seed        = 42
)

# --- 3. Inspect the posterior summaries -------------------------------------
res$post_mean_mus        # posterior mean latent performance scores
res$post_mean_mubeta     # posterior mean self-evaluation bias
res$post_median_Varbeta  # posterior median of the bias variance

# Aggregated ranking: higher latent score => better rank (rank 1 = best)
agg_rank <- rank(-res$post_mean_mus)
print(agg_rank)

# --- 4. Accuracy against the known ground truth -----------------------------
cat("Top-1 accuracy:", calculate_top1(agg_rank, dat$true_rank), "\n")
cat("Top-3 accuracy:", calculate_top3(agg_rank, dat$true_rank), "\n")
