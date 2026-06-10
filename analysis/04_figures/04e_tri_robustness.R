load("RData/02_tri_muBeta1.RData")

load("RData/02_tri_muBeta05.RData")
res_05 <- results
res_05$center <- 0.5

load("RData/02_tri_muBeta15.RData")
res_15 <- results
res_15$center <- 1.5

load("RData/02_tri_muBeta5.RData")
res_5 <- results
res_5$center <- 5

library(tidyverse)
load("RData/02_compare_1000_seed42_42.RData")
res_uni <- results_multi%>% 
  filter(muBeta==1, rho==0.55, sigmaBeta==sqrt(0.5))%>%
  select(c(1:6, 21, 36))
res_uni$center <- "uni"

tri_res <- rbind.data.frame(res_05, res_1, res_15, res_5, res_uni)
tri_res$mu_scale =tri_res$center 

library(ggplot2)
library(dplyr)

# Settings
sigma_beta <- sqrt(0.5)
mu_scales <- c(0.5, 1, 1.5, 5)

x_grid <- seq(-10, 10, length.out = 3000)

# Default prior: beta_j ~ N(mu_beta = 1, sigma_beta^2 = 0.5)
normal_df <- data.frame(
  beta = x_grid,
  density = dnorm(x_grid, mean = 1, sd = sigma_beta),
  distribution = "Default normal"
)

# Tri-modal beta distributions
tri_df <- do.call(rbind, lapply(mu_scales, function(mu_beta){
  data.frame(
    beta = x_grid,
    density =
      0.25 * dnorm(x_grid, mean = -mu_beta, sd = sigma_beta) +
      0.25 * dnorm(x_grid, mean = 0,        sd = sigma_beta) +
      0.50 * dnorm(x_grid, mean =  mu_beta, sd = sigma_beta),
    distribution = paste0("Tri-modal, ", mu_beta)
  )
}))

plot_df <- rbind(normal_df, tri_df)

plot_df$distribution <- factor(
  plot_df$distribution,
  levels = c(
    "Tri-modal, 0.5",
    "Tri-modal, 1",
    "Tri-modal, 1.5",
    "Tri-modal, 5",
    "Default normal"
  )
)

cols <- c(
  "Tri-modal, 0.5" = "#6BAED6",
  "Tri-modal, 1"   = "#4A90E2",
  "Tri-modal, 1.5" = "#7B6FD6",
  "Tri-modal, 5"   = "#5E3C99",
  "Default normal" = "black"
)

labs_legend <- c(
  "Tri-modal, 0.5" = expression(paste("Tri-modal, ", eta, "=0.5")),
  "Tri-modal, 1"   = expression(paste("Tri-modal, ", eta, "=1")),
  "Tri-modal, 1.5" = expression(paste("Tri-modal, ", eta, "=1.5")),
  "Tri-modal, 5"   = expression(paste("Tri-modal, ", eta, "=5")),
  "Default normal" = expression(paste("N(", eta, "=1, ", sigma[beta]^2, "=0.5)"))
)

p_density  =ggplot(plot_df, aes(x = beta, y = density, color = distribution)) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(
    values = cols,
    labels = labs_legend
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 0.6)) +
  labs(
    x = expression(beta[j]),
    y = "Density",
    color = expression(paste("Distribution of ", beta[j]))
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.74),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(data.table)
library(ggplot2)
library(patchwork)

setDT(tri_res)

# ============================================================
# 1. Settings
# ============================================================

sigma_beta <- sqrt(0.5)
mu_scales <- c(0.5, 1, 1.5, 5)

# Replace these with your actual metric column names
metric_cols <- c(
  "BayeSRank",
  "top1_BayeSRank",
  "top3_BayeSRank"
)

metric_labels <- c(
  BayeSRank = "Spearman correlation (%)",
  top1_BayeSRank = "Top-1 coverage (%)",
  top3_BayeSRank = "Top-3 coverage (%)"
)

# ============================================================
# 2. Panel 1: density curves
# ============================================================

x_grid <- seq(-10, 10, length.out = 3000)

normal_df <- data.frame(
  beta = x_grid,
  density = dnorm(x_grid, mean = 1, sd = sigma_beta),
  distribution = "Default normal"
)

tri_df <- do.call(rbind, lapply(mu_scales, function(mu_beta) {
  data.frame(
    beta = x_grid,
    density =
      0.25 * dnorm(x_grid, mean = -mu_beta, sd = sigma_beta) +
      0.25 * dnorm(x_grid, mean = 0,        sd = sigma_beta) +
      0.50 * dnorm(x_grid, mean =  mu_beta, sd = sigma_beta),
    distribution = paste0("Tri-modal, ", mu_beta)
  )
}))

plot_df <- rbind(normal_df, tri_df)

plot_df$distribution <- factor(
  plot_df$distribution,
  levels = c(
    "Tri-modal, 0.5",
    "Tri-modal, 1",
    "Tri-modal, 1.5",
    "Tri-modal, 5",
    "Default normal"
  )
)

density_cols <- c(
  "Tri-modal, 0.5" = "#6BAED6",
  "Tri-modal, 1"   = "#4A90E2",
  "Tri-modal, 1.5" = "#7B6FD6",
  "Tri-modal, 5"   = "#5E3C99",
  "Default normal" = "black"
)

density_labs <- c(
  "Tri-modal, 0.5" = expression(paste("Tri-modal, ", eta, "=0.5")),
  "Tri-modal, 1"   = expression(paste("Tri-modal, ", eta, "=1")),
  "Tri-modal, 1.5" = expression(paste("Tri-modal, ", eta, "=1.5")),
  "Tri-modal, 5"   = expression(paste("Tri-modal, ", eta, "=5")),
  "Default normal" = expression(paste("N(", mu[beta], "=1, ", sigma[beta]^2, "=0.5)"))
)

p_density <- ggplot(plot_df, aes(x = beta, y = density, color = distribution)) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(values = density_cols, labels = density_labs) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 0.6)) +
  labs(
    x = expression(beta[j]),
    y = "Density",
    color = expression(paste("Distribution of ", beta[j]))
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = c(0.2, 0.65),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_density
# ============================================================
# 3. Prepare long data for metric panels
#    center values: "0.5", "5", and "uni"
# ============================================================

tri_res[, center := as.character(center)]

tri_long <- melt(
  tri_res,
  id.vars = c("n_teams", "center"),
  measure.vars = metric_cols,
  variable.name = "metric",
  value.name = "value"
)

setDT(tri_long)

tri_long[, metric := metric_labels[as.character(metric)]]
tri_long[, value := 100 * value]

tri_long[, dist_lab := fifelse(
  center == "uni",
  "Default normal",
  paste0("Tri-modal, eta = ", center)
)]

tri_long[, dist_lab := factor(
  dist_lab,
  levels = c(
    "Tri-modal, eta = 0.5",
    "Tri-modal, eta = 5",
    "Default normal"
  )
)]

metric_cols_plot <- c(
  "Tri-modal, eta = 0.5" = "#6BAED6",
  "Tri-modal, eta = 5"   = "#5E3C99",
  "Default normal"             = "black"
)


# Important: average over simulation replicates
tri_summ <- tri_long[
  ,
  .(
    mean_value = mean(value, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(.N)
  ),
  by = .(n_teams, center, dist_lab, metric)
]
# ============================================================
# 4. Metric panels
# ============================================================

make_metric_panel <- function(dat, metric_name, ylab) {
  ggplot(
    dat[metric == metric_name],
    aes(x = n_teams, y = mean_value, color = dist_lab, group = dist_lab)
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.1) +
    scale_color_manual(values = metric_cols_plot, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      labels = function(x) paste0(x, "%")
    ) +
    scale_x_continuous(breaks = sort(unique(dat$n_teams))) +
    labs(x = expression(n), y = ylab, color = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p_spear <- make_metric_panel(
  tri_summ,
  "Spearman correlation (%)",
  "Spearman correlation (%)"
)

p_top1 <- make_metric_panel(
  tri_summ,
  "Top-1 coverage (%)",
  "Top-1 coverage (%)"
)

p_top3 <- make_metric_panel(
  tri_summ,
  "Top-3 coverage (%)",
  "Top-3 coverage (%)"
)

# ============================================================
# 5. Combine panels
# ============================================================

p_tri_modal <- p_density + p_spear + p_top1 + p_top3 +
  plot_layout(widths = c(1.5, 1, 1, 1))

p_tri_modal

ggsave(
  "robust_tri.png",
  plot = p_tri_modal,
  width = 16,
  height = 3.5,
  dpi = 600
)

