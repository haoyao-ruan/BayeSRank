library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

library(tidyverse)
load("RData/02_compare_1000_seed42_42.rdata")
results_100 <- results_multi %>% dplyr::select(c(1:6, 21, 36))%>%
  filter(rho==0.55, muBeta==1, sigmaBeta==sqrt(0.5))%>%
  mutate(distribution="N(0,1)")

load("RData/SA/02_t3_100.RData")
res_t3 <- results %>% mutate(distribution="t-3")

load("RData/SA/02_t5_100.RData")
res_t5 <- results %>% mutate(distribution="t-5")

load("RData/SA/02_t15_100.RData")
res_t15 <- results %>% mutate(distribution="t-15")

load("RData/SA/02_t25_100.RData")
res_t25 <- results %>% mutate(distribution="t-25")


all_res_t <- rbind.data.frame(res_t3,res_t5,
                             res_t15, res_t25,res_default)

# -------------------------------------------
# 1. Define and order distributions
# -------------------------------------------
levels_t <- c("t-3", "t-5", "t-15", "t-25", "N(0,1)")

# -------------------------------------------
# 2. Summarize results
# -------------------------------------------
res_t_summary <- all_res_t %>%
  group_by(n_teams, distribution) %>%
  summarise(
    mean_corr = mean(BayeSRank, na.rm = TRUE),
    mean_top1 = mean(top1_BayeSRank, na.rm = TRUE),
    mean_top3 = mean(top3_BayeSRank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    distribution = factor(distribution, levels = levels_t),
    n_teams = factor(n_teams, levels = sort(unique(n_teams))),
    mean_corr = mean_corr * 100,
    mean_top1 = mean_top1 * 100,
    mean_top3 = mean_top3 * 100
  )

# -------------------------------------------
# 3. Define custom color palette (emerald → lime → pear → basil → black)
# -------------------------------------------
# Approximate hex codes for natural green tones:
color_vals <- c("#046307", "#7AC74F", "#C3D825", "#895129", "#000000")
names(color_vals) <- levels_t

# Thinner line weights for non-black curves
line_sizes <- c(0.9, 0.9, 0.9, 0.9, 1.2)
names(line_sizes) <- levels_t

# -------------------------------------------
# 4. Common theme
# -------------------------------------------
theme_custom <- theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------------------------------------------
# 5. Function for performance panels
# -------------------------------------------
plot_panel_t <- function(yvar, ylabel) {
  ggplot(res_t_summary, aes(
    x = n_teams, y = !!sym(yvar),
    group = distribution, color = distribution, linewidth = distribution
  )) +
    geom_line() +
    geom_point(size = 2.3) +
    scale_color_manual(values = color_vals) +
    scale_linewidth_manual(values = line_sizes, guide = "none") +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    labs(x = "n", y = ylabel) +
    theme_custom
}

p2 <- plot_panel_t("mean_corr", "Spearman correlation (%)")
p3 <- plot_panel_t("mean_top1", "Top-1 coverage (%)")
p4 <- plot_panel_t("mean_top3", "Top-3 coverage (%)")

# -------------------------------------------
# 6. First panel: distribution shapes
# -------------------------------------------
curve_data_t <- data.frame()
for (dist in levels_t) {
  if (dist == "N(0,1)") {
    x <- seq(-6, 6, length.out = 400)
    density <- dnorm(x)
  } else {
    df_val <- as.numeric(gsub("t-", "", dist))
    x <- seq(-6, 6, length.out = 400)
    # Standardize to unit variance
    density <- dt(x * sqrt((df_val - 2) / df_val), df = df_val) / sqrt(df_val / (df_val - 2))
  }
  curve_data_t <- rbind(
    curve_data_t,
    data.frame(x = x, density = density, distribution = dist)
  )
}
curve_data_t$distribution <- factor(curve_data_t$distribution, levels = levels_t)

# p1 <- ggplot(curve_data_t, aes(
#   x = x, y = density, color = distribution, linewidth = distribution
# )) +
#   geom_line() +
#   scale_color_manual(values = color_vals) +
#   scale_linewidth_manual(values = line_sizes) +
#   labs(x = "x", y = "Density", 
#        title = NULL) +
#   theme_custom +
#   theme(legend.position = "none")

# Move legend inside the first panel (distribution shapes)
p1 <- ggplot(curve_data_t, aes(
  x = x, y = density, color = distribution, linewidth = distribution
)) +
  geom_line() +
  scale_color_manual(values = color_vals) +
  scale_linewidth_manual(values = line_sizes, guide = "none") +
  labs(
    x = expression(mu[i]),
    y = "Density", title = NULL, 
    color = "Distribution-df") +
  theme_custom +
  theme(
    legend.position = c(0.8, 0.73),   # adjust coordinates (x, y) inside panel
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

# Remove legends from other panels
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

# -------------------------------------------
# 7. Combine panels (1×4 layout, one legend only)
# -------------------------------------------
final_plot_t <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 4,widths = c(1.5, 1, 1, 1))

final_plot_t


ggsave(
  "robust_t.png",
  plot = final_plot_t,
  width = 16,
  height = 3.5,
  dpi = 600
)

