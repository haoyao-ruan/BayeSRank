library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

load("RData/02_kappa1.rdata")
res_kappa1 <- results
res_kappa1$distribution <- "G(1,1)=Exp(1) – 2.00"

load("RData/02_kappa5.rdata")
res_kappa5 <- results
res_kappa5$distribution <- "G(5, √5) – 0.89"

load("RData/02_kappa10.rdata")
res_kappa10 <- results
res_kappa10$distribution <- "G(10, √10) – 0.63"

load("RData/02_kappa20.rdata")
res_kappa20 <- results
res_kappa20$distribution <- "G(20, √20) – 0.45"

library(tidyverse)
load("RData/02_compare_1000_seed42_42.rdata")
res_default <- results_multi %>% select(c(1:6, 21, 36))%>%
  filter(rho==0.55, muBeta==1, sigmaBeta==sqrt(0.5))%>%
  mutate(kappa=NA,
         distribution="N(0,1) – 0.00")

all_res_gamma<- rbind.data.frame(res_kappa1,res_kappa5,
                                 res_kappa10, res_kappa20,
                                 res_default)


# -------------------------------------------
# 1. Define and order distributions (descending skewness)
# -------------------------------------------
levels_gamma <- c("G(1,1)=Exp(1) – 2.00",
                  "G(5, √5) – 0.89",
                  "G(10, √10) – 0.63",
                  "G(20, √20) – 0.45",
                  "N(0,1) – 0.00")

# Set factor levels in the desired descending order of skewness
all_res_gamma$distribution <- factor(
  all_res_gamma$distribution,
  levels = levels_gamma
)
# -------------------------------------------
# 2. Summarize results
# -------------------------------------------
res_gamma_summary <- all_res_gamma %>%
  group_by(n_teams, distribution) %>%
  summarise(
    mean_corr = mean(BayeSRank, na.rm = TRUE),
    mean_top1 = mean(top1_BayeSRank, na.rm = TRUE),
    mean_top3 = mean(top3_BayeSRank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    distribution = factor(distribution, levels = levels_gamma),
    n_teams = factor(n_teams, levels = sort(unique(n_teams))),
    mean_corr = mean_corr * 100,
    mean_top1 = mean_top1 * 100,
    mean_top3 = mean_top3 * 100
  )

# -------------------------------------------
# 3. Define color gradient (chilli → salmon → amber → grey → black)
# -------------------------------------------
color_vals <- c("#A40E11", "#F27272", "#FFC300", "#888888", "#000000")
names(color_vals) <- levels_gamma

# Thinner non-black lines, thicker black
line_sizes <- c(0.9, 0.9, 0.9, 0.9, 1.3)
names(line_sizes) <- levels_gamma

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
plot_panel_gamma <- function(yvar, ylabel) {
  ggplot(res_gamma_summary, aes(
    x = n_teams, y = !!sym(yvar),
    group = distribution, color = distribution, 
    linewidth = distribution
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
    labs(x = "n", y = ylabel, color = "Distribution-skewness") +
    theme_custom
}

p2 <- plot_panel_gamma("mean_corr", "Spearman correlation (%)")
p3 <- plot_panel_gamma("mean_top1", "Top-1 coverage (%)")
p4 <- plot_panel_gamma("mean_top3", "Top-3 coverage (%)")

# -------------------------------------------
# 6. First panel: distribution shapes
# -------------------------------------------
dist_info <- data.frame(
  label = levels_gamma,
  shape = c(1, 5, 10, 20, NA),
  rate = c(1, sqrt(5), sqrt(10), sqrt(20), NA)
)

curve_data_gamma <- data.frame()
for (i in 1:nrow(dist_info)) {
  dist <- dist_info$label[i]
  if (is.na(dist_info$shape[i])) {
    # Normal benchmark
    x <- seq(-5, 15, length.out = 400)
    density <- dnorm(x)
  } else {
    shape <- dist_info$shape[i]
    rate <- dist_info$rate[i]
    x <- seq(0, 15, length.out = 400)
    density <- dgamma(x, shape = shape, rate = rate)
  }
  curve_data_gamma <- rbind(curve_data_gamma,
                            data.frame(x = x, density = density, distribution = dist))
}

curve_data_gamma$distribution <- factor(curve_data_gamma$distribution, levels = levels_gamma)

p1 <- ggplot(curve_data_gamma, aes(
  x = x, y = density, color = distribution, linewidth = distribution
)) +
  geom_line() +
  scale_color_manual(values = color_vals) +
  scale_linewidth_manual(values = line_sizes, guide = "none") +
  labs(x = "x", y = "Density", title = "Distribution shapes") +
  theme_custom +
  theme(legend.position = "none")

# -------------------------------------------
# 7. Combine panels (1×4 layout, single shared legend)
# -------------------------------------------
# final_plot_gamma <- p1 + p2 + p3 + p4 +
#   plot_layout(ncol = 4, guides = "collect") &
#   theme(legend.position = "right")
# 
# final_plot_gamma
# 
# # Remove unwanted Gamma(0.5, 0.71) and Gamma(2, 1.41)
# all_res_gamma <- all_res_gamma %>%
#   filter(!distribution %in% c("G(0.5,0.71)- 2.83", "G(2,1.41)- 1.41"))

# Recode labels into polished display form
# all_res_gamma$distribution <- recode(all_res_gamma$distribution,
#                                      "G(1,1)- 2.00"   = "G(1,1)=Exp(1) – 2.00",
#                                      "G(5,2.24)- 0.89" = "G(5, √5) – 0.89",
#                                      "G(10,3.16)- 0.63" = "G(10, √10) – 0.63",
#                                      "G(20,4.47)- 0.45" = "G(20, √20) – 0.45",
#                                      "N(0,1)- 0.00"     = "N(0,1) – 0.00"
# )


# Ruby → Salmon → Amber → Grey → Black
color_vals <- c("#E0115F", "#F27272", "#FFC300", "#888888", "#000000")
names(color_vals) <- levels_gamma

plot_panel_gamma <- function(yvar, ylabel) {
  ggplot(res_gamma_summary, aes(
    x = n_teams, y = !!sym(yvar),
    group = distribution, color = distribution, linewidth = distribution
  )) +
    geom_line() +
    geom_point(size = 2.3) +
    scale_color_manual(values = color_vals) +
    scale_linewidth_manual(values = line_sizes) +
    scale_y_continuous(
      limits = c(0, 100),
      oob = scales::squish,
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%")
    ) +
    labs(x = "n", y = ylabel, color = "Distribution-skewness") +
    theme_custom
}

# Move legend inside the first panel (distribution shapes)
p1 <- ggplot(curve_data_gamma, aes(
  x = x, y = density, color = distribution, linewidth = distribution
)) +
  geom_line() +
  scale_color_manual(values = color_vals) +
  scale_linewidth_manual(values = line_sizes, guide = "none") +
  labs(
    x = expression(mu[i]),
    y = "Density", title = NULL, 
    color = "Distribution-skewness") +
  theme_custom +
  theme(
    legend.position = c(0.75, 0.75),   # adjust coordinates (x, y) inside panel
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

# Combine panels (1×4 layout)
final_plot_gamma <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 4,widths = c(1.5, 1, 1, 1))

final_plot_gamma

ggsave(
  "robust_gamma.png",
  plot = final_plot_gamma,
  width = 16,
  height = 3.5,
  dpi = 600
)

