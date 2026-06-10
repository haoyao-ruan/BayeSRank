
load("RData/02_compareBiGER_1000_seed42_42.RData")
res_biger <- results

colnames(res_biger)<-gsub("BiGER",
  "BayePeer", colnames(res_biger))


load("RData/02_compare_1000_seed42_42.RData")
colnames(results_multi)<-gsub("BayeSRank","BayeSPeer", colnames(results_multi))

res_bayespeer <- results_multi%>% select(1:6, 21, 36)
res <- merge(res_biger, res_bayespeer)

res$muBeta <- sprintf("%.1f", as.numeric(res$muBeta))


library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

setDT(res)

base_rho   <- 0.55
base_mu    <- "1.0"
base_sigma <- sqrt(0.5)
base_n     <- 15

metric_cols <- c(
  "BayeSPeer", "BayePeer",
  "top1_BayeSPeer", "top1_BayePeer",
  "top3_BayeSPeer", "top3_BayePeer"
)

long_results <- melt(
  res,
  id.vars = c("n_teams", "rho", "muBeta", "sigmaBeta", "sim_id"),
  measure.vars = metric_cols,
  variable.name = "raw_method",
  value.name = "value"
)

plot_sum <- long_results %>%
  mutate(
    Metric = case_when(
      grepl("^top1_", raw_method) ~ "Top-1\nCoverage Rate (%)",
      grepl("^top3_", raw_method) ~ "Top-3\nCoverage Rate (%)",
      TRUE ~ "Spearman\nCorrelation (%)"
    ),
    Method = case_when(
      grepl("BayeSPeer", raw_method) ~ "BayeSPeer",
      grepl("BayePeer", raw_method) ~ "BayePeer"
    )
  ) %>%
  group_by(n_teams, rho, muBeta, sigmaBeta, Metric, Method) %>%
  summarise(
    mean_value = 100 * mean(value, na.rm = TRUE),
    se = 100 * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

plot_df <- bind_rows(
  plot_sum %>%
    filter(rho == base_rho, muBeta == base_mu, abs(sigmaBeta - base_sigma) < 1e-8) %>%
    mutate(Section = "n", Group = as.character(n_teams)),
  
  plot_sum %>%
    filter(n_teams == base_n, muBeta == base_mu, abs(sigmaBeta - base_sigma) < 1e-8) %>%
    mutate(Section = "ρ", Group = as.character(rho)),
  
  plot_sum %>%
    filter(n_teams == base_n, rho == base_rho, abs(sigmaBeta - base_sigma) < 1e-8) %>%
    mutate(Section = "μᵦ", Group = as.character(muBeta)),
  
  plot_sum %>%
    filter(n_teams == base_n, rho == base_rho, muBeta == base_mu) %>%
    mutate(
      Section = "σᵦ",
      Group = case_when(
        abs(sigmaBeta - sqrt(0.25)) < 1e-8 ~ "√0.25",
        abs(sigmaBeta - sqrt(0.5))  < 1e-8 ~ "√0.5",
        abs(sigmaBeta - 1)          < 1e-8 ~ "√1",
        abs(sigmaBeta - sqrt(2))    < 1e-8 ~ "√2",
        abs(sigmaBeta - 2)          < 1e-8 ~ "√4",
        TRUE ~ as.character(round(sigmaBeta, 3))
      )
    )
) %>%
  mutate(
    Section = factor(Section, levels = c("n", "ρ", "μᵦ", "σᵦ")),
    Method = factor(Method, levels = c("BayePeer", "BayeSPeer")),
    Metric = factor(
      Metric,
      levels = c(
        "Spearman\nCorrelation (%)",
        "Top-1\nCoverage Rate (%)",
        "Top-3\nCoverage Rate (%)"
      )
    ),
    Group = case_when(
      Section == "n"  ~ factor(Group, levels = c("5", "10", "15", "20", "25", "50")),
      Section == "ρ"  ~ factor(Group, levels = c("0.15", "0.35", "0.55", "0.75", "0.95")),
      Section == "μᵦ" ~ factor(Group, levels = c("-1.5","-1.0", "-0.5", "0.0", "0.5", "1.0", "1.5", "5.0")),
      Section == "σᵦ" ~ factor(Group, levels = c("√0.25", "√0.5", "√1", "√2", "√4")),
      TRUE ~ factor(Group)
    ),
    FillKey = interaction(Method, Metric, sep = "_")
  )

make_panel <- function(metric_name, ylab, colors) {
  
  ggplot(plot_df %>% filter(Metric == metric_name),
         aes(x = Group, y = mean_value, fill = Method)) +
    
    geom_col(
      position = position_dodge(0.75),
      width = 0.65,
      color = "black",
      linewidth = 0.35
    ) +
    
    geom_errorbar(
      aes(ymin = mean_value - 1.96*se, ymax = mean_value + 1.96*se),
      position = position_dodge(0.75),
      width = 0.18
    ) +
    
    geom_text(
      aes(label = round(mean_value,1)),
      position = position_dodge(0.95),
      vjust = -0.5,
      size = 2.8
    ) +
    
    facet_grid(
      . ~ Section,
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_fill_manual(
      values = colors,
      name = "Method"
    ) +
    
    coord_cartesian(ylim = c(0, 105), clip = "off") +
    
    labs(x = NULL, y = ylab) +
    
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

# Row 1 (blue-gray)
p1 <- make_panel(
  "Spearman\nCorrelation (%)",
  "Spearman\nCorrelation (%)",
  c("BayePeer" = "#9EADB4", "BayeSPeer" = "#2F4D5F")
)

# Row 2 (red)
p2 <- make_panel(
  "Top-1\nCoverage Rate (%)",
  "Top-1\nCoverage Rate (%)",
  c("BayePeer" = "#E8A092", "BayeSPeer" = "#D95745")
)

# Row 3 (orange)
p3 <- make_panel(
  "Top-3\nCoverage Rate (%)",
  "Top-3\nCoverage Rate (%)",
  c("BayePeer" = "#F4C79A", "BayeSPeer" = "#e68d45")
)

library(patchwork)

final_plot <- p1 / p2 / p3 +
  plot_layout(guides = "keep")  # ← CRUCIAL

final_plot

ggsave(
  "barplot_3metrics.png",
  plot = final_plot,
  width = 15,
  height = 6,
  dpi = 600
)
