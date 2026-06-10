load("RData/02_compareBiGER_1000_seed42_42.RData")
res_biger <- results

colnames(res_biger) <- gsub("BiGER", "BayePeer", colnames(res_biger))

load("RData/02_compare_1000_seed42_42.RData")
colnames(results_multi) <- gsub("BayeSRank", "BayeSPeer", colnames(results_multi))

res_bayespeer <- results_multi %>% dplyr::select(1:6, 21, 36)

res_nonmix <- merge(res_biger, res_bayespeer)

load("RData/02_compare_mix.RData")
res_mix <- results%>% dplyr::select(1:6, 11,20,25, 34,39)
colnames(res_mix) <- gsub("BiGER", "BayePeer", colnames(res_mix))
colnames(res_mix) <- gsub("BayeSRank", "BayeSPeer", colnames(res_mix))

res <- rbind.data.frame(res_nonmix, res_mix)%>% 
  mutate(sigma2Beta=sigmaBeta^2,
         rho=ifelse(rho=="U(0.15, 0.95)", "U(.15,.95)", rho))

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
  id.vars = c("n_teams", "rho", "muBeta", "sigma2Beta","sigmaBeta", "sim_id"),
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
  group_by(n_teams, rho, muBeta, sigma2Beta,sigmaBeta, Metric, Method) %>%
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
    mutate(Section = "rho", Group = as.character(rho)),
  
  plot_sum %>%
    filter(n_teams == base_n, rho == base_rho, abs(sigmaBeta - base_sigma) < 1e-8) %>%
    mutate(Section = "muBeta", Group = as.character(muBeta)),
  
  plot_sum %>%
    filter(n_teams == base_n, rho == base_rho, muBeta == base_mu) %>%
    mutate(
      Section = "sigma2Beta",
      sigma2Beta = sigmaBeta^2,
      Group = sprintf("%.2f", sigma2Beta)
    )
) %>%
  mutate(
    Section = factor(
      Section,
      levels = c("n", "rho", "muBeta", "sigma2Beta")
    ),
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
      Section == "n" ~ factor(
        Group,
        levels = c("5", "10", "15", "20", "25", "50")
      ),
      Section == "rho" ~ factor(
        Group,
        levels = c("0.15", "0.35", "0.55", "0.75", "0.95","U(.15,.95)")
      ),
      Section == "muBeta" ~ factor(
        Group,
        levels = c("-1.5", "-1.0", "-0.5", "0.0", "0.5", "1.0", "1.5", "5.0")
      ),
      Section == "sigma2Beta" ~ factor(
        Group,
        levels = c("0.25", "0.50", "1.00", "2.00", "4.00")
      ),
      TRUE ~ factor(Group)
    )
 )

# section_labeller <- as_labeller(
#   c(
#     n = "n",
#     rho = "rho",
#     mu_beta = "mu[beta]",
#     sigma_beta2 = "sigma[beta]^2"
#   ),
#   label_parsed
# )


facet_labeller <- as_labeller(c(
  "n"          = "**n**",
  "rho"        = "**ρ**",
  "muBeta"     = "**μ**<sub>**β**</sub>",
  "sigma2Beta" = "**σ**<sub>**β**</sub><sup>**2**</sup>"
))

make_panel <- function(metric_name, ylab, colors, ymin) {
  
  panel_data <- plot_df %>%
    filter(Metric == metric_name) %>%
    mutate(
      x_num = as.numeric(Group), #how to fix this because we also use scale_xcontinuous
      dodge_offset = ifelse(Method == "BayePeer", -0.2, 0.2),
      xmin = x_num + dodge_offset - 0.2,
      xmax = x_num + dodge_offset + 0.15
    )
  
  top_labels <- panel_data %>%
    distinct(Section, Group, x_num) %>%
    mutate(y = 112) # <- moved higher
  
  ggplot(panel_data) +
    
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = mean_value,
        fill = Method
      ),
      color = "black",
      linewidth = 0.25
    ) +
    
    geom_errorbar(
      aes(
        x = x_num + dodge_offset,
        ymin = mean_value - 1.96 * se,
        ymax = mean_value + 1.96 * se,
        group = Method
      ),
      width = 0.16,
      linewidth = 0.2
    ) +
    
    geom_text(
      aes(
        x = x_num + dodge_offset +
          ifelse(Method == "BayePeer", -0.04, 0.04),
        y = mean_value +
          ifelse(Method == "BayePeer", 2.5, 5.5), # stagger labels
        label = round(mean_value, 1)
      ),
      size = 3.5#,
     # fontface = "bold"
    ) +
    
    geom_text(
      data = top_labels,
      aes(
        x = x_num,
        y = y,
        label = Group
      ),
      inherit.aes = FALSE,
      size = 4.5,
      angle = 45,
      hjust = 0,    # left-align since text fans upward-right
      vjust = 0
    )+
    
    facet_grid(
      . ~ Section,
      scales = "free_x",
      space = "free_x",
      labeller = facet_labeller
    ) +
    
    scale_fill_manual(
      values = colors,
      name = "Method"
    ) +
    
    scale_x_continuous(
      breaks = sort(unique(panel_data$x_num)),
      labels = NULL,
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    
    coord_cartesian(
      ylim = c(ymin, 108),  # was 112 — grid now ends lower
      clip = "off"
    ) +
    
    labs(x = NULL, y = ylab) +
    
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_markdown(
        size = 18,
        margin = margin(b = 33)
      ),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 14), 
      panel.spacing.x = unit(0.5, "lines"),
      
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0) 
    )
}
(p1 <- make_panel(
  "Spearman\nCorrelation (%)",
  "Spearman\nCorrelation (%)",
  c("BayePeer" = "#9EADB4", "BayeSPeer" = "#2F4D5F"),
  ymin = 40
))

p2 <- make_panel(
  "Top-1\nCoverage Rate (%)",
  "Top-1\nCoverage Rate (%)",
  c("BayePeer" = "#E8A092", "BayeSPeer" = "#D95745"),
  ymin = 20
)
  
p3 <- make_panel(
  "Top-3\nCoverage Rate (%)",
  "Top-3\nCoverage Rate (%)",
  c("BayePeer" = "#F4C79A", "BayeSPeer" = "#e68d45"),
  ymin = 40
)

final_plot <- p1 / p2 / p3 +
  plot_layout(guides = "keep")

final_plot

ggsave(
  "barplot_3metrics.png",
  plot = final_plot,
  width = 16.9,
  height = 8.6,
  dpi = 600
)
