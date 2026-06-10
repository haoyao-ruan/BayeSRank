library(ggplot2)
library(dplyr)
library(stringr)

my_method <- "BayeSRank"

best_methods <- c("BayeSRank","BIRRA","Stuart","MC3","CEMC.k","Mean")

# Unique methods
#df_expanded<-readRDS("RData/03_processed_top1_seed2222.rds")
df_expanded<-readRDS("RData/03_processed_corr_1000_seed42_42.rds")

#df <-df_expanded %>% filter(str_detect(method, "selfrm|BayeSRank") )

df <-df_expanded %>% filter(!str_detect(method, "selfrm|BayeRank")) 
  #filter(method %in% best_methods)

methods <- unique(df$method)

my_method <- "BayeSRank"
other_methods <- setdiff(methods, my_method)


width = 0   # disables horizontal lines at top/bottom
theme(panel.border = element_rect(color = "grey40", fill = NA))

order_methods <- rev(unique(df$method[!str_detect(df$method, "selfrm")]))
df$method = factor(df$method,
                               levels=order_methods)

df_main  <- df %>% filter(method != "BayeSRank")
df_bayes <- df %>% filter(method == "BayeSRank")

# Colors
method_colors <- c(
  "BayeSRank" = "gray30",
  "BIRRA"     = "red",#"#999999",
  
  "RRA"       = "#F781BF",  # soft magenta
  "Stuart"    = "#4DAF4A",
  "BayeRank"  = "tomato3",
  "min"       = "#984EA3",  # forest green

  "CEMC.k"    = "orange",  # gold/yellow
  "CEMC.s"    = "deepskyblue4",  # bright green
  
  "MC1"       = "firebrick",
  "MC2"       = "olivedrab4",  # turquoise
  "MC3"       = "darkviolet",  # hot pink
  
  "Mean"      = "#377EB8",
  "Median"    = "gold",  # deep violet
  "Geomean"   = "mediumpurple4",  # vivid blue
  "L2norm"    = "salmon" # purple
  
)
# method_colors <- c(
#   "BIRRA.selfrm"     = "#999999",
# 
#   "CEMC.k.selfrm"    = "orange",  # gold/yellow
#   "CEMC.s.selfrm"    = "#4DAF4A",  # bright green
#   "Geomean.selfrm"   = "#377EB8",  # vivid blue
#   "L2norm.selfrm"    = "salmon",  # purple
# 
#   "MC2.selfrm"       = "firebrick",
#   "MC3.selfrm"       = "darkviolet",  # hot pink
#   "Mean.selfrm"      = "gold",  # deep violet
#   "Median.selfrm"    = "mediumpurple4",
#   "RRA.selfrm"       = "#F781BF",  # soft magenta
#   "Stuart.selfrm"    = "#984EA3",  # forest green
#   "BayeRank.selfrm"  = "tomato3",
#   "min.selfrm"       = "deepskyblue4",
#   "MC1.selfrm"       = "olivedrab4",  # turquoise
#   "BayeSRank" = "gray30"
# )

df$method <- factor(df$method,
                             levels = rev(levels(df$method)))

# df$method <- factor(df$method, 
#                     levels = c("BayeSRank", "BIRRA", "Stuart", "MC3", "CEMC.k","Mean"))

ggplot(df, aes(x = Group, y = mean_value, group = method)) +
  # geom_errorbar(
  #   aes(
  #     ymin = mean_value - 1.96 * se,
  #     ymax = mean_value + 1.96 * se,
  #     color = method,
  #     alpha = ifelse(method == my_method, 0.8, 0.6)
  #   ),
  #   width = 0,
  #   position = position_jitter(width = 0),
  #   show.legend = FALSE
  # ) +
  
  scale_x_discrete(position = "top") +
  geom_line(
    aes(
      color = method,
      size = ifelse(method == my_method, 1.8, 1.2),
      alpha = ifelse(method == my_method, 0.8, 0.7)
    ),
    position = position_jitter(width = 0.02, height = 0)
  ) +
  
  geom_point(
    aes(
      color = method,
      shape = method,
      alpha = ifelse(method == my_method, 0.8, 0.6)
    ),
    size = 4,
    position = position_jitter(width = 0.01, height = 0)
  ) +
  geom_line(
    data = df_bayes,
    aes(x = Group, y = mean_value, group = method),
    color = "black",
    alpha=0.8,
    linewidth = 1.5,
    show.legend = FALSE,
    position = position_jitter(width = 0.01, height = 0)
  ) +
  geom_point(
    data = df_bayes,
    aes(x = Group, y = mean_value, group = method),
    color = "black",
    alpha=0.8,
    size = 4,
    shape = 5,
    show.legend = FALSE,
    position = position_jitter(width = 0.01, height = 0)
  ) +
 
  
  # Facets and styling
  facet_grid(~Section, scales = "free_x", space = "free_x") +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = 0:13) +
  scale_size_identity() +
  scale_alpha_identity() +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "gray40", fill = NA),  # panel borders
    strip.text.x = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.placement = "outside",  # Keeps facet labels above plot
    axis.ticks.length = unit(5, "pt"),  # Increase tick spacing
    panel.spacing = unit(1.2, "lines"),  # Add space between facets if needed
  ) +
  labs(
    title = NULL,
    x = NULL,
    y ="Spearman Correlation (%)",
    #y = "Top 1 Coverage rate (%)",
    color = "Method",
    shape = "Method"
   # subtitle = "Robustness checking: true performance ~ folded Normal distribution"
  ) +ylim(c(30,100))
# 1645 x 432


#-------------- performance envelope ---------------
library(dplyr)
library(ggplot2)

# ungroup first so old grouping doesn't interfere
df2 <- df %>% ungroup()

# choose baselines only
baseline_df <- df2 %>%
  filter(method != "BayeSRank")

# BayeSRank only
bayesrank_df <- df2 %>%
  filter(method == "BayeSRank")

# performance envelope: min/max across baselines at each x-position
envelope_df <- baseline_df %>%
  group_by(Section, Group) %>%
  summarise(
    ymin = min(mean_value, na.rm = TRUE),
    ymax = max(mean_value, na.rm = TRUE),
    .groups = "drop"
  )

x_map <- df2 %>%
  distinct(Section, Group) %>%
  group_by(Section) %>%
  arrange(Section, Group, .by_group = TRUE) %>%
  mutate(x = row_number()) %>%
  ungroup()
envelope_df  <- envelope_df  %>% left_join(x_map, by = c("Section", "Group"))
bayesrank_df <- bayesrank_df %>% left_join(x_map, by = c("Section", "Group"))
baseline_df  <- baseline_df  %>% left_join(x_map, by = c("Section", "Group"))

ggplot() +
  geom_ribbon(
    data = envelope_df,
    aes(x = x, ymin = ymin, ymax = ymax, group = 1),
    fill = "grey90",
    alpha = 0.5
  ) +
  geom_line(
    data = bayesrank_df,
    aes(x = x, y = mean_value, group = 1),
    linewidth = 1.2,
    color = "black"
  ) +
  geom_point(
    data = bayesrank_df,
    aes(x = x, y = mean_value),
    size = 2.2,
    color = "black"
  ) +
  facet_wrap(~ Section, scales = "free_x", nrow = 1) +
  scale_x_continuous(
    breaks = x_map$x,
    labels = x_map$Group
  ) +
  labs(
    x = NULL,
    y = "Performance",
    title = "BayeSRank versus baseline performance envelope"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )



df3 <- df2 %>%
  mutate(category = case_when(
    method %in% c("Mean", "Median", "L2norm", "GeoMean") ~ "Borda",
    method %in% c("MC1", "MC2", "MC3") ~ "Markov",
    method %in% c("BIRRA", "Stuart", "RRA", "CEMC.k", "CEMC.s") ~ "Optimization",
    method %in% c("BayeSRank") ~ "Bayesian",
    TRUE ~ "Other"
  ))
top_methods <- df3 %>%
  filter(method != "BayeSRank") %>%
  group_by(category, method) %>%
  summarise(avg_perf = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(order_by = avg_perf, n = 1, with_ties = FALSE)

top_df <- df3 %>%
  semi_join(top_methods, by = c("category", "method")) %>%
  left_join(x_map, by = c("Section", "Group"))

bayesrank_df2 <- df3 %>%
  filter(method == "BayeSRank") %>%
  left_join(x_map, by = c("Section", "Group"))

ggplot() +
  geom_ribbon(
    data = envelope_df,
    aes(x = x, ymin = ymin, ymax = ymax, group = 1),
    fill = "grey85",
    alpha = 0.7
  ) +
  geom_line(
    data = top_df,
    aes(x = x, y = mean_value, color = method, group = method),
    linewidth = 0.9
  ) +
  geom_point(
    data = top_df,
    aes(x = x, y = mean_value, color = method),
    size = 1.8
  ) +
  geom_line(
    data = bayesrank_df2,
    aes(x = x, y = mean_value, group = 1),
    linewidth = 1.4,
    color = "black"
  ) +
  geom_point(
    data = bayesrank_df2,
    aes(x = x, y = mean_value),
    size = 2.4,
    color = "black"
  ) +
  facet_wrap(~ Section, scales = "free_x", nrow = 1) +
  scale_x_continuous(
    breaks = x_map$x,
    labels = x_map$Group
  ) +
  labs(
    x = NULL,
    y = "Performance",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside", 
    legend.position = "bottom"
  )
