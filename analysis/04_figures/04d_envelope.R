# 1.	Identify the Range: At each data point (e.g., each n or rho or mu), calculate the minimum and maximum values across all baseline methods (BIRRA, Mean, MC3, etc.).
# 
# 2.	Create the Shaded Area: Fill the area between that minimum and maximum with a light, neutral color (like light gray or a pale blue). This represents the "envelope" of standard performance.
# 3.	Highlight BayeSRank and other top methods in each category : Plot the BayeSRank and other best-in-category methods line as a thick, bold, black and other brightly colored line on top of that shaded area.
library(ggplot2)
library(dplyr)
library(stringr)

my_method <- "BayeSRank"

best_methods <- c("BayeSRank","BIRRA","Stuart","MC3","CEMC.k","Mean")

# Unique methods
#df_expanded<-readRDS("RData/03_processed_top1_seed2222.rds")
df_expanded<-readRDS("RData/03_processed_top3_seed1111.rds")

#df <-df_expanded %>% filter(str_detect(method, "selfrm|BayeSRank") )

df <-df_expanded %>% filter(!str_detect(method, "selfrm|BayeRank")) %>%
  filter(#method %in% best_methods,
    n_teams !=50)

method_colors <- c(
  "BayeSRank" = "gray30",
  "BIRRA"     = "red",
  "RRA"       = "#F781BF",
  "Stuart"    = "#4DAF4A",
  "BayeRank"  = "tomato3",
  "min"       = "#984EA3",
  "CEMC.k"    = "orange",
  "CEMC.s"    = "deepskyblue4",
  "MC1"       = "firebrick",
  "MC2"       = "olivedrab4",
  "MC3"       = "darkviolet",
  "Mean"      = "#377EB8",
  "Median"    = "gold",
  "Geomean"   = "mediumpurple4",
  "L2norm"    = "salmon"
)

plot_base <- df %>%
  ungroup() %>%
  mutate(
    Section = factor(Section, levels=c("n","ρ","μᵦ","σᵦ")),
    Group   = as.character(Group),
    method  = as.character(method)
  )


n_order     <- c("5", "10", "15", "20", "25")
rho_order   <- c("0.15", "0.35", "0.55", "0.75", "0.95")
mu_order    <- c("-1.5","-0.5","0.0","0.5","1.0","1.5","5.0")
sigma_order <- c("√0.25", "√0.5", "√1", "√2", "√4")

x_map <- bind_rows(
  tibble(Section = "n",   Group = n_order,     x = seq_along(n_order)),
  tibble(Section = "ρ",   Group = rho_order,   x = seq_along(rho_order)),
  tibble(Section = "μᵦ",  Group = mu_order,    x = seq_along(mu_order)),
  tibble(Section = "σᵦ",  Group = sigma_order, x = seq_along(sigma_order))
)

plot_base <- plot_base %>%
  left_join(x_map, by = c("Section", "Group"))
plot_base %>% filter(is.na(x)) %>% distinct(Section, Group)


highlight_df <- plot_base %>%
  filter(method %in% best_methods)

baseline_df <- plot_base %>%
  filter(!method %in% best_methods)

envelope_df <- baseline_df %>%
  group_by(Section, x) %>%
  summarise(
    ymin = min(mean_value, na.rm = TRUE),
    ymax = max(mean_value, na.rm = TRUE),
    .groups = "drop"
  )


section_labs <- c(
  "n"  = "n",
  "ρ"  = "ρ",
  "μᵦ" = "μᵦ",
  "σᵦ" = "σᵦ"
)

axis_labels <- x_map %>%
  arrange(Section, x)

ggplot() +
  # performance envelope for non-highlighted methods
  geom_ribbon(
    data = envelope_df,
    aes(x = x, ymin = ymin, ymax = ymax, group = 1),
    fill = "grey85",
    alpha = 0.55
  ) +
  
  # highlighted benchmark methods
  geom_line(
    data = highlight_df %>% filter(method != "BayeSRank"),
    aes(x = x, y = mean_value, color = method, group = method),
    linewidth = 1.0,
    alpha = 0.85
  ) +
  geom_point(
    data = highlight_df %>% filter(method != "BayeSRank"),
    aes(x = x, y = mean_value, color = method),
    size = 2,
    alpha = 0.9
  ) +
  
  # BayeSRank on top
  geom_line(
    data = highlight_df %>% filter(method == "BayeSRank"),
    aes(x = x, y = mean_value, group = method),
    color = method_colors["BayeSRank"],
    linewidth = 1.8
  ) +
  geom_point(
    data = highlight_df %>% filter(method == "BayeSRank"),
    aes(x = x, y = mean_value), shape=2,
    color = method_colors["BayeSRank"],
    size = 2
  ) +
  
  facet_wrap(~ Section, scales = "free_x", nrow = 1,
             labeller = as_labeller(section_labs)) +
  
  scale_color_manual(
    values = method_colors[best_methods[best_methods != "BayeSRank"]],
    breaks = best_methods[best_methods != "BayeSRank"]
  ) +
  
  scale_x_continuous(
    breaks = 1:8,
    labels = function(z) z
  ) +
  
  labs(
    x = NULL,
    y = "Performance",
    color = NULL
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey82", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 18),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )

plot_base2 <- plot_base %>%
  mutate(
    Group_plot = case_when(
      Section == "n"  ~ factor(Group, levels = n_order),
      Section == "ρ"  ~ factor(Group, levels = rho_order),
      Section == "μᵦ" ~ factor(Group, levels = mu_order),
      Section == "σᵦ" ~ factor(Group, levels = sigma_order)
    )#,
   # Group_plot = factor(Group_plot, levels=c("n", "ρ","μᵦ","σᵦ"))
  )

highlight_df <- plot_base2 %>%
  filter(method %in% best_methods)

baseline_df <- plot_base2 %>%
  filter(!method %in% best_methods)

envelope_df <- baseline_df %>%
  group_by(Section, Group_plot) %>%
  summarise(
    ymin = min(mean_value, na.rm = TRUE),
    ymax = max(mean_value, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  geom_ribbon(
    data = envelope_df,
    aes(x = Group_plot, ymin = ymin, ymax = ymax, group = 1),
    fill = "grey85",
    alpha = 0.55
  ) +
  
  geom_line(
    data = highlight_df %>% filter(method != "BayeSRank"),
    aes(x = Group_plot, y = mean_value, color = method, group = method),
    linewidth = 1.0,
    alpha = 0.85
  ) +
  geom_point(
    data = highlight_df %>% filter(method != "BayeSRank"),
    aes(x = Group_plot, y = mean_value, color = method), shape=2,
    size = 2,
    alpha = 0.9
  ) +
  
  geom_line(
    data = highlight_df %>% filter(method == "BayeSRank"),
    aes(x = Group_plot, y = mean_value, group = method),
    color = method_colors["BayeSRank"],
    linewidth = 1.5
  ) +
  geom_point(
    data = highlight_df %>% filter(method == "BayeSRank"),
    aes(x = Group_plot, y = mean_value),
    color = method_colors["BayeSRank"],
    size = 2
  ) +
  #scale_x_discrete(position = "top") +
  facet_wrap(~ Section, scales = "free_x", nrow = 1,
             labeller = as_labeller(section_labs)) +
  
  scale_color_manual(
    values = method_colors[best_methods[best_methods != "BayeSRank"]],
    breaks = best_methods[best_methods != "BayeSRank"]
  ) +
  
  labs(
    x = NULL,
    y = "Spearman Correlation",
    color = NULL
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
   # panel.grid.minor = element_blank(),
   # panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey84", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 18),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )

