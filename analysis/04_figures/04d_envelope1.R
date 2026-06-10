library(dplyr)
library(ggplot2)
library(forcats)
library(tidyverse)

df_expanded<-readRDS("RData/03_processed_top3_1000_seed42_42.rds")

df <-df_expanded %>%
   filter(!str_detect(method, "selfrm|BayePeer"))

# df <-df_expanded %>%
#    filter(str_detect(method, "selfrm|BayeSPeer") )
#   
# -------------------------------------------------
# 1. methods to highlight
# -------------------------------------------------
best_methods <- c("BayeSPeer", "BIRRA", "Stuart", "MC3", "CEMC.k", "Mean")

# best_methods <- c("BayeSPeer", "Stuart.selfrm", "MC3.selfrm", 
#                    "CEMC.k.selfrm", "Mean.selfrm", "BIRRA.selfrm")
# -------------------------------------------------
# 2. enforce correct order of x-values within each section
# -------------------------------------------------
df2 <- df %>%
  mutate(
    Section = factor(Section, levels = c("n", "rho", "muBeta", "sigma2Beta")),
    
    Group = case_when(
      Section == "n"   ~ factor(as.character(Group),
                                levels = c("5", "10", "15", "20", "25")),
      Section == "muBeta"  ~ factor(as.character(Group),
                                levels = c("-1.5", "-0.5", "0.0", "0.5", "1.0", "1.5", "5.0")),
      Section == "rho"   ~ factor(as.character(Group),
                                levels = c("0.15", "0.35", "0.55", "0.75", "0.95","U(.15,.95)")),
      Section == "sigma2Beta"  ~ factor(as.character(Group),
                                levels = c("0.25", "0.50", "1.00", "2.00", "4.00")),
      TRUE ~ as.factor(Group)
    ),
    method = factor(method, levels = c(
      "BayeSPeer",  # 👑 first
      "Stuart",
      "MC3",
      "CEMC.k",
      "Mean",
      "BIRRA"
    ))
    # method = factor(method, levels = c(
    #   "BayeSPeer",  # 👑 first
    #   "Stuart.selfrm",
    #   "MC3.selfrm",
    #   "CEMC.k.selfrm",
    #   "Mean.selfrm",
    #   "BIRRA.selfrm"
    # ))
  )

envelope_df <- df2 %>%
  filter(!method %in% best_methods) %>%
  group_by(Section, Group) %>%
  summarise(
    ymin = min(mean_value, na.rm = TRUE),
    ymax = max(mean_value, na.rm = TRUE),
    .groups = "drop"
  )


highlight_df <- df2 %>%
  filter(method %in% best_methods) %>%
  mutate(method = factor(method, levels = best_methods))

method_colors <- c(
  "BayeSPeer" = "black",
  "Stuart"    = "#4DAF4A",
  "Mean"       = "#7B4CC2",
  "MC3"    = "#F28E2B",
  "BIRRA"      = "#3778BF",
  "CEMC.k"     = "red"
)

method_linewidths <- c(
  "BayeSPeer" = 1.1,
  "Stuart"    = 0.9,
  "MC3"       = 0.9,
  "CEMC.k"    = 0.9,
  "Mean"      = 0.9,
  "BIRRA"     = 0.9
)

method_linewidths <- c(
  "BayeSPeer" = 1.1,
  "Stuart"    = 0.9,
  "MC3"       = 0.9,
  "CEMC.k"    = 0.9,
  "Mean"      = 0.9,
  "BIRRA"     = 0.9
)


method_shapes <- c(
  "BayeSPeer" = 1,  # circle
  "Stuart"    = 0,  # square
  "MC3"       = 2,  # triangle up
  "CEMC.k"    = 3,  # triangle down
  "Mean"      = 4,  # diamond
  "BIRRA"     = 5   # reuse if needed
)


# method_colors <- c(
#   "BayeSPeer" = "black",
#   "Stuart.selfrm"    = "#4DAF4A",
#   "Mean.selfrm"       = "#7B4CC2",
#   "MC3.selfrm"    = "#F28E2B",
#   "BIRRA.selfrm"      = "#3778BF",
#   "CEMC.k.selfrm"     = "red"
# )
# 
# 
# method_linewidths <- c(
#   "BayeSPeer" = 1.1,
#   "Stuart.selfrm"    = 0.9,
#   "MC3.selfrm"       = 0.9,
#   "CEMC.k.selfrm"    = 0.9,
#   "Mean.selfrm"      = 0.9,
#   "BIRRA.selfrm"     = 0.9
# )
# 
# 
# method_shapes <- c(
#   "BayeSPeer" = 1,  # circle
#   "Stuart.selfrm"    = 0,  # square
#   "MC3.selfrm"       = 2,  # triangle up
#   "CEMC.k.selfrm"    = 3,  # triangle down
#   "Mean.selfrm"      = 4,  # diamond
#   "BIRRA.selfrm"     = 5   # reuse if needed
# )

library(ggtext)
facet_labeller <- as_labeller(c(
  "n"          = "**n**",
  "rho"        = "**ρ**",
  "muBeta"     = "**μ**<sub>**β**</sub>",
  "sigma2Beta" = "**σ**<sub>**β**</sub><sup>**2**</sup>"
))

p <- ggplot() +
  # 🌫️ Envelope first (very back)
  geom_ribbon(
    data = envelope_df,
    aes(x = Group, ymin = ymin, ymax = ymax, group = 1),
    fill = "grey85",
    alpha = 0.8
  ) +
  
  # 🎨 Other methods (middle layer)
  geom_line(
    data = highlight_df %>% filter(method != "BayeSPeer"),
    aes(x = Group, y = mean_value, 
        color = method, 
        group = method),
    linewidth = 0.9, 
  ) +
  
  geom_point(
    data = highlight_df %>% filter(method != "BayeSPeer"),
    aes(x = Group, y = mean_value, 
        color = method, 
        shape = method,
        group = method),
    size = 2.5,
    stroke = 1,        # border thickness
    fill = NA,
    position = position_jitter(width = 0.15, height = 0)
  )+
  
  # 🖤 BayeSPeer LAST (top layer)
  geom_line(
    data = highlight_df %>% filter(method == "BayeSPeer"),
    aes(x = Group, y = mean_value, 
        color = method, 
        group = method),
    linewidth = 1.1
  ) +
  
    geom_point(
      data = highlight_df %>% filter(method == "BayeSPeer"),
      aes(x = Group, y = mean_value,
          color = method, 
          shape = method,
          group = method),
      size = 2.5,
      stroke = 1,
      fill = NA
    )+
  
  facet_wrap(~ Section, 
             scales = "free_x", 
             nrow = 1,
             strip.position = "top",
             labeller = facet_labeller) +
  
  scale_color_manual(values = method_colors,
                     breaks = c("BayeSPeer", "MC3","BIRRA",
                                "CEMC.k", "Stuart",  "Mean")) +
                     # breaks = c("BayeSPeer", "BIRRA.selfrm", "Stuart.selfrm",
                     #           "MC3.selfrm", "CEMC.k.selfrm",  "Mean.selfrm"))+
    scale_shape_manual(values = method_shapes,
                       breaks = c("BayeSPeer", "MC3","BIRRA",
                                  "CEMC.k", "Stuart",  "Mean"))+
                     # breaks = c("BayeSPeer", "BIRRA.selfrm", "Stuart.selfrm",
                     #            "MC3.selfrm", "CEMC.k.selfrm",  "Mean.selfrm"))+
  labs(
    x = NULL,
    y = "Top 3 Coverage rate (%)",#"Spearman Correlation",
    color = NULL,
    linetype = NULL
  ) +
  
  theme_minimal(base_size = 16) +
  ylim(c(0,100))+
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 18, face = "bold"),
    strip.text.x = element_markdown(size = 18, margin = margin(b = 0)),
    strip.background = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = c(0.63, 0.77),
    legend.text  = element_text(size = 12),  # method names
    legend.title = element_text(size = 12),   # "Method"
    legend.justification = c(0, 0),
    legend.key.height = unit(0.9, "lines"),
    legend.spacing.y = unit(0.5, "lines"),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    panel.grid.minor = element_blank()
  )+
  labs(color = "Method", shape = "Method")+
  scale_x_discrete(position = "top") +
  
  theme(
  #  axis.title.x = element_blank(),
  # axis.text.x.top = element_text(size = 12),
    panel.spacing.x = unit(1.5, "lines"),   # widens spacing between panels
    panel.border = element_rect(color = "gray0", fill = NA, linewidth = 0.5), # optional border around each panel
    
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12,angle = 60, vjust=-0.01, hjust = 0),
    axis.text.y = element_text(size = 12) 
  )+guides(
    color = guide_legend(nrow = 2),
    shape = guide_legend(nrow = 2)
  )
p
ggsave(
  "envelope_top3.png",
  plot = p,
  width = 13.5,
  height = 5.5,
  dpi = 600
)


