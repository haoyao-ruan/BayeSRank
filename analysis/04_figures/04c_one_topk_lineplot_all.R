library(ggplot2)
library(dplyr)
library(stringr)

my_method <- "BayeSPeer"

# Unique methods
df_expanded<-readRDS("RData/03_processed_top1_1000_seed42_42.rds")

df <-df_expanded %>% filter(!str_detect(method, "selfrm|BayePeer") )

#choose if selfrm
df <-df_expanded  %>% filter(str_detect(method, "selfrm|BayeSPeer") )

#methods <- unique(df$method)

methods_to_compare= c("BayeSPeer",  #"BayePeer",
                      "BIRRA", "BIRRA.selfrm",
                      
                      "RRA", "RRA.selfrm", # distributino based
                      "Stuart", "Stuart.selfrm", # summary stats based
                      "min", "min.selfrm",
                      
                      "MC1", "MC1.selfrm",
                      "MC2","MC2.selfrm",
                      "MC3", "MC3.selfrm", # MC based
                      
                      
                      "CEMC.s", "CEMC.s.selfrm",
                      "CEMC.k","CEMC.k.selfrm", # optimization based
                      
                      
                      "Mean", "Mean.selfrm",
                      "Geomean","Geomean.selfrm",
                      "L2norm", "L2norm.selfrm",
                      "Median", "Median.selfrm")#20-27

my_method <- "BayeSPeer"
other_methods <- setdiff(methods_to_compare, my_method)
other_methods <- factor(other_methods,
                        levels =  methods_to_compare)

width = 0   # disables horizontal lines at top/bottom
theme(panel.border = element_rect(color = "grey40", fill = NA))


# order_methods <- rev(unique(df$method[!str_detect(df$method, "selfrm")]))
# df$method = factor(df$method,
#                                levels=order_methods)

df_main  <- df %>% filter(method != "BayeSPeer")
df_bayes <- df %>% filter(method == "BayeSPeer")

# Colors
# method_colors <- c(
#   "BayeSPeer" = "gray30",
#   "BayePeer"  = "tomato3",
#   "BIRRA"     = "#3778BF",#"#999999",
# 
#   "RRA"       = "#F781BF",  # soft magenta
#   "Stuart"    = "#4DAF4A",  # bright green
#   "min"       = "deepskyblue4",
# 
#   "MC1"       = "#984EA3",  # turquoise
#   "MC2"       = "firebrick",
#   "MC3"       = "#F28E2B", #"darkviolet",
#   "CEMC.k"    = "red",  # gold/yellow
#   "CEMC.s"    = "olivedrab4",  # forest green
# 
#   "Mean"      = "#7B4CC2",  # vivid blue
#   "Geomean"   = "gold",  # deep violet
#   "L2norm"    = "salmon",  # purple
#   "Median"    = "mediumpurple4"
# 
# )
# 
# method_shapes <- c(
#   "BayeSPeer" = 1,
#   "Stuart"    = 0,
#   "MC3"       = 2,
#   "CEMC.k"    = 3,
#   "Mean"      = 4,
#   "BIRRA"     = 5,
# 
#   "BayePeer"  = 6,
# 
#   "RRA"       = 7,
#   "min"       = 8,
# 
#   "MC1"       = 9,
#   "MC2"       = 10,
#   "CEMC.s"    = 11,
# 
#   "Geomean"   = 12,
#   "L2norm"    = 13,
#   "Median"    = 14
# )
method_colors <- c(
  "BayeSPeer" = "gray30",
  "BayePeer"  = "tomato3",
  "BIRRA.selfrm"     = "#3778BF",#"#999999",

  "RRA.selfrm"       = "#F781BF",  # soft magenta
  "Stuart.selfrm"    = "#4DAF4A",  # bright green
  "min.selfrm"       = "deepskyblue4",

  "MC1.selfrm"       = "#984EA3",  # turquoise
  "MC2.selfrm"       = "firebrick",
  "MC3.selfrm"       = "#F28E2B", #"darkviolet",
  "CEMC.k.selfrm"    = "red",  # gold/yellow
  "CEMC.s.selfrm"    = "olivedrab4",  # forest green

  "Mean.selfrm"      = "#7B4CC2",  # vivid blue
  "Geomean.selfrm"   = "gold",  # deep violet
  "L2norm.selfrm"    = "salmon",  # purple
  "Median.selfrm"    = "mediumpurple4"

)

method_shapes <- c(
  "BayeSPeer" = 1,
  "Stuart.selfrm"    = 0,
  "MC3.selfrm"       = 2,
  "CEMC.k.selfrm"    = 3,
  "Mean.selfrm"      = 4,
  "BIRRA.selfrm"     = 5,

  "BayePeer.selfrm"  = 6,

  "RRA.selfrm"       = 7,
  "min.selfrm"       = 8,

  "MC1.selfrm"       = 9,
  "MC2.selfrm"       = 10,
  "CEMC.s.selfrm"    = 11,

  "Geomean.selfrm"   = 12,
  "L2norm.selfrm"    = 13,
  "Median.selfrm"    = 14
)


df$method <- factor(df$method, 
                             levels =  methods_to_compare)



library(ggtext)
facet_labeller <- as_labeller(c(
  "n"          = "**n**",
  "rho"        = "**ρ**",
  "muBeta"     = "**μ**<sub>**β**</sub>",
  "sigma2Beta" = "**σ**<sub>**β**</sub><sup>**2**</sup>"
))


p <- ggplot(df, aes(x = Group, y = mean_value, group = method)) +
  
  scale_x_discrete(position = "top") +
  geom_line(
    aes(
      color = method,
      size = ifelse(method == my_method, 1.8, 1.2),
      alpha = ifelse(method == my_method, 0.8, 0.7)
    ),
    position = position_jitter(width = 0.01, height = 0)
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
    shape = 22,
    show.legend = FALSE,
    position = position_jitter(width = 0.01, height = 0)
  ) +
 
  
  # Facets and styling
  facet_grid(~Section, scales = "free_x", space = "free_x",
             labeller = facet_labeller) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_size_identity() +
  scale_alpha_identity() +
  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "gray40", fill = NA),  # panel borders
    strip.text.x = element_markdown(size = 18, margin = margin(b = 0)),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12,angle = 45, vjust=-0.01, hjust = 0),
    axis.text.y = element_text(size = 12), 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.placement = "outside",  # Keeps facet labels above plot
    axis.ticks.length = unit(5, "pt"),  # Increase tick spacing
    panel.spacing = unit(1.2, "lines"),  # Add space between facets if needed
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Top 1 Coverage rate (%)",#"Spearman Correlation",#"
    color = "Method",
    shape = "Method"
   # subtitle = "Robustness checking: true performance ~ folded Normal distribution"
  ) +ylim(c(1,99))
# 1645 x 432

p
ggsave(
  "lineplot_top1.selfrm.png",
  plot = p,
  width = 13.5,
  height = 5,
  dpi = 600
)
