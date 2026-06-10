#-------------
library(dplyr)
library(readr)
library(tidyverse)
library(scales)  # for rescale

df_expanded <-readRDS("RData/03_processed_corr_1000_seed42_42.rds")
df_expanded <- df_expanded %>% filter(method !="BayePeer")


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

n_methods <- length(levels(df_expanded$method))
y_lines <- c(seq(2.5, n_methods - 0.5, by = 2))


df_expanded$method = factor(df_expanded$method, 
                            levels = rev(methods_to_compare))
df <- df_expanded %>%
  mutate(rho=as.factor(rho))%>%
  group_by(Group, Section) %>%
  mutate(
    rank_in_column = rank(mean_value, ties.method = "last")  # rank 1 = best
  ) %>%
  ungroup()


library(ggtext)


# Step 1: Rank methods within each Group + Section (lowest rank = best)
df <- df %>% group_by(Group, Section) %>%
  mutate(
    rank_in_column = rank(mean_value, ties.method = "first")  # 1 = best
  ) %>%
  ungroup()

# Step 2: Number of unique ranks (e.g. 14)
n_ranks <- length(unique(df$rank_in_column))

# Step 3: Plot

facet_labeller <- as_labeller(c(
  "n"          = "**n**",
  "rho"        = "**ρ**",
  "muBeta"     = "**μ**<sub>**β**</sub>",
  "sigma2Beta" = "**σ**<sub>**β**</sub><sup>**2**</sup>"
))

section_labeller <- as_labeller(
  c(
    n = "n",
    rho = "rho",
    muBeta = "mu[beta]",
    sigma2Beta = "sigma[beta]^2"
  ),
  label_parsed
)

p <- ggplot(df%>%rename, aes(x = Group, y = method)) +
  geom_point(
    aes(size = mean_value, fill = rank_in_column),
    shape = 21,
    color = "white", stroke = 0.4
  ) +
  
  # Your custom color scale (maps to rank)
  scale_fill_gradientn(
    colours = c("deepskyblue4", 
                "cadetblue3",
                "lightblue",
                "beige",
                "darksalmon",
                "tomato3",
                "firebrick"),
    values = rescale(seq(1, n_ranks, length.out = 7)),  # spreads palette evenly over rank range
    name = "Performance\nin Setting",
    breaks = c(1, ceiling(n_ranks/2), n_ranks),
    labels = c("Worst", "Middle", "Best")
  ) +
  
  scale_size(
    range = c(0.2, 7),
    name = "Spearman\nCorrelation",
    breaks = c(30, 50, 70, 90, 95),
    labels = c("0.3","0.5","0.7", "0.9","0.95")
  )  +
 # geom_text(aes(label = label), size = 3, color = "black") +
  
  guides(
    size = guide_legend(
      override.aes = list(shape = 21,
                          color = "black", 
                          fill = "white", 
                          stroke = 0.1)
    )
  ) +
  facet_grid(~Section, scales = "free_x", space = "free_x",
             labeller = facet_labeller) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(
   # title = "Simulation Results: Spearman Correlation\n(Color = Relative Performance in Column, Size = Actual Correlation Value)",
   title = NULL,
     x = NULL, y = "Method"
  ) +
  theme(
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.5, "lines"),   # widens spacing between panels
    panel.border = element_rect(color = "gray0", fill = NA, linewidth = 0.5), # optional border around each panel

    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 12),

    strip.text.x = element_markdown(size = 18, margin = margin(b = 0)),
    axis.text.x = element_text(size = 12,angle = 60, vjust=-0.01, hjust = 0)
)+
  geom_hline(yintercept = y_lines, 
                 color = "gray60", linewidth = 0.2)
    
p
ggsave(
  "heatmap_corr_circle.png",
  plot = p,
  width = 13.5,
  height = 7.5,
  dpi = 600
)



# summary stats
# Step 3: Plot
ggplot(df%>%rename, aes(x = Group, y = method)) +
  geom_point(
    aes(size = mean_value),
    shape = 22,
    color = "white", stroke = 0.4
  )  +
  geom_text(aes(label = label), size = 3, color = "black") +
  
  guides(
    size = guide_legend(
      override.aes = list(shape = 22,
                          color = "black", 
                          fill = "white", 
                          stroke = 0.1)
    )
  ) +
  
  facet_grid(~Section, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  labs(
    # title = "Simulation Results: Spearman Correlation\n(Color = Relative Performance in Column, Size = Actual Correlation Value)",
    title = NULL,
    x = NULL, y = "Method"
  ) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 14, face = "bold"),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.5, "lines"),   # widens spacing between panels
    panel.border = element_rect(color = "gray0", fill = NA, linewidth = 0.5), # optional border around each panel
    
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )+geom_hline(yintercept = y_lines, 
               color = "gray60", linewidth = 0.2)

