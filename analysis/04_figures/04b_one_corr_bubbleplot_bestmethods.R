#-------------
library(dplyr)
library(readr)
library(tidyverse)
library(scales)  # for rescale

df_expanded <-readRDS("RData/03_processed_corr_seed42_42.rds")
df_expanded <- df_expanded %>% filter(method !="BayePeer")

methods_to_compare= c("BayeSPeer",  "BayePeer",
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

df_expanded$method = factor(df_expanded$method, 
                            levels = methods_to_compare)

df_expanded <- df_expanded %>%
  ungroup() %>%
  mutate(
    category = case_when(
      
      # --- Bayesian ---
      str_detect(method, "Baye|BIRRA") ~ "Bayesian",
      
      # --- Distribution-based (RRA family) ---
      str_detect(method, "min|RRA|Stuart") ~ "Distribution-based",
      
      # --- Optimization-based ---
      str_detect(method,"CEMC") ~ "Optimization-based",
      
      # --- Markov chain-based ---
      str_detect(method,"MC") ~ "Markov chain–based",
      
      # --- Non-optimization-based (heuristics like Borda) ---
      str_detect(method,"Mean|Median|Geo|L2") ~ 
        "Non-optimization-based",
      
      TRUE ~ "Other"
    )
  )

df_expanded %>% 
  group_by(method, category) %>%
  summarise(mean_mean_corr = mean(mean_value)) %>%
  group_by(category) %>% 
  mutate(best_corr=max(mean_mean_corr),
         best_method=method[mean_mean_corr==best_corr]) -> pick_methods

(c(as.character(unique(pick_methods$best_method)),"BIRRA" )-> best_methods)

n_methods <- length(best_methods)
y_lines <- c( seq(1.5, n_methods-1.5, by = 1))


df <- df_expanded%>%filter(method %in% best_methods, n_teams!=50)%>%
  group_by(Group, Section) %>%
  mutate(
    color_scaled = rescale(mean_value, to = c(0, 1))
  ) %>%
  ungroup()

ggplot(df, 
       aes(x = Group, y = method)) +
  geom_point(aes(size = mean_value, fill = color_scaled),
             shape = 22,
             color = "white",      # lighter border than black
             stroke = 0.3) +
  
  guides(
    size = guide_legend(
      override.aes = list(shape = 22, color = "black",
                          fill = "white", stroke = 0.5)
    )
  )+
  scale_size(
             range = c(0.5, 7.5),
             name = "Spearman\nCorrelation",  # <- customize this label as needed
             breaks= c( 0.5, 0.75, 0.85, 
                        0.9, 0.95, 0.98),
             labels = c("0.50","0.75","0.85",
                        "0.90", "0.95","0.98")
         ) +
  scale_fill_gradientn(
           colors = c("deepskyblue4", 
                      "cadetblue3",
                      "lightblue",
                      "beige", #"salmon1",
                      "darksalmon",
                      "tomato3",
                      "firebrick"#,"firebrick4"
                      ),
           values = rescale(c(#0,0.25,
             0.5, 0.75, 0.85, 
             0.9, 0.95, 0.98)),  # scaled between 0 and 1
           breaks = c(50, 75,85, 90, 95, 98), 
           name = "Spearman\nCorrelation",
           guide = guide_colorbar(barwidth = 1, barheight = 10)
         ) +
  facet_grid(~Section, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    title = "Simulation Results: Spearman Correlation",
    x = NULL, y = "Method"
  ) +
  scale_x_discrete(position = "top")+
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 14, face = "bold"),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(1.5, "lines"),   # widens spacing between panels
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),  # optional border around each panel
    
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )+geom_hline(yintercept = y_lines, 
               color = "gray80", linewidth = 0.3)





##########

df <- df_expanded%>%filter(method %in% best_methods)%>%
  group_by(Group, Section) %>%
  mutate(
    rank_in_column = rank(mean_value, ties.method = "last")  # rank 1 = best
  ) %>%
  ungroup()

df$method <- factor(df$method, 
                    levels = c("BayeSPeer", "BIRRA", "Stuart", "MC3", "CEMC.k","Mean"))

# Step 1: Rank methods within each Group + Section (lowest rank = best)
df <- df %>%
  group_by(Group, Section) %>%
  mutate(
    rank_in_column = rank(mean_value, ties.method = "first")  # 1 = best
  ) %>%
  ungroup()

# Step 2: Number of unique ranks (e.g. 14)
n_ranks <- length(unique(df$rank_in_column))

# Step 3: Plot
ggplot(df%>%rename, aes(x = Group, y = method)) +
  geom_point(
    aes(size = mean_value, fill = rank_in_column),
    shape = 22,
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
    range = c(0.5, 7.5),
    name = "Spearman\nCorrelation",
    breaks = c(30, 50, 70, 90, 95),
    labels = c("0.3","0.5","0.7", "0.9","0.95")
  )  +
 # geom_text(aes(label = label), size = 3, color = "black") +
  
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
    

