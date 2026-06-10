#---------------------------
library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)


load("RData/02_compare_mix_mc1.RData")
res_mc1 <- results%>% select(-c(BayeSRank, top1_BayeSRank, top3_BayeSRank))

load("RData/02_compare_mix.RData")
res_mix <- merge(results, res_mc1, by = c("n_teams","rho","muBeta","sigmaBeta","sim_id"))

load("RData/02_compare_1000_seed42_42.RData")
results_multi <- rbind.data.frame(results_multi, res_mix)


load("RData/02_compare_selfrm_mix.RData")
res_selfrm <- results_selfrm

load("RData/02_compare_rm_1000_seed42_42.RData")
results_selfrm <- rbind.data.frame(res_selfrm, results_selfrm)

library(data.table)

combined_res <- merge(results_multi, results_selfrm,
                      by = c("n_teams","rho","muBeta","sigmaBeta","sim_id"))%>% 
  as.data.table()

colnames(combined_res)<-gsub("BiGER","BayePeer", colnames(combined_res))
colnames(combined_res)<-gsub("BayeSRank","BayeSPeer", colnames(combined_res))


# Combine data and filter for
methods= c("BayeSPeer",  #"BayePeer",
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
                      
corcolumns <- colnames(combined_res)[colnames(combined_res) %in% methods]

top1columns <- colnames(combined_res)[str_detect(colnames(combined_res),"top1")]

top3columns <- colnames(combined_res)[str_detect(colnames(combined_res),"top3")]


combined_res_long <- combined_res %>% 
  dplyr::select(c('n_teams', 'rho', 'muBeta', 'sigmaBeta',
                  'sim_id', 
                  top3columns))%>%
  
  pivot_longer(
    cols = -c(n_teams, rho, muBeta, sigmaBeta,sim_id), #, post_mean_mubeta, post_median_Varbeta),  
    names_to = "method",  
    values_to = "value"
  ) %>%
  mutate(muBeta = formatC(muBeta, format = "f", digits = 1),
       #  rho = formatC(rho, format = "f", digits = 2),
         n_teams=as.factor(n_teams),
         rho=as.factor(ifelse(rho=="U(0.15, 0.95)", "U(.15,.95)",rho)),
         muBeta=factor(muBeta,
                       levels=c("-1.5","-0.5","0.0","0.5","1.0","1.5","5.0","10.0")),
         sigma2Beta=factor(sigmaBeta,
                          levels = c(0.5, 0.707106781186548, 1, 1.4142135623731, 2),
                          labels = c("0.25", "0.50", "1.00", "2.00", "4.00")),
         
         #method = factor(method, levels = rev(methods))
         method = factor(substr(method,6,100), levels = rev(methods))
  ) 


# Compute mean and standard error correctly
combined_summ_long <- combined_res_long %>%
  group_by(n_teams, rho, muBeta, sigma2Beta, method) %>%
  summarise(
    mean_value = mean(value, na.rm=TRUE) * 100,  # Convert mean to percentage
   # median_value = median(value, na.rm=TRUE) * 100,  
    se = (sd(value, na.rm=TRUE) / sqrt(n())) * 100  # Convert SE to percentage
  ) 

#########################
df<- combined_summ_long

# Define default values
default_n_teams <- "15"
default_rho <- "0.55"
default_muBeta <- "1.0"
default_sigma2Beta <- "0.50"

# Create 'Group' and 'Section' columns to categorize parameters correctly
df <- combined_summ_long %>%mutate(
  Group = case_when(
    n_teams == default_n_teams & rho == default_rho & muBeta == default_muBeta & sigma2Beta == default_sigma2Beta ~"default", 
    rho == default_rho & muBeta == default_muBeta & sigma2Beta == default_sigma2Beta~ as.character(n_teams),
    n_teams == default_n_teams & muBeta == default_muBeta & sigma2Beta == default_sigma2Beta~ as.character(rho),
    n_teams == default_n_teams & rho == default_rho & sigma2Beta == default_sigma2Beta ~ as.character(muBeta),
    n_teams == default_n_teams & rho == default_rho & muBeta == default_muBeta~ as.character(sigma2Beta),
    TRUE ~ NA_character_
  )) 


df_default <- df %>%
  filter(n_teams == default_n_teams, 
         rho == default_rho, 
         muBeta == default_muBeta,
         sigma2Beta == default_sigma2Beta) 

df_default_n <- df_default%>% mutate(Group="15", Section="n")
df_default_rho <- df_default%>% mutate(Group="0.55", Section="rho")
df_default_muBeta <- df_default%>% mutate(Group="1.0", Section="muBeta")
df_default_sigma2Beta <- df_default%>% mutate(Group="0.50", Section="sigma2Beta")

df_expanded <-df %>% filter(Group !="default")%>%
  mutate(
    Section = case_when(
      Group %in% unique(n_teams) ~ "n",
      Group %in% unique(rho) ~ "rho",
      Group %in% unique(muBeta) ~ "muBeta",
      Group %in% unique(sigma2Beta) ~ "sigma2Beta",
    TRUE ~ NA_character_
    )
    )%>%bind_rows(df_default_n, 
                  df_default_rho,
                  df_default_muBeta,
                  df_default_sigma2Beta) %>%  # Add back the original data
  mutate(
    label = paste0(formatC(mean_value, format = "f", digits = 1), 
                   " (", formatC(se, format = "f", digits = 1), ")")
  ) %>%
  mutate(Group= factor(Group, levels=c(levels(n_teams),
                                       levels(rho),
                                       levels(muBeta),
                                       levels(sigma2Beta))),
         Section=factor(Section, levels=c("n","rho","muBeta","sigma2Beta")))


saveRDS(df_expanded,"RData/03_processed_top3_1000_seed42_42.rds")
#---------------------------
# Testing testing 
ggplot(df_expanded, 
       aes(x = Group, y = method, fill = mean_value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3, color = "black") +  # Updated text labels
  # scale_fill_gradient2(low = "deepskyblue4", mid = "white", high = "tomato3", midpoint = 50) +  # Adjusted midpoint for percentages
  scale_fill_gradientn(
    colors = c(#"royalblue4",
      "deepskyblue4", 
               "cadetblue3",
      #"lightblue", 
               "beige", 
               "salmon1","tomato3"), 
    values = scales::rescale(c(0, 50, 75, 90,95)),  # Adjust breakpoints
    breaks = c(50, 75, 95),  # Legend ticks
    labels = c("50%","75%", #"90-95%",
               "95%")  # Custom labels
  ) +
  facet_grid(~Section, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    title = "testing",
    x = NULL, y = "Method", fill = "Spearman\nCorrelation"
  ) +
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.text.y = element_text(size = 14, face = "bold"),
    strip.placement = "outside",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

  )




