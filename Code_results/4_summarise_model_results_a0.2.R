###################################
# title: Results from Output without rewiring and same different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 05-2024
###################################

# Load Packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)

#Load data ####
land_pv<- read.csv("./Data/land_pv.csv") %>% 
  column_to_rownames(var = "X")

host_info <- read.csv("./Results/host_info.csv") %>% 
  column_to_rownames(var = "X")

prst_info <- read.csv("./Results/prst_info.csv") %>% 
  column_to_rownames(var = "X")

net.metric_mean <- read.csv("./Results/net.metric_mean.csv") %>% column_to_rownames(var = "X")

interaction_df <- read.csv("./Results/interaction_df.csv") %>% column_to_rownames(var = "X")

m_p_trait_match<- read.csv("./Results/m_p_traitMatch.csv") %>% column_to_rownames(var = "X")

# Model results   ####
model_results <- data.frame(forest_cover=double(),
                            p_value=double(),
                            patch_no=integer(),
                            guild= vector(),
                            species = vector(),
                            z= double())


for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    
    l <- read.csv(paste0("./Output_v2b/metanetwork_fc",f,"_p",p,"_er1_ec0.5_cr1_cc0.6_mr0.5_mc0.7_phi0.5_alpha0.2_eps5.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=as.vector(l$patch_no),
                    guild=as.vector(l$guild),
                    species = as.vector(l$species),
                    z= as.vector(l$z))
    
    model_results <- rbind(model_results, l)
    
  }
}

# prepare species data 
host_df <- host_info %>%
  rownames_to_column(var= "species") %>% 
  mutate(guild = "resources")

prst_df <- prst_info %>% 
  mutate(guild = "consumers",
         ZoonoticStatus = case_when(
           ZoonoticStatus == 0 ~ "Non-zoonotic",
           ZoonoticStatus == 1 ~ "Zoonotic"))

############### Plot landscapes ###############
land_pv$Land_use <- factor(land_pv$Land_use , levels = c("Core", "Edge", "Disturbed"))

Land_plot=ggplot(data=land_pv, aes(x=X, y=Y, fill=as.factor(Land_use))) +
  geom_tile() +
  facet_grid(p_value~forest_cover) +
  coord_fixed(ratio=1) +
  scale_x_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_fill_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank())
Land_plot

ggsave("./Figures/Land_plot.png", width = 15, height = 15, units = "cm",  dpi = 600)

###################### Plot host pf ################
pf_plot<- host_df %>%
  ggplot(aes(x = pf)) +
  geom_density(aes(fill = as.factor(HostOrder), color = as.factor(HostOrder)), alpha = 0.3) +
  labs(x = "Occurrence frequency in forest areas")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
pf_plot

ggsave("./Figures/Host_of.png", width = 15, height = 12, units = "cm",  dpi = 600)

#Summarizing Model Data  ####
## Calculate parasite extinction rate data ####
extinction <-  model_results %>%
  filter(guild == "consumers") %>% 
  inner_join(prst_df) %>% 
  select(forest_cover, p_value, Name, species) %>% 
  mutate(Name = as.factor(Name)) %>%  
  group_by(forest_cover, p_value) %>%
  reframe(Prst_richness= as.numeric(n_distinct(Name)), Ext_rate = (103 - Prst_richness)/103) %>% 
  mutate(Landscape = paste(forest_cover, p_value, sep = "_")) %>% 
  mutate(frag_value = 1 -p_value)

## Summarise host data ####

p_res_h_fc <- land_pv %>%
  select(forest_cover, p_value, patch_no, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count) %>% 
  full_join(model_results) %>%
  filter(guild == "resources") %>% 
  mutate(species = as.character(species)) %>% 
  inner_join(host_df) %>% 
  select(forest_cover, p_value, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count, Name, z, pf) %>% 
  dplyr::mutate(Name = as.factor(Name)) %>%  
  group_by(forest_cover, p_value, Land_use, core_patch_count, edge_patch_count, disturbed_patch_count, Name, pf) %>%
  dplyr::summarise(n= n()) %>% 
  dplyr::mutate(Abd_type= case_when(
    Land_use == "Core" ~ n/core_patch_count,
    Land_use == "Edge" ~ n/edge_patch_count,
    Land_use == "Disturbed" ~ n/disturbed_patch_count)) %>% 
  ungroup() %>% 
  dplyr::select(- c(core_patch_count, edge_patch_count, disturbed_patch_count)) %>% 
  dplyr::rename(host_species = Name) %>% 
  mutate_if(is.character, as.factor)

h_abd_df <- land_pv %>% 
  select(- c(X, Y)) %>% 
  dplyr::rename(patch_id = patch_no) %>% 
  mutate(patch_id = as.character(patch_id)) %>% 
  right_join(interaction_df) %>% 
  unite("Interaction", parasite_species, host_species, sep= "/.") %>% 
  group_by(forest_cover,p_value, Land_use, core_patch_count, edge_patch_count,  disturbed_patch_count, Interaction, value) %>% 
  dplyr::summarise(Count = sum(value), .groups = "drop") %>% 
  separate(Interaction, into = c("parasite_species", "host_species"), sep= "/.") %>% 
  inner_join(p_res_h_fc)  %>% 
  select(forest_cover,p_value, Land_use, parasite_species, Abd_type) %>% 
  group_by(forest_cover,p_value, Land_use, parasite_species) %>% 
  dplyr::summarise(Host_abd = mean(Abd_type), Host_abd_sd = sd(Abd_type)) %>%
  mutate(CoefVar_host = Host_abd_sd/Host_abd) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.factor) 

str(h_abd_df)

## Parasite and Host attributes ####
prst_net <- prst_df %>% 
  select(Name, ZoonoticStatus, Degree, ParGroup) %>% 
  dplyr::rename(parasite_species = Name)

host_net<- host_df %>% 
  dplyr::rename(host_species = Name)

## Calculate parasite relative degree ####

nh_df <- land_pv %>% 
  select(- c(X, Y)) %>% 
  dplyr::rename(patch_id = patch_no) %>% 
  mutate(patch_id = as.character(patch_id)) %>% 
  right_join(interaction_df) %>%
  select(forest_cover,p_value, Land_use, patch_id, parasite_species) %>%    
  group_by(forest_cover,p_value, Land_use, patch_id, parasite_species) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  inner_join(prst_net) %>% 
  mutate(Degree_rel = (Count/Degree)) %>% 
  group_by(forest_cover, p_value, Land_use, parasite_species, ZoonoticStatus, ParGroup, Degree) %>%
  dplyr::summarise(Mean_count = mean(Count), SD_rel= sd(Degree_rel), Mean_rel = mean(Degree_rel), .groups = "drop") %>% 
  mutate(zs_dif = Mean_rel/SD_rel)

str(nh_df)

## Gather outcomes in one dataframe ####

outcomes_model <- nh_df %>%
  select(- Mean_count) %>% 
  inner_join(h_abd_df) %>% 
  inner_join(net.metric_mean) %>% 
  mutate_if(is.character, as.factor) %>% 
  inner_join(m_p_trait_match) %>% 
  select(-c(A_mtheta, A_nopartners, Count)) %>% 
  inner_join(extinction) %>% 
  inner_join(betadiv.b_true)
str(outcomes_model)

# Analyse outcomes distribution
outcomes_model$Mean_rel_sc <- as.numeric(scale(outcomes_model$Mean_rel))
outcomes_model$Amatch_sc <- as.numeric(scale(outcomes_model$A_match))
outcomes_model$Ext_rate_sc <- as.numeric(scale(outcomes_model$Ext_rate))
outcomes_model$Host_abd_sc <- as.numeric(scale(outcomes_model$Host_abd))

# PCA axis 
pca_data <- outcomes_model %>%
  select(Mean_size, Mean_connectance, Mean_NODF, Mean_mod) %>%
  scale()  
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
biplot(pca_result)
axis_pca <- data.frame(pca_result[["rotation"]])

outcomes_model$PCA_axis1 <- pca_result$x[, 1] #connectance positive, size and modularity negative

# save data
write.csv(outcomes_model, "./Results/outcomes_model.csv")
