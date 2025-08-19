###################################
# title: Beta-diversity analysis from Output without rewiring and different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 03-2025
###################################

# Load packages ####

pacman::p_load(car, igraph, plotly, data.table, betapart, vegan, reshape2, ggplot2, betalink, knitr, plyr, bipartite, igraph, tidyverse, janitor, stringr, dplyr, hrbrthemes, viridis, ggpubr, ggsci, float, ggpmisc, lmPerm, performance, cowplot, sjPlot)

# Load data ####

land_pv<- read.csv("./Data/land_pv.csv") %>% 
  column_to_rownames(var = "X")

interaction_df <- read.csv("./Results/interaction_df.csv") %>% column_to_rownames(var = "X")

# Work on dataframe ####

land_pv_clean <- land_pv %>%
  select(-c(X, Y)) %>%
  rename(patch_id = patch_no) %>%
  mutate(patch_id = as.character(patch_id))

interaction_clean <- interaction_df %>%
  mutate(patch_id = as.character(patch_id))

# Sample local networks by habitat type ####
set.seed(123)

# prepare dataframe

patch_count <- land_pv %>% 
  select("forest_cover", "p_value", "core_patch_count", "edge_patch_count", "disturbed_patch_count") %>%
  dplyr::rename(Core = core_patch_count, Edge =edge_patch_count, Disturbed=disturbed_patch_count) %>% 
  gather("count", "value", 3:5)

# Sample networks

forest_vals <- c(0.1, 0.3, 0.5, 0.7)
p_vals <- c(0.01, 0.25, 0.5)
land_use_vals <- c("Core", "Edge", "Disturbed")

metanetwork_list <- list()

  for (f in forest_vals) {
  for (p in p_vals) {
    
    list <- list()
    
    for (h in land_use_vals) {
      
      patch_count_i <- patch_count %>% 
        filter(forest_cover == f,
               p_value == p,
               count == h) %>%
        dplyr::select(count, value) %>% 
        dplyr::summarise(m = max(4, min(round(value * 0.1)))) %>% 
        pull(m)
        
      patch_filtered <- land_pv_clean  %>% 
      filter(forest_cover == f,
             p_value == p,
             Land_use == h) %>% 
        distinct(patch_id) %>% 
      slice_sample(n = patch_count_i) %>% 
        pull(patch_id)
      

    for (j in patch_filtered) {
      
      var_name <- paste0("networks_", p,"_", f,"_",h,"_",j)
      
      n_filtered_data <- land_pv_clean %>% 
        right_join(interaction_clean) %>% 
        filter(forest_cover == f,
               p_value == p,
               Land_use == h, 
               patch_id == j)
      
      n_filtered_data <- n_filtered_data %>% 
        group_by(host_species, parasite_species) %>% 
        dplyr::mutate(Count = n()) %>% 
        ungroup()
      
      reg_web = n_filtered_data %>%  
        select(host_species, parasite_species, Count) %>% 
        #filter(Count > 1L) %>% 
        pivot_wider(names_from= parasite_species, values_from= Count, values_fill = 0) %>% 
        remove_rownames %>% 
        column_to_rownames(var="host_species")
    
      list[[var_name]] <- as.matrix(reg_web)
      
    }
      
    }
    
    var_list <- paste0("land_", p,"_", f)
    metanetwork_list[[var_list]] <- list
  }
}

# Beta indexes ####
## Betalink metric ####

beta_result<-list()

for (f in forest_vals) {
  for (p in p_vals) {
    
    var_name<-paste0("beta_", p,"_", f)
    
    meta<- metanetwork_list[[paste0("land_", p,"_", f)]]
    
    networks <- betalink::prepare_networks(meta, directed = FALSE)
    
    betadiversity <- network_betadiversity(networks)
    
    beta_result[[var_name]]<- data.frame(betadiversity) 
    
    print(var_name)
    
  }
  }


## Bipartite metric ####

# Beta with sampled local networks

beta_result_bi<-list()

for (f in forest_vals) {
  for (p in p_vals) {
    
    var_name<-paste0("beta_", p,"_", f)
    
    meta<- metanetwork_list[[paste0("land_", p,"_", f)]]
    
    network_betalink <- bipartite::webs2array(meta)
    
    betadiv.b <- betalinkr_multi(network_betalink, partition.st=T, binary=TRUE)
    
    beta_result_bi[[var_name]]<- data.frame(betadiv.b) 
    
    print(var_name)
    
  }
}

# Work on results ####
## results with sampled local networks ####
str(beta_result_bi[["beta_0.01_0.1"]])

df_beta<- names(beta_result_bi) 
beta_df<- data.frame()

for (i in df_beta) {
  
test<- beta_result_bi[[i]]

betadiv_df <- test %>% 
  separate(i, into = c("net_i", "p_value_i", "forest_cover_i", "Land_use_i", "patch_id_i"), sep= "_") %>% 
  separate(j, into = c("net_i", "p_value_j", "forest_cover_j", "Land_use_j", "patch_id_j"), sep= "_") %>% 
  select(-c("net_i", S, OS, ST.lh)) %>% 
  mutate(same_habitat = ifelse(Land_use_i == Land_use_j, TRUE, FALSE)) %>%
  unite("p_value", p_value_i, p_value_j, sep= "_") %>%
  unite("forest_cover", forest_cover_i, forest_cover_j, sep= "_")  %>% 
  unite("join_hab", Land_use_i, Land_use_j, sep= "_") %>% 
  mutate(p_value = case_when(
    p_value == "0.01_0.01" ~ "0.01",
    p_value == "0.25_0.25" ~ "0.25",
    p_value == "0.5_0.5" ~ "0.5",
    TRUE ~ p_value),
    forest_cover = case_when(
    forest_cover == "0.1_0.1" ~ "0.1",
    forest_cover == "0.3_0.3" ~ "0.3",
    forest_cover == "0.5_0.5" ~ "0.5",
    forest_cover == "0.7_0.7" ~ "0.7",
    TRUE ~ forest_cover)) %>% 
    mutate(forest_cover = as.double(forest_cover), p_value = as.double(p_value)) %>%
    mutate(Landscape = paste(forest_cover, p_value, sep = "_"),
           frag_value = 1 - p_value) %>%
    select(forest_cover, frag_value, Landscape, join_hab, same_habitat, WN, ST, ST.l, ST.h)

beta_df <- rbind(beta_df, betadiv_df)

}
  
  # save data
  write.csv(betadiv_df, "./Results/beta_values_random.csv")
  
