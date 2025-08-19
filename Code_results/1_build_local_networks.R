##################################
# date: 09-2024
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
##################################

# Load Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(bipartite)
library(test)

# Landscape data  ####
landscapes <- data.frame(forest_cover=double(),
                         p_value=double(),
                         patch_no=integer(),
                         patch_state=integer())

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("./Data/landscapes/forestcover_",f,"_pvalue_",p,"_edge.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l))
    landscapes <- rbind(landscapes, l)
  }
}

landscapes <- landscapes %>%
  dplyr::rename(patch_state="V1") %>%
  left_join(., data.frame(patch_no=1:2500, expand.grid(Y=50:1, X=1:50)))

land_pv = landscapes %>% 
  #select(- c(X, Y)) %>% 
  mutate(Land_use = 
           case_when(patch_state  == 5 ~ "Core",
                     patch_state  == 4 ~ "Edge",
                     patch_state  == 3 ~ "Edge",
                     patch_state  == 2 ~ "Edge",
                     patch_state  == 1 ~ "Edge",
                     patch_state  ==  0 ~ "Disturbed",
                     TRUE ~ NA)) %>% 
  select(-patch_state) %>% 
  mutate(Fragmentation =  case_when(p_value  == 0.01 ~ "High",
                                    p_value  == 0.25 ~ "Medium",
                                    p_value  == 0.5 ~ "Low"))

str(land_pv)

grouped <- land_pv %>%
  select(p_value, forest_cover, Land_use, patch_no) %>% 
  group_by(p_value, forest_cover, Land_use) %>%
  dplyr::summarise(count = n(), .groups = 'drop')

# Creating new columns for forest and non-forest patch counts
grouped <- grouped %>%
  mutate(
    core_patch_count = ifelse(Land_use == 'Core', count, 0),
    edge_patch_count = ifelse(Land_use == 'Edge', count, 0),
    disturbed_patch_count = ifelse(Land_use == 'Disturbed', count, 0)
  )

# Summing up the patch counts for each forest_cover and fragmentation_level
final_df <- grouped %>%
  group_by(p_value, forest_cover) %>%
  dplyr::summarise(
    core_patch_count = sum(core_patch_count),
    edge_patch_count = sum(edge_patch_count),
    disturbed_patch_count = sum(disturbed_patch_count),
    .groups = 'drop'
  )

land_pv<- inner_join(land_pv, final_df)

write.csv(land_pv, "./Data/land_pv.csv")

# Load Model data #####

host_info <- read.csv("./Results/host_info.csv") %>% column_to_rownames(var = "X")
prst_info <- read.csv("./Results/prst_info.csv") %>% column_to_rownames(var = "X")
Minc <- read.csv("./Data/M_inc.csv") %>% column_to_rownames(var = "X") %>% 
  rename('juquitiba-like_virus_strain_on576' = 'juquitiba.like_virus_strain_on576')

host_df <- host_info %>%
  rownames_to_column(var= "species") %>% 
  mutate(guild = "resources")

prst_df <- prst_info %>% 
  mutate(guild = "consumers",
         ZoonoticStatus = case_when(
           ZoonoticStatus == 0 ~ "Non-zoonotic",
           ZoonoticStatus == 1 ~ "Zoonotic"))

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

# Separate data   ####
str(model_results)

model_p <- model_results %>%
  filter(guild == "consumers") %>%
  full_join(land_pv) %>%
  left_join(prst_df) %>% 
  select(forest_cover,p_value, patch_no, Land_use, Name, species, z, ZoonoticStatus)

model_h <- model_results %>%
  filter(guild == "resources") %>% 
  full_join(land_pv) %>%
  mutate(species = as.character(species)) %>% 
  inner_join(host_df) %>% 
  select(forest_cover,p_value, patch_no, Land_use, Name, z, pf)

unique(model_h$Name)

# Build Local networks  ####
# Get the names of host and parasite species
Minc <- as.matrix(Minc)

host_species <- rownames(Minc)
parasite_species <- colnames(Minc)

for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){

    result_hosts <-  model_h %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f) 
    
    result_parasites <-model_p %>% 
      filter(p_value == p) %>% 
      filter(forest_cover == f) 

  # Number of sites
  num_sites <- unique(result_parasites$patch_no)

  # List to store the interaction matrices for each site
  var_name <- paste0("p", p, "_lnet_", f)
  
  site_list <- vector("list", length(num_sites))
  
  names(site_list) <- paste0("p_", num_sites)

  # Loop over each site
    for (site in 1:length(site_list)) {
  
  # Get the presence data for this site
    host_present <- result_hosts[result_hosts$patch_no == site, ]
    parasite_present <- result_parasites[result_parasites$patch_no == site, ]
 
    current_incidence <- as.matrix(Minc[host_present$Name, parasite_present$Name])
    current_incidence <- as.matrix(current_incidence[,colSums(current_incidence)>0, drop = FALSE])
    current_incidence <- as.matrix(current_incidence[rowSums(current_incidence)>0, , drop = FALSE])

  # Store the site-specific interaction matrix in the list
  
    site_list[[paste0("p_", site)]] <- current_incidence

  }
  
  assign(var_name, site_list)
  
    }
  }

  for(p in c(0.01,0.25,0.5)){
    for(f in c(0.1,0.3,0.5,0.7)){
    
    var_name <- paste0("p", p,"_lnet_", f)
    
    filtered_matrices <- lapply(get(paste0("p", p,"_lnet_", f)), function(mat) {
      if (nrow(mat) > 2 && ncol(mat) > 2) {
        return(mat)
      } else {
        return(NULL)
      }
    })
    
    # Remove NULL elements (i.e., matrices that were too small)
    filtered_matrices <- Filter(Negate(is.null), filtered_matrices)
    assign(var_name, filtered_matrices)
    
      }
    }

save(p0.01_lnet_0.1,file="./Data/p0.01_lnet_0.1.RData")
save(p0.01_lnet_0.3,file="./Data/p0.01_lnet_0.3.RData")
save(p0.01_lnet_0.5,file="./Data/p0.01_lnet_0.5.RData")
save(p0.01_lnet_0.7,file="./Data/p0.01_lnet_0.7.RData")

save(p0.25_lnet_0.1,file="./Data/p0.25_lnet_0.1.RData")
save(p0.25_lnet_0.3,file="./Data/p0.25_lnet_0.3.RData")
save(p0.25_lnet_0.5,file="./Data/p0.25_lnet_0.5.RData")
save(p0.25_lnet_0.7,file="./Data/p0.25_lnet_0.7.RData")

save(p0.5_lnet_0.1,file="./Data/p0.5_lnet_0.1.RData")
save(p0.5_lnet_0.3,file="./Data/p0.5_lnet_0.3.RData")
save(p0.5_lnet_0.5,file="./Data/p0.5_lnet_0.5.RData")
save(p0.5_lnet_0.7,file="./Data/p0.5_lnet_0.7.RData")

## Build a list for all local interactions   ####
interaction_list <- list()
idx <- 1  # Index for tracking list position

for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){
    
      local_net <- get(paste0("p", p,"_lnet_", f))
      
      unique_patch <- names(local_net)
      
      # Loop through each patch matrix
      for(patch in unique_patch) {
        
        # Convert the interaction matrix to a data frame and add metadata
        df <- as.data.frame(as.table(local_net[[patch]])) %>% 
          mutate(patch_id = sub(".*_(.*)", "\\1", patch), 
                 forest_cover = f, 
                 p_value = p) %>% 
          dplyr::rename(host_species = Var1, 
                        parasite_species = Var2, 
                        value = Freq) %>% 
          filter(value == 1)
        
        # Store the data frame in the list
        interaction_list[[idx]] <- df
        idx <- idx + 1
      }
        
    }
  }


# Combine all data frames into one
interaction_df <- dplyr::bind_rows(interaction_list)

# Save data frames

write.csv(interaction_df, "./Results/interaction_df.csv")



