###################################
# date: 09-2024
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

# Load Packages
library(data.table)
  library(readr)
  library(dplyr)
  library(purrr)

######## Load data ########

  networkName = "metanetwork"
  # Read network incidence matrix
  M_inc <- as.matrix(fread(file.path("./Data", networkName, "Minc.csv"), header = FALSE))
  # Read theta values
  coev_h <- fread(file.path("./Data", networkName, "theta_r.csv"), header = FALSE)
  coev_p <- fread(file.path("./Data", networkName, "theta_c.csv"), header = FALSE)
  
  # Measuring trait-matching
  for(p in c(0.01,0.25,0.5)){
    for(f in c(0.1,0.3,0.5,0.7)){
      
  # Number of host and parasite species
  n_h <- nrow(M_inc)
  n_p <- ncol(M_inc)
  
       # Import simulation results
      result_hosts <-  model_results %>% 
        filter(p_value == p) %>% 
        filter(forest_cover == f) %>% 
        filter(guild== "resources")
      
      result_parasites <- model_results %>% 
        filter(p_value == p) %>% 
        filter(forest_cover == f) %>% 
        filter(guild== "consumers")
      

  theta_values <- data.frame(
    guild = c(rep("hosts", n_h), rep("parasites", n_p)),
    species = c(seq_len(n_h), seq_len(n_p)),
    theta = c(coev_h$V1, coev_p$V1)
  )
  
  # Create empty dataframes
  matching_network <- data.frame(
    patch_no = rep(seq_len(50 * 50)),
    match = rep(NA_real_),
    match_theta = rep(NA_real_),
    match_hosts = rep(NA_real_),
    match_parasites = rep(NA_real_)
  )
  
  matching_species <- rbind(
    data.frame(
      patch_no = rep(seq_len(50 * 50)),
      species = rep(seq_len(n_h),  each = 50 * 50),
      guild = rep("hosts", each = 50 * 50 * n_h),
      match = rep(NA_real_),
      match_theta = rep(NA_real_),
      no_partners = rep(NA_integer_)
    ),
    data.frame(
      patch_no = rep(seq_len(50 * 50)),
      species = rep(seq_len(n_p), each = 50 * 50),
      guild = rep("parasites", each = 50 * 50 * n_p),
      match = rep(NA_real_),
      match_theta = rep(NA_real_),
      no_partners = rep(NA_integer_)
    )
  )
  
 
      
      # Patches with data
      patches <- unique(c(result_hosts$patch_no, result_parasites$patch_no))
      
      for (patch in patches) {
        # Filter data for current patch
        current_hosts <- result_hosts[result_hosts$patch_no == patch,]
        current_hosts <- current_hosts[order(current_hosts$species),]
        current_parasites <- result_parasites[result_parasites$patch_no == patch,]
        current_parasites <- current_parasites[order(current_parasites$species),]
        
        # Hosts and parasites present
        host_sp <- current_hosts$species
        parasite_sp <- current_parasites$species
        n_hosts <- length(host_sp)
        n_parasites <- length(parasite_sp)
        
        # Define the final trait value of the species after coevolution
        z <- c(current_hosts$z, current_parasites$z)
        
        theta <- c(
          theta_values[theta_values$guild == "hosts" & theta_values$species %in% host_sp, "theta"],
          theta_values[theta_values$guild == "parasites" & theta_values$species %in% parasite_sp, "theta"]
        )
        
        # Define the sensitivity parameter
        alpha <- 0.2
        
        # Current incidence matrix
        current_incidence <- as.matrix(M_inc[host_sp, parasite_sp])
        
        # Number of partners of each species
        host_part <- rowSums(current_incidence)
        parasite_part <- colSums(current_incidence)
        
        # Store results in dataframe
        matching_species[matching_species$patch_no == patch & matching_species$guild == "hosts" & matching_species$species %in% host_sp, "no_partners"] <- host_part
        matching_species[matching_species$patch_no == patch & matching_species$guild == "parasites" & matching_species$species %in% parasite_sp, "no_partners"] <- parasite_part
        
        # Calculate trait matching if at least one interaction present
        if (sum(current_incidence) != 0) {
          # Convert to adjacency matrix
          current_adj <- matrix(NA_integer_, nrow = n_hosts + n_parasites, ncol = n_hosts + n_parasites)
          current_adj[1:n_hosts, 1:n_hosts] <- 0
          current_adj[(n_hosts + 1):(n_hosts + n_parasites), (n_hosts + 1):(n_hosts + n_parasites)] <- 0
          current_adj[1:n_hosts, (n_hosts + 1):(n_hosts + n_parasites)] <- current_incidence
          current_adj[(n_hosts + 1):(n_hosts + n_parasites), 1:n_hosts] <- t(current_incidence)
          
          # Calculate the degree of trait matching between species of different guilds
          matching_diff <- outer(z, z, FUN = function(x, y) abs(x - y))
          matching <- exp(-alpha * matching_diff^2)
          
          # Guild trait matching
          hosts_matching <- (sum(matching[1:n_hosts, 1:n_hosts]) - n_hosts) / (n_hosts * n_hosts - n_hosts)
          parasites_matching <- (sum(matching[(n_hosts + 1):(n_hosts + n_parasites), (n_hosts + 1):(n_hosts + n_parasites)]) - n_parasites) / (n_parasites * n_parasites - n_parasites)
          
          # Assign NA to matching between non-interacting species
          matching[is.na(current_adj)] <- NA
          
          # Calculate the mean trait matching value for the network
          network_matching <- mean(matching, na.rm = TRUE)
          
          # Compute the degree of trait matching for each species
          species_matching <- rowMeans(matching, na.rm = TRUE)
          
          # Store results in dataframe
          matching_network[matching_network$patch_no == patch, "match"] <- network_matching
          
          matching_network[matching_network$patch_no == patch, "match_hosts"] <- hosts_matching
          
          matching_network[matching_network$patch_no == patch, "match_parasites"] <- parasites_matching
          
          matching_species[matching_species$patch_no == patch & matching_species$guild == "hosts" & matching_species$species %in% host_sp, "match"] <- species_matching[1:n_hosts]
          
          matching_species[matching_species$patch_no == patch & matching_species$guild == "parasites" & matching_species$species %in% parasite_sp, "match"] <- species_matching[(n_hosts + 1):length(species_matching)]
          
        }
        
        # Calculate the degree of trait matching with theta
        matching_theta <- exp(-alpha * (theta - z)^2)
        
        # Calculate the mean trait matching value for the network
        network_matching_theta <- mean(matching_theta, na.rm = TRUE)
        
        # Store results in dataframe
        matching_network[matching_network$patch_no == patch, "match_theta"] <- network_matching_theta
        
        matching_species[matching_species$patch_no == patch & matching_species$guild == "hosts" & matching_species$species %in% host_sp, "match_theta"] <- matching_theta[1:n_hosts]
        
        matching_species[matching_species$patch_no == patch & matching_species$guild == "parasites" & matching_species$species %in% parasite_sp, "match_theta"] <- matching_theta[(n_hosts + 1):length(matching_theta)]
  
      }
  
      var_name <- paste0("p", p,"TraitMatch_", f)
      
      assign(var_name, list(matching_network = matching_network, matching_species = matching_species))
      
        } 
  }

  
# Organize trait matching data

trait_match = data.frame(NULL)
  
  for(p in c(0.01,0.25,0.5)){
    for(f in c(0.1,0.3,0.5,0.7)){
  
      land_f <- land_pv %>% 
        select(-c(Y, X, core_patch_count, edge_patch_count, disturbed_patch_count)) %>% 
        filter(p_value == p,
               forest_cover == f)
      
  TraitMatch <- get(paste0("p",p ,"TraitMatch_", f))[["matching_species"]] %>%
  filter(match != is.na(match)) %>% 
    mutate(p_value = p,
           forest_cover = f)
  
  trait_match = rbind(trait_match, TraitMatch)
  
    }
  }

# summarise trait matching data

trait_match_mean = data.frame(NULL)

for(p in c(0.01,0.25,0.5)){
  for(f in c(0.1,0.3,0.5,0.7)){
    land_f <- land_pv %>% 
      select(-c(Y, X, core_patch_count, edge_patch_count, disturbed_patch_count)) %>% 
      filter(p_value == p,
             forest_cover == f)
    
    TraitMatch <- get(paste0("p",p ,"TraitMatch_", f))[["matching_species"]] %>%
      filter(match != is.na(match)) %>% 
      inner_join(land_f) %>% 
      dplyr::group_by(Land_use, species, guild) %>% 
      dplyr::summarise(A_match = mean(match),SD_match= sd(match), A_mtheta = mean(match_theta), A_nopartners= mean(no_partners), .groups = "drop") %>% 
      dplyr::mutate(p_value = p,
             forest_cover = f,
             Coef_matchr = SD_match/A_match)
    
    trait_match_mean = rbind(trait_match_mean, TraitMatch)
    
  }
}

# save data
write.csv(trait_match_mean, "./Results/Mean_traitMatch.csv")
write.csv(trait_match, "./Results/Res_traitMatch.csv")

# filter parasite data
prst_df1 <- prst_df %>% 
  mutate(guild = "parasites")
str(prst_df1)

m_p_trait_match <- trait_match_mean %>% 
  filter(guild == "parasites") %>%
  inner_join(prst_df1) %>% 
  dplyr::select(- c(species, ParStatus, guild)) %>% 
  dplyr::rename("parasite_species"=Name)

write.csv(m_p_trait_match, "./Results/m_p_traitMatch.csv")
