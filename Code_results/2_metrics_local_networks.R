###################################
# date: 09-2024
# title: Local Networks topology of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# Project: Towards Dilution Effect
###################################

# Load Packages
library(betalink)
library(bipartite)
library(iNEXT)
library(vegan)
library(ggplot2)
library(network)
library(sna)
library(igraph)
library(ade4)
library(intergraph)
library(tibble)
library(dplyr)
library(tidyverse)
library(plyr)
library(reshape2)
library(cowplot)
library(ggpubr)

# Load data ####

load("./Data/p0.01_lnet_0.1.RData")
load("./Data/p0.01_lnet_0.3.RData")
load("./Data/p0.01_lnet_0.5.RData")
load("./Data/p0.01_lnet_0.7.RData")
load("./Data/p0.25_lnet_0.1.RData")
load("./Data/p0.25_lnet_0.3.RData")
load("./Data/p0.25_lnet_0.5.RData")
load("./Data/p0.25_lnet_0.7.RData")
load("./Data/p0.5_lnet_0.1.RData")
load("./Data/p0.5_lnet_0.3.RData")
load("./Data/p0.5_lnet_0.5.RData")
load("./Data/p0.5_lnet_0.7.RData")

#Topology Measures ####

Net_metric <- data.frame() 

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){


# Include meta_web in the list of local webs

result_matrices <- get(paste0("p", p,"_lnet_", f))

##Connectance ####

connect_web <- sapply(result_matrices, networklevel, index="connectance")

connect_web <- as.data.frame(connect_web)
connect_web <- connect_web %>%
  rownames_to_column() %>% 
  separate(rowname, into = c("Local", "Metric"), sep = "\\.") %>% 
  select(-Metric)


## Make graph object ####
local_graphs <- list()

for (i in unique(names(result_matrices))) {
  
 net_graph = graph_from_biadjacency_matrix(result_matrices[[i]], weighted = NULL)
 local_graphs[[i]]<- net_graph
 
}

## Make connected graphs ####
net_connect <- list()

for (i in unique(names(local_graphs))) {
  
  net_1<-(local_graphs[[i]])
  
  # Find the components of the graph
  components <- components(net_1)
  
  # Extract component membership and sizes
  membership <- components$membership
  sizes <- components$csize
  
  # Print components
  print(sizes)
  
  if(is_connected(net_1) == FALSE){
    # Connect the components by adding an edge
    # Find representative nodes from each component
    component_ids <- unique(membership)
    
    # Connect each pair of components by adding an edge between them
    for (l in 1:(length(component_ids) - 1)) {
      # Find a node from the current component
      node_from <- which(membership == component_ids[l])[1]
      # Find a node from the next component
      node_to <- which(membership == component_ids[l + 1])[1]
      
      # Add an edge between these nodes
      net_1 <- add_edges(net_1, c(node_from, node_to))
    }
    
  }
  # Check if the graph is now connected
  print(is_connected(net_1))  # Should return TRUE
  
  net_connect[[i]]<-net_1
  
}

## Size ####
network_sizes <- lapply(local_graphs, vcount)
network_edge_counts <- lapply(local_graphs, ecount)

head(network_sizes)
network_sizes = unlist(network_sizes)

network_sizes <- as.data.frame(network_sizes)
network_sizes <- network_sizes %>% 
  rownames_to_column() %>% 
  dplyr::rename("Local"= rowname)

## Modularity #### 
mod.groups <- sapply(net_connect, cluster_fast_greedy, weights = NULL)
like.m.groups <- sapply(mod.groups, modularity)

like.m.groups <- data.frame(Modularity = like.m.groups)
like.m.groups <- like.m.groups %>% 
  rownames_to_column() %>% 
  dplyr::rename("Local"= rowname)

## Nestedness ####

nest.groups <- sapply(result_matrices, nestednodf)
nest.groups= as_tibble(nest.groups)
nest.g = nest.groups[3,c(1:length(result_matrices))]  %>% unnest(cols = c(1:length(result_matrices)))

nest.g<- as.data.frame(t(nest.g))
colnames(nest.g)= c("N.collums", "N.rows", "NODF")

nest.g <- nest.g %>% 
  mutate("Local"= as.character(unique(names(local_graphs))))

## Creating Data frame with all metrics for all networks #### 

groups.metric = full_join(network_sizes, nest.g, by = "Local")
groups.metric = full_join(groups.metric, like.m.groups, by = "Local")
groups.metric = full_join(groups.metric, connect_web, by = "Local")

colnames(groups.metric) = c("Local","Size", "N.collums", "N.rows", "NODF", "Modularity", "Connectance")

groups.metric= groups.metric %>% 
  select(Local, Size, NODF, Modularity, Connectance) %>% 
  mutate(p_value = p,
         forest_cover = f)

Net_metric <- rbind(Net_metric, groups.metric)

  }
}

net.metric <- Net_metric %>% 
  separate(Local, c(NA, "patch_no")) %>% 
  mutate(patch_no = as.integer(patch_no)) %>% 
  mutate(NODF_st = NODF/100) %>% 
  inner_join(land_pv) %>% 
  select(patch_no, forest_cover, p_value, Land_use, Size, NODF_st, Modularity, Connectance)

write.csv(net.metric, "./Results/net_metric.csv")


# Summarise metrics ####

net.metric_mean <- net.metric %>%
             select(patch_no, forest_cover, p_value, Land_use, Size, NODF_st, Modularity, Connectance) %>% 
  dplyr::group_by(forest_cover, p_value, Land_use) %>% 
  dplyr::summarise(Mean_connectance = mean(Connectance, na.rm = TRUE),
                   sd_con = sd(Connectance),
            Mean_size = mean(Size, na.rm = TRUE),
            sd_size= sd(Size),
            Mean_NODF = mean(NODF_st, na.rm = TRUE),
            sd_NODF = sd(NODF_st),
            Mean_mod = mean(Modularity, na.rm = TRUE),
            sd_mod= sd(Modularity),
            .groups = "drop") %>% 
  mutate(zc_con = Mean_connectance/sd_con, zc_size = Mean_size/sd_size, zc_NODF = Mean_NODF/sd_NODF, zc_mod = Mean_mod/sd_mod)


write.csv(net.metric_mean, "./Results/net.metric_mean.csv")

# Make Regional networks ####

for(f in c(0.1,0.3,0.5,0.7)){
  for(p in c(0.01,0.25,0.5)){
   
    n_filtered_data <- interaction_df %>% 
      filter(forest_cover == f) %>% 
      filter(p_value == p)

n_filtered_data <- n_filtered_data %>% 
  group_by(host_species, parasite_species) %>% 
  dplyr::mutate(Count = n()) %>% 
  ungroup()

reg_web = n_filtered_data %>%  
  select(host_species, parasite_species, Count) %>% 
  group_by(host_species, parasite_species) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  filter(n > 1L) %>% 
  pivot_wider(names_from= parasite_species, values_from= n)

reg_web <- replace(reg_web, is.na(reg_web), 0)
reg_web <- replace(reg_web, is.na(reg_web), 0)
reg_web <- reg_web %>% remove_rownames %>% column_to_rownames(var="host_species")
reg_web <- mutate_all(reg_web, ~ replace(., . > 1, 1))
 var_name <- paste0("p", p,"reg_", f)

 assign(var_name, reg_web)

  }
}

# Make Density plots for each metric for local webs #### 
Net_metric_c <- Net_metric %>% 
  gather("Metric", "value", 2:5) 
str(Net_metric_c)

meta_size <- vcount(net_graph)

meta_nodf <- data.frame(nestednodf(Minc)[["statistic"]])[3,]

meta_mod <- modularity(cluster_fast_greedy(net_graph))

meta_con <- networklevel(Minc, index="connectance")

NODF_plot<- Net_metric_c %>%
  filter(Metric == "NODF") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = as.factor(forest_cover), color = as.factor(forest_cover)), alpha = 0.3) +
  geom_vline(xintercept = meta_nodf, linetype = "dashed", color = "orange", size = 1) +
  labs(x = "NODF - nestedness")+
  facet_wrap(~ as.factor(p_value), scales ="free")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
NODF_plot

Mod_plot<- Net_metric_c %>%
  filter(Metric == "Modularity") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = as.factor(forest_cover), color = as.factor(forest_cover)), alpha = 0.3) +
  geom_vline(xintercept = meta_mod, linetype = "dashed", color = "orange", size = 1) +
  labs(x = "Modularity")+
  facet_wrap(~ as.factor(p_value), scales ="free")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Mod_plot

Con_plot<- Net_metric_c %>%
  filter(Metric == "Connectance") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = as.factor(forest_cover), color = as.factor(forest_cover)), alpha = 0.3) +
  geom_vline(xintercept = meta_con, linetype = "dashed", color = "orange", size = 1) +
  labs(x = "Connectance")+
  facet_wrap(~ as.factor(p_value), scales ="free")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Con_plot

Size_plot<- Net_metric_c %>%
  filter(Local != 2501) %>% 
  filter(Metric == "Size") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = as.factor(forest_cover), color = as.factor(forest_cover)), alpha = 0.3) +
  geom_vline(xintercept = meta_size, linetype = "dashed", color = "orange", size = 1) +
  labs(x = "Size")+
  facet_wrap(~ as.factor(p_value), scales ="free")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Size_plot


plot.list <- list(Con_plot, Size_plot, NODF_plot, Mod_plot)

all.plots<-ggarrange(plotlist = plot.list, ncol = 1,
                   nrow = 4, common.legend = TRUE, legend = "bottom")

png(filename="./Figures/FIG.densityPlots.png", width=25, height=25, units="cm", res=600,  bg = "transparent")

all.plots

dev.off()

# Plot Metaweb using igraph ####
meta_matrix=as.matrix(Minc)
net_graph = graph_from_biadjacency_matrix(Minc, weighted = T)

host_vert<- rownames(Minc)
para_vert<- colnames(Minc)

host_color<- "green"
parasites_color<- "#DFAE82"

# Set parasites vertices shapes
V(net_graph)[names(V(net_graph)) %in% host_vert]$shape <- "circle"

V(net_graph)[names(V(net_graph)) %in% para_vert]$shape <- "square"

# Set species vertices color
V(net_graph)[names(V(net_graph)) %in% para_vert]$color <- parasites_color
V(net_graph)[names(V(net_graph)) %in% host_vert]$color <- host_color

# Find nodes with highest degrees
top_nodes <- names(which(igraph::degree(net_graph) >= 8))

# Set vertex.label for top nodes
V(net_graph)$label <- ifelse(names(V(net_graph)) %in% top_nodes, names(V(net_graph)), NA)  

png(filename="./FIG.metanetwork.png", width=35, height=25, units="cm", res=600,  bg = "transparent")

plot(net_graph,
     vertex.label= V(net_graph)$label,
     vertex.label.cex = 0.8, # Set the label size
     vertex.label.color = "", # Set the label color
     vertex.size=0.8*igraph::degree(net_graph),
     vertex.shape = V(net_graph)$shape,
     edge.width=1.4,
     edge.color="gray50",
     layout= layout_nicely)

dev.off()
