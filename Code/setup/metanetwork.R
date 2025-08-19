# Building Model Inputs

######## Load Packages ######## 
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

######## Load data ########

metalist<-read.csv("C:/Users/anapa/OneDrive/Desktop/Nuvem/Projetos/PDJ-Fiocruz/Dados_Host-parasite-interactions/Dados_Gabi/S5_Host_Parasite_Brazil.csv", sep=";")
Sp_category<-read.csv("C:/Users/anapa/OneDrive/Desktop/Nuvem/Projetos/PDJ-Fiocruz/Dados_Salve/Generated_data/Mammals_category.csv")

############################### Metanetwork data ##################################
str(metalist)
unique(metalist$HostOrder)

meta_parasite = metalist %>%
  filter(HostOrder %in% c("Rodentia","Carnivora", "Didelphimorphia", "Primates", "Cingulata", "Lagomorpha", "Pilosa", "Perissodactyla")) %>% 
  filter(Biome == "atlantic_forest") %>% 
  filter(!ParStatus %in% c("family", "genus")) %>% 
  filter(ParGroup %in% c("Bacteria", "DNA_virus",  "RNA_virus", "Virus")) %>% 
  filter(!ParasiteCurrentName %in% c("not_found",  "not_reported")) 
str(meta_parasite$HostMach_COL_SiBBr)

meta_web = meta_parasite %>% 
  select(HostMach_COL_SiBBr, ParasiteCurrentName) %>% 
  group_by(HostMach_COL_SiBBr, ParasiteCurrentName) %>%
  dplyr::summarise(Occur= n(), .groups = 'drop') %>%
  pivot_wider(names_from= ParasiteCurrentName, values_from= Occur)

meta_web <- meta_web %>% remove_rownames %>% column_to_rownames(var="HostMach_COL_SiBBr")

meta_web <- replace(meta_web, is.na(meta_web), 0)
meta_web <- mutate_all(meta_web, ~ replace(., . > 1, 1))
meta_web=meta_web[,colSums(meta_web)>0]
meta_web=meta_web[rowSums(meta_web)>0,]

Sp_category = Sp_category %>% 
  mutate(Nome_Cient = case_when(Nome_Cient ==  "alouatta_guariba_guariba" ~ "alouatta_guariba",
                                .default = as.character(Nome_Cient)))

host_class <- Sp_category[Sp_category$Nome_Cient %in% rownames(meta_web), ]

sp_notin <- meta_web[!rownames(meta_web) %in% Sp_category$Nome_Cient, ]

meta_web_f <- meta_web[rownames(meta_web) %in% host_class$Nome_Cient, ]
meta_web_f=meta_web_f[,colSums(meta_web_f)>0]
meta_web_f=meta_web_f[rowSums(meta_web_f)>0,]

all(rownames(meta_web_f) == host_class$Nome_Cient)
reorder = match(host_class$Nome_Cient, rownames(meta_web_f))
meta_web_f= meta_web_f[reorder,]

host_class = replace(host_class, is.na(host_class), 0)

host_info = data.frame(species= 1:51, Name = rownames(meta_web_f))

prst_info = data.frame(species= 1:103, Name = colnames(meta_web_f))

meta_pf <- meta_parasite %>% 
  select(ParasiteCurrentName, ZoonoticStatus, ParStatus, ParGroup) %>% 
  group_by(ParasiteCurrentName, ZoonoticStatus, ParStatus, ParGroup) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::rename(Name = ParasiteCurrentName)

prst_info <- prst_info %>% 
  inner_join(meta_pf)

meta_hf <- meta_parasite %>% 
  select(HostMach_COL_SiBBr, HostOrder) %>% 
  group_by(HostMach_COL_SiBBr, HostOrder) %>% 
  dplyr::summarise(Count = n()) %>% 
  dplyr::rename(Name = HostMach_COL_SiBBr)

host_type <- host_class %>% 
  select(Nome_Cient, Category) %>% 
  dplyr::rename(Name = Nome_Cient, type = Category) %>% 
  remove_rownames() %>% 
  inner_join(meta_hf) %>% 
  column_to_rownames(var = "Name")

host_df <- host_class %>% 
  select(Nome_Cient, Forest) %>% 
  dplyr::rename(Name = Nome_Cient, pf = Forest) %>% 
  remove_rownames() %>% 
  inner_join(meta_hf) %>% 
  column_to_rownames(var = "Name")

Minc= as.matrix(meta_web_f)
 write.table(Minc, "./Data/Minc.csv", row.names = FALSE, col.names = FALSE)
 write.csv(host_type, "./Data/host_class.csv")
 write.table(host_df, "./Data/host_forest_occurrence.csv", row.names = FALSE, col.names = FALSE)

 ######################## Measure species degree #####################################
 net.asso=as.matrix(meta_web_f)
 net_graph = graph_from_biadjacency_matrix(net.asso, weighted = T)
 
 Species_degree = data.frame(Degree = degree(net_graph)) %>% 
   rownames_to_column(var = "Name")
 
 prst_info <- prst_info %>% 
   inner_join(Species_degree)
 
 host_info <- host_df %>%
   rownames_to_column(var = "Name") %>% 
   inner_join(Species_degree)
 
 write.csv(host_info, "./Results/host_info.csv")
 write.csv(prst_info, "./Results/prst_info.csv")
 write.csv(Minc, "./Data/M_inc.csv")
 
 ############## Plot metanetwork #############
 
 host_vert<- rownames(net.asso)
 para_vert<- colnames(net.asso)
 
 host_color<- "orange"
 parasites_color<- "black"
 
 # Set parasites vertices shapes
 V(net_graph)[names(V(net_graph)) %in% host_vert]$shape <- "circle"
 
 V(net_graph)[names(V(net_graph)) %in% para_vert]$shape <- "square"
 
 # Set species vertices color
 V(net_graph)[names(V(net_graph)) %in% para_vert]$color <- parasites_color
 V(net_graph)[names(V(net_graph)) %in% host_vert]$color <- host_color
 
 # tkid <- tkplot(net_graph)
 # l <- tkplot.getcoords(tkid)
 # tk_close(tkid, window.close = T)
 
 l<-layout_with_dh(net_graph)
 layout=layout_in_circle(net_graph, order = sorted_vertices)
 
 png(filename="./Figures/FIG.metanetwork2.png", width=35, height=25, units="cm", res=600,  bg = "transparent")
 
 plot(net_graph,vertex.label=NA,
      vertex.size=0.6*igraph::degree(net_graph),
      vertex.shape = V(net_graph)$shape,
      edge.width=1,
      edge.color="gray50",
      layout= layout_nicely)
 
 dev.off()
 

