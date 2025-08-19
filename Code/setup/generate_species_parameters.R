# GENERATE SPECIES PARAMETERS

# this script generates theta parameters for the coevolution model

# load packages
library(dplyr)

# set working directory
setwd("~/HL_coevolution_host_parasite")

# specify network to download
networkName = "metanetwork"

# specify mean theta value for disturbance-adapted and forest specialist hosts
theta_d = 10
theta_f = 20

# import network incidence matrix
Minc = read.table(paste0("Data/",networkName,"/Minc.csv"))

# OPTION 1 - hosts classified as "forest-specialist" or "disturbance-adapted" ----

# specify number of replicas
reps = 10

# number of resources
n_r = nrow(Minc)

# dataframe with host habitat preference and mean theta value
host_types = data.frame(host=1:n_r, 
                        type=read.table(paste0("Data/",networkName,"/host_types.csv"))$V1) %>%
  mutate(theta_mean=ifelse(type=="disturbance-adapted", theta_d, theta_f))

# loop through replicas and draw theta values
for(r in 1:reps){
  
  # set seed to replica number
  set.seed(r)
  
  # sample theta values for hosts
  theta_r = rnorm(n_r, mean=host_types$theta_mean, sd=1)
  
  # parasite theta = mean theta of its hosts
  theta_c = theta_r * Minc
  theta_c[theta_c==0] = NA
  theta_c = as.vector(colMeans(theta_c, na.rm=TRUE))
  
  # write out
  write.table(theta_r, paste0("Data/",networkName,"/theta_r_r",r,".csv"), row.names=FALSE, col.names=FALSE)
  write.table(theta_c, paste0("Data/",networkName,"/theta_c_r",r,".csv"), row.names=FALSE, col.names=FALSE)
}

# OPTION 2 - hosts assigned probability of occurrence in forest ----

# host forest occurrence probability
pf = read.table(paste0("Data/",networkName,"/host_forest_occurrence.csv"))$V1

# calculate host theta as weighted average of theta_f and theta_d
theta_r = pf*theta_f + (1-pf)*theta_d

# parasite theta = mean theta of its hosts
theta_c = theta_r * Minc
theta_c[theta_c==0] = NA
theta_c = as.vector(colMeans(theta_c, na.rm=TRUE))

# write out
write.table(theta_r, paste0("Data/",networkName,"/theta_r.csv"), row.names=FALSE, col.names=FALSE)
write.table(theta_c, paste0("Data/",networkName,"/theta_c.csv"), row.names=FALSE, col.names=FALSE)
