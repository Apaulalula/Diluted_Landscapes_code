# install NLMR package
library(devtools)
devtools::install_github("ropensci/NLMR")

library(NLMR)
library(landscapemetrics)
library(raster)
library(ggplot2)
library(tidyr)
library(dplyr)

# set working directory
setwd("~/HL_coevolution_host_parasite")

# GENERATE LANDSCAPES ----

# function for generating landscapes
# INPUTS:
# n - grid size (discretisation in x and y direction)
# f - forest cover [0,1]
# p_value - proportion of elements randomly selected to form clusters [0.01 ~ random, 0.5 - clustered]

generate_landscapes <- function(n, f, p_value){

  # generate landscape
  l <- NLMR::nlm_randomcluster(nrow=n, ncol=n, resolution=1,
                               p=p_value, ai=c((1-f),f),
                               neighbourhood=4, rescale=TRUE)
    
  # calculate forest cover
  fc <- round(sum(values(l))/length(values(l)), digits=2)

  # initialise iteration count
  it_count <- 1
    
  # regenerate landscape until specified forest cover
  while(fc!=f & it_count<=100){
      
    # regenerate landscape
    l <- NLMR::nlm_randomcluster(nrow=n, ncol=n, resolution=1,
                                 p=p_value, ai=c((1-f),f),
                                 neighbourhood=4, rescale=TRUE)
    # calculate forest cover
    fc <- round(sum(values(l))/length(values(l)), digits=2)

    # update iteration count
    it_count <- it_count + 1
  }
    
  # if not converged on forest cover, print warning
  if(fc!=f) {print("not converged")}
  
  # write out results
  write.table(values(l), paste0("Data/landscapes/forestcover_",f,"_pvalue_",p_value,".csv", sep=""), row.names=FALSE, col.names=FALSE)
}

# run function
generate_landscapes(n=50, f=0.1, p_value=0.01)
generate_landscapes(n=50, f=0.3, p_value=0.01)
generate_landscapes(n=50, f=0.5, p_value=0.01)
generate_landscapes(n=50, f=0.7, p_value=0.01)
generate_landscapes(n=50, f=0.1, p_value=0.25)
generate_landscapes(n=50, f=0.3, p_value=0.25)
generate_landscapes(n=50, f=0.5, p_value=0.25)
generate_landscapes(n=50, f=0.7, p_value=0.25)
generate_landscapes(n=50, f=0.1, p_value=0.5)
generate_landscapes(n=50, f=0.3, p_value=0.5)
generate_landscapes(n=50, f=0.5, p_value=0.5)
generate_landscapes(n=50, f=0.7, p_value=0.5)


# 100% forest landscape
n <- 50
l <- rep(1,n*n)
write.table(l, paste0("Data/landscapes/forestcover_1_pvalue_0.csv", sep=""), row.names=FALSE, col.names=FALSE)


# PLOTS ----

landscapes <- data.frame(forest_cover=double(),
                         p_value=double(),
                         patch_no=integer(),
                         patch_state=integer())
for(f in c(0.1,0.3,0.5)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("Data/landscapes/forestcover_",f,"_pvalue_",p,".csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l))
    landscapes <- rbind(landscapes, l)
  }
}
landscapes <- landscapes %>%
  rename(patch_state="V1") %>%
  left_join(., data.frame(patch_no=1:2500, expand.grid(Y=50:1, X=1:50)))

ggplot(data=landscapes, aes(x=X, y=Y, fill=as.factor(patch_state))) + 
  geom_tile() +
  facet_grid(forest_cover~p_value) +
  coord_fixed(ratio=1) +
  scale_x_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_fill_manual(values=c("#FFBF61", "#7D891D"), name=NULL)


landscapes_edge <- data.frame(forest_cover=double(),
                         p_value=double(),
                         patch_no=integer(),
                         patch_state=integer())
for(f in c(0.1,0.3,0.5)){
  for(p in c(0.01,0.25,0.5)){
    l <- read.table(paste0("Data/landscapes/forestcover_",f,"_pvalue_",p,"_edge.csv"), quote="\"", comment.char="")
    l <- data.frame(forest_cover=f,
                    p_value=p,
                    patch_no=1:2500,
                    patch_state=as.vector(l))
    landscapes_edge <- rbind(landscapes_edge, l)
  }
}
landscapes_edge <- landscapes_edge %>%
  rename(patch_state="V1") %>%
  left_join(., data.frame(patch_no=1:2500, expand.grid(Y=50:1, X=1:50)))

ggplot(data=landscapes_edge, aes(x=X, y=Y, fill=as.factor(patch_state))) + 
  geom_tile() +
  facet_grid(forest_cover~p_value) +
  coord_fixed(ratio=1) +
  scale_x_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(1-0.5,50+0.5), expand=c(0,0)) +
  labs(fill=NULL) +
  scale_fill_brewer(palette="BrBG", name=NULL)

