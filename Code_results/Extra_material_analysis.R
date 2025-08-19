###################################
# title: Results from Output without rewiring and same different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 10 - 2024
###################################

# Load Packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(betareg)
library(mgcv)
library(marginaleffects)
library(MuMIn)
library(car)
library(piecewiseSEM)
library(glmmTMB)
library(lavaan)
library(semPlot)

# Load data

p_res_h_fc <- read.csv("./Results/p_res_h_fc.csv")

################# Organize host data ###########
str(p_res_h_fc)

land_use_dummies <- model.matrix(~ Land_use - 1, data = p_res_h_fc)
p_res_h_fc <- cbind(p_res_h_fc, land_use_dummies)

################# Analyse host global abundance ###########
# See data distribution
hist(p_res_h_fc$Abd_type)
# look at dataframe
str(p_res_h_fc)

# scale  data
p_res_h_fc$Abd_type_sc <- as.numeric(scale(p_res_h_fc$Abd_type))

#see data distribution

hist(p_res_h_fc$Abd_type_sc)

# Make PCA axis for network structure
pca_data <- p_res_h_fc %>%
  select(Mean_size, Mean_connectance, Mean_NODF, Mean_mod) %>%
  scale()  
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
biplot(pca_result)
axis_pca_host <- data.frame(pca_result[["rotation"]])

# Include axis in the dataframe
p_res_h_fc$PCA_axis <- pca_result$x[, 1]
p_res_h_fc$PCA_axis2 <- pca_result$x[, 2]

### Build linear model #####
m_fit <- glmmTMB(Abd_type_sc ~ forest_cover +FragmentationLow + Land_use +PCA_axis + PCA_axis2, data =p_res_h_fc, family= gaussian)

# Analyse model results
summary(m_fit)

# Analyse model residuals
DHARMa::testResiduals(m_fit)

# See best model configuration
options(na.action = "na.fail")
dredge(m_fit, m.lim = c(0,5))

 ### Final model with select variables

f_fit <- glm(Abd_type_sc ~ forest_cover+ FragmentationHigh + Land_useCore + Land_useDisturbed + PCA_axis2, data =p_res_h_fc, family= gaussian)

# Analyse model results
summary(f_fit)

# Analyse inflation factor
vif(f_fit)

#test model residuals
DHARMa::testResiduals(f_fit)


plot_predictions(f_fit, condition = c('PCA_axis2'), vcov = TRUE,
                 type = 'link')

# Build linear model for pca axis
fit <- lmer(PCA_axis ~ forest_cover + Fragmentation + Land_use+ (1|Landscape), data =p_res_h_fc, REML=FALSE)
# Analyse model results
summary(fit)
drop1(fit)

#test fit#test model residuals
DHARMa::testResiduals(fit)

# See best model configuration
options(na.action = "na.fail")
dredge(fit, m.lim = c(0,5))

########### Work on Structural Equation model for host global abundance ########
abd_psem <- psem(
  glm(Abd_type_sc ~ Land_useCore + Land_useDisturbed + PCA_axis, data =p_res_h_fc, family= gaussian),
  glm(PCA_axis ~ forest_cover +FragmentationHigh + Land_useCore+ Land_useDisturbed , data =p_res_h_fc, family= gaussian)
)

# Analyse psem
LLchisq(abd_psem) # good fit
fisherC(abd_psem)

piecewiseSEM::coefs(abd_psem)
a_h_psem<- summary(abd_psem)
a_h_psem
plot(abd_psem)

# Make lavan model for comparison
model <- '
  Abd_type_sc ~ Land_useCore + Land_useDisturbed + PCA_axis
  PCA_axis ~ forest_cover +FragmentationHigh +Land_useCore+ Land_useDisturbed
'

fit_sem <- sem(model, data = p_res_h_fc)
varTable(fit_sem)
summary(fit_sem)

semPaths(fit_sem, what = "std", layout = "tree",edge.label.cex = 2.5, intercepts = FALSE, residuals = FALSE)


#### Host abd plt ####

############# Parasite global abundance plots ####
Abd_plot1 <-p_res_h_fc %>% 
  select(-p_value) %>% 
  ggplot(aes(x= PCA_axis, y= Abd_type_sc)) + 
  geom_point(alpha= 0.1, aes(color= Land_use))+
  geom_smooth(method = lm, se = TRUE, aes(color = "grey50"))+
  #scale_color_manual(values=c("Zoonotic" = '#FF6F59', "Non-zoonotic" = '#254441'), name=NULL) +
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Habitat type")+
  labs( x = "Network structure PC1", y= "Host relative abundance (scaled)")+
  theme_bw()+
  theme(legend.position= "none", #c(0.8,0.15),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Abd_plot1

str(p_res_h_fc)
p_res_h_fc$Land_use <- factor(p_res_h_fc$Land_use , levels = c("Core", "Edge", "Disturbed"))

Abd_plot <- p_res_h_fc %>% 
  select(-p_value) %>% 
  group_by(Land_use) %>%
  dplyr::summarise(mean_value = mean(Abd_type_sc),
                   se_value = sd(Abd_type_sc) / sqrt(n()), .groups = "drop") %>% 
  ggplot(aes(x= as.factor(Land_use), y= mean_value), position = position_dodge(0.3)) + 
  geom_pointrange(aes(ymin = mean_value - se_value, ymax = mean_value + se_value, color = Land_use)) +
  scale_color_manual(values=c("Core" = '#7D891D', "Edge" = '#636b24', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs( x = "Patch land use", y= "Relative Abundance (scaled)")+
  theme_bw()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Abd_plot
