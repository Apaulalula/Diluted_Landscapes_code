###################################
# title: Results from Output without rewiring and same different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 09 - 2024
###################################

# Load Packages
pacman::p_load(tidyverse, ggplot2, ggpubr, rstatix, betareg, mgcv, lme4, marginaleffects, semPlot,piecewiseSEM,lavaanPlot,semPlot, GGally, glmmTMB, MuMIn, car, cowplot, ggeffects, ggpredict, performance)

# Load data #########################

prev_HT_df<- read.csv("./Results/prev_HT_df.csv")
outcomes_model<- read.csv("./Results/outcomes_model.csv")

# Work on data set ################################### 

# Verify unique levels of categorical variables
outcomes_model$parasite_species <- factor(outcomes_model$parasite_species)
outcomes_model$Landscape <- factor(outcomes_model$Landscape)

levels(outcomes_model$parasite_species)

## Other outcomes data frame ####
str(outcomes_model)

outcomes_model$ZoonoticStatus_dummy <- ifelse(outcomes_model$ZoonoticStatus == "Zoonotic", 1, 0)

land_use_dummies <- model.matrix(~ Land_use - 1, data = outcomes_model)
outcomes_model <- cbind(outcomes_model, land_use_dummies)

group_dummies <- model.matrix(~ ParGroup - 1, data = outcomes_model)
outcomes_model <- cbind(outcomes_model, group_dummies)

# Model analysis ####
## Fit and analyse extinction model ####
ext = glm(Ext_rate ~as.factor(forest_cover)+as.factor(frag_value), data = extinction)
summary(ext)
DHARMa::testResiduals(ext)
vif(ext)

tab_model(ext,
          show.se = TRUE,
          show.stat = TRUE,
          show.p = TRUE,
          show.ci = FALSE,  # ou TRUE se quiser intervalo de confiança
          transform = NULL, # ou "exp" se for modelo log-link
          dv.labels = "Parasite extinction rate",
          title = "GLM observed results for parasite extinction rate relationship with landscape factors",
          file = "tabela_gam_ext.doc")

## Fit and analyse network structure models ####
# Fit a GAM
gam_fit <- gam(PCA_axis1 ~ as.factor(forest_cover)+ as.factor(frag_value) + Land_use + s(Landscape, bs = "re"), 
               data = outcomes_model, method = "REML")

# Run dredge
options(na.action = "na.fail")
dredge(gam_fit)

plot(gam_fit, pages = 1)# Visualize the random effect
gam.check(gam_fit)# Check model diagnostics

# analyse results
gam_pc1<-summary(gam_fit)
r_squared <- r2(gam_fit)

tab_model(gam_fit,
          show.se = TRUE,
          show.stat = TRUE,
          show.p = TRUE,
          show.ci = FALSE,  # ou TRUE se quiser intervalo de confiança
          transform = NULL, # ou "exp" se for modelo log-link
          dv.labels = "Network structure - PC1",
          title = "GAM observed results for network structure relationship with landscape factors",
          file = "tabela_gam_pc1.doc")

## Relative degree and trait matching models ###########
### Analyse outcomes distribution ###########
hist(outcomes_model$Mean_rel_sc) # degree difference 
hist(outcomes_model$Amatch_sc) # species trait matching
hist(outcomes_model$Ext_rate_sc) # parasite global extinction rate
hist(outcomes_model$Host_abd) #host mean abundance

### Construct global models  ###########
str(outcomes_model)

##### Relative Degree ####
m_fit <- lmer(Mean_rel_sc ~ forest_cover*frag_value +Land_useDisturbed*ZoonoticStatus_dummy+ Land_useCore + Host_abd + ParGroupBacteria + ParGroupDNA_virus+PCA_axis1 + Ext_rate_sc+ (1|parasite_species), data =outcomes_model)

##### Species trait matching ####
am_fit <- lmer(Amatch_sc ~ forest_cover*frag_value +Land_useDisturbed*ZoonoticStatus_dummy+ Land_useCore + ParGroupBacteria + ParGroupDNA_virus+ Host_abd+ PCA_axis1 + Ext_rate_sc+ (1|parasite_species), data =outcomes_model)

### Analyse global models ####
DHARMa::testResiduals(m_fit)# relative degree 
DHARMa::testResiduals(am_fit)# species trait matching

#### Remove outlier's for trait matching model ####
cooksD <- cooks.distance(am_fit) # Calculate Cook's Distance for the GLM

plot(cooksD, type = "h", main = "Cook's Distance", ylab = "Cook's D") # Plot Cook's Distance

abline(h =30 / nrow(outcomes_model), col = "red")  # Rule of thumb threshold

influential <- which(cooksD > 30 / nrow(outcomes_model)) # Identify influential points

df_outliers <- outcomes_model[-influential, ] # Remove influential points

### Global model without outliers ####
##### Relative Degree ####
m_fit <- lmer(Mean_rel_sc ~ forest_cover*frag_value +Land_useDisturbed*ZoonoticStatus_dummy+ Land_useCore + Host_abd +ParGroupBacteria + ParGroupDNA_virus+PCA_axis1 + Ext_rate_sc+ (1|parasite_species), data =df_outliers)

##### Species trait matching ####
am_fit <- lmer(Amatch_sc ~ forest_cover*frag_value +Land_useDisturbed*ZoonoticStatus_dummy+ Land_useCore + ParGroupBacteria + ParGroupDNA_virus +PCA_axis1 + Ext_rate_sc+ (1|parasite_species), data =df_outliers)

##### Analyse models fit ####
DHARMa::testResiduals(m_fit)# relative degree  
DHARMa::testResiduals(am_fit)# species trait matching

##### Select best model for each outcome ####
options(na.action = "na.fail")

# Relative degree 
dredge(m_fit, m.lim = c(0,5))
drop1(m_fit, test = "Chisq")

# species trait matching
dredge(am_fit, m.lim = c(0,5))
drop1(am_fit, test = "Chisq")

### Final model configuration for each outcome ####
#### Relative degree ####
m_fit <- lmer(Mean_rel_sc ~ Land_useDisturbed + Land_useCore + Host_abd+ ZoonoticStatus_dummy + Ext_rate_sc + (1|parasite_species),REML=FALSE, data =df_outliers)

#### Species trait matching ####
am_fit <- lmer(Amatch_sc ~ forest_cover + Land_useCore+ Land_useDisturbed +PCA_axis1 +Host_abd + (1|parasite_species), REML=FALSE, data =df_outliers)

### Analyse final models ####
#### Relative degree ####
summary(m_fit)
vif(m_fit)
drop1(m_fit)

#### species trait matching ####
summary(am_fit)
vif(am_fit)
drop1(am_fit)

### Structural Equation model for the different outcomes #######
str(df_outliers)

#### Relative degree ####
m_fit <- lmer(Mean_rel_sc ~ Land_useDisturbed + Land_useCore + Host_abd+ ZoonoticStatus_dummy + Ext_rate_sc + (1|parasite_species),REML=FALSE, data =df_outliers)

df_psem <- psem(
  lmer(Mean_rel_sc ~ Land_useDisturbed + Land_useCore + Host_abd+ ZoonoticStatus_dummy + Ext_rate_sc + (1|parasite_species), data =df_outliers),
  lmer(Host_abd ~ Ext_rate_sc + forest_cover + Land_useCore +Land_useDisturbed + (1|parasite_species), data = df_outliers),
  glm(Ext_rate_sc ~forest_cover+frag_value, data = df_outliers)
)

#### PSEM results summary ####
d_psem<-summary(df_psem)
d_psem$AIC
d_psem$Cstat
d_psem$coefficients
d_psem$R2

df_d_psem<- data.frame(d_psem$coefficients, Aic= d_psem$AIC, d_psem$Cstat)
df_d_r2 <- data.frame(r2=d_psem$R2)

write.csv(df_d_psem, "./Results/df_d_psem.csv")
write.csv(df_d_r2, "./Results/d_r2_psem.csv")

### Species trait matching ####
am_fit <- lmer(Amatch_sc ~ forest_cover + Land_useCore+ Land_useDisturbed +PCA_axis1 +Host_abd + (1|parasite_species),REML=FALSE, data =df_outliers)

tm_psem <- psem(
  lmer(Amatch_sc ~  forest_cover + Land_useCore +Land_useDisturbed + Host_abd + PCA_axis1 + (1|parasite_species), data =df_outliers),
  lmer(Host_abd ~ PCA_axis1 + Land_useCore +Land_useDisturbed +frag_value + (1|parasite_species), data = df_outliers),
  lmer(PCA_axis1 ~ forest_cover + frag_value + Land_useCore + Land_useDisturbed+ (1|Landscape), data = df_outliers)
)

#### PSEM results summary ####
t_psem<- summary(tm_psem)
t_psem$AIC
t_psem$Cstat
t_psem$coefficients
t_psem$R2

tm_df_psem<- data.frame(t_psem$coefficients, Aic= t_psem$AIC, t_psem$Cstat)
tm_r2 <- data.frame(r2=t_psem$R2)

write.csv(tm_df_psem, "./Results/tm_df_psem.csv")
write.csv(tm_r2, "./Results/tm_r2_psem.csv")
