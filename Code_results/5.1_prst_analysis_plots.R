###################################
# title: Results from Output without rewiring and same different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 09 - 2024
###################################

# Load Packages ####
pacman::p_load(tidyverse, ggplot2, ggpubr, rstatix, betareg, mgcv, lme4, marginaleffects, semPlot,piecewiseSEM,lavaanPlot,semPlot, GGally, glmmTMB, MuMIn, car, cowplot, ggeffects, ggpredict, performance)

# Extinction plot ####
Ex_plot <- plot_predictions(ext, condition = c('forest_cover', 'frag_value'),vcov = TRUE,
                            type = 'response')+ 
  scale_color_manual(values=c("0.99" = '#C2A878', "0.75" = '#6E633D', "0.5"= '#355834'), name=NULL) +
  #scale_color_manual(values=c("0.5" = '#7D891D', "0.25" = '#636b24', "0.01" = '#FFBF61'), name="Patch land use")+
  labs( x = "Proportion of forest cover", y= "Parasite extinction rate")+
  theme_bw()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
Ex_plot

# Network structure plots ####

PC1_plot2 <- plot_predictions(gam_fit, condition = c('forest_cover', 'Land_use'),vcov = TRUE,
                              type = 'response')+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  ylim(-3, 4)+
  labs( x = "Proportion of forest cover", y= "Network structure PC1")+
  theme_bw()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")

PC1_plot1<- plot_predictions(gam_fit, condition = c('forest_cover', 'frag_value'),vcov = TRUE,
                             type = 'link')+ 
  scale_color_manual(values=c("0.99" = '#C2A878', "0.75" = '#6E633D', "0.5"= '#355834'), name=NULL) +
  ylim(-2, 4)+
  labs( x = "Proportion of forest cover", y= "Network structure PC1")+
  theme_bw()+
  theme(legend.position= "none",#c(0.8,0.9),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")


# Combine and save plot ####
final_plot<- plot_grid(PC1_plot2, PC1_plot1, Ex_plot, ncol = 3, labels = c("A", "B", "C"))

# Display the final plot
print(final_plot)

png(filename="./Figures/PC_ns.png", width=30, height=10, units="cm", res=600)
plot(final_plot)
dev.off()

# Outcomes Plots #######################################
## Relative degree plots #####
### Plot relation between relative degree,zoonotic status and extinction rate  ####
dg_plot1 <-df_outliers %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Ext_rate, y= Mean_rel_sc, color= ZoonoticStatus)) + 
  geom_point(alpha= 0.1)+
  geom_smooth(method = lm, se = TRUE)+
  scale_color_manual(values=c("Zoonotic" = '#FF6F59', "Non-zoonotic" = '#254441'), name=NULL) +
  labs( x = "Parasite extinction rate", y= "Relative degree (scaled)")+
  theme_bw()+
  theme(legend.position="none", #c(0.8,0.15),
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal")
dg_plot1

png(filename="./Figures/DG_plot1.png", width=15, height=7, units="cm", res=600)
plot(dg_plot1)
dev.off()

### Plot relation between relative degree, land use and host abundance  ####
dg_plot <-df_outliers %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Host_abd, y= Mean_rel_sc)) + 
  geom_point(alpha= 0.5, aes(color= as.factor(Land_use)))+
  geom_smooth(method = lm, se = TRUE, color= "grey50")+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs( x = "Host abundance", y= "Relative degree (scaled)")+
  theme_bw()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank())
dg_plot

png(filename="./Figures/DG_plot.png", width=15, height=7, units="cm", res=600)
plot(dg_plot)
dev.off()

## Trait matching plots ####
### Plot relation between trait matching, forest cover and network structure ####
TM_plot1 <-df_outliers %>% 
  select(-p_value) %>% 
  ggplot(aes(x= PCA_axis1, y= Amatch_sc)) + 
  geom_point(alpha= 0.15, aes(color= as.factor(Land_use)))+
  geom_smooth(method = lm, se = TRUE, color= "grey50")+
  #scale_color_continuous_sequential(palette = "Greens", h2 = 40)+
 # scale_color_manual(values=c("0.1" = '#CFE0BC', "0.3" = '#7FA653',"0.5" = '#63783D',"0.7" = '#234D21'), name=NULL)+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs( x = "Network structure PC1", y= "Trait matching (scaled)")+
  theme_bw()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank())
TM_plot1

png(filename="./Figures/tm_plot1.png", width=10, height=10, units="cm", res=600)
plot(TM_plot1)
dev.off()

### Plot relation between trait matching, habitat type and host abundance ####
TM_plot2 <-df_outliers %>% 
  select(-p_value) %>% 
  ggplot(aes(x= Host_abd, y= Amatch_sc)) + 
  geom_point(alpha= 0.15, aes(color= as.factor(Land_use)))+
  geom_smooth(method = "lm", se = TRUE, aes(color= Land_use))+
  #scale_color_continuous_sequential(palette = "Greens", h2 = 40)+
  scale_color_manual(values=c("Core" = '#273B09', "Edge" = '#7B904B', "Disturbed" = '#FFBF61'), name="Patch land use")+
  labs( x = "Host relative abundance", y= "Trait matching (scaled)")+
  theme_bw()+
  theme(legend.position="none",
        #legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_blank())
TM_plot2

png(filename="./Figures/tm_plot2.png", width=11, height=10, units="cm", res=600)
plot(TM_plot2)
dev.off()


