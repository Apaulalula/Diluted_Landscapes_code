###################################
# title: Beta-diversity analysis from Output without rewiring and different strength of co-evolutionary selection for hosts and parasites
# Author: Ana Paula Lula Costa
# Parameters used: alpha: 0.2; eps: 5; mr: 0.5, mc: 0.7, phi: 0.5
# date: 03-2025
###################################

# Load packages ####

pacman::p_load(car, igraph, plotly, data.table, betapart, vegan, reshape2, ggplot2, betalink, knitr, plyr, bipartite, igraph, tidyverse, janitor, stringr, dplyr, hrbrthemes, viridis, ggpubr, ggsci, float, ggpmisc, lmPerm, performance, cowplot, sjPlot)

# Load data ###

betadiv_df<- read.csv(betadiv_df, "./Results/beta_values_random.csv")

# Work on data set ####
str(beta_df)

beta_df_true<- beta_df %>% 
  filter(same_habitat == TRUE)

hist((beta_df_true$ST))

beta_df_true$frag_value <- factor(beta_df_true$frag_value , levels = c("0.99", "0.75", "0.5"))
beta_df_true$Landscape <- factor(beta_df_true$Landscape)
beta_df_true$join_hab <- factor(beta_df_true$join_hab)
beta_df_true$forest_cover <- factor(beta_df_true$forest_cover)
str(beta_df_true)

beta_df_false<- beta_df %>% 
  filter(same_habitat == FALSE)

beta_df_false$frag_value <- factor(beta_df_false$frag_value , levels = c("0.99", "0.75", "0.5"))
beta_df_false$Landscape <- factor(beta_df_false$Landscape)
beta_df_false$join_hab <- factor(beta_df_false$join_hab)
beta_df_false$forest_cover <- factor(beta_df_false$forest_cover)
str(beta_df_false)

# St model within habitat types ####

a<-gam(ST ~ forest_cover*join_hab + frag_value + s(Landscape, bs = "re"), data =beta_df_true)
sum_a <- summary(a)

gam.check(a)

gam_a<- data.frame(sum_a$p.table, random= sum_a$s.table, r2= sum_a$r.sq)

write.csv(gam_a, "./Results/gam_st_with.csv")

tab_model(a,
          show.se = TRUE,
          show.stat = TRUE,
          show.p = TRUE,
          show.ci = FALSE,  # ou TRUE se quiser intervalo de confiança
          transform = NULL, # ou "exp" se for modelo log-link
          dv.labels = "ST index",
          title = "GAM observed results of within habitats interaction dissimilarity based on turnover (Bst)",
          file = "tabela_gam_stw.doc")

# St model between habitat types ####
## Try different models ####
f<-gam(ST ~ forest_cover+ frag_value + join_hab+ s(Landscape, bs = "re"), data =beta_df_false)
performance::r2(f)
summary(f)

f_inter <- gam(ST ~ forest_cover * join_hab + frag_value + s(Landscape, bs = "re"),
               data = beta_df_false)
summary(f_inter)

f_inter2 <- gam(ST ~ forest_cover + frag_value * join_hab + s(Landscape, bs = "re"), data = beta_df_false)

f_full <- gam(ST ~ forest_cover * join_hab + frag_value * join_hab + s(Landscape, bs = "re"), data = beta_df_false)

# Compare models
anova(f, f_inter, f_inter2,f_full, test = "Chisq")

# Run dredge
options(na.action = "na.fail")
dredge_results <- dredge(f)

### model result ####
plot(f_full, pages = 1)# Visualize the random effect
gam.check(f_full)

sum_f <- summary(f_full)

gam_f<- data.frame(sum_f$p.table, random= sum_f$s.table, r2= sum_f$r.sq)

write.csv(gam_f, "./Results/gam_st_bet.csv")

tab_model(f_full,
          show.se = TRUE,
          show.stat = TRUE,
          show.p = TRUE,
          show.ci = FALSE,  # ou TRUE se quiser intervalo de confiança
          transform = NULL, # ou "exp" se for modelo log-link
          dv.labels = "ST index",
          title = "Resultados do modelo GAM",
          file = "tabela_gam.doc")

# ST plots ####
#within habitats
preds_df_a <- as.data.frame(ggpredict(a, terms = c("forest_cover", "join_hab")))
str(preds_df_a)

a_plot<- ggplot(preds_df_a, aes(x = x, y = predicted, color = factor(group))) +
  geom_point(size = 3, alpha= 0.8) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, alpha= 0.5) +
  scale_y_continuous(position = "left", limits = c(0.3, 0.8))+
  scale_color_manual(
  values= c("Core_Core" = '#273B09', "Edge_Edge" = '#7B904B', "Disturbed_Disturbed" = '#FFBF61'),
  labels = c("Core_Core" = "Forest core", "Edge_Edge" = "Forest edge", "Disturbed_Disturbed"= "Disturbed habitat"), name = "Habitat type")+
  theme_bw()+
  theme(title = element_text(size= 14),
        axis.text = element_text(size = 12),
        legend.position= "none",
        legend.background = element_rect(color="black", linetype="solid"),
        legend.direction = "vertical",
        legend.box = "vertical")+
  labs(title= "A - Within habitat dissimilarity",y = "ST index", x = "Proportion of forest cover")
a_plot

# between habitats plots

preds_df_f <- as.data.frame(ggpredict(f_full, terms = c("forest_cover", "frag_value", "join_hab")))

f_plot <- ggplot(preds_df_f, aes(x = x, y = predicted, color = factor(group))) +
  geom_point(size = 3, alpha= 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, alpha= 0.5) +
  scale_y_continuous(position = "right", limits = c(0.3, 0.8))+
  scale_color_manual(values = c("0.99" = '#C2A878', "0.75" = '#6E633D', "0.5" = '#355834'), name = NULL) +
  facet_wrap(~facet, labeller = as_labeller(c(
    "Core_Disturbed" = "Forest core - disturbed",
    "Core_Edge" = "Forest core - edge",
    "Edge_Disturbed" = "Forest edge - disturbed"
  ))) +
  #ylim(0.3,0.8)+
  labs(title = "B - Between habitat dissimilarity",
       x = "Proportion of forest cover",
       y = "ST index") +
  theme_bw() +
  theme(title = element_text(size= 14),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(color = "black", size = 13),
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y =element_blank()
  )
f_plot


## Merge plots ####
a_plot_c <- cowplot::ggdraw() +
  cowplot::draw_plot(a_plot,  0, 0, 1, .94)  # a
a_plot_c

grid <- ggarrange(a_plot_c, f_plot, nrow = 1, ncol = 2, widths = c(1, 3))

grid

png(filename="./Figures/beta_plot_grid.png", width=40, height=12, units="cm", res=600)
plot(grid)
dev.off()
