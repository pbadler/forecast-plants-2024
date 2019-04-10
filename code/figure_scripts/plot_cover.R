rm(list = ls())

library(tidyverse)
library(egg)
library(lme4)

# call from removal_analysis_wrapper.r

load('code/figure_scripts/my_plotting_theme.Rdata')

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

cover <- readRDS('data/temp_data/all_cover.RDS')

cover$spp <- factor( cover$spp, labels = c('Artemisia', 'Hesperostipa', 'Poa', 'Pseudoroegneria'))

#calculate treatment means by year
cover_per_treatment <- 
  cover %>% 
  group_by( Treatment, year, spp ) %>% 
  summarise( cover = mean(cover, na.rm = T), n = n())

cover_per_treatment %>% 
  ggplot(aes( x = year, y = cover , color = Treatment)) + 
  geom_line() + 
  facet_wrap(~spp )

# calculate year-to-year log changes
cover_change <- 
  cover %>%
  arrange(Treatment, QuadName, spp, year) %>% 
  group_by( Treatment, QuadName, spp ) %>% 
  mutate( log_change = log( cover / lag(cover, 1))) %>% 
  mutate( log_change = ifelse( !is.finite(log_change), NA, log_change) )

# calculate deviations from pretreatment year
cover_diff <- 
  cover %>% 
  filter( year > 2005 ) %>% 
  group_by( Treatment, QuadName, spp ) %>% 
  mutate( cover_diff = cover - cover[year == 2011] ) %>%
  arrange( Treatment, QuadName, spp, year)

mean_cover_diff <- 
  cover_diff %>% 
  group_by( Treatment, spp, year) %>% 
  summarise(mean_diff = mean(cover_diff,na.rm = T))

mean_log_cover_diff <- 
  cover_diff %>% 
  mutate( log_diff = log(cover/cover[year == 2011 ])) %>%
  mutate( log_diff = ifelse(!is.finite(log_diff), NA, log_diff)) %>% 
  group_by( Treatment, spp, year) %>% 
  summarise( mean_log_diff = mean(log_diff, na.rm = T))

# statistical tests ####################################################
dARTR <- 
  cover_change %>% 
  filter( spp == 'Artemisia', year > 2010) 

mARTR <- lmer(log_change ~ Treatment + (1|year), data = dARTR)
summary( mARTR )

dHECO <- 
  cover_change %>% 
  filter( spp == 'Hesperostipa', year > 2010) 

mHECO <- lmer(log_change ~ Treatment + (1|year), data = dHECO)
summary( mHECO )

dPOSE <- 
  cover_change %>% 
  filter( spp == 'Poa', year > 2010) 

mPOSE <- lmer(log_change ~ Treatment + (1|year), data = dPOSE)
summary( mPOSE )

dPSSP <- 
  cover_change %>% 
  filter( spp == 'Pseudoroegneria', year > 2010) 

mPSSP <- lmer(log_change ~ Treatment + (1|year), data = dPSSP)
summary( mPSSP )

statsOutput <- paste0( "output/", species, "_lmer_pgr_stats_table.text")
output <- list(mARTR,mHECO,mPOSE,mPSSP)

for(i in 1:length(species)){ 
  texreg::texreg(output[[i]], ci.force=TRUE,caption="Cover change models",
                 caption.above=TRUE,file=statsOutput[i])
} 


# figures ########################################################################

p1 <- 
  cover_per_treatment %>% 
  filter( year > 2005) %>%
  ggplot( aes( x = year, y = cover, color = Treatment)) + 
  geom_line() + 
  facet_wrap( ~ spp, nrow = 1 ) + 
  ylab( 'Mean cover (%)' ) + 
  scale_color_manual(values = my_colors[2:4]) + 
  scale_x_continuous(breaks = c(2008, 2012, 2016)) + 
  my_theme + 
  theme( strip.text = element_text(face = 'italic'), 
         axis.title.x = element_blank())

p1.ARTR <- 
  p1 %+% 
  (p1$data %>% filter( spp == 'Artemisia')) +  
  ylim(c(0, 20) ) + 
  theme(legend.position = c(0.4, 0.1), 
        legend.key.width = unit(2, 'line'))

p1.HECO <- 
  p1 %+% 
  (p1$data %>% filter( spp == 'Hesperostipa')) +  
  ylim(c(0, 3.5)) + 
  guides(color = F) + 
  theme(axis.title.y = element_blank())

p1.POSE <- 
  p1 %+% 
  (p1$data %>% filter( spp == 'Poa')) +  
  ylim(c(0, 3.5)) + 
  guides(color = F) + 
  theme(axis.title.y = element_blank())

p1.PSSP <- 
  p1 %+% 
  (p1$data %>% filter( spp == 'Pseudoroegneria')) +
  ylim(c(0, 3.5)) + 
  guides(color = F) + 
  theme(axis.title.y = element_blank())


#1. Average cover treatment and year

p1.all_cover <- 
  ggarrange(p1.ARTR + theme(legend.position = c(0.4, 0.2)), 
          p1.HECO, 
          p1.POSE, 
          p1.PSSP, nrow = 1, 
          padding = unit(1,'line'), 
          bottom = 'Year')

ggsave('figures/treatment_trends_cover_new.png', 
       p1.all_cover, 
       height = 4, 
       width = 8, 
       units = 'in', dpi = 400)

#2. Log cover change by treatment and year 

mean_cov_change <- 
  cover_change %>% 
  group_by( Treatment, year, spp ) %>% 
  summarise( mean_log_change = mean(log_change, na.rm = T)) %>% 
  filter( year > 2006)  

log_change_plot <- 
  mean_cov_change %>% 
  ggplot( aes( x = year, y = mean_log_change, color = Treatment)) + 
    geom_line() +
    geom_point() + 
    geom_hline(aes( yintercept = 0 ), linetype = 2) + 
    facet_wrap( ~ spp, nrow = 1 ) + 
    ylab( 'Avg Log Change in Cover' ) + 
    scale_color_manual(values = my_colors[2:4]) + 
    scale_x_continuous(breaks = c(2008, 2012, 2016)) + 
    my_theme + 
    theme( strip.text = element_text(face = 'italic'), 
           axis.title.x = element_blank(), 
           legend.position = c(0.1, 0.2)) + 
    ylim(-2, 2)

ggsave('figures/treatment_trends_logChange_new.png', 
       log_change_plot,
       height = 4, 
       width = 8, 
       units = 'in', dpi = 400)

# 
#2. Average cover deviation (w.r.t. pretreatment year)

cover_dev_plot <- 
  mean_cover_diff %>% 
  ungroup() %>% 
  ggplot( aes( x = year, y = mean_diff, color = Treatment ))  + 
  geom_line() +
  geom_point() + 
  geom_hline(aes( yintercept = 0 ), linetype = 2) + 
  facet_wrap( ~ spp, nrow = 1 ) + 
  ylab( 'Percentage point deviation from 2011 cover' ) + 
  scale_color_manual(values = my_colors[2:4]) + 
  scale_x_continuous(breaks = c(2008, 2012, 2016)) + 
  my_theme + 
  theme( strip.text = element_text(face = 'italic'), 
         axis.title.x = element_blank(), 
         legend.position = c(0.92, 0.8))


ggsave(filename = 'figures/cover_deviation_new.png',
        cover_dev_plot, 
        height = 4, 
        width = 8, 
        units = 'in', dpi = 400)

# 3. Change in cover bar chart 

dat <- 
  cover_change %>% 
  ungroup() %>% 
  group_by( QuadName, spp, Treatment ) %>% 
  summarise( lc = log( cover[year == 2016]/cover[year == 2011] ) ) %>% 
  mutate( lc = ifelse(!is.finite(lc), NA, lc))

lc.ARTR <- lm(dat =  subset( dat, spp == 'Artemisia'), lc ~ Treatment )
lc.HECO <- lm(dat =  subset( dat, spp == 'Hesperostipa'), lc ~ Treatment )
lc.POSE <- lm(dat =  subset( dat, spp == 'Poa'), lc ~ Treatment )
lc.PSSP <- lm(dat =  subset( dat, spp == 'Pseudoroegneria'), lc ~ Treatment )

summary(lc.ARTR)
summary(lc.HECO)
summary(lc.POSE)
summary(lc.PSSP)

library(xtable)
xt1 <- xtable(lc.ARTR, 
              caption = 'Treatment effects on log cover change for \textit{Artemisia} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changeARTR')
xt2 <- xtable(lc.HECO, caption = 'Treatment effects on log cover change for \textit{Hesperostipa} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changeHECO')
xt3 <- xtable(lc.POSE, caption = 'Treatment effects on log cover change for \textit{Poa} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changePOSE')
xt4 <- xtable(lc.PSSP, caption = 'Treatment effects on log cover change for \textit{Pseudoroegneria} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changePSSP')

print(xt1, 'manuscript/ARTR_cover_change.tex', type = 'latex', caption.placement ="top")
print(xt2, 'manuscript/HECO_cover_change.tex', type = 'latex', caption.placement ="top")
print(xt3, 'manuscript/POSE_cover_change.tex', type = 'latex', caption.placement ="top")
print(xt4, 'manuscript/PSSP_cover_change.tex', type = 'latex', caption.placement ="top")

artr <- data.frame( summary(lc.ARTR)$coefficients)
heco <- data.frame(summary(lc.HECO)$coefficients)
pose <- data.frame(summary(lc.POSE)$coefficients)
pssp <- data.frame(summary(lc.PSSP)$coefficients)

table_df <- rbind( 
  c('ARTR', '', '', '' ), 
  colnames(artr), 
  round(artr, 2), 
  c('HECO', '', '', ''), 
  colnames( heco), 
  round( heco, 2), 
  c('POSE', '', '', ''), 
  colnames( pose), 
  round( pose, 2 ) , 
  c('PSSP', '', '', ''), 
  colnames( pssp), 
  round( pssp, 2)) 

overall_xtable <- xtable( table_df , 
                          caption = 'Treatment effects on log cover change for each species from 2011 to 2016. Intercept gives control effects.', 
                          label = 'table:coverChange')

print(overall_xtable, 'manuscript/overall_cover_change.tex', type = 'latex', caption.placement = 'top',  digits=c(2,2,2,2,2))

# plot average cover chag

start_to_finish_plot <- 
  dat %>% 
  ggplot(aes( x = Treatment, y = lc, color = Treatment )) + 
  geom_boxplot(show.legend = F) + 
  geom_point()  + 
  geom_hline(aes(yintercept = 0 ), linetype = 2) + 
  ylab( 'Log change in cover 2011 to 2016')  + 
  xlab ( '' ) + 
  facet_grid( . ~ spp) + 
  scale_color_manual( values = my_colors[ 2:4]) + 
  my_theme + 
  theme(strip.text = element_text(face = 'italic'))


ggsave('figures/start_to_finish_cover_change_new.png', 
       start_to_finish_plot,  
       height = 4, 
       width = 8, 
       units = 'in', dpi = 400)
       



