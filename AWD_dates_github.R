#Code from Abrahms et al. 2022 "Long-term, climate driven phenological shift in a tropical large carnivore"

#######################
# Import libraries ####
#######################

library(tidyverse)
library(lme4) #for linear mixed models
library(gamm4) #for GAMMS
library(MuMIn) #for AIC model selection

#################
# Load Data ####
#################
  
dates <- read.csv("AWD_dates_github.csv")
  
#############################################
# Have birthing dates changed over time? ####
#############################################
  
lm <- lmer(WhelpYday~Year + (1|Pack) + (1|Year), data=dates)
 
################################################################################
# Have birthing dates changed with respect to the annual temperature nadir? ####
################################################################################

lm <- lmer(WhelpDaysBeforeNadir~Year + (1|Pack), data=dates)

#################################################################################
# Have temperatures on birthing dates and denning periods changed over time? ####
#################################################################################

lm <- lmer(max_whelping_temp~Year + (1|Pack) + (1|Year), data=dates)

lm <- lmer(mean_max_denning_temp~Year + (1|Pack) + (1|Year), data=dates)

######################################
# What is driving later breeding? ####
######################################

fit1 <- gamm4(WhelpYday ~ Year, random = ~(1|Pack) , data=dates)
fit2 <- gamm4(WhelpYday ~ Year, random = ~(1|Pack) + (1|Year), data=dates)
fit3 <- gamm4(WhelpYday ~ s(yearly_temp_max) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit4 <- gamm4(WhelpYday ~ s(yearly_temp_mean) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit5 <- gamm4(WhelpYday ~ s(Sep_temp_mean) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit6 <- gamm4(WhelpYday ~ s(March_temp_mean) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit7 <- gamm4(WhelpYday ~ s(Annual_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit8 <- gamm4(WhelpYday ~ s(prev_Sep_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit9 <- gamm4(WhelpYday ~ s(March_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit10 <- gamm4(WhelpYday ~ s(May_rain) + Year , random = ~(1|Pack) + (1|Year), data=dates)
fit11 <- gamm4(WhelpYday ~ s(Sep_temp_max) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit12 <- gamm4(WhelpYday ~ s(March_temp_max) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit13 <- gamm4(WhelpYday ~ s(May_temp_max) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit14 <- gamm4(WhelpYday ~ s(Sep_temp_mean) + s(prev_Sep_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit15 <- gamm4(WhelpYday ~ s(May_temp_mean) + s(May_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit16 <- gamm4(WhelpYday ~ s(March_temp_mean) + s(March_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit17 <- gamm4(WhelpYday ~ s(yearly_temp_max) + s(Annual_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit18 <- gamm4(WhelpYday ~ s(yearly_temp_mean) + s(Annual_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit19 <- gamm4(WhelpYday ~ s(prev_yearly_temp_max) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit20 <- gamm4(WhelpYday ~ s(prev_yearly_temp_mean) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit21 <- gamm4(WhelpYday ~ s(prev_yearly_temp_max) + s(prev_Annual_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)
fit22 <- gamm4(WhelpYday ~ s(prev_yearly_temp_mean) + s(prev_Annual_rain) + Year, random = ~(1|Pack) + (1|Year), data=dates)

table <- data.frame(formula = c(Reduce(paste, deparse(fit1$gam$formula)),Reduce(paste, deparse(fit2$gam$formula)),Reduce(paste, deparse(fit3$gam$formula)),Reduce(paste, deparse(fit4$gam$formula)),Reduce(paste, deparse(fit5$gam$formula)),
                                Reduce(paste, deparse(fit6$gam$formula)),Reduce(paste, deparse(fit7$gam$formula)),Reduce(paste, deparse(fit8$gam$formula)),Reduce(paste, deparse(fit9$gam$formula)),Reduce(paste, deparse(fit10$gam$formula)),
                                Reduce(paste, deparse(fit11$gam$formula)),Reduce(paste, deparse(fit12$gam$formula)),Reduce(paste, deparse(fit13$gam$formula)),Reduce(paste, deparse(fit14$gam$formula)),Reduce(paste, deparse(fit15$gam$formula)),
                                Reduce(paste, deparse(fit16$gam$formula)),Reduce(paste, deparse(fit17$gam$formula)),Reduce(paste, deparse(fit18$gam$formula)),Reduce(paste, deparse(fit19$gam$formula)),Reduce(paste, deparse(fit20$gam$formula)),
                                Reduce(paste, deparse(fit21$gam$formula)),Reduce(paste, deparse(fit22$gam$formula))), 
                    AICc=c(AICc(fit1$mer), AICc(fit2$mer), AICc(fit3$mer), AICc(fit4$mer), AICc(fit5$mer), AICc(fit6$mer),AICc(fit7$mer),  AICc(fit8$mer),AICc(fit9$mer),AICc(fit10$mer),
                           AICc(fit11$mer),    AICc(fit12$mer),   AICc(fit13$mer),AICc(fit14$mer), AICc(fit15$mer),AICc(fit16$mer), AICc(fit17$mer), AICc(fit18$mer), AICc(fit19$mer), AICc(fit20$mer), 
                           AICc(fit21$mer), AICc(fit22$mer)))
(table <- table %>% arrange(AICc) %>% mutate(AICc=round(AICc,2), deltaAIC = round(AICc-min(AICc),2)))

############################################
# Does temperature affect litter size?  ####
############################################

fit1 <- gamm4(littersize ~ s(mean_max_conception_temp) + s(mean_max_denning_temp) +s(max_whelping_temp) + (Year) + (PackSize), random = ~(1|Pack) + (1|Year), data=dates, family="poisson")
fit2 <- gamm4(littersize ~ s(mean_conception_temp) + s(mean_denning_temp) +s(max_whelping_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit3 <- gamm4(littersize ~ s(mean_max_denning_temp) +s(max_whelping_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit4 <- gamm4(littersize ~ s(mean_max_conception_temp) + s(mean_max_denning_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit5 <- gamm4(littersize ~ s(mean_max_denning_temp) +s(max_whelping_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit6 <- gamm4(littersize ~ s(mean_max_denning_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit7 <- gamm4(littersize ~ s(mean_max_conception_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit8 <- gamm4(littersize ~ s(max_whelping_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit9 <- gamm4(littersize ~ s(Rain_DenningMonths) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit10 <- gamm4(littersize ~ s(Rain_ConceptionMonth) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit11 <- gamm4(littersize ~ s(mean_max_denning_temp) + s(Rain_DenningMonths) + (Year) + (PackSize), random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")
fit12 <- gamm4(littersize ~ s(mean_denning_temp) + (Year) + (PackSize) , random = ~ (1|Pack)+ (1|Year), data=dates, family="poisson")


table <- data.frame(formula = c(Reduce(paste, deparse(fit1$gam$formula)),Reduce(paste, deparse(fit2$gam$formula)),Reduce(paste, deparse(fit3$gam$formula)),Reduce(paste, deparse(fit4$gam$formula)),Reduce(paste, deparse(fit5$gam$formula)),
                                Reduce(paste, deparse(fit6$gam$formula)),Reduce(paste, deparse(fit7$gam$formula)),Reduce(paste, deparse(fit8$gam$formula)),Reduce(paste, deparse(fit9$gam$formula)),Reduce(paste, deparse(fit10$gam$formula)),
                                Reduce(paste, deparse(fit11$gam$formula)),Reduce(paste, deparse(fit12$gam$formula))), 
                    AICc=c(AICc(fit1$mer), AICc(fit2$mer), AICc(fit3$mer), AICc(fit4$mer), AICc(fit5$mer), AICc(fit6$mer),AICc(fit7$mer),  AICc(fit8$mer),AICc(fit9$mer),AICc(fit10$mer),
                           AICc(fit11$mer),    AICc(fit12$mer)))
(table <- table %>% arrange(AICc) %>% mutate(AICc=round(AICc,2), deltaAIC = round(AICc-min(AICc),2)))


#####################################################
# Does temperature affect yearling recruitment?  ####
#####################################################

fit1 <- gamm4(yrlsurv ~ (littersize) + s(mean_max_postdenning_temp) + (Year) + (PackSize), random=~ (1|Pack) + (1|Year), data=dates, family="poisson")
fit2 <- gamm4(yrlsurv ~ (littersize) + s(Rain_PostDenningMonths) + (Year) + (PackSize), random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit3 <- gamm4(yrlsurv ~ (littersize) + s(mean_max_postdenning_temp) + s(Rain_PostDenningMonths) + (Year) + (PackSize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit4 <- gamm4(yrlsurv ~ (littersize) + s(mean_max_postdenning_temp) + (Year) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit5 <- gamm4(yrlsurv ~ (littersize) + s(Rain_PostDenningMonths) + (Year) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit6 <- gamm4(yrlsurv ~ (littersize) + s(mean_max_postdenning_temp) + s(Rain_PostDenningMonths) + (Year) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit7 <- gamm4(yrlsurv ~ (littersize) +  (Year) + (PackSize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit8 <- gamm4(yrlsurv ~ (littersize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit9 <- gamm4(yrlsurv ~ (littersize) + (Year) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit10 <- gamm4(yrlsurv ~ (littersize) + (PackSize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit11 <- gamm4(yrlsurv ~ (littersize) + s(mean_max_postdenning_temp) + (PackSize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")
fit12 <- gamm4(yrlsurv ~ (littersize) + s(Rain_PostDenningMonths) + (PackSize) , random=~  (1|Pack)+ (1|Year), data=dates, family="poisson")

table <- data.frame(formula = c(Reduce(paste, deparse(fit1$gam$formula)),Reduce(paste, deparse(fit2$gam$formula)),Reduce(paste, deparse(fit3$gam$formula)),Reduce(paste, deparse(fit4$gam$formula)),Reduce(paste, deparse(fit5$gam$formula)),
                                Reduce(paste, deparse(fit6$gam$formula)),Reduce(paste, deparse(fit7$gam$formula)),Reduce(paste, deparse(fit8$gam$formula)),Reduce(paste, deparse(fit9$gam$formula)),Reduce(paste, deparse(fit10$gam$formula)),
                                Reduce(paste, deparse(fit11$gam$formula)),Reduce(paste, deparse(fit12$gam$formula))), 
                    AICc=c(AICc(fit1$mer), AICc(fit2$mer), AICc(fit3$mer), AICc(fit4$mer), AICc(fit5$mer), AICc(fit6$mer),AICc(fit7$mer),  AICc(fit8$mer),AICc(fit9$mer),AICc(fit10$mer),
                           AICc(fit11$mer),    AICc(fit12$mer)))
(table <- table %>% arrange(AICc) %>% mutate(AICc=round(AICc,2), deltaAIC = round(AICc-min(AICc),2)))



 



