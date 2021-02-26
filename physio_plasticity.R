setwd("/Users/melissadrown/Documents/School/RSMAS/Research/Summer 2020/Writing/Physio_Plasticity_Manuscript/Scripts and Data")
data_wam <- read.csv("2019WAM_backgroun_corrected_flat.csv")
data_ctmax <- read.csv("OCNJ_s19_ctmax_2.csv")
data_cam <- read.csv("mass_cor_OCNJ_s19_cam_2.csv")
data_all <- read.csv("1_DLC_MKD_Physiol_trait_Oct2020_JMP.csv")

library(dplyr)
library(ggplot2)
library(lme4)
library(broom)
library(MASS)
library(raster)
library(car)
library(lubridate)
library(tidyr)

# look at data files
head(data_wam)
head(data_ctmax)
head(data_cam)

# function to calculate SE and mean summary table
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Date ecognized as a chronological date in all the files.
data_cam1 <- data_cam
data_cam <- subset(data_cam, data_cam$MR_pmol_s!="NA")
data_ctmax <- subset(data_ctmax, data_ctmax$ctmax_date!="4/15/19")

data_wam$wam_date <- as.POSIXct(as.character(data_wam$wam_date), format="%m/%d/%y")
data_ctmax$ctmax_date <- mdy(data_ctmax$ctmax_date)
data_cam$cam_date <- mdy(data_cam$cam_date)
data_cam1$cam_date <- mdy(data_cam1$cam_date)

# make acclimation order a factor
data_wam$accl_order <- as.factor(data_wam$accl_order)
data_ctmax$accl_order <- as.factor(data_ctmax$accl_order)

# make temp a factor
data_wam$temp <- as.factor(data_wam$temp)
data_ctmax$temp <- as.factor(data_ctmax$temp)
data_cam$temp <- as.factor(data_cam$temp)
data_cam1$temp <- as.factor(data_cam1$temp)

#######################################################################
# Whole Animal Metabolic Rate
# Linear model
wam_model <- data_wam %>% do(model = stepAIC(aov(log10(tenth_perc_MO2_mg_hr)~ temp + log10(mass_kg) + sex + accl_order + pop, data=.), scope = . ~ .^2, direction='both'))
# view best fit model
wam_model$model

# build model of best fit
wam_model_final <- aov(formula = log10(tenth_perc_MO2_mg_hr) ~ temp + log10(mass_kg) + 
      accl_order + pop + temp:accl_order + temp:log10(mass_kg), 
    data = data_wam)

# view ANOVA table
anova(wam_model_final)

# w/% var explained
af1 <- anova(wam_model_final)
af1ss <- af1$"Sum Sq"
print(cbind(af1,PctExp=af1ss/sum(af1ss)*100))

# Fig S1
wam.mass.model.temp <- summary(lm(log10(data_wam$tenth_perc_MO2_mg_hr)~ log10(data_wam$mass_kg)))
data_wam$residuals_log_mass2 <- wam.mass.model.temp$residuals

data_wam$temp <- factor(data_wam$temp, levels=c("12", "28"))

data.wam.temp <- summarySE(data_wam, measurevar="residuals_log_mass2", groupvars=c("temp","pop","accl_order"))

pd <- position_dodge(0.2) # move them .05 to the left and right

data.wam.temp
wam.temp <- ggplot(data.wam.temp, aes(x=temp, y=residuals_log_mass2, col=pop, shape=accl_order)) + 
  geom_errorbar(aes(ymin=residuals_log_mass2-se, ymax=residuals_log_mass2+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_color_manual(values=c("blue", "red", "purple"))+
  theme_bw()
wam.temp
#ggsave("combined_wam.eps", wam.temp, width=6, height=6)

#####
# plasticity in WAM
wam_cold <- subset(data_wam, data_wam$temp=="12")
wam_hot <- subset(data_wam, data_wam$temp=="28")

wam_all_wide <- merge(wam_cold, wam_hot, by="FishID")

# variable.x = 12°C, variable.y = 28°C
wam_diff <- wam_all_wide %>% group_by("FishID") %>%
  mutate(diff_wam_btw_accl=(tenth_perc_MO2_mg_hr.y-tenth_perc_MO2_mg_hr.x))

wam_diff <- wam_all_wide %>% group_by("FishID") %>%
  mutate(log2_ratio_accl=(log2(tenth_perc_MO2_mg_hr.y/tenth_perc_MO2_mg_hr.x)))

# log2 transformed ratio =1 is a 2x change in MR
# view log2 distribution and upper and lower 10%
quantile(wam_diff$log2_ratio_accl, c(0.1, 0.9))
mean(wam_diff$log2_ratio_accl) 

# group by level of plasticity
wam_diff <- mutate(wam_diff, plas_group = ifelse(log2_ratio_accl <0.3831368, "low_10th",
                                     ifelse(log2_ratio_accl > 1.7999633 , "high_10th", "middle")))

# are 12°C and 28°C metabolic rate correlated? No significant correlation p = 0.45, R2 = 0.010
summary(lm(wam_diff$tenth_perc_MO2_mg_hr.y~wam_diff$tenth_perc_MO2_mg_hr.x))

# split by level of plasticity
fits_plas_wam <- lmList(tenth_perc_MO2_mg_hr.y ~ tenth_perc_MO2_mg_hr.x | plas_group, data=wam_diff)
fits_plas_wam
summary(fits_plas_wam)
sapply(fits_plas_wam,function(x) summary(x)$r.squared)
# plots
ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.x, y=tenth_perc_MO2_mg_hr.y, col=plas_group)) + geom_point() + geom_smooth(method=lm) + 
  labs(x="12°C Metabolic Rate", y="28°C Metabolic Rate") + 
  theme_bw()

# plasticity in WAM
# 28°C p < 2e-16 ***
summary(lm(wam_diff$diff_wam_btw_accl~wam_diff$tenth_perc_MO2_mg_hr.y))

# 12°C p = 0.000553 ***
summary(lm(wam_diff$diff_wam_btw_accl~wam_diff$tenth_perc_MO2_mg_hr.x))

#####
# Figure 2
wam_12_plasticity <- ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.x, y=log2_ratio_accl)) + geom_point(aes(col=plas_group)) + 
  geom_smooth(method=lm, color="black", alpha=1) + 
  labs(x="12°C Whole Animal Metabolic Rate", y="Plasticity in Whole Animal Metabolic Rate (log2 ratio)") + 
  annotate("text", x=2, y = 5.5, label="p < 0.0001***, R2 = 0.19") +
  theme_bw() +
  theme(legend.position=c(3.2, 4))
wam_12_plasticity

wam_28_plasticity <- ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.y, y=log2_ratio_accl)) + geom_point(aes(col=plas_group)) + 
  geom_smooth(method=lm, color="black", alpha=1) + 
  labs(x="28°C Whole Animal Metabolic Rate", y="Plasticity in Whole Animal Metabolic Rate (log2 ratio)") + 
  annotate("text", x=4, y = 5.5, label="p < 0.0001***, R2 = 0.72") + 
  theme_bw()  +
  theme(legend.position=c(6, 4))
wam_28_plasticity

panelc <- ggplot(wam_diff, aes(x=log2_ratio_accl)) + 
  geom_histogram(bins=10, fill="lightgrey", color="black") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(log2_ratio_accl)),
             color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=0.3831368),
             color="blue", linetype="dashed",size=1)+
  geom_vline(aes(xintercept=1.7999633),
             color="red",linetype="dashed", size=1)+
  annotate("text", x=2.5, y = 17, label="mean=1.07") + 
  theme_bw()

paneld <- ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.x, y=tenth_perc_MO2_mg_hr.y, col=plas_group)) + geom_point() + geom_smooth(method=lm, alpha=1) + 
  labs(x="12°C Metabolic Rate", y="28°C Metabolic Rate") + 
  theme_bw() +
  theme(legend.position="none")
library(ggpubr)

wam_plas_fig <- ggarrange(wam_12_plasticity, wam_28_plasticity, panelc, paneld, ncol=2, nrow=2)
wam_plas_fig

#ggsave("wam_plas_fig.eps", wam_plas_fig, width=10, height=8)
#######################################################################
# Critical Thermal Maximum
# Linear model
ctmax_model <- data_ctmax %>%  do(model = stepAIC(aov(ctmax~ temp + weight_wam + sex + accl_order + pop, data=.), scope = . ~ .^2, direction='both'))
# view best fit model 
ctmax_model$model

# build model of best fit
ctmax_model_final <- aov(formula = ctmax ~ temp + sex + accl_order + temp:accl_order + 
                           + pop, data = data_ctmax)
# view ANOVA table
anova(ctmax_model_final)

# w/% var explained
af2 <- anova(ctmax_model_final)
af2ss <- af2$"Sum Sq"
print(cbind(af2,PctExp=af2ss/sum(af2ss)*100))

# Fig S2
data_ctmax$temp <- factor(data_ctmax$temp, levels=c("12", "28"))

data.ctmax.temp <- summarySE(data_ctmax, measurevar="ctmax", groupvars=c("temp","pop","accl_order"))

pd <- position_dodge(0.2) # move them .05 to the left and right


data.ctmax.temp
ctmax.temp <- ggplot(data.ctmax.temp, aes(x=temp, y=ctmax, col=pop, shape=accl_order)) + 
  geom_errorbar(aes(ymin=ctmax-se, ymax=ctmax+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_color_manual(values=c("blue","red", "purple"))+
  theme_bw() + facet_grid(rows = vars(temp), scale="free")
ctmax.temp

#ggsave("combined_ctmax.eps", ctmax.temp, width=6, height=6)

#####
# plasticity in CTmax
ctmax_cold <- subset(data_ctmax, data_ctmax$temp=="12")
ctmax_hot <- subset(data_ctmax, data_ctmax$temp=="28")

# variable.x = 12°C, variable.y =28°C
ctmax_all_wide <- merge(ctmax_cold, ctmax_hot, by="FishID", all = T)

# variable.x = 12°C, variable.y = 28°C
ctmax_diff <- ctmax_all_wide %>% group_by(FishID) %>%
  mutate(ctmax_log2ratio=log2(ctmax.y/ctmax.x))

# R2= 0.09, p =0.00384**
summary(lm(ctmax_all_wide$ctmax.y~ctmax_all_wide$ctmax.x))

ggplot(ctmax_all_wide, aes(x=ctmax.x, y=ctmax.y)) + geom_point() + geom_smooth(method=lm) + 
  labs(x="12°C CTmax", y="28°C CTmax") + theme_bw()
#####

# Figure 3
# higher 12°C = lower plasticity
ctmax_12_plasticity <- ggplot(ctmax_diff, aes(x=ctmax.x, y=ctmax_log2ratio)) + geom_point() + 
  geom_smooth(method=lm, alpha=1) + 
  labs(x="12°C CTmax", y="Plasticity in CTmax (log2 ratio)") + 
  annotate("text", x=37, y = 0.32, label="p < 0.001***, R2 = 0.93") + 
  theme_bw()
ctmax_12_plasticity

# p < 2e-16 ***, R2= 0.9012
summary(lm(ctmax_diff$ctmax_log2ratio~ctmax_diff$ctmax.x))

# horizontal lines, no relationship between 28°C and plasticity
ctmax_28_plasticity <- ggplot(ctmax_diff, aes(x=ctmax.y, y=ctmax_log2ratio)) + geom_point() + 
  geom_smooth(method=lm, alpha=1) + 
  labs(x="28°C CTmax", y="Plasticity in CTmax") + 
  annotate("text", x=42.8, y = 0.4, label="p = 0.71, R2 = -0.009") + 
  theme_bw()
ctmax_28_plasticity

# p =  0.9293, R2 = 9.194e-05
summary(lm(ctmax_diff$ctmax_log2ratio~ctmax_diff$ctmax.y))

#ggsave("ctmax_12_plasticity.eps", ctmax_12_plasticity)
#ggsave("ctmax_28_plasticity.eps", ctmax_28_plasticity)

#######################################################################
# Cardiac Metabolic Rate
# Linear model
cam_model <- data_cam %>% do(model = stepAIC(aov(MR_pmol_s ~ temp + cam_mass + substrate + sex + pop, data=.), scope = . ~ .^2, direction = 'both'))
# view best fit model
cam_model$model

# build model of best fit
cam_model_final <- aov(MR_pmol_s~temp + cam_mass + substrate + sex + 
                         temp:substrate + temp:cam_mass + cam_mass:substrate + cam_mass:sex + 
                         temp:sex, data = data_cam)

# view ANOVA table
anova(cam_model_final)

# w/% var explained
af3 <- anova(cam_model_final)
af3ss <- af3$"Sum Sq"
print(cbind(af3,PctExp=af3ss/sum(af3ss)*100))

# Box plot of CaM with no residuals
cb_data <- summarySE(data_cam_hearts, measurevar="MR_pmol_s", groupvars=c("temp","substrate"))

cb_boxes <- ggplot(cb_data, aes(y=MR_pmol_s, 
                                x=substrate:temp, col=temp))  + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(aes(ymin=MR_pmol_s-se, ymax=MR_pmol_s+se), width=.2) +
  theme_bw() + scale_color_manual(values=c("blue", "orange"))+
  labs(y="Cardiac Metabolic Rate (pmol/sec)", x="Substrate") +
  annotate("text", x=1, y=47, label="a") +
  annotate("text", x=2, y=53, label="a") +
  annotate("text", x=3, y=42, label="b") +
  annotate("text", x=4, y=36, label="b") +
  annotate("text", x=5, y=32, label="b") +
  annotate("text", x=6, y=32, label="b") +
  annotate("text", x=7, y=32, label="b") +
  annotate("text", x=8, y=21, label="c") +
  labs(y="Cardiac Metabolic Rate (pmol/sec)", x="Substrate")

cb_boxes
#ggsave("no_mass_cam_boxplot.eps", cb_boxes, width=10, height=6)

# Model without any mass covariate
cam_no_mass <- data_cam %>% do(model = stepAIC(aov(MR_pmol_s ~ temp + substrate + sex + pop, data=.), scope = . ~ .^2, direction = 'both'))
# view best fit model
cam_no_mass$model
cam_no_mass_mod <- aov(formula = MR_pmol_s ~ temp + substrate + pop + temp:substrate, data = data_cam)

# no change in temp:substrate within a substrate
# changes relationship among substrates

af10 <- anova(cam_no_mass_mod)
af10ss <- af10$"Sum Sq"
print(cbind(af10,PctExp=af10ss/sum(af10ss)*100))

anova(cam_no_mass_mod)
TukeyHSD(cam_no_mass_mod, which="temp:substrate")

# Figure 4
data_cam_mass <- lm(data_cam$MR_pmol_s~data_cam$cam_mass)
data_cam$MR_res <- data_cam_mass$residuals
data_cam$substrate <- factor(data_cam$substrate, levels=c("Glucose", "FA", "LKA","INH"))

cam_boxes_data <- summarySE(data_cam, measurevar="MR_res", groupvars=c("temp","substrate"))

cam_boxes3 <- ggplot(cam_boxes_data, aes(y=MR_res, 
                                         x=substrate:temp, col=temp))  + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(aes(ymin=MR_res-se, ymax=MR_res+se), width=.2) +
  theme_bw() + scale_color_manual(values=c("blue", "orange")) +
  annotate("text", x=1, y=16, label="a") +
  annotate("text", x=2, y=22, label="a") +
  annotate("text", x=3, y=8, label="b") +
  annotate("text", x=4, y=4, label="b") +
  annotate("text", x=5, y=0, label="bc") +
  annotate("text", x=6, y=0, label="bc") +
  annotate("text", x=7, y=-2, label="c") +
  annotate("text", x=8, y=-12, label="d") +
  labs(y="Cardiac Metabolic Rate (pmol/sec)", x="Substrate")
cam_boxes3

#ggsave("cam_box_new.eps", cam_boxes3, height=6, width=8)

# CaM Model with heart mass instead of body mass
cam_model_h <- data_cam %>% do(model = stepAIC(aov(MR_pmol_s ~ temp + heart_mass + substrate + sex + pop, data=.), scope = . ~ .^2, direction = 'both'))
# view best fit model
cam_model_h$model
cam_model_h_final <- aov(formula = MR_pmol_s ~ temp + heart_mass + substrate + sex + pop + temp:heart_mass + temp:substrate + heart_mass:sex + heart_mass:substrate + heart_mass:pop + temp:sex, data = data_cam)

# no change in temp:substrate within a substrate
# changes relationship among substrates

af4 <- anova(cam_model_h_final)
af4ss <- af4$"Sum Sq"
print(cbind(af4,PctExp=af4ss/sum(af4ss)*100))

anova(cam_model_final)
anova(cam_model_h_final)

TukeyHSD(cam_model_h_final, which="temp:substrate")
TukeyHSD(cam_model_final, which="temp:substrate")

# second order terms for heart mass residual model
# Due to outlier? Removing one outlier removes signif of heart_mass:substrate and heart_mass:sex
data_cam2 <- subset(data_cam, data_cam$heart_mass<0.04)

# CaM Model with heart mass instead of body mass
cam_model_h2 <- data_cam2 %>% do(model = stepAIC(aov(MR_pmol_s ~ temp + heart_mass + substrate + sex + pop, data=.), scope = . ~ .^2, direction = 'both'))
# view best fit model
cam_model_h2$model

cam_model_h2 <-  aov(formula = MR_pmol_s ~ temp + heart_mass + substrate + sex + 
                       pop + temp:substrate + temp:heart_mass + heart_mass:pop + 
                       heart_mass:sex + temp:sex, data = data_cam2)
anova(cam_model_h2)
TukeyHSD(cam_model_h2, which="temp:substrate")

# MR ~heart_mass:pop
ggplot(data_cam2, aes(x=heart_mass, y=MR_res, col=pop))+
  geom_point()+
  geom_smooth(method="lm") +
  facet_wrap(~temp) +
  scale_color_manual(values=c("blue","red","purple"))+
  labs(y="Cardiac Metabolic Rate Heart Mass Residuals (pmol/sec)", x="Heart Mass (g)")+
  theme_bw()

# MR ~temp:heart_mass, see Supple. Fig 2. slope of MR vs heart mass differs by temp

## Figure S4 (CaM Boxplot with heart mass residuals)
data_cam_hearts <- subset(data_cam, data_cam$heart_mass!="NA")
data_cam_mass2 <- lm(data_cam_hearts$MR_pmol_s~data_cam_hearts$heart_mass)
data_cam_hearts$MR_res_heart <- data_cam_mass2$residuals
data_cam_hearts$substrate <- factor(data_cam_hearts$substrate, levels=c("Glucose", "FA", "LKA","INH"))

cam_boxes_data_hearts <- summarySE(data_cam_hearts, measurevar="MR_res_heart", groupvars=c("temp","substrate"))

cam_boxes4 <- ggplot(cam_boxes_data_hearts, aes(y=MR_res_heart, 
                                         x=substrate:temp, col=temp))  + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(aes(ymin=MR_res_heart-se, ymax=MR_res_heart+se), width=.2) +
  theme_bw() + scale_color_manual(values=c("blue", "orange"))+
  labs(y="Cardiac Metabolic Rate heart mass residuals (pmol/sec)", x="Substrate") +
  annotate("text", x=1, y=12, label="a") +
  annotate("text", x=2, y=22, label="a") +
  annotate("text", x=3, y=8, label="b") +
  annotate("text", x=4, y=6, label="bc") +
  annotate("text", x=5, y=-2, label="c") +
  annotate("text", x=6, y=1, label="c") +
  annotate("text", x=7, y=-2, label="c") +
  annotate("text", x=8, y=-9, label="d") +
  labs(y="Cardiac Metabolic Rate Heart Mass Residuals (pmol/sec)", x="Substrate")
cam_boxes4

#ggsave("cam_box_heart.eps", cam_boxes4, height=6, width=8)

# Figure 5
# Heart mass and body mass by acclimation temp
t.test(data_cam$heart_mass~data_cam$temp)
t.test(data_cam$cam_mass~data_cam$temp)

scd <- summarySE(data_cam, measurevar="heart_mass", groupvars="temp", na.rm=TRUE)

cam_heartmass <- ggplot(scd, aes(x=temp, y=heart_mass)) + 
  geom_errorbar(aes(ymin=heart_mass-se, ymax=heart_mass+se), width=.1) +
  geom_point(col=c("blue", "orange"), size=5) + 
  annotate("text", x=2.3, y= 0.0175, label="p < 0.0001") +
  ylab("Ventricular Mass (g)") +
  xlab("Acclimation Temperature (°C)") +
  theme_bw()
cam_heartmass
scd2 <- summarySE(data_cam, measurevar="cam_mass", groupvars="temp", na.rm=TRUE)

cam_bodymass <- ggplot(scd2, aes(x=temp, y=cam_mass)) + 
  geom_errorbar(aes(ymin=cam_mass-se, ymax=cam_mass+se), width=.1) +
  geom_point(col=c("blue", "orange"), size=5) + 
  annotate("text", x=2.3, y=12.35, label="p=0.43") +
  ylab("Body Mass (g)") +
  xlab("Acclimation Temperature (°C)") +
  theme_bw()

cam_heartmass
cam_bodymass

# calculate heart/body mass ratio
hb_ratio <- data_cam %>% group_by(FishID) %>%
  mutate(hb_ratio=heart_mass/cam_mass)

hb_ratio <- subset(hb_ratio, hb_ratio!="NA")
hb_rat_dat <- summarySE(hb_ratio, measurevar="hb_ratio", groupvars=c("temp"))

hb_rat_dat

hb_ratio_plot <- ggplot(data=hb_rat_dat, aes(x=temp, y=hb_ratio)) +  geom_errorbar(aes(ymin=hb_ratio-se, ymax=hb_ratio+se), width=0.2)+ geom_point(color=c("blue","orange"),size=2.5)+theme_bw() +
  labs(y="Heart Mass to Body Mass Ratio", x="Acclimation Temperature (°C)")

#ggsave("hb_ratio_plot.eps", hb_ratio_plot, width=6, height=6)
#ggsave("cam_bodymass.eps", cam_bodymass, height=6, width=6)
#ggsave("cam_heartmass.eps", cam_heartmass, height=6, width=6)

# Figure S3
hm_bm_temp <- ggplot(data=data_cam, aes(x=cam_mass, y=heart_mass, col=temp)) +
  geom_point() +
  geom_smooth(method="lm", alpha=1)+
  scale_color_manual(values=c("blue", "orange"))+
  labs(y="Heart Mass (g)", x="Body Mass (g)")+
  theme_bw()
hm_bm_temp
#ggsave("hm_bm_temp.eps",  hm_bm_temp, width=6, height=6)

lm(data=data_cam, heart_mass~cam_mass+temp+cam_mass:temp)

fits <- lmList(heart_mass ~ cam_mass | temp, data=data_cam)
summary(fits)


####################################################################3
# Figure S6 -- Partial Correlations
## p.cor data from JMP
# 12°C
mat_cor_12 <- read.csv("Partial_correlations_12.csv")
mat_cor_12_ps <- read.csv("Partial_correlations_12_ps.csv")

rownames(mat_cor_12) <- mat_cor_12[,1]
mat_cor_12 <- mat_cor_12[,2:7]

rownames(mat_cor_12_ps) <- mat_cor_12_ps[,1]
mat_cor_12_ps <- mat_cor_12_ps[,2:7]

#28°C 
mat_cor_28 <- read.csv("Partial_correlations_28.csv")
mat_cor_28_ps <- read.csv("Partial_correlations_28_ps.csv")

rownames(mat_cor_28) <- mat_cor_28[,1]
mat_cor_28 <- mat_cor_28[,2:7]

rownames(mat_cor_28_ps) <- mat_cor_28_ps[,1]
mat_cor_28_ps <- mat_cor_28_ps[,2:7]

library(corrplot)

corrplot(as.matrix(mat_cor_12), type = "upper", p.mat=as.matrix(mat_cor_12_ps), sig.level=0, insig="p-value", tl.col = "black", tl.srt = 45, method="color", diag=FALSE)

corrplot(as.matrix(mat_cor_28), type = "upper", p.mat=as.matrix(mat_cor_28_ps), sig.level=0, insig="p-value", tl.col = "black", tl.srt = 45, method="color", diag=FALSE)

## LKA_CaM has an interaction between heart mass and temp, no other substrate does.
summary(lm(data_all$GLU_CamMass_Residuals~data_all$heart_mass*data_all$Acc_Temp))
summary(lm(data_all$LKA__CamMass_Residuals~data_all$heart_mass*data_all$Acc_Temp))
summary(lm(data_all$FA_CamMass_Residuals~data_all$heart_mass*data_all$Acc_Temp))
summary(lm(data_all$END_CamMass_Residuals~data_all$heart_mass*data_all$Acc_Temp))

######################################################################
#### Mean and SD for all traits
summary(aov(CtMax_.C~Acc_Temp*heart_mass, data=data_all))
summary(aov(Log10_WAM_mg_hr~Acc_Temp*heart_mass, data=data_all))

## Means
data_cam_hot <- subset(data_cam, data_cam$temp=="28")
data_cam_cold <- subset(data_cam, data_cam$temp=="12")

data_cam_hot_glu <- subset(data_cam_hot, data_cam_hot$substrate=="Glucose")
data_cam_hot_FA <- subset(data_cam_hot, data_cam_hot$substrate=="FA")
data_cam_hot_LKA <- subset(data_cam_hot, data_cam_hot$substrate=="LKA")
data_cam_hot_end <- subset(data_cam_hot, data_cam_hot$substrate=="INH")

data_cam_cold_glu <- subset(data_cam_cold, data_cam_cold$substrate=="Glucose")
data_cam_cold_FA <- subset(data_cam_cold, data_cam_cold$substrate=="FA")
data_cam_cold_LKA <- subset(data_cam_cold, data_cam_cold$substrate=="LKA")
data_cam_cold_end <- subset(data_cam_cold, data_cam_cold$substrate=="INH")

mean(wam_hot$tenth_perc_MO2_mg_hr)
mean(wam_cold$tenth_perc_MO2_mg_hr)
mean(wam_hot$residuals_log_mass2)
mean(wam_cold$residuals_log_mass2)

mean(ctmax_hot$ctmax)
mean(ctmax_cold$ctmax)

mean(data_cam_hot_glu$MR_pmol_s)
mean(data_cam_hot_FA$MR_pmol_s)
mean(data_cam_hot_LKA$MR_pmol_s)
mean(data_cam_hot_end$MR_pmol_s)

mean(data_cam_cold_glu$MR_pmol_s)
mean(data_cam_cold_FA$MR_pmol_s)
mean(data_cam_cold_LKA$MR_pmol_s)
mean(data_cam_cold_end$MR_pmol_s)

sd(wam_hot$tenth_perc_MO2_mg_hr)
sd(wam_cold$tenth_perc_MO2_mg_hr)
sd(wam_hot$residuals_log_mass2)
sd(wam_cold$residuals_log_mass2)

sd(ctmax_hot$ctmax)
sd(ctmax_cold$ctmax)

sd(data_cam_hot_glu$MR_pmol_s)
sd(data_cam_hot_FA$MR_pmol_s)
sd(data_cam_hot_LKA$MR_pmol_s)
sd(data_cam_hot_end$MR_pmol_s)

sd(data_cam_cold_glu$MR_pmol_s)
sd(data_cam_cold_FA$MR_pmol_s)
sd(data_cam_cold_LKA$MR_pmol_s)
sd(data_cam_cold_end$MR_pmol_s)

mean(data_ctmax$weight_wam)

#################################################################
### Plasticity Models ###
head(wam_diff)

## WAM
# Linear model
wam_plas_mod <- wam_diff %>% do(model = stepAIC(aov(log2_ratio_accl~ mass_kg.y + sex.y + accl_order.y + pop.y, data=.), scope = . ~ .^2, direction='both'))

# view best fit model
wam_plas_mod$model

wam_plas_mod_final <-  aov(formula = log2_ratio_accl ~ mass_kg.y + sex.y + accl_order.y + 
                             pop.y + mass_kg.y:pop.y,data = wam_diff)

af_wamplas <- anova(wam_plas_mod_final)
af_wamplas_ss <- af_wamplas$"Sum Sq"
print(cbind(af_wamplas,PctExp=af_wamplas_ss/sum(af_wamplas_ss)*100))

ggplot(wam_diff, aes(x=mass_kg.y, y=log2_ratio_accl, col=pop.y)) +
  geom_point()+
  geom_smooth(method="lm") +
  scale_color_manual(values=c("blue","red","purple"))+
  theme_bw()

## CTmax
# Linear model
# model using all factors 
ctmax_plas_mod_final <-  aov(formula =ctmax_log2ratio  ~ weight_wam.y + sex.y + accl_order.y + 
                             pop.y + weight_wam.y:pop.y + weight_wam.y:sex.y + weight_wam.y:accl_order.y + sex.y:accl_order.y + sex.y:pop.y + accl_order.y:pop.y, data = ctmax_diff)

# find best fit
stepAIC(ctmax_plas_mod_final, scope=.~.^2, direction='both')

ctmax_plas_mod_final <- aov(formula = ctmax_log2ratio ~ weight_wam.y + sex.y + accl_order.y + 
                              pop.y + weight_wam.y:pop.y + weight_wam.y:sex.y + sex.y:pop.y + 
                              accl_order.y:pop.y, data = ctmax_diff)

af_ctmaxplas <- anova(ctmax_plas_mod_final)
af_ctmaxplas_ss <- af_ctmaxplas$"Sum Sq"
print(cbind(af_ctmaxplas,PctExp=af_ctmaxplas_ss/sum(af_ctmaxplas_ss)*100))

ggplot(ctmax_diff, aes(x=accl_order.y, y=ctmax_log2ratio)) +
  geom_boxplot()+
  scale_color_manual(values=c("blue","red","purple"))+
  labs(y="Plasticity in CTmax (°C)", x="Acclimation Order")+
  theme_bw()

TukeyHSD(ctmax_plas_mod_final, which="accl_order.y")