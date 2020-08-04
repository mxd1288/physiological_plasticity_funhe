setwd("/Users/melissadrown/Documents/School/RSMAS/Research/Summer 2020/Writing/Phenotyping/Spring 2019 Phenotyping/Data for Phenotypes")
data_wam <- read.csv("2019WAM_backgroun_corrected_flat.csv")
data_ctmax <- read.csv("OCNJ_s19_ctmax_2.csv")
data_cam <- read.csv("mass_cor_OCNJ_s19_cam_2.csv")

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

# Fix the Date to be recognized as a chronological date in all the files.
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

# Fig 2
wam.mass.model.temp <- summary(lm(log10(data_wam$tenth_perc_MO2_mg_hr)~ log10(data_wam$mass_kg)))
data_wam$residuals_log_mass2 <- wam.mass.model.temp$residuals

data_wam$temp <- factor(data_wam$temp, levels=c("12", "28"))

data.wam.temp <- summarySE(data_wam, measurevar="residuals_log_mass2", groupvars=c("temp","pop"))

pd <- position_dodge(0.2) # move them .05 to the left and right

data.wam.temp
wam.temp <- ggplot(data.wam.temp, aes(x=temp, y=residuals_log_mass2, col=pop)) + 
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

# are 12°C and 28°C metabolic rate correlated? No significant correlation p = 0.45, R2 = 0.010
summary(lm(wam_all_wide$tenth_perc_MO2_mg_hr.y~wam_all_wide$tenth_perc_MO2_mg_hr.x))

ggplot(wam_all_wide, aes(x=tenth_perc_MO2_mg_hr.x, y=tenth_perc_MO2_mg_hr.y)) + geom_point() + geom_smooth(method=lm) + 
  labs(x="12°C Metabolic Rate", y="28°C Metabolic Rate") + theme_bw()

# plasticity in WAM
# 28°C p < 2e-16 ***
summary(lm(wam_diff$diff_wam_btw_accl~wam_diff$tenth_perc_MO2_mg_hr.y))

# 12°C p = 0.000553 ***
summary(lm(wam_diff$diff_wam_btw_accl~wam_diff$tenth_perc_MO2_mg_hr.x))

#####
# Supple. Fig 1
wam_12_plasticity <- ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.x, y=diff_wam_btw_accl)) + geom_point() + 
  geom_smooth(method=lm) + 
  labs(x="12°C Whole Animal Metabolic Rate", y="Plasticity in Whole Animal Metabolic Rate") + 
  annotate("text", x=3.1, y = 5.5, label="p < 0.0001***, R2 = 0.19") + 
  theme_bw()

wam_28_plasticity <- ggplot(wam_diff, aes(x=tenth_perc_MO2_mg_hr.y, y=diff_wam_btw_accl)) + geom_point() + 
  geom_smooth(method=lm) + 
  labs(x="28°C Whole Animal Metabolic Rate", y="Plasticity in Whole Animal Metabolic Rate") + 
  annotate("text", x=6, y = 5.5, label="p < 0.0001***, R2 = 0.72") + 
  theme_bw()

wam_plas_fig <- grid.arrange(wam_12_plasticity, wam_28_plasticity, nrow=1)

#ggsave("wam_plas_fig.eps", wam_plas_fig, width=10, height=6)
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

# Fig 3
data_ctmax$temp <- factor(data_ctmax$temp, levels=c("12", "28"))

data.ctmax.temp <- summarySE(data_ctmax, measurevar="ctmax", groupvars=c("temp","pop"))

pd <- position_dodge(0.2) # move them .05 to the left and right

data.ctmax.temp
ctmax.temp <- ggplot(data.ctmax.temp, aes(x=temp, y=ctmax, col=pop)) + 
  geom_errorbar(aes(ymin=ctmax-se, ymax=ctmax+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_color_manual(values=c("blue", "red", "purple"))+
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
  mutate(diff_ctmax_btw_accl=ctmax.y-ctmax.x)

# R2= 0.09, p =0.00384**
summary(lm(ctmax_all_wide$ctmax.y~ctmax_all_wide$ctmax.x))

ggplot(ctmax_all_wide, aes(x=ctmax.x, y=ctmax.y)) + geom_point() + geom_smooth(method=lm) + 
  labs(x="12°C CTmax", y="28°C CTmax") + theme_bw()
#####

# Fig 4
# higher 12°C = lower plasticity
ctmax_12_plasticity <- ggplot(ctmax_diff, aes(x=ctmax.x, y=diff_ctmax_btw_accl)) + geom_point() + 
  geom_smooth(method=lm) + 
  labs(x="12°C CTmax", y="Plasticity in CTmax") + 
  annotate("text", x=34.25, y = 8.75, label="p < 0.001***, R2 = 0.90") + 
  theme_bw()
ctmax_12_plasticity

# p < 2e-16 ***, R2= 0.9012
summary(lm(ctmax_diff$diff_ctmax_btw_accl~ctmax_diff$ctmax.x))

# horizontal lines, no relationship between 28°C and plasticity
ctmax_28_plasticity <- ggplot(ctmax_diff, aes(x=ctmax.y, y=diff_ctmax_btw_accl)) + geom_point() + 
  geom_smooth(method=lm) + 
  labs(x="28°C CTmax", y="Plasticity in CTmax") + 
  annotate("text", x=42, y = 8.75, label="p = 0.93, R2 < 0.0001") + 
  theme_bw()
ctmax_28_plasticity

# p =  0.9293, R2 = 9.194e-05
summary(lm(ctmax_diff$diff_ctmax_btw_accl~ctmax_diff$ctmax.y))

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

# Fig 5
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

#Fig 6
# Heart mass and body mass by acclimation temp
t.test(data_cam$heart_mass~data_cam$temp)
t.test(data_cam$cam_mass~data_cam$temp)

scd <- summarySE(data_cam, measurevar="heart_mass", groupvars="temp", na.rm=TRUE)

cam_heartmass <- ggplot(scd, aes(x=temp, y=heart_mass)) + 
  geom_errorbar(aes(ymin=heart_mass-se, ymax=heart_mass+se), width=.1) +
  geom_point(col=c("blue", "orange"), size=3) + 
  annotate("text", x=2.3, y= 0.0175, label="p < 0.0001") +
  ylab("Ventricular Mass (g)") +
  xlab("Acclimation Temperature (°C)") +
  theme_bw()

scd2 <- summarySE(data_cam, measurevar="cam_mass", groupvars="temp", na.rm=TRUE)

cam_bodymass <- ggplot(scd2, aes(x=temp, y=cam_mass)) + 
  geom_errorbar(aes(ymin=cam_mass-se, ymax=cam_mass+se), width=.1) +
  geom_point(col=c("blue", "orange"), size=3) + 
  annotate("text", x=2.3, y=12.35, label="p=0.43") +
  ylab("Body Mass (g)") +
  xlab("Acclimation Temperature (°C)") +
  theme_bw()

cam_heartmass
cam_bodymass
#ggsave("cam_bodymass.eps", cam_bodymass, height=6, width=6)
#ggsave("cam_heartmass.eps", cam_heartmass, height=6, width=6)

# Supple. Fig 2
data_cam_cold <- subset(data_cam, data_cam$temp=="12")
data_cam_hot <- subset(data_cam, data_cam$temp=="28")

summary(lm(data=data_cam_cold, heart_mass~cam_mass))

cam_body_heart_12 <- ggplot(data_cam_cold, aes(x=cam_mass, y=heart_mass), fill="temp") + 
  geom_point(color="blue") +
  geom_smooth(method=lm, color="blue") + 
  annotate("text", x=8.5, y=0.051, label="R2=0.18, p<0.0001") +
  labs(x="Body Mass (g)", y="Heart Mass (g)") + 
  theme_bw()

cam_body_heart_12

summary(lm(data=data_cam_hot, heart_mass~cam_mass))

cam_body_heart_28 <- ggplot(data_cam_hot, aes(x=cam_mass, y=heart_mass), fill="temp") + 
  geom_point(color="red") +
  geom_smooth(method=lm, color="red") + 
  annotate("text", x=8.5, y=0.051, label="R2=0.64, p<0.0001") +
  labs(x="Body Mass (g)", y="Heart Mass (g)") + 
  theme_bw()

cam_body_heart_28

cam_body_heart <- grid.arrange(cam_body_heart_12, cam_body_heart_28, nrow=1)

#ggsave("cam_heart_bodymass_fig.eps", cam_body_heart, width=10, height=6)
######################################################################
# Correlation among phenotypes
# set up data
data_cam1$MR_pmol_s[is.na(data_cam1$MR_pmol_s)] <- 0
data_cam_model <- lm(data_cam1$MR_pmol_s~data_cam1$cam_mass)
data_cam1$MR_res <- data_cam_model$residuals

data_cam_sub <- data_cam1[c(5,7,8,9,10,14,17)]
wide_cam <- data_cam_sub %>%
  spread(substrate, MR_res)
wide_cam2<-na.omit(wide_cam)
head(wide_cam2, 24)

all_data_wide <- merge(data_ctmax, data_wam, by="FishID")
all_data_wide <- merge(all_data_wide, wide_cam2, by="FishID")
head(all_data_wide)

all_data_hot <- subset(all_data_wide, all_data_wide$temp.x=="28")
all_data_cold <- subset(all_data_wide,all_data_wide$temp.x=="12")

# Between temps 
# WAM
summary(lm(wam_all_wide$tenth_perc_MO2_mg_hr.y~wam_all_wide$tenth_perc_MO2_mg_hr.x))

# CTmax
summary(lm(ctmax_all_wide$ctmax.y~ctmax_all_wide$ctmax.x))

# Correlation plots - within temps
library(corrplot)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
str(all_data_hot)

wide_all_corr1 <- cor(all_data_hot[,c(13,34,39:42)])
p.mat <- cor.mtest(wide_all_corr1)
head(p.mat[, 1:6])
M1 <- corrplot(wide_all_corr1, addCoef.col = "black", p.mat = p.mat, sig.level = 0.01, insig = "blank", type="upper", order="hclust")
colnames(M1) <- c("CTmax", "Metabolic Rate", "Fatty Acids", "Glucose", "Endogenous", "LKA")
rownames(M1) <- c("CTmax", "Metabolic Rate", "Fatty Acids", "Glucose", "Endogenous", "LKA")

correlations_28 <- corrplot(M1, addCoef.col = "black", p.mat = p.mat, sig.level = 0.05, insig = "blank", type="upper", order="hclust")

wide_all_corr2 <- cor(all_data_cold[,c(13,34,39:42)])
p.mat2 <- cor.mtest(wide_all_corr2)
head(p.mat2[, 1:6])
M2 <- corrplot(wide_all_corr2, addCoef.col = "black", p.mat = p.mat2, sig.level = 0.01, insig = "blank", type="upper", order="hclust")
colnames(M2) <- c("CTmax", "Metabolic Rate", "Fatty Acids", "Glucose", "Endogenous", "LKA")
rownames(M2) <- c("CTmax", "Metabolic Rate", "Fatty Acids", "Glucose", "Endogenous", "LKA")

correlations_12 <- corrplot(M2, addCoef.col = "black", p.mat = p.mat2, sig.level = 0.05, insig = "blank", type="upper", order="hclust")
