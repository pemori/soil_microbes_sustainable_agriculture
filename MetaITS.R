  # ---
  # title: "Soil diversity 16s"
  # authors: "Pedro Mondaca - Francisco E. Fonturbel"
  # date: "07/nov/2022"
  # ---
  
  ## A meta-analysis on agricultural management effects on soil diversity
  #  This meta-analysis assesses the effects of different types of agricultural management on microbial diversity of the soil using 16s marker.
  # 
  setwd("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023")
  
  #Loading data file
  #Loading packages required
  library(meta)
  library(metafor)
  library(readxl)
  library(dplyr)
  dat<-read_xlsx("MetaITS.xlsx", sheet = "Sheet1")
  
  dat$pair<-as.character(dat$pair) 
  dat$paper_code<-as.character(dat$paper_code)
  dat$ID<-as.character(dat$ID)
  str(dat)
  #data $ sample_duration <- paste (data $ sample, data $ duration, sep = " - ")
  # dat<-data %>% 
  #   filter(Hedges_g > -10 & Hedges_g < 10) %>% 
  #   filter(Sc !=0 &  Se !=0 & Xc !=0 & Xe !=0 & Hedges_g != "#DIV/0!")
  # dat$pair<-as.character(dat$pair) 
  # str(dat)
  #SECTION 1: ALPHA DIVERSITY ####

#----------R_family-----------

R_family<-dat %>%
  dplyr::filter(response=="R_family")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = R_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)
#No Significant REM
#No Significant heterogeneity

dmetar::find.outliers(overall.mod)
#No outliers

#There is no **significant overall effect** on soil microbial diversity. There is no **significant heterogeneity** in the model.
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

R_family <- filter(R_family, pair != "21" & pair != "24" & pair !=  "19")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc,
                      data = R_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")

summary.meta(overall.mod)
#ns RME and heterogeneity

#### (2) Baujat plot 

#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")

R_family <- filter(R_family, pair != "1" & pair != "26")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc,
                      data = R_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)
#ns RME and heterogeneity

#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pub-bias.html


#### (3) Correlation

#We will examine the correlation between the effect size and sample size using Kendall's approach.
cor.test(overall.mod$TE, R_family$Variance2, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, R_family$Variance2)
#There is no correlation

#### (4) Trim & Fill 
trimfill(overall.mod)
summary(overall.mod)

# FOREST PLOT de FF
forest.meta(overall.mod, 
            sortvar = TE,
            leftlabs = c("management", "g", "SE"))

#FOREST PLOT alternativo
forest.meta(overall.mod, layout = "RevMan5")

#pendiente plotear solo el overall effect

#PRUEBA DE ORCHARD PLOT
# 
# install.packages('pacman')
# rm(list = ls())
# devtools::install_github("daniel1noble/orchaRd", force = TRUE)
# pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans,
#               ape, phytools, flextable)

library(orchaRd)

R_family
# We need to calculate the effect sizes, in this case d
R_family <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                  n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = R_family)

R_fam_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                     data = R_family)
summary(R_fam_MA)

model_results <- mod_results(R_fam_MA, mod = "1", at = NULL, data = R_family, group = "paper_code")
model_results

#ORCHARD PLOT
library(ggplot2)
orchaRd::orchard_plot(R_fam_MA, mod = "1", data = R_family, group = "paper_code",
                      xlab = "Standardised mean difference", 
                      fill=FALSE, 
                      transfm = "none")

orchard_pot(R_fam_MA, mod = "1", data = R_family, group = "paper_code",
                      xlab = "Standardised mean difference", 
                      fill=FALSE, 
                      transfm = "none") + ylim(-2,4)

#600 x 220

orchaRd::caterpillars(R_fam_MA, mod = "1", data = R_family, group = "paper_code", 
                             xlab = "Standardised mean difference",
                      xlim(c(0,2)))
                          
#forest(R_fam_MA, xlim=c(0,2))

#Ahora con management

R_fam_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                   data = R_family)
summary(R_fam_MA)

model_results <- mod_results(R_fam_MA, mod = "management", at = NULL, data = R_family, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(R_fam_MA, mod = "management", data = R_family, group = "paper_code", xlab = "Standardised mean difference",
                      transfm = "none") + ylim(-2,4)


# 600 x 350

#ojo que difieren los resultados con meta

## Examining agric management effect ####

manage.mod<-update.meta(overall.mod,
                        subgroup = management)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

# USE OF PESTICIDES HAVE SIGNIFICANT EFFECT ON FUNGAL RICHNESS

#Duration of the experiment
manage.mod<-update.meta(overall.mod,
                        subgroup = duration)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))
#NO EFFECT OF DURATION

#Effect of climate
manage.mod<-update.meta(overall.mod,
                        subgroup = climate)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))
#No EFFECT of CLIMATE

#Crop
manage.mod<-update.meta(overall.mod,
                        subgroup = plant_man)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))
#No crop effect

#----------H_family-----------

H_family<-dat %>%
  dplyr::filter(response=="H_family")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = H_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)
#No significant REM
#No significant heterogeneity


#Interpretaci?n de heterogeneidad
#Output: "Quantifying heterogeneity"
#tau^2: If the confidence interval does not contain zero indicates that some between-study
#heterogeneity exists in the data
#The value of tau (e.g. tau=0.5) means that the true effect sizes have an estimated 
#standard deviation of 0.5, expressed on the scale of the effect size metric (Hedges'g in this case)

#Q is compared with K-1 degrees of freedom. If it is significant, we say heterogeneity 
#is high, but we should not base our assessment on the Q alone.

#install.packages("remotes")
#remotes::install_github("MathiasHarrer/dmetar")

library(dmetar)

dmetar::find.outliers(overall.mod)
#No outliers

#There is no **significant overall effect** on soil microbial diversity. Also, there is a **significant heterogeneity** in the model (which is expected in this kind of models)
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot
funnel.meta(overall.mod, studlab = TRUE)

  cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

  H_family <- filter(H_family, pair != "1" & pair != "21")
  
  overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                        data = H_family,
                        studlab = pair,
                        fixed = FALSE,
                        random = TRUE,
                        method.tau = "REML",
                        hakn = TRUE,
                        prediction = TRUE,
                        sm = "SMD",
                        method.smd = "Hedges")
  
  summary.meta(overall.mod)
#ns REM and Q

#### (2) Baujat plot 
  
#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")  

#H_family <- filter(H_family, pair != "26" )

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = H_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")

summary.meta(overall.mod)

    
#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pub-bias.html

#### (3) Correlation

#We will examine the correlation between the effect size and sample size using Kendall's approach.
cor.test(overall.mod$TE, overall.mod$seTE, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, overall.mod$seTE)

#The correlation consider "d" in the calculation which lead to bias when among-study variance are large.
#Thus, we will examine the correlation using equation 4 of Hamman et al., 2018, as an alternative method of the within study variance that avoid using the estimate d in its calculation (proposed by Hedges 1982)
cor.test(overall.mod$TE, H_family$Variance2, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, H_family$Variance2)

# No significant correlation

#### (4) Trim & Fill 
trimfill(overall.mod)

summary(overall.mod)

#REM is not significant in both cases.

#Well, well, well... it seems that our meta-analysis is robust and our results are reliable. Now that we have checked it, we'll examine our moderator variables to disentangle the complexity of these responses:

#
forest.meta(overall.mod, 
            sortvar = TE,
            leftlabs = c("Pair", "g", "SE"))

#Print 3000 x 3000 if it's necessary

#orchard plot

# We need to calculate the effect sizes, in this case d
H_family <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                   n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = H_family)

H_fam_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                   data = H_family)
summary(H_fam_MA)

model_results <- mod_results(H_fam_MA, mod = "1", at = NULL, data = H_family, group = "paper_code")
model_results

#ORCHARD PLOT
library(ggplot2)
orchaRd::orchard_plot(H_fam_MA, mod = "1", data = H_family, group = "paper_code",
                      xlab = "Standardised mean difference", 
                      fill=FALSE, 
                      transfm = "none")

orchard_pot(H_fam_MA, mod = "1", data = H_family, group = "paper_code",
            xlab = "Standardised mean difference", 
            fill=FALSE, 
            transfm = "none") + ylim(-2,4)

#600 x 220



orchaRd::caterpillars(H_fam_MA, mod = "1", data = H_family, group = "paper_code", 
                      xlab = "Standardised mean difference",
                      xlim(c(0,2)))

#forest(H_fam_MA, xlim=c(0,2))

#Ahora con management

H_fam_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                   data = H_family)
summary(H_fam_MA)

model_results <- mod_results(H_fam_MA, mod = "management", at = NULL, data = H_family, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(H_fam_MA, mod = "management", data = H_family, group = "paper_code", xlab = "Standardised mean difference",
            transfm = "none") + ylim(-2,4)


# 600 x 350

## Examining agric management effect ####

manage.mod<-update.meta(overall.mod,
                        subgroup = management)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


#type of sample and duration of the experiment
manage.mod<-update.meta(overall.mod,
                        subgroup = duration)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

#nothing significant

manage.mod<-update.meta(overall.mod,
                        subgroup = plant_man)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

#nada

manage.mod<-update.meta(overall.mod,
                        subgroup = climate)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))
#No effect of climate


#----------J_family----------
  J_family<-dat %>%
  dplyr::filter(response=="J_family")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = J_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)
#No significant REM
#Significant heterogeneity

dmetar::find.outliers(overall.mod)
#NO OUTLIERS
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour


summary.meta(overall.mod)


#### (2) Baujat plot 

#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")

J_family <- filter(J_family, pair != "26")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = J_family,
                      studlab = pair,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")

summary.meta(overall.mod)

#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero


#### (3) Correlation

#We will examine the correlation between the effect size and sample size using Kendall's approach.
cor.test(overall.mod$TE, J_family$Variance2, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, J_family$Variance2)
#There is no correlation

#### (4) Trim & Fill 
trimfill(overall.mod)

summary(overall.mod)
#After adding 4 studies, REM does not change

#
forest.meta(overall.mod, 
            sortvar = TE,
            leftlabs = c("Pair", "g", "SE"))

# Print 3000 x 3000


## Examining agric management effect ####

manage.mod<-update.meta(overall.mod,
                        subgroup = management)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

# Evennes was not modified by agric management

#type of sample and duration of the experiment
manage.mod<-update.meta(overall.mod,
                        subgroup = sample_duration)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

#SAMPLE_DURATION WAS NOT CRITICAL

manage.mod<-update.meta(overall.mod,
                        subgroup = climate)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

#CLIMATE WAS NOT CRITICAL

manage.mod<-update.meta(overall.mod,
                        subgroup = duration)
summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))

#NO EFFECT



#SECTION 2: BETA DIVERSITY ####

#Cambia la composición?

abund<-dat %>%
  dplyr::filter(cat_response == "abund")

#   Let's see it using a random-effects model:

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = abund,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

manage.mod<-update.meta(overall.mod,
                        subgroup = response)

summarymeta<-summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


# INTENTO DE GRAFICAR
#1
summarymeta <- as.data.frame(summarymeta)

library(meta)
library(ggplot2)

# Create forest plot
forest_plot <- ggplot(summarymeta, aes(x = TE, y = studlab, ymin = lower, ymax = upper, labels = abund$management)) +
  geom_pointrange(shape = 21, size = 0.1, fill = "darkgreen") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Study", y = "Hedges'D") +
  theme_bw()

forest_plot

#2
# Create a data frame
df <- data.frame(
  Category = c("A", "B"),
  Estimate = c(2.5, 3.2),
  Lower = c(1.2, 2.1),
  Upper = c(3.8, 4.3),
  vi = c(3.14, 5),
  sei = c(3.14,5)
)

# Generate forest plot
forest(
  x = summarymeta$TE,
  lower = summarymeta$lower,
  upper = summarymeta$upper,
  vi= 5,
  xlab = "Study",
  alim = c(0.5, 4.5),
  at = c(1, 2),
  labels = abund$management,
  col = "black",
  cex = 0.7,
  lwd = 2,
  pch = 16
)

