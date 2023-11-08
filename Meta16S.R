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
dat<-read_xlsx("Meta16s.xlsx", sheet = "trimmed")

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)
dat<-dat %>%
  filter(Hedges_g > -10 & Hedges_g < 10)
#   filter(Sc !=0 &  Se !=0 & Xc !=0 & Xe !=0 & Hedges_g != "#DIV/0!")

dat $ sample_duration <- paste (dat $ sample, dat $ duration, sep = " - ")
dat

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
#Significant REM
#Significant heterogeneity


dmetar::find.outliers(overall.mod)
R_family <- filter(R_family, pair != "2" & pair != "129")

#Second model after outliers remotion
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
#Significant REM
#No significant heterogeneity

dmetar::find.outliers(overall.mod)

#There is no **significant overall effect** on soil microbial diversity. There is **significant heterogeneity** in the model.
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

R_family <- filter(R_family, pair != "120" & pair != "108" & pair !="106")

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
# significant REM, ns heterogeneity

#### (2) Baujat plot 

#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
# y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
# plot(y, "baujat")
# 
# R_family <- filter(R_family, pair != "91")
# 
# overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
#                       data = R_family,
#                       studlab = pair,
#                       fixed = FALSE,
#                       random = TRUE,
#                       method.tau = "REML",
#                       hakn = TRUE,
#                       prediction = TRUE,
#                       sm = "SMD",
#                       method.smd = "Hedges")
# 
# summary.meta(overall.mod)
# 
# #significant REM

#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pub-bias.html


#### (3) Correlation

#We will examine the correlation between the effect size and sample size using Kendall's approach.
cor.test(overall.mod$TE, R_family$Variance2, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, R_family$Variance2)
#There is a correlation


#### (4) Trim & Fill 
trimfill(overall.mod)

summary(overall.mod)
#After adding 14 studies, REM is still not significant

#
forest.meta(overall.mod, 
            sortvar = TE,
            leftlabs = c("Pair", "g", "SE"))

# Print 3000 x 3000

#Sobre SMD: https://wviechtb.github.io/metafor/reference/escalc.html
# We need to calculate the effect sizes, in this case d
R_family <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                   n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = R_family)

#Results Test for Heterogeneity: Q(df = 85) = 121.0843, p-val = 0.0062

R_fam_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                   data = R_family)
summary(R_fam_MA)

library(orchaRd)
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
            fill=TRUE, 
            transfm = "none") + ylim(-2,4)

#600 x 220

#Ahora con management

R_fam_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                   data = R_family)
summary(R_fam_MA)
#Test for Residual Heterogeneity: QM(df = 3) = 13.2815, p-val = 0.0041


model_results <- mod_results(R_fam_MA, mod = "management", at = NULL, data = R_family, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(R_fam_MA, mod = "management", data = R_family, group = "paper_code", xlab = "Standardised mean difference",
            transfm = "none") + ylim(-2,4)


# 600 x 350

# ## Examining agric management effect ####
# 
# manage.mod<-update.meta(overall.mod,
#                         subgroup = management)
# summary.meta(manage.mod)
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# # USE OF FERTILIZERS AND PESTICIDES HAVE SIGNIFICANT EFFECT ON BACTERIAL RICHNESS
# 
# #type of sample and duration of the experiment
# manage.mod<-update.meta(overall.mod,
#                         subgroup = sample_duration)
# summary.meta(manage.mod)
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# #SST WAS SIGNIFICANT POSSIBLY DUE TO THE EFFECT OF PESTICIDES
# 
# manage.mod<-update.meta(overall.mod,
#                         subgroup = climate)
# summary.meta(manage.mod)
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# #POSITIVE RESPONSES WERE SHOWN IN WARM FARMLANDS
# 
# manage.mod<-update.meta(overall.mod,
#                         subgroup = duration)
# summary.meta(manage.mod)
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# #ONLY SHORT TERM
# 
# manage.mod<-update.meta(overall.mod,
#                         subgroup = plant_man)
# summary.meta(manage.mod)
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# #ONLY SHORT TERM
# 




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
#Significant heterogeneity


#Interpretación de heterogeneidad
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
H_family <- filter(H_family, pair != "1" & pair != "2" & pair != "108" & pair != "125" )

#Second model after outliers remotion
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
#Significant heterogeneity

dmetar::find.outliers(overall.mod)
H_family <- filter(H_family, pair != "130" & pair !="115")

#Second model after outliers remotion
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
#No Significant heterogeneity

dmetar::find.outliers(overall.mod)
#No outliers detected



#There is no **significant overall effect** on soil microbial diversity. 
#Also, there is a **significant heterogeneity** in the model (which is expected in this kind of models)
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:



#### (1) Funnel plot
funnel.meta(overall.mod, studlab = TRUE)

  cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

  H_family <- filter(H_family, pair != "61" )
  
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
  
# ns heterogenity nor REM
  
#### (2) Baujat plot 
  
#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")  

#H_family <- filter(H_family, pair != "91" & pair != "89" & pair !="3")

# overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
#                       data = H_family,
#                       studlab = pair,
#                       fixed = FALSE,
#                       random = TRUE,
#                       method.tau = "REML",
#                       hakn = TRUE,
#                       prediction = TRUE,
#                       sm = "SMD",
#                       method.smd = "Hedges")
# 
# summary.meta(overall.mod)
#ns REM nor Q 
    
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

# We need to calculate the effect sizes, in this case d
H_family <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                   n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = H_family)

H_fam_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                   data = H_family)
summary(H_fam_MA)
#Results Test for Heterogeneity: Q(df = 82) = 133.4456, p-val = 0.0003

library(orchaRd)
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
            fill=TRUE, 
            transfm = "none") + ylim(-2,4)

#600 x 220

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


orchard_plot(H_fam_MA, mod = "management", data = H_family, group = "paper_code", xlab = "Standardised mean difference",
            transfm = "none") + ylim(-2,4)

#hasta aquí llego



## Examining agricultural management effect ####

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

# NO MANAGEMENT HAD SIGNIFICANT EFFECT

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

#RLT IS SIGNIFICANT BUT ONLY 4 OBSERVATIONS

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

#NO CLIMATE TYPE HAS EFFECT ON THE RESPONSE

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
J_family <- filter(J_family, pair != "1" & pair != "2" & pair != "19" & pair != "67" & pair !="115")

#Second model after outliers remotion
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
#No significant heterogeneity

dmetar::find.outliers(overall.mod)

dmetar::find.outliers(overall.mod)
J_family <- filter(J_family, pair != "12")

#Second model after outliers remotion
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
#No significant REM (marginal)
#No significant heterogeneity

dmetar::find.outliers(overall.mod)

#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

J_family <- filter(J_family, pair != "61" & pair != "6")

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


#### (2) Baujat plot 

#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")

J_family <- filter(J_family, pair != "89" & pair !="87" & pair !="3" & pair !="91")

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

#Cambia la composici?n?

# abund<-data %>%
#   dplyr::filter(cat_response == "abund")

abund<-dat %>%
  dplyr::filter(cat_response == "abund")

abund[abund == 0] <- 0.0000000001
min(abund$Xc)
min(abund$Xe)
dim(abund)

abund2<-abund %>% 
  mutate(oli = case_when(Xe == 0.0000000001 & Xc == 0.0000000001 ~ 0))

dim(abund2)

abund3<-abund2 %>%
  filter(Xe != 0 & Xc != 0) 

dim(abund3)

Chequear



overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = abund2,
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

summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


abund_cyan<-abund2 %>%
  dplyr::filter(response == "Cyanobacteria")
abund_cyan2<-dat %>%
  dplyr::filter(response == "Cyanobacteria")
abund_cyan3<-abund %>%
  dplyr::filter(response == "Cyanobacteria")

Descubrir xq se me van los valores de Cyanobacteria al usar el filter!

dim(abund_cyan)
dim(abund_cyan2)
dim(abund_cyan3)
