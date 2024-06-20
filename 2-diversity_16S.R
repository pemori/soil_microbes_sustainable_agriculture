# ---
# title: "Soil diversity 16s"
# authors: "Pedro Mondaca - Francisco E. Fonturbel"
# date: "07/nov/2022"
# ---

setwd("E:\\2022\\1-Meta-análisis microbioma suelos agrícolas")

#Loading data file
#Loading packages required
library(meta)
library(metafor)
library(readxl)
library(dplyr)
data<-read_xlsx("Meta16s.xlsx", sheet = "all")

#data $ sample_duration <- paste (data $ sample, data $ duration, sep = " - ")
dat<-data %>% 
  filter(Hedges_g > -10 & Hedges_g < 10)
dat<-filter(dat, Hedges_g != "#DIV/0!")
dat

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
str(dat)
#SECTION 1: ALPHA DIVERSITY ####


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

library(dmetar)

dmetar::find.outliers(overall.mod)
H_family <- filter(H_family, pair != "1" & pair != "2" & pair != "12" & pair != "108" & pair != "125" & pair != "130")

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
H_family <- filter(H_family, pair != "115")

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

#No outliers detected

#There is no **significant overall effect** on soil microbial diversity. 
#Also, there is a **significant heterogeneity** in the model (which is expected in this kind of models)
#Before examining the results using the `moderator variables` defined, we must assess the robustness of our data, using three tests:

#### (1) Funnel plot
funnel.meta(overall.mod, studlab = TRUE)

  cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

  H_family <- filter(H_family, pair != "19" & pair != "61" )
  
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
  
#heterogeneity marginally no significant
  
#### (2) Baujat plot 
  
#The Baujat plot allows us to visually explore heterogeneity in our meta-analysis, visualizing the individual contribution of each case study on overall heterogeneity and the influence in the result.
y <- dmetar::InfluenceAnalysis(overall.mod, random=TRUE)
plot(y, "baujat")  


#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero 

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
R_family <- filter(R_family, pair != "2" & pair != "49" & pair != "89" & pair != "129")

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
# ns REM nor heterogeneity


#Is there a significant bias on this? let's test it;
metabias(overall.mod, method.bias = "Pustejovsky")
# No bias, intercept is close to zero #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pub-bias.html


#### (3) Correlation

#We will examine the correlation between the effect size and sample size using Kendall's approach.
cor.test(overall.mod$TE, R_family$Variance2, method="kendall", continuity = TRUE, alternative = "two.sided")
plot(overall.mod$TE, R_family$Variance2)

#### (4) Trim & Fill 
trimfill(overall.mod)

summary(overall.mod)
#After adding 15 studies, REM is still not significant

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

# USE OF FERTILIZERS AND PESTICIDES HAVE SIGNIFICANT EFFECT ON BACTERIAL RICHNESS