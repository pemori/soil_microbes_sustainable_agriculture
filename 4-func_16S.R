setwd("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023")

#Loading data file
#Loading packages required
library(meta)
#library(metafor)
library(readxl)
library(dplyr)
library(metafor)
library(orchaRd)
dat<-read_xlsx("Meta16s.xlsx", sheet = "trimmed")
dat<-dat %>%
  filter(Hedges_g > -10 & Hedges_g < 10)
# dat<-filter(dat, Hedges_g != "#DIV/0!")
str(dat)
dat

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)


#functions according to abundance#####

funcA<-dat %>%
  dplyr::filter(cat_response == "funcA") 

funcab<-funcA

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

funcab <- filter(funcab, ID != "3287" & ID != "1139" & ID != "1075" & ID!="1079" & ID != "1479" & ID != "2807" & ID != "599" & ID!="3747" & ID != "3737" & ID != "3393" & ID != "2687" & ID!="2685" & ID != "1821" & ID != "269" & ID != "3219" & ID!="4033" & ID != "1711" & ID!="3795" & ID!="5009")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)
#REM significant, Q no, pero obvio

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

dmetar::find.outliers(overall.mod)
#no outlier

manage.mod<-update.meta(overall.mod,
                        subgroup = response)


summary.meta(manage.mod)

# We need to calculate the effect sizes, in this case d
funcab <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab)

#Ahora con responses
funcab_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                         data = funcab)
summary(funcab_MA)
#QM(df = 12) = 16.8813, p-val = 0.1541

model_results <- mod_results(funcab_MA, mod = "response", at = NULL, data = funcab, group = "paper_code")
model_results

#increase ureolysis, nitrification
#estimate 0.15, 0.21

#ORCHARD PLOT
library(ggplot2)
orchard_plot(funcab_MA, mod = "response", data = funcab, group = "paper_code", legend.pos = "top.right", xlab = "Standardised mean difference", angle = 0,
                      transfm = "none") + ylim(c(-2.5,5))

# 550 x 500


#fert
funcab_fert<-funcab %>%
  dplyr::filter(management == "F")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab_fert,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

dmetar::find.outliers(overall.mod)


manage.mod<-update.meta(overall.mod,
                        subgroup = response)


# We need to calculate the effect sizes, in this case d
funcab_fert <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                   n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab_fert)

#Ahora con responses

funcab_fert_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                   data = funcab_fert)
summary(funcab_fert_MA)
#QM(df = 12) = 8.3721, p-val = 0.7554

model_results <- mod_results(funcab_fert_MA, mod = "response", at = NULL, data = funcab_fert, group = "paper_code")
model_results

#increase chitinolysis, nitrification
# estimate = 0.30 y 0.35

#ORCHARD PLOT
orchard_plot(funcab_fert_MA, mod = "response", data = funcab_fert, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
                      transfm = "none") + ylim(c(-2.5,2.5))

# 500 x 550

funcab_till<-funcab %>%
  dplyr::filter(management == "T")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab_till,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

dmetar::find.outliers(overall.mod)


manage.mod<-update.meta(overall.mod,
                        subgroup = response)


# We need to calculate the effect sizes, in this case d
funcab_till <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab_till)

#Ahora con responses

funcab_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                         data = funcab_till)
summary(funcab_till_MA)

model_results <- mod_results(funcab_till_MA, mod = "response", at = NULL, data = funcab_till, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_plot(funcab_till_MA, mod = "response", data = funcab_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-2.5,2.5))

#decrease cellulo y xilano
#-0.36 y -0.33
# 500x550

#carbon
funcab_carb<-funcab %>%
  dplyr::filter(management == "C")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab_carb,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

dmetar::find.outliers(overall.mod)


manage.mod<-update.meta(overall.mod,
                        subgroup = response)


# We need to calculate the effect sizes, in this case d
funcab_carb <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab_carb)

#Ahora con responses

funcab_carb_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                         data = funcab_carb)
summary(funcab_carb_MA)

model_results <- mod_results(funcab_carb_MA, mod = "response", at = NULL, data = funcab_carb, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_plot(funcab_carb_MA, mod = "response", data = funcab_carb, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-2.5,2.5))

#no changes
#  400 x 550

#pest
funcab_pest<-funcab %>%
  dplyr::filter(management == "P")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcab_pest,
                      studlab = ID,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

funnel.meta(overall.mod, studlab = TRUE)

cc <- funnel(overall.mod, common = TRUE, studlab = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour

dmetar::find.outliers(overall.mod)


manage.mod<-update.meta(overall.mod,
                        subgroup = response)


# We need to calculate the effect sizes, in this case d
funcab_pest <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab_pest)

#Ahora con responses

funcab_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                         data = funcab_pest)
summary(funcab_pest_MA)

model_results <- mod_results(funcab_pest_MA, mod = "response", at = NULL, data = funcab_pest, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_plot(funcab_pest_MA, mod = "response", data = funcab_pest, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-2.5,2.5))

#no changes

