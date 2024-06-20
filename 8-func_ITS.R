# 16s fert= QM(df = 12) = 8.3721, p-val = 0.7554
# 16s tilla= QM(df = 12) = 14.1775, p-val = 0.2895
# 16 carb= QM(df = 12) = 5.8331, p-val = 0.9243
# 16 pest= QM(df = 12) = 5.5220, p-val = 0.9382
# Its= QM(df = 20) = 21.5332, p-val = 0.3664
# Its fert= QM(df = 20) = 21.0254, p-val = 0.3956
# Its till= QM(df = 20) = 15.7998, p-val = 0.7290
# Its carb= QM(df = 20) = 11.6462, p-val = 0.9277
# Its pest= QM(df = 20) = 12.5483, p-val = 0.8959




setwd("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023")

#Loading data file
#Loading packages required
library(meta)
#library(metafor)
library(readxl)
library(dplyr)
library(metafor)
library(orchaRd)
dat<-read_xlsx("MetaITS.xlsx", sheet = "Sheet1")

#data $ sample_duration <- paste (data $ sample, data $ duration, sep = " - ")
 dat<-data %>% 
   filter(Hedges_g > -10 & Hedges_g < 10)
 dat<-filter(dat, Hedges_g != "#DIV/0!")
str(dat)
dat

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)

#func según Richness####
funcR<-dat %>%
  dplyr::filter(cat_response == "funcR")
  #dplyr::filter(response != "R_epiphyte" & response !="R_nectar_tap_saprotroph" & response !="R_pollen_saprotroph" & response !="R_unspecified_saprotroph" & response != "R_unspecified_symbiotroph")




# We need to calculate the effect sizes, in this case d
funcR <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcR)

#Ahora con responses

funcR_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                   data = funcR)
summary(funcR_MA)

library(orchaRd)
model_results <- orchaRd::mod_results(funcR_MA, mod = "response", at = NULL, data = funcR, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcR_MA, mod = "response", data = funcR, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

#Orden
# R_sooty_mold
# R_arthropod_parasite
# R_wood_saprotroph
# R_root_endophyte
# R_ectomycorrhizal
# R_protistan_parasite
# R_arthropod_associated
# R_litter_saprotroph
# R_foliar_endophyte
# R_soil_saprotroph
# R_dung_saprotroph
# R_mycoparasite
# R_plant_pathogen
# R_root_associated
# R_arbuscular_mycorrhizal
# R_fungal_decomposer
#FERTILIZATION CONTRAST #####

# Para fertilization
funcR_fert<-funcR %>%
  dplyr::filter(management == "F")

funcR_fert_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = funcR_fert)
summary(funcR_fert_MA)

model_results <- mod_results(funcR_fert_MA, mod = "response", at = NULL, data = funcR, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcR_fert_MA, mod = "response", data = funcR, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

# Para carbon
funcR_carb<-funcR %>%
  dplyr::filter(management == "C")

funcR_carb_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = funcR_carb)
summary(funcR_carb_MA)

model_results <- mod_results(funcR_carb_MA, mod = "response", at = NULL, data = funcR, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcR_carb_MA, mod = "response", data = funcR, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

# Para tillage
funcR_till<-funcR %>%
  dplyr::filter(management == "T")

funcR_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = funcR_till)
summary(funcR_till_MA)

model_results <- mod_results(funcR_till_MA, mod = "response", at = NULL, data = funcR, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcR_till_MA, mod = "response", data = funcR, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")



# Para pesticides
funcR_pest<-funcR %>%
  dplyr::filter(management == "P")

funcR_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = funcR_pest)
summary(funcR_pest_MA)

model_results <- mod_results(funcR_pest_MA, mod = "response", at = NULL, data = funcR, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcR_pest_MA, mod = "response", data = funcR, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")


# NINGUNA DIFERENCIA.. PROBARÉ CON ABUND

#func según abundancia#####

funcA<-dat %>%
  dplyr::filter(cat_response == "funcA") 

# funcA[funcA == 0] <- 0.0000000001
# min(funcA$Xc)
# min(funcA$Xe)
# dim(funcA)
# funcA2<-funcA %>%
#   mutate(CW = case_when(Xe == 0.0000000001 & Xc == 0.0000000001 ~ "R"))
# 
# funcA2$CW <- replace(funcA2$CW,is.na(funcA2$CW),0)
# 
# funcA3<-funcA2 %>%
#   filter(CW != "R")
# 
# funcab<-funcA3
funcab<-funcA

#Ver funcab primero...
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

funcab <- filter(funcab, ID != "365" & ID != "59" & ID !="2339" & ID !="2375")

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

dmetar::find.outliers(overall.mod)


manage.mod<-update.meta(overall.mod,
                        subgroup = response)


summary.meta(manage.mod)


#abund_plant_pathogen is decreased. Increase unspecified symbiotroph

# We need to calculate the effect sizes, in this case d
funcab <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = funcab)

#Ahora con responses

funcab_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                         data = funcab)
summary(funcab_MA)
#decrease plant pathogen

model_results <- mod_results(funcab_MA, mod = "response", at = NULL, data = funcab, group = "paper_code")
model_results

#ORCHARD PLOT
library(ggplot2)
orchard_plot(funcab_MA, mod = "response", data = funcab, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
                      transfm = "none") + ylim(c(-2.5,5))

#increase of unspecified symbiotroph
#estimate = 0.30
# 500 x 750


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

model_results <- mod_results(funcab_fert_MA, mod = "response", at = NULL, data = funcab_fert, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(funcab_fert_MA, mod = "response", data = funcab_fert, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
                      transfm = "none") + ylim(c(-2.5,2.5))

#decrease of protistan parasite group
# 400 x 645

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
orchard_pot(funcab_till_MA, mod = "response", data = funcab_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-2.5,2.5))
# 400 x 645
#no changes


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
orchard_pot(funcab_carb_MA, mod = "response", data = funcab_carb, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-2.5,2.5))

#no changes
#  400 x 645

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
orchard_pot(funcab_pest_MA, mod = "response", data = funcab_pest, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
              transfm = "none") + ylim(c(-2.5,2.5))

#no changes
#  400 x 645


























































# CARBON MANAGEMENT ####


funcab_carb<-funcab %>%
  dplyr::filter(management == "C")
# #func[func == 0] <- 0.0000000001
# min(func$Xc)
# min(func$Xe)
# dim(func)
# func2<-func %>% 
#   mutate(CW = case_when(Xe == 0.0000000001 & Xc == 0.0000000001 ~ "R"))
# 
# func2$CW <- replace(func2$CW,is.na(func2$CW),0)
# 
# func3<-func2 %>%
#   filter(CW != "R") 

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

funcab_carb <- filter(funcab_carb, ID != "3264" & ID != "3459")

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

summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


# ninguna función cambia!
# q me dice func de shannon en vez de richness?

funcA<-dat %>%
  dplyr::filter(cat_response == "funcA")


funcA[funcA == 0] <- 0.0000000001
min(funcA$Xc)
min(funcA$Xe)
dim(funcA)
funcA2<-funcA %>%
  mutate(CW = case_when(Xe == 0.0000000001 & Xc == 0.0000000001 ~ "R"))

funcA2$CW <- replace(funcA2$CW,is.na(funcA2$CW),0)

funcA3<-funcA2 %>%
  filter(CW != "R")

funcA_carb<-funcA3 %>%
  dplyr::filter(management == "C")


overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcA_carb,
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

#Se reducirían los root endophyte

# TILLAGE ####

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

funcab_till <- filter(funcab_till, ID != "1755" & ID != "3630" & ID != "183" & ID != "1740" & ID != "3558" & ID != "1953")

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

summary.meta(manage.mod)
forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


# cambia arbuscular mycorrhizal y arthropod associated
# q me dice func de shannon en vez de richness?

funcA_till<-funcA3 %>%
  dplyr::filter(management == "T")


overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcA_till,
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

#ningún cambio

#PESTICIDE CONTRAST #####

funcR_pest<-funcR %>%
  dplyr::filter(management == "P")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcR_pest,
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

# increase arthropod parasite and soil saprotroph

# q me dice func de shannon en vez de richness?

funcA_pest<-funcA3 %>%
  dplyr::filter(management == "P")


overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcA_pest,
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

funcA_pest <- filter(funcA_pest, ID != "2339" & ID != "2375")

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = funcA_pest,
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

#increase nectar tap saprotroph and unspecified symbiotroph pero decrease wood saprotroph


