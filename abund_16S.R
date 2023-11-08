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
dat<-filter(dat, Hedges_g != "#DIV/0!")


#SECTION 2: BETA DIVERSITY ####

#Cambia la composición?

abund<-dat %>%
  dplyr::filter(cat_response == "abundance")
   # dplyr::filter(response != "Classiculomycetes" & response != "Collemopsidiomycetes" & response != "Endogonomycetes" &
                # response != "Entorrhizomycetes" & response != "Exobasidiomycetes" & response != "Harpellomycetes" &
                #   response != "Laboulbeniomycetes" & response !=  "Lichinomycetes" & response != "Monoblepharidomycetes" &
                #    response != "Mucoromycetes" & response != "Pneumocystidomycetes" & response != "Wallemiomycetes ")


abundA<-abund
abundA[abundA == 0] <- 0.0000000001
min(abundA$Xc)
min(abundA$Xe)
dim(abundA)
abundA2<-abundA %>%
  mutate(CW = case_when(Xe == 0.0000000001 & Xc == 0.0000000001 ~ "R"))

abundA2$CW <- replace(abundA2$CW,is.na(abundA2$CW),0)

abundA3<-abundA2 %>%
  filter(CW != "R")

#compA<-abundA3
compA<-abund

#based on richness####
#   Let's see it using a random-effects model:

# overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
#                       data = abund,
#                       fixed = FALSE,
#                       random = TRUE,
#                       method.tau = "REML",
#                       hakn = TRUE,
#                       prediction = TRUE,
#                       sm = "SMD",
#                       method.smd = "Hedges")
# summary.meta(overall.mod)
# 
# manage.mod<-update.meta(overall.mod,
#                         subgroup = response)
# 
# summary.meta(manage.mod)
# 
# forest.meta(manage.mod,
#             random = TRUE,
#             layout = "subgroup",
#             overall = FALSE,
#             prediction = FALSE,
#             hetstat = FALSE,
#             xlim=c(-1,2.5))
# 
# 
# # We need to calculate the effect sizes, in this case d
# abund <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
#                      n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = abund)
# 
# #Ahora con responses
# 
# abund_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
#                         data = abund)
# summary(abund_MA)
# 
# library(orchaRd)
# model_results <- orchaRd::mod_results(abund_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
# model_results
# 
# abund_MA 
# 
# # x <- c("Leotiomycetes","Pezizomycetes","Agaricomycetes","Tremellomycetes","Glomeromycetes","Dothideomycetes","Saccharomycetes","Sordariomycetes","Pucciniomycetes","Ustilaginomycetes","Chytridiomycetes","Wallemiomycetes","Eurotiomycetes")
# # x
# 
# #ORCHARD PLOT
# orchard_pot(abund_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
#                       transfm = "none")
# 
# #Print
# # Para fertilization
# abund_fert<-abund %>%
#   dplyr::filter(management == "F")
# 
# abund_fert_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
#                    data = abund_fert)
# summary(abund_fert_MA)
# 
# model_results <- mod_results(abund_fert_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
# model_results
# 
# #ORCHARD PLOT
# orchard_pot(abund_fert_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
#             transfm = "none")
# 
# #Lower Tremellomycetes
# 
# # Para carbon
# abund_carb<-abund %>%
#   dplyr::filter(management == "C")
# 
# abund_carb_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
#                         data = abund_carb)
# summary(abund_carb_MA)
# 
# model_results <- mod_results(abund_carb_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
# model_results
# 
# # Lower Glomeromycetes
# 
# #ORCHARD PLOT
# orchard_pot(abund_carb_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
#             transfm = "none")
# 
# #ME FALTA ESTANDARIZAR LA SMD A 4
# 
# # Para tillage
# abund_till<-abund %>%
#   dplyr::filter(management == "T")
# 
# abund_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
#                         data = abund_till)
# summary(abund_till_MA)
# 
# model_results <- mod_results(abund_till_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
# model_results
# 
# # Lower Glomeromycetes
# 
# #ORCHARD PLOT
# orchard_plot(abund_till_MA, mod = "response", data = abund_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
#             transfm = "none")
# 
# # Para pesticide
# abund_pest<-abund %>%
#   dplyr::filter(management == "P")
# 
# abund_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
#                         data = abund_pest)
# summary(abund_pest_MA)
# 
# model_results <- mod_results(abund_pest_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
# model_results
# 
# 
# 
# #ORCHARD PLOT
# orchard_pot(abund_pest_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
#              transfm = "none")
# 
# 
# #No se puede graficar para T y P # si se puede al transformar data
###### Este queda #####
#Ahora para abundance, instead richness # el otro no lo modifiqué, al menos no por management.. solo de aquí en adelante

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = compA,
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


# We need to calculate the effect sizes, in this case d
compA <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = compA)

#Ahora con responses

compA_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                   data = compA)
summary(compA_MA)
#Test of Moderators: QM(df = 15) = 57.4497, p-val < .0001

library(orchaRd)
model_results <- orchaRd::mod_results(compA_MA, mod = "response", at = NULL, data = compA, group = "paper_code")
model_results
 

#ORCHARD PLOT
orchard_pot(compA_MA, mod = "response", data = compA, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,10)) 
#printed 500 x 480
#lower gemma y chloroflexi, higher spiro, fibro y bacteroide
#estimate = -0.45; ±95% CIs = -0.65, -0.23 and -0.35; ±95% CIs = -0.55, -0.14
#estimate = 0.29; ±95% CIs = 0.025, 0.55; 0.31; ±95% CIs = 0.049, 0.566; 0.32; ±95% CIs = 0.11, 0.52

# Para fertilization
compA_fert<-compA %>%
  dplyr::filter(management == "F")

compA_fert_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_fert)
summary(compA_fert_MA)

model_results <- mod_results(compA_fert_MA, mod = "response", at = NULL, data = compA, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(compA_fert_MA, mod = "response", data = compA_fert, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,5)) 


#Lower Gemma and chloro
#higher firmi, cyano, spiro, fibro and bactero
#estimate -0.73 y -0.56; 0.48, 0.72, 0.54, 0.63 y 0.59

#printed 400 x 522

# Para carbon
compA_carb<-compA %>%
  dplyr::filter(management == "C")

compA_carb_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_carb)
summary(compA_carb_MA)

model_results <- mod_results(compA_carb_MA, mod = "response", at = NULL, data = compA_carb, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(compA_carb_MA, mod = "response", data = compA_carb, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,5)) 

# higher bacteroid fibro y spiro, lower plancto
# estimate=0.30;0.51 y 0.48. plancto=-0.39
#printed 400 x 522

# Para till
compA_till<-compA %>%
  dplyr::filter(management == "T")

compA_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_till)
summary(compA_till_MA)

model_results <- mod_results(compA_till_MA, mod = "response", at = NULL, data = compA_till, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(compA_till_MA, mod = "response", data = compA_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,5))

# Lower acidobacteria, proteobacteria, chloroflexi y gemmatimonadetes. higher actino
# -0.53, -0.53, -0.56,-0.72. 0.73
#printed 400 x 522


# Para pest
compA_pest<-compA %>%
  dplyr::filter(management == "P")

compA_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_pest)
summary(compA_pest_MA)

model_results <- mod_results(compA_pest_MA, mod = "response", at = NULL, data = compA_pest, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(compA_pest_MA, mod = "response", data = compA_pest, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,5))

# ns
#printed 400 x 522



# resultador por phylum
# Gemmatimonadetes
compA_Gemm<-compA %>%
  dplyr::filter(response == "16-Gemmatimonadetes")

compA_Gemm_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                        data = compA_Gemm)
compA_Gemm_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                        data = compA_Gemm)
summary(compA_Gemm_MA)

library(orchaRd)
model_results <- mod_results(compA_Gemm_MA, mod = "management", at = NULL, data = compA, group = "paper_code")
model_results

# #
# Test for Residual Heterogeneity:
#   QE(df = 82) = 158.0397, p-val < .0001
# 
# Test of Moderators (coefficients 2:4):
#   QM(df = 3) = 4.9307, p-val = 0.1769


