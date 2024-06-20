setwd("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023")

#Loading data file
#Loading packages required
library(meta)
library(metafor)
library(readxl)
library(dplyr)
dat<-read_xlsx("MetaITS.xlsx", sheet = "Sheet1")

dat<-dat %>% 
  filter(Hedges_g > -10 & Hedges_g < 10)
dat<-filter(dat, Hedges_g != "#DIV/0!")

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)


#SECTION 2: BETA DIVERSITY ####

abund<-dat %>%
  dplyr::filter(cat_response == "abund") %>%
    dplyr::filter(response != "Classiculomycetes" & response != "Collemopsidiomycetes" & response != "Endogonomycetes" &
                response != "Entorrhizomycetes" & response != "Exobasidiomycetes" & response != "Harpellomycetes" &
                  response != "Laboulbeniomycetes" & response !=  "Lichinomycetes" & response != "Monoblepharidomycetes" &
                   response != "Mucoromycetes" & response != "Pneumocystidomycetes" & response != "Wallemiomycetes ")

compA<-abund
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

summary.meta(manage.mod)

forest.meta(manage.mod,
            random = TRUE,
            layout = "subgroup",
            overall = FALSE,
            prediction = FALSE,
            hetstat = FALSE,
            xlim=c(-1,2.5))


# We need to calculate the effect sizes, in this case d
abund <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                     n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = abund)

#Ahora con responses

abund_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = abund)
summary(abund_MA)

library(orchaRd)
model_results <- orchaRd::mod_results(abund_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
model_results

abund_MA 

# x <- c("Leotiomycetes","Pezizomycetes","Agaricomycetes","Tremellomycetes","Glomeromycetes","Dothideomycetes","Saccharomycetes","Sordariomycetes","Pucciniomycetes","Ustilaginomycetes","Chytridiomycetes","Wallemiomycetes","Eurotiomycetes")
# x

#ORCHARD PLOT
orchard_pot(abund_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
                      transfm = "none")

#Print
# Para fertilization
abund_fert<-abund %>%
  dplyr::filter(management == "F")

abund_fert_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                   data = abund_fert)
summary(abund_fert_MA)

model_results <- mod_results(abund_fert_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(abund_fert_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

#Lower Tremellomycetes
#estimate= -0.70

# Para carbon
abund_carb<-abund %>%
  dplyr::filter(management == "C")

abund_carb_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = abund_carb)
summary(abund_carb_MA)

model_results <- mod_results(abund_carb_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
model_results

# Lower Glomeromycetes an LEotyiomicetes. Higher Sordariomycetes
# -0.55, -0.67. 0.49

#ORCHARD PLOT
orchard_pot(abund_carb_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

#ME FALTA ESTANDARIZAR LA SMD A 4

# Para tillage
abund_till<-compA %>%
  dplyr::filter(management == "T")

abund_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = abund_till)
summary(abund_till_MA)

model_results <- mod_results(abund_till_MA, mod = "response", at = NULL, data = compA, group = "paper_code")
model_results

# Lower Glomeromycetes

#ORCHARD PLOT
orchard_plot(abund_till_MA, mod = "response", data = abund_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")

# Para pesticide
abund_pest<-abund %>%
  dplyr::filter(management == "P")

abund_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = abund_pest)
summary(abund_pest_MA)

model_results <- mod_results(abund_pest_MA, mod = "response", at = NULL, data = abund, group = "paper_code")
model_results



#ORCHARD PLOT
orchard_pot(abund_pest_MA, mod = "response", data = abund, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
             transfm = "none")


#No se puede graficar para T y P

#Ahora para abundance, instead richness

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

library(orchaRd)
model_results <- orchaRd::mod_results(compA_MA, mod = "response", at = NULL, data = compA, group = "paper_code")
model_results

compA_MA 

# x <- c("Leotiomycetes","Pezizomycetes","Agaricomycetes","Tremellomycetes","Glomeromycetes","Dothideomycetes","Saccharomycetes","Sordariomycetes","Pucciniomycetes","Ustilaginomycetes","Chytridiomycetes","Wallemiomycetes","Eurotiomycetes")
# x

#ORCHARD PLOT
orchard_pot(compA_MA, mod = "response", data = compA, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,9)) 
#printed 500 x 450
#nothing significant

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


#Lower Tremellomycetes

#printed 400 x 500

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

# Lower Glomeromycetes, lower leotiomycetes, higher sordariomycetes

# Para till
compA_till<-abund %>%
  dplyr::filter(management == "T")

compA_till_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_till)
summary(compA_till_MA)

model_results <- mod_results(compA_till_MA, mod = "response", at = NULL, data = compA_till, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_pot(compA_till_MA, mod = "response", data = compA_till, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none") + ylim(c(-5,5))

# Higher Leotiomycetes
#estimate= 0.62

# Para pest
compA_pest<-abund %>%
  dplyr::filter(management == "P")

compA_pest_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~response ,random = list(~1 | ID),
                        data = compA_pest)
summary(compA_pest_MA)

model_results <- mod_results(compA_pest_MA, mod = "response", at = NULL, data = compA_pest, group = "paper_code")
model_results

#ORCHARD PLOT
orchard_plot(compA_pest_MA, mod = "response", data = compA_pest, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")# + ylim(c(-5,5))

#Lower Dothideomycete
#higher chytridiomycetes and sordariomycetes
#-0.79, 0.70 y 0.63
