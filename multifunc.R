#Boxplot
library(readxl)
setwd("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\Nuevos")
library(ggplot2)
data_16<-read_xlsx("MF_16s.xlsx", sheet = "MF")
head(data_16)
dat<-data_16
#Loading data file
#Loading packages required
library(meta)
library(metafor)
library(readxl)
library(dplyr)


# dat<-dat %>% 
#   filter(Hedges_g > -10 & Hedges_g < 10)
# dat<-filter(dat, Hedges_g != "#DIV/0!")

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)


#SECTION 4: MF ####


#   Let's see it using a random-effects model:

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = dat,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

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


# We need to calculate the effect sizes, in this case d
dat <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
                n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = dat)

dat_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                   data = dat)
summary(dat_MA)
model_results <- mod_results(dat_MA, mod = "1", at = NULL, data = dat, group = "paper_code")
model_results

orchard_pot(dat_MA, mod = "1", data = dat, group = "paper_code",
            xlab = "Standardised mean difference", 
            fill=TRUE, 
            transfm = "none") + ylim(-2,2)

#600 x 220

library(orchaRd)

#Ahora con responses

dat_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                   data = dat)
model_results <- orchaRd::mod_results(dat_MA, mod = "management", at = NULL, data = dat, group = "paper_code")
model_results
summary(dat_MA)

orchard_pot(dat_MA, mod = "management", data = dat, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")+ ylim(-2,2)

# 600 x 350

#ahora fungi

data_ITS<-read_xlsx("MF_ITS.xlsx", sheet = "MF")
head(data_ITS)
dat<-data_ITS
#Loading data file
#Loading packages required

# dat<-dat %>% 
#   filter(Hedges_g > -10 & Hedges_g < 10)
# dat<-filter(dat, Hedges_g != "#DIV/0!")

dat$pair<-as.character(dat$pair) 
dat$paper_code<-as.character(dat$paper_code)
dat$ID<-as.character(dat$ID)
str(dat)


#SECTION 4: MF ####


#   Let's see it using a random-effects model:

overall.mod<-metacont(Ne, Xe, Se, Nc, Xc, Sc, 
                      data = dat,
                      fixed = FALSE,
                      random = TRUE,
                      method.tau = "REML",
                      hakn = TRUE,
                      prediction = TRUE,
                      sm = "SMD",
                      method.smd = "Hedges")
summary.meta(overall.mod)

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


# We need to calculate the effect sizes, in this case d
dat <- escalc(measure = "SMD", n1i = Ne, sd1i = Se, m1i = Xe,
              n2i = Nc, sd2i = Sc, m2i = Xc, var.names = c("SMD", "vSMD"), data = dat)
dat_MA <- rma.mv(yi = SMD, V = vSMD, random = list(~1 | ID),
                   data = dat)
summary(dat_MA)

model_results <- orchaRd::mod_results(dat_MA, mod = "1", at = NULL, data = dat, group = "paper_code")
model_results

orchard_pot(dat_MA, mod = "1", data = dat, group = "paper_code",
            xlab = "Standardised mean difference", 
            fill=TRUE, 
            transfm = "none") + ylim(-2,2)

#600 x 220

library(orchaRd)
model_results <- orchaRd::mod_results(dat_MA, mod = "management", at = NULL, data = dat, group = "paper_code")
model_results
#Ahora con responses

dat_MA <- rma.mv(yi = SMD, V = vSMD, mods = ~management ,random = list(~1 | ID),
                 data = dat)
summary(dat_MA)

orchard_pot(dat_MA, mod = "management", data = dat, group = "paper_code", xlab = "Standardised mean difference", angle = 0,
            transfm = "none")+ ylim(-2,2)

# 600 x 350


