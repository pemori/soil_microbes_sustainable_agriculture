# Analysis of 16s OTU or ASV datasets
# Author: Pedro Mondaca (pedro.mondaca@pucv.cl)

library(vegan)
library(readxl)
library(tidyverse)
data <- read_excel("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\16S\\80_16S.xlsx", 
                   sheet = "Hoja1", na = "0")
data[is.na(data)] <- 0
data<-data[-c(1)] 

nrow(data) #should equals 563
df0 <- data
colSums(df0) # must equals 100!

df0 %>% 
  mutate(across(everything(), function(x) {
    x*100/sum(x)
  }))
colSums(df0)

df <- as.data.frame(t(df0))

pm <- read_excel("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\pm_16s.xlsx", 
                 sheet = "Final", col_types = c("text","text","text","text","text", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                 "numeric"))

colnames(df)<-pm$family #family to factor

#alpha diversity
H_family <- diversity(df, index="shannon") #Shannon index using vegan package
R_family <- specnumber(df) # Richness

tab_family <- data.frame(H_family,R_family)
t_family <- t(tab_family)

#community analysis
df_phylum0 <- aggregate(df0, by = list(pm$phylum), FUN = sum)
df_phylum <- as.data.frame(t(df_phylum0[-c(1)]))
colnames(df_phylum)<-df_phylum0$Group.1

#Phyla abundance
#Abditibacteriota<-df_phylum$Abditibacteriota
Acidobacteria<-df_phylum$Acidobacteria
Actinobacteria<-df_phylum$Actinobacteria
#Aquificae<-df_phylum$Aquificae
Armatimonadetes<-df_phylum$Armatimonadetes
Bacteroidetes<-df_phylum$Bacteroidetes
#Balneolaeota<-df_phylum$Balneolaeota
#C.Cryosericota<-df_phylum$C.Cryosericota
#C.Krumholzibacteriota<-df_phylum$C.Krumholzibacteriota
#Caldiserica<-df_phylum$Caldiserica
#Calditrichaeota<-df_phylum$Calditrichaeota
Chlamydiae<-df_phylum$Chlamydiae
#Chlorobi<-df_phylum$Chlorobi
Chloroflexi<-df_phylum$Chloroflexi
#Chrysiogenetes<-df_phylum$Chrysiogenetes
#Coprothermobacterota<-df_phylum$Coprothermobacterota
Cyanobacteria<-df_phylum$Cyanobacteria
#Deferribacteres<-df_phylum$Deferribacteres
#DeinococcusThermus<-df_phylum$DeinococcusThermus
#Dictyoglomi<-df_phylum$Dictyoglomi
#Elusimicrobia<-df_phylum$Elusimicrobia
Fibrobacteres<-df_phylum$Fibrobacteres
Firmicutes<-df_phylum$Firmicutes
#Fusobacteria<-df_phylum$Fusobacteria
Gemmatimonadetes<-df_phylum$Gemmatimonadetes
#Ignavibacteriae<-df_phylum$Ignavibacteriae
#iritimatiellaeota<-df_phylum$Kiritimatiellaeota
#Lentispharae<-df_phylum$Lentispharae
Nitrospinae<-df_phylum$Nitrospinae
Planctomycetes<-df_phylum$Planctomycetes
Proteobacteria<-df_phylum$Proteobacteria
#Rhodothermaeota<-df_phylum$Rhodothermaeota
Spirochaetes<-df_phylum$Spirochaetes
#Synergistetes<-df_phylum$Synergistetes
Tenericutes<-df_phylum$Tenericutes
#Thermodesulfobacteria<-df_phylum$Thermodesulfobacteria
#Thermotogae<-df_phylum$Thermotogae
Verrucomicrobia<-df_phylum$Verrucomicrobia

phyll<-cbind(Acidobacteria, Actinobacteria, Armatimonadetes, Bacteroidetes, Chlamydiae,Chloroflexi,Cyanobacteria,Fibrobacteres,Firmicutes,Gemmatimonadetes,Nitrospinae,Planctomycetes,Proteobacteria,Spirochaetes,Tenericutes,Verrucomicrobia)
phyll<-t(phyll)


#potential functions
nitrification<-df*pm$nitrification
denitrification<-df*pm$denitrification
methanol_oxidation<-df*pm$methanol_oxidation
nitrogen_fixation<-df*pm$nitrogen_fixation
chitinolysis<-df*pm$chitinolysis
ligninolysis<-df*pm$ligninolysis
respiration_of_sulfur_compounds<-df*pm$respiration_of_sulfur_compounds
xylanolysis<-df*pm$xylanolysis
methanotrophy<-df*pm$methanotrophy						
cellulolysis<-df*pm$cellulolysis
ureolysis<-df*pm$ureolysis
plant_pathogen<-df*pm$plant_pathogen
predatory_or_exoparasitic<-df*pm$predatory_or_exoparasitic

#diversity of functions
#nitrification
Abund_nitrification <-rowSums(nitrification)
R_nitrification <- specnumber(nitrification) 
tab_nitrification <- data.frame(Abund_nitrification,R_nitrification)
t_nitrification <- t(tab_nitrification)

#denitrification
Abund_denitrification <-rowSums(denitrification)
R_denitrification <- specnumber(denitrification) 
tab_denitrification <- data.frame(Abund_denitrification,R_denitrification)
t_denitrification <- t(tab_denitrification)

#methanol_oxidation
Abund_methanol_oxidation <-rowSums(methanol_oxidation)
R_methanol_oxidation <- specnumber(methanol_oxidation) 
tab_methanol_oxidation <- data.frame(Abund_methanol_oxidation,R_methanol_oxidation)
t_methanol_oxidation <- t(tab_methanol_oxidation)

#nitrogen_fixation
Abund_nitrogen_fixation <-rowSums(nitrogen_fixation)
R_nitrogen_fixation <- specnumber(nitrogen_fixation) 
tab_nitrogen_fixation <- data.frame(Abund_nitrogen_fixation,R_nitrogen_fixation)
t_nitrogen_fixation <- t(tab_nitrogen_fixation)

#chitinolysis
Abund_chitinolysis <-rowSums(chitinolysis)
R_chitinolysis <- specnumber(chitinolysis) 
tab_chitinolysis <- data.frame(Abund_chitinolysis,R_chitinolysis)
t_chitinolysis <- t(tab_chitinolysis)

#ligninolysis
Abund_ligninolysis <-rowSums(ligninolysis)
R_ligninolysis <- specnumber(ligninolysis) 
tab_ligninolysis <- data.frame(Abund_ligninolysis,R_ligninolysis)
t_ligninolysis <- t(tab_ligninolysis)

#respiration_of_sulfur_compounds
Abund_respiration_of_sulfur_compounds <-rowSums(respiration_of_sulfur_compounds)
R_respiration_of_sulfur_compounds <- specnumber(respiration_of_sulfur_compounds) 
tab_respiration_of_sulfur_compounds <- data.frame(Abund_respiration_of_sulfur_compounds,R_respiration_of_sulfur_compounds)
t_respiration_of_sulfur_compounds <- t(tab_respiration_of_sulfur_compounds)

#xylanolysis
Abund_xylanolysis <-rowSums(xylanolysis)
R_xylanolysis <- specnumber(xylanolysis) 
tab_xylanolysis <- data.frame(Abund_xylanolysis,R_xylanolysis)
t_xylanolysis <- t(tab_xylanolysis)

#methanotrophy	
Abund_methanotrophy	<-rowSums(methanotrophy)
R_methanotrophy	<- specnumber(methanotrophy) 
tab_methanotrophy <- data.frame(Abund_methanotrophy,R_methanotrophy)
t_methanotrophy	<- t(tab_methanotrophy)

#cellulolysis
Abund_cellulolysis <-rowSums(cellulolysis)
R_cellulolysis <- specnumber(cellulolysis) 
tab_cellulolysis <- data.frame(Abund_cellulolysis,R_cellulolysis)
t_cellulolysis <- t(tab_cellulolysis)

#ureolysis
Abund_ureolysis <-rowSums(ureolysis)
R_ureolysis <- specnumber(ureolysis) 
tab_ureolysis <- data.frame(Abund_ureolysis,R_ureolysis)
t_ureolysis <- t(tab_ureolysis)

#plant_pathogen
Abund_plant_pathogen <-rowSums(plant_pathogen)
R_plant_pathogen <- specnumber(plant_pathogen) 
tab_plant_pathogen <- data.frame(Abund_plant_pathogen,R_plant_pathogen)
t_plant_pathogen <- t(tab_plant_pathogen)

#predatory_or_exoparasitic
Abund_predatory_or_exoparasitic <-rowSums(predatory_or_exoparasitic)
R_predatory_or_exoparasitic <- specnumber(predatory_or_exoparasitic) 
tab_predatory_or_exoparasitic <- data.frame(Abund_predatory_or_exoparasitic,R_predatory_or_exoparasitic)
t_predatory_or_exoparasitic <- t(tab_predatory_or_exoparasitic)

todo<-rbind(t_phylum,t_nitrification,t_denitrification,t_methanol_oxidation,t_nitrogen_fixation,t_chitinolysis,t_ligninolysis,t_respiration_of_sulfur_compounds,t_xylanolysis,t_methanotrophy,t_cellulolysis,t_ureolysis,t_plant_pathogen,t_predatory_or_exoparasitic)
todo<-data.matrix(todi)
todo<-rbind(todi,phyll)

storage.mode(todo) <-"numeric"
as.data.frame(todo)

library(xlsx)
write.xlsx(todo,"D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\Analizados-16s\\Ready-80-16s.xlsx",sheetName = "diver")

