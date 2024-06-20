# Analysis of ITS OTU or ASV datasets
# Author: Pedro Mondaca (pedro.mondaca@pucv.cl)

library(vegan)
library(readxl)
library(tidyverse)
data <- read_excel("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\ITS\\744_ITS.xlsx",
                                      sheet = "Hoja1", na = "0")
data[is.na(data)] <- 0
dat<-data[-c(1)] 

nrow(data) #should equals 831
df0 <- dat
colSums(df0) # should equals 100 

df0 %>% 
  mutate(across(everything(), function(x) {
    x*100/sum(x)
  }))
colSums(df0)

class(data$family_planilla_maestra)
data$family_planilla_maestra <- as.factor(data$family_planilla_maestra)
class(data$family_planilla_maestra)
df <- as.data.frame(t(df0))
df

pm<- read_excel("D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\pm_fung.xlsx", 
                 sheet = "final", col_types = c("text","text","text","text", "text", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric","numeric"))

colnames(df)<-pm$family #family to factor

#alpha diversity
H_family <- diversity(df, index="shannon") #calculate Shannon index using vegan package
R_family <- specnumber(df) # Richness

#community analysis
df_class0 <- aggregate(df0, by = list(pm$class), FUN = sum)
df_class <- as.data.frame(t(df_class0[-c(1)]))
colnames(df_class)<-df_class0$Group.1

#Ver abundancia de class por tratamiento
Agaricomycetes<-df_class$Agaricomycetes
Agaricostilbomycetes<-df_class$Agaricostilbomycetes
Archaeorhizomycetes<-df_class$Archaeorhizomycetes
Arthoniomycetes<-df_class$Arthoniomycetes
Atractiellomycetes<-df_class$Atractiellomycetes
Bartheletiomycetes<-df_class$Bartheletiomycetes
Basidiobolomycetes<-df_class$Basidiobolomycetes
Blastocladiomycetes<-df_class$Blastocladiomycetes
Candelariomycetes<-df_class$Candelariomycetes
Chytridiomycetes<-df_class$Chytridiomycetes
Classiculomycetes<-df_class$Classiculomycetes
Collemopsidiomycetes<-df_class$Collemopsidiomycetes
Coniocybomycetes<-df_class$Coniocybomycetes
Cryptomycocolacomycetes<-df_class$Cryptomycocolacomycetes
Cystobasidiomycetes<-df_class$Cystobasidiomycetes
Dacrymycetes<-df_class$Dacrymycetes
Dimargaritomycetes<-df_class$Dimargaritomycetes
Dothideomycetes<-df_class$Dothideomycetes
Endogonomycetes<-df_class$Endogonomycetes
Entomophthoromycetes<-df_class$Entomophthoromycetes
Entorrhizomycetes<-df_class$Entorrhizomycetes
Eurotiomycetes<-df_class$Eurotiomycetes
Exobasidiomycetes<-df_class$Exobasidiomycetes
Geminibasidiomycetes<-df_class$Geminibasidiomycetes
Geoglossomycetes<-df_class$Geoglossomycetes
Glomeromycetes<-df_class$Glomeromycetes
Harpellomycetes<-df_class$Harpellomycetes
Kickxellomycetes<-df_class$Kickxellomycetes
Laboulbeniomycetes<-df_class$Laboulbeniomycetes
Lecanoromycetes<-df_class$Lecanoromycetes
Leotiomycetes<-df_class$Leotiomycetes
Lichinomycetes<-df_class$Lichinomycetes
Malasseziomycetes<-df_class$Malasseziomycetes
Microbotryomycetes<-df_class$Microbotryomycetes
Mixiomycetes<-df_class$Mixiomycetes
Moniliellomycetes<-df_class$Moniliellomycetes
Monoblepharidomycetes<-df_class$Monoblepharidomycetes
Mortierellomycetes<-df_class$Mortierellomycetes
Mucoromycetes<-df_class$Mucoromycetes
Neocallimastigomycetes<-df_class$Neocallimastigomycetes
Neolectomycetes<-df_class$Neolectomycetes
Neozygitomycetes<-df_class$Neozygitomycetes
Orbiliomycetes<-df_class$Orbiliomycetes
Pezizomycetes<-df_class$Pezizomycetes
Pneumocystidomycetes<-df_class$Pneumocystidomycetes
Pucciniomycetes<-df_class$Pucciniomycetes
Saccharomycetes<-df_class$Saccharomycetes
Sanchytriomycetes<-df_class$Sanchytriomycetes
Sareomycetes<-df_class$Sareomycetes
Schizosaccharomycetes<-df_class$Schizosaccharomycetes
Sordariomycetes<-df_class$Sordariomycetes
Spiculogloeomycetes<-df_class$Spiculogloeomycetes
Taphrinomycetes<-df_class$Taphrinomycetes
Tremellomycetes<-df_class$Tremellomycetes
Tritirachiomycetes<-df_class$Tritirachiomycetes
Umbelopsidomycetes<-df_class$Umbelopsidomycetes
Ustilaginomycetes<-df_class$Ustilaginomycetes
Wallemiomycetes<-df_class$Wallemiomycetes
Xylobotryomycetes<-df_class$Xylobotryomycetes
Xylonomycetes<-df_class$Xylonomycetes
Zoopagomycetes<-df_class$Zoopagomycetes

class<-cbind(Agaricomycetes, Agaricostilbomycetes, Archaeorhizomycetes, Arthoniomycetes, Atractiellomycetes, Bartheletiomycetes, Basidiobolomycetes, Blastocladiomycetes, Candelariomycetes, Chytridiomycetes, Classiculomycetes, Collemopsidiomycetes, Coniocybomycetes, Cryptomycocolacomycetes, Cystobasidiomycetes, Dacrymycetes, Dimargaritomycetes, Dothideomycetes, Endogonomycetes, Entomophthoromycetes, Entorrhizomycetes, Eurotiomycetes, Exobasidiomycetes, Geminibasidiomycetes, Geoglossomycetes, Glomeromycetes, Harpellomycetes, Kickxellomycetes, Laboulbeniomycetes, Lecanoromycetes, Leotiomycetes, Lichinomycetes, Malasseziomycetes, Microbotryomycetes, Mixiomycetes, Moniliellomycetes, Monoblepharidomycetes, Mortierellomycetes, Mucoromycetes, Neocallimastigomycetes, Neolectomycetes, Neozygitomycetes, Orbiliomycetes, Pezizomycetes, Pneumocystidomycetes, Pucciniomycetes, Saccharomycetes, Sanchytriomycetes, Sareomycetes, Schizosaccharomycetes, Sordariomycetes, Spiculogloeomycetes, Taphrinomycetes, Tremellomycetes, Tritirachiomycetes, Umbelopsidomycetes, Ustilaginomycetes, Wallemiomycetes, Xylobotryomycetes, Xylonomycetes, Zoopagomycetes)
#class<-cbind(Agaricomycetes,  Chytridiomycetes, Classiculomycetes, Collemopsidiomycetes,  Dothideomycetes, Endogonomycetes, Entorrhizomycetes, Eurotiomycetes, Exobasidiomycetes, Glomeromycetes, Harpellomycetes, Laboulbeniomycetes, Leotiomycetes, Lichinomycetes, Monoblepharidomycetes, Mucoromycetes, Pezizomycetes, Pneumocystidomycetes, Pucciniomycetes, Saccharomycetes, Sordariomycetes, Tremellomycetes, Ustilaginomycetes, Wallemiomycetes)
class_abund<-t(class)

#potential functions
arbuscular_mycorrhizal<-df*pm$arbuscular_mycorrhizal
arthropod_parasite<-df*pm$arthropod_parasite
arthropod_associated<-df*pm$arthropod_associated
dung_saprotroph<-df*pm$dung_saprotroph
ectomycorhizal<-df*pm$ectomycorhizal
epiphyte<-df*pm$epiphyte
foliar_endophyte<-df*pm$foliar_endophyte
fungal_decomposer<-df*pm$fungal_decomposer
litter_saprotroph <-df*pm$litter_saprotroph
mycoparasite <-df*pm$mycoparasite
nectar_tap_saprotroph<-df*pm$nectar_tap_saprotroph
plant_pathogen <-df*pm$plant_pathogen
pollen_saprotroph<-df*pm$pollen_saprotroph
protistan_parasite <-df*pm$protistan_parasite
root_endophyte<-df*pm$root_endophyte
soil_saprotroph <-df*pm$soil_saprotroph
sooty_mold <-df*pm$sooty_mold
root_associated<-df*pm$root_associated
unspecified_saprotroph<-df*pm$unspecified_saprotroph
unspecified_symbiotroph<-df*pm$unspecified_symbiotroph
wood_saprotroph<-df*pm$wood_saprotroph

#alpha diversity df genes
#arbuscular_mycorrhizal
H_arbuscular_mycorrhizal <- diversity(arbuscular_mycorrhizal, index="shannon") #calculate Shannon index using vegan package
Abund_arbuscular_mycorrhizal <-rowSums(arbuscular_mycorrhizal)
R_arbuscular_mycorrhizal <- specnumber(arbuscular_mycorrhizal) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_arbuscular_mycorrhizal <- data.frame(H_arbuscular_mycorrhizal,Abund_arbuscular_mycorrhizal,R_arbuscular_mycorrhizal)
t_arbuscular_mycorrhizal <- t(tab_arbuscular_mycorrhizal)

#arthropod_parasite
H_arthropod_parasite <- diversity(arthropod_parasite, index="shannon") #calculate Shannon index using vegan package
Abund_arthropod_parasite <-rowSums(arthropod_parasite)
R_arthropod_parasite <- specnumber(arthropod_parasite) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_arthropod_parasite <- data.frame(H_arthropod_parasite,Abund_arthropod_parasite,R_arthropod_parasite)
t_arthropod_parasite <- t(tab_arthropod_parasite)

#arthropod_associated
H_arthropod_associated <- diversity(arthropod_associated, index="shannon") #calculate Shannon index using vegan package
Abund_arthropod_associated <-rowSums(arthropod_associated)
R_arthropod_associated <- specnumber(arthropod_associated) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_arthropod_associated <- data.frame(H_arthropod_associated,Abund_arthropod_associated,R_arthropod_associated)
t_arthropod_associated <- t(tab_arthropod_associated)

#dung_saprotroph
H_dung_saprotroph <- diversity(dung_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_dung_saprotroph <-rowSums(dung_saprotroph)
R_dung_saprotroph <- specnumber(dung_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_dung_saprotroph <- data.frame(H_dung_saprotroph,Abund_dung_saprotroph,R_dung_saprotroph)
t_dung_saprotroph <- t(tab_dung_saprotroph)

#ectomycorhizal
H_ectomycorhizal <- diversity(ectomycorhizal, index="shannon") #calculate Shannon index using vegan package
Abund_ectomycorhizal <-rowSums(ectomycorhizal)
R_ectomycorhizal <- specnumber(ectomycorhizal) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_ectomycorhizal <- data.frame(H_ectomycorhizal,Abund_ectomycorhizal,R_ectomycorhizal)
t_ectomycorhizal <- t(tab_ectomycorhizal)

#epiphyte
H_epiphyte <- diversity(epiphyte, index="shannon") #calculate Shannon index using vegan package
Abund_epiphyte <-rowSums(epiphyte)
R_epiphyte <- specnumber(epiphyte) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_epiphyte <- data.frame(H_epiphyte,Abund_epiphyte,R_epiphyte)
t_epiphyte <- t(tab_epiphyte)

#foliar_endophyte
H_foliar_endophyte <- diversity(foliar_endophyte, index="shannon") #calculate Shannon index using vegan package
Abund_foliar_endophyte <-rowSums(foliar_endophyte)
R_foliar_endophyte <- specnumber(foliar_endophyte) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_foliar_endophyte <- data.frame(H_foliar_endophyte,Abund_foliar_endophyte,R_foliar_endophyte)
t_foliar_endophyte <- t(tab_foliar_endophyte)

#fungal_decomposer
H_fungal_decomposer <- diversity(fungal_decomposer, index="shannon") #calculate Shannon index using vegan package
Abund_fungal_decomposer <-rowSums(fungal_decomposer)
R_fungal_decomposer <- specnumber(fungal_decomposer) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_fungal_decomposer <- data.frame(H_fungal_decomposer,Abund_fungal_decomposer,R_fungal_decomposer)
t_fungal_decomposer <- t(tab_fungal_decomposer)

#litter_saprotroph
H_litter_saprotroph <- diversity(litter_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_litter_saprotroph <-rowSums(litter_saprotroph)
R_litter_saprotroph <- specnumber(litter_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_litter_saprotroph <- data.frame(H_litter_saprotroph,Abund_litter_saprotroph,R_litter_saprotroph)
t_litter_saprotroph <- t(tab_litter_saprotroph)

#mycoparasite
H_mycoparasite <- diversity(mycoparasite, index="shannon") #calculate Shannon index using vegan package
Abund_mycoparasite <-rowSums(mycoparasite)
R_mycoparasite <- specnumber(mycoparasite) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_mycoparasite <- data.frame(H_mycoparasite,Abund_mycoparasite,R_mycoparasite)
t_mycoparasite <- t(tab_mycoparasite)

#nectar_tap_saprotroph
H_nectar_tap_saprotroph <- diversity(nectar_tap_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_nectar_tap_saprotroph <-rowSums(nectar_tap_saprotroph)
R_nectar_tap_saprotroph <- specnumber(nectar_tap_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_nectar_tap_saprotroph <- data.frame(H_nectar_tap_saprotroph,Abund_nectar_tap_saprotroph,R_nectar_tap_saprotroph)
t_nectar_tap_saprotroph <- t(tab_nectar_tap_saprotroph)

#plant_pathogen
H_plant_pathogen <- diversity(plant_pathogen, index="shannon") #calculate Shannon index using vegan package
Abund_plant_pathogen <-rowSums(plant_pathogen)
R_plant_pathogen <- specnumber(plant_pathogen) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_plant_pathogen <- data.frame(H_plant_pathogen,Abund_plant_pathogen,R_plant_pathogen)
t_plant_pathogen <- t(tab_plant_pathogen)

#pollen_saprotroph
H_pollen_saprotroph <- diversity(pollen_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_pollen_saprotroph <-rowSums(pollen_saprotroph)
R_pollen_saprotroph <- specnumber(pollen_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_pollen_saprotroph <- data.frame(H_pollen_saprotroph,Abund_pollen_saprotroph,R_pollen_saprotroph)
t_pollen_saprotroph <- t(tab_pollen_saprotroph)

#protistan_parasite
H_protistan_parasite <- diversity(protistan_parasite, index="shannon") #calculate Shannon index using vegan package
Abund_protistan_parasite <-rowSums(protistan_parasite)
R_protistan_parasite <- specnumber(protistan_parasite) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_protistan_parasite <- data.frame(H_protistan_parasite,Abund_protistan_parasite,R_protistan_parasite)
t_protistan_parasite <- t(tab_protistan_parasite)

#root_endophyte
H_root_endophyte <- diversity(root_endophyte, index="shannon") #calculate Shannon index using vegan package
Abund_root_endophyte <-rowSums(root_endophyte)
R_root_endophyte <- specnumber(root_endophyte) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_root_endophyte <- data.frame(H_root_endophyte,Abund_root_endophyte,R_root_endophyte)
t_root_endophyte <- t(tab_root_endophyte)

#soil_saprotroph
H_soil_saprotroph <- diversity(soil_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_soil_saprotroph <-rowSums(soil_saprotroph)
R_soil_saprotroph <- specnumber(soil_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_soil_saprotroph <- data.frame(H_soil_saprotroph,Abund_soil_saprotroph,R_soil_saprotroph)
t_soil_saprotroph <- t(tab_soil_saprotroph)

#sooty_mold
H_sooty_mold <- diversity(sooty_mold, index="shannon") #calculate Shannon index using vegan package
Abund_sooty_mold <-rowSums(sooty_mold)
R_sooty_mold <- specnumber(sooty_mold) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_sooty_mold <- data.frame(H_sooty_mold,Abund_sooty_mold,R_sooty_mold)
t_sooty_mold <- t(tab_sooty_mold)

#root_associated
H_root_associated <- diversity(root_associated, index="shannon") #calculate Shannon index using vegan package
Abund_root_associated <-rowSums(root_associated)
R_root_associated <- specnumber(root_associated) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_root_associated <- data.frame(H_root_associated,Abund_root_associated,R_root_associated)
t_root_associated <- t(tab_root_associated)

#unspecified_saprotroph
H_unspecified_saprotroph <- diversity(unspecified_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_unspecified_saprotroph <-rowSums(unspecified_saprotroph)
R_unspecified_saprotroph <- specnumber(unspecified_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_unspecified_saprotroph <- data.frame(H_unspecified_saprotroph,Abund_unspecified_saprotroph,R_unspecified_saprotroph)
t_unspecified_saprotroph <- t(tab_unspecified_saprotroph)

#unspecified_symbiotroph
H_unspecified_symbiotroph <- diversity(unspecified_symbiotroph, index="shannon") #calculate Shannon index using vegan package
Abund_unspecified_symbiotroph <-rowSums(unspecified_symbiotroph)
R_unspecified_symbiotroph <- specnumber(unspecified_symbiotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_unspecified_symbiotroph <- data.frame(H_unspecified_symbiotroph,Abund_unspecified_symbiotroph,R_unspecified_symbiotroph)
t_unspecified_symbiotroph <- t(tab_unspecified_symbiotroph)

#wood_saprotroph
H_wood_saprotroph <- diversity(wood_saprotroph, index="shannon") #calculate Shannon index using vegan package
Abund_wood_saprotroph <-rowSums(wood_saprotroph)
R_wood_saprotroph <- specnumber(wood_saprotroph) # Richness #Saqu? J ya que al estar dividido por H=0 enmuchos casos, resultaba en un valor infinito
tab_wood_saprotroph <- data.frame(H_wood_saprotroph,Abund_wood_saprotroph,R_wood_saprotroph)
t_wood_saprotroph <- t(tab_wood_saprotroph)


traits<-rbind(t_arbuscular_mycorrhizal,t_arthropod_parasite,t_arthropod_associated,
            t_dung_saprotroph,t_ectomycorhizal,t_epiphyte,t_foliar_endophyte,
            t_fungal_decomposer, t_litter_saprotroph, t_mycoparasite, t_nectar_tap_saprotroph, 
            t_plant_pathogen, t_pollen_saprotroph, t_protistan_parasite, t_root_endophyte,
            t_soil_saprotroph, t_sooty_mold, t_root_associated, t_unspecified_saprotroph, 
            t_unspecified_symbiotroph, t_wood_saprotroph)


#Final step
colSums(df)
S
todi<-rbind(t_class,t_order,t_family,S)
todi<-data.matrix(todi)
todi
todi<-rbind(todi,class_abund,traits)
traits
storage.mode(todi) <-"numeric"
as.data.frame(todi)

library(xlsx)
write.xlsx(todi,"D:\\Backup\\2022\\1-Meta-análisis microbioma suelos agrícolas\\2023\\Analizados-ITS\\Ready-744-ITS-2023-Prueba.xlsx",sheetName = "diver")
