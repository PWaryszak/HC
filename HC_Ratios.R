#Look at Nap/PAH_Total ratio:
#As per: Wang, Zhendi, and Merv F Fingas. 2003. “Development of Oil Hydrocarbon Fingerprinting and Identification Techniques.” 
#Marine Pollution Bulletin 47 (9): 423–52.
#https://doi.org/https://doi.org/10.1016/S0025-326X(03)00215-7.

#LOAD packages and DATA:========
# Start the code
rm(list = ls(all = TRUE))

library(tidyverse)
library(Hmisc)

hc <- read.csv("hc1.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek2.csv") # read in C&N data, comes from Maria's CN-analysis of 9 replicate cores

#Compute mean dry.fraction + mean C_percent:=======
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars( C_percent, Dry_Fraction), mean)#mean values as hc samples were more coarse than sc

#Join hc and sc_high data======
hc_joined <- left_join(hc,sc_high, by = "SampleID_hc") %>% #join data with mean C_percent and  Dry_Fraction
  mutate (ratio1 = Naphthalene/Total_PAH,
           ratio2 = Total_PAH/ Total_HC)

#Estimate PAH and TPH per dry mass by dividing content in wet sample (off hc) by mean dry fraction (off sc_high):
hc_dry <- hc_joined %>%
  #Total = sum of all TPH-s, Total_TPH is sum of TPH in dry sediments
  mutate_at(vars(Naphthalene:Total_HC), function(x, na.rm=T)(x / hc_joined$Dry_Fraction ))%>% #estimate mg/kg per dry mass
  mutate(Total_TPH = rowSums(hc_joined[,(which(names(hc_joined)=='TPH_C_06to09'):which(names(hc_joined)=='TPH_C_37UP'))])) %>%#sum up all TPH-s to double check if it is the same as Total_HC (it is!)
  
  #PAH has to be normalized ot 0.2-10% C_percent. If outside 0.2-10% limit divide by min or max if below or above.
  #WEB: https://www.waterquality.gov.au/anz-guidelines/guideline-values/default/sediment-quality-toxicants
  mutate(OC_norm = ifelse(C_percent < 0.2, 0.2,  ifelse(C_percent > 10, 10, hc_joined$C_percent))) %>%
  
  #ABBREVIATE names to standard abbreviations and normalize all PAHs to OC_norm:
  mutate(Nap = Naphthalene /OC_norm,
         Met = X2.Methylnaphthalene /OC_norm  ,
         Acy = Acenaphthylene       /OC_norm  ,
         Ace = Acenaphthene         /OC_norm  ,
         Ant = Anthracene           /OC_norm  ,
         Flu = Fluorene             /OC_norm  ,
         Flt = Fluoranthene         /OC_norm  ,
         Phe = Phenanthrene         /OC_norm ,
         Pyr = Pyrene               /OC_norm  ,
         BaA = Benzaanthracene     /OC_norm  ,
         Chr = Chrysene            /OC_norm  ,
         BbF = Benzobfluoranthene  /OC_norm  ,
         BkF = Benzokfluoranthene   /OC_norm  ,
         BaP = Benzoapyrene         /OC_norm  ,
         IcdP = Indeno1.2.3c.dpyrene/OC_norm  ,
         DahA = Dibenza.hanthracene /OC_norm  ,
         BghiP = Benzoghiperylene   /OC_norm  ,
         Total_PAH_norm = Total_PAH /OC_norm,
         TPH_06to09  = TPH_C_06to09  ,
         TPH_10to14 =  TPH_C_10to14  ,
         TPH_15to28 = TPH_C_15to28   ,
         TPH_29to36 = TPH_C_29to36   ,
         TPH_37toUP = TPH_C_37UP ) %>%
  mutate(Nap_PAH = Nap/Total_PAH_norm,
         PAH_TPH = Total_PAH_norm/ Total_TPH,
         FltPyr = Flt/Pyr,
         FluPyr = Flu/Pyr,
         PheAnt = Phe/Ant)

#EXPLORE Values:
summary(hc_dry$FluPyr) #0.062500  but we used the wron PAH here. It should be flt! as per paper below
range(hc_dry$FluPyr) #LOW, if below 0.4 petrogening source
#WEB: https://www.sciencedirect.com/science/article/abs/pii/S0025326X02001170

#Correct ratio:
summary(hc_dry$FltPyr) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1538  0.5000  0.7500  0.6905  0.9062  1.0000 
range(hc_dry$FltPyr)#0.1538462 1.0000000
ggplot(hc_dry, aes(y=HC_Depth_Range, x=FltPyr)) + geom_point() +coord_flip()


summary(hc_dry$PheAnt) #Phenanthrene:Anthracene
###Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
#0.2667  1.0000  2.0000  2.4015  4.0000  8.0000 
ggplot(hc_dry, aes(y=HC_Depth_Range, x=PheAnt)) + geom_point() +coord_flip()
#This plot can go into suppl





#ANOVAS:
summary(lm(Nap_PAH~HC_Depth_Range, data = hc_dry))
summary(lm(PAH_TPH~HC_Depth_Range, data = hc_dry))

#HC in Wet soil:
summary(lm(ratio1~HC_Depth_Range, data = hc_joined))
summary(lm(ratio2~HC_Depth_Range, data = hc_joined))

