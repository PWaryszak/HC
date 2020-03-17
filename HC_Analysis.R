#HC_Analysis:
library(tidyverse)
library(sjPlot)
library(sjmisc)

#Load Data:====
hc <- read.csv("hc.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek.csv") # read in C&N data, comes from Maria's CN-analysis of 9 replicate cores

#Compute mean dry.fraction + mean C_percent:=======
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars(CompactionIn, CompactionOut, C_percent, Dry_Fraction),mean)#mean values as hc samples were more coarse than sc

#Join hc and sc_high data======
hc_dry <- left_join(hc,sc_high, by = "SampleID_hc") #join data with Dry Fraction

#Estimate PAH and TPH per dry mass by dividing content in wet sample (off hc) by mean dry fraction (off sc_high):
hc_dry <- hc_dry %>%
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x / hc_dry$Dry_Fraction )) #estimate mg/kg per dry mass

#Normalize PAG & TPH per Total Organic Carbon:=====
#Normalized to enable comparison with toxicity guideline values.
hc_dry_norm <- hc_dry %>%
  mutate(Nap = Naphthalene / C_percent,
         Met = X2.Methylnaphthalene /C_percent,
         Acy = Acenaphthylene / C_percent,
         Ace = Acenaphthene / C_percent,
         Ant = Anthracene / C_percent,
         Flu = Fluorene / C_percent,
         Flt = Fluoranthene / C_percent,
         Phe = Phenanthrene/C_percent,
         Pyr = Pyrene / C_percent,
         BaA = Benzaanthracene / C_percent,
         Chr = Chrysene / C_percent,
         BbF = Benzobfluoranthene / C_percent,
         BkF = Benzokfluoranthene / C_percent,
         BaP = Benzoapyrene / C_percent,
         IcdP = Indeno1.2.3c.dpyrene / C_percent,
         DahA = Dibenza.hanthracene / C_percent,
         BghiP = Benzoghiperylene / C_percent,
         TPH_06to09  =TPH_C6.9 / C_percent,
         TPH_10to14 = TPH_C10.14 / C_percent,
         TPH_15to28 = TPH_C15.28 / C_percent ,
         TPH_29to36 = TPH_C29.36 / C_percent ,
         TPH_36toMore = TPH_C.36/ C_percent)

hc_dry_norm$TotPAH_norm <- rowSums(hc_dry_norm[,56:72]) #sum up all normalized  PAH 
hc_dry_norm$TotTPH_norm <- rowSums(hc_dry_norm[,73:77]) #sum up all normalized  TPH
hc_dry_norm$TotALL_norm <- rowSums(hc_dry_norm[,56:77]) #sum up all normalized  Hydrocarbons

View(hc_dry_norm)
#Example:If sediment contains 5 mg/kg of total PAHs and 0.55% OC, then
#1% OC normalised concentration = 5/0.55 = 9.1 mg/kg of total PAHs (1% OC). 

#ANOVA on total normalized TPH & PAH in relation to mangrove arrival time (below 20cm):======
#Below 20cm depth as it was the affected area in all cores
#Top 20cm was the mangrove organic matter:
hist(log(hc_dry_norm$TotTPH_norm))
tph_model <- lm(log(TotALL_norm)~SiteYear, data = hc_dry_norm[hc_dry_norm$HC_DepthFrom.cm > 20,])
summary(tph_model)
tab_model(tph_model)
###################Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        5.81060    0.86409   6.725 9.71e-06 ***
#SiteYearStony1996 -0.02105    1.49664  -0.014   0.9890    
#SiteYearStony2006  3.02188    1.39330   2.169   0.0478 *  

#Correlation of TPH ~ PAH (entire core):========
hc_corr_pah <- hc_dry_norm %>%
  select(SampleID,TotPAH_norm) %>%
  gather(Type, pah_content, -SampleID) #Prep long format for stack-plotting.

hc_corr_tph <- hc_dry_norm %>%
  select(SampleID,TotALL_norm ) %>%
  gather(Type, tph_content, -SampleID) #Prep long format for stack-plotting.

hc_corr <- cbind(hc_corr_pah,hc_corr_tph)

cor.test(hc_corr$pah_content,hc_corr$tph_content, method =  "pearson")
#t = 10.2, df = 41, p-value =  8.654e-13
# cor =  0.85
