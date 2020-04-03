#HC_Analysis:
library(tidyverse)
library(sjPlot)
library(sjmisc)

#Load Data:====
hc <- read.csv("hc1.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek2.csv") # read in C&N data, comes from Maria's CN-analysis of 9 replicate cores

#Compute mean dry.fraction + mean C_percent:=======
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars( C_percent, Dry_Fraction), mean)#mean values as hc samples were more coarse than sc
View(sc_high)
#Join hc and sc_high data======
hc_dry <- left_join(hc,sc_high, by = "SampleID_hc") #join data with mean C_percent and  Dry_Fraction

#Estimate PAH and TPH per dry mass by dividing content in wet sample (off hc) by mean dry fraction (off sc_high):
hc_dry <- hc_dry %>%
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x / hc_dry$Dry_Fraction )) #estimate mg/kg per dry mass

View(hc_dry)
#Normalize PAG & TPH per Total Organic Carbon:=====
#Normalized to enable comparison with toxicity guideline values.
#Example:If sediment contains 5 mg/kg of total PAHs and 0.55% OC, then
#1% OC normalised concentration = 5/0.55 = 9.1 mg/kg of total PAHs (1% OC). 

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

View(hc_dry_norm)
hc_dry_norm_total <- hc_dry_norm %>% #Computing Total content of TPH and PAH in samples:
mutate(TotPAH_norm = rowSums(hc_dry_norm[,(which(names(hc_dry_norm)=='Nap'):which(names(hc_dry_norm)=='BghiP'))])) %>% #sum up all normalized PAH 
mutate(TotTPH_norm = rowSums(hc_dry_norm[,(which(names(hc_dry_norm)=='TPH_06to09'):which(names(hc_dry_norm)=='TPH_36toMore'))])) #sum up all normalized PAH 

#ANOVA on total normalized TPH & PAH in relation to mangrove arrival time (below 20cm):======
#Below 20cm depth as it was the affected area in all cores
#Top 20cm was the mangrove organic matter:
hc_dry_norm_total_20to80 <- hc_dry_norm_total[hc_dry_norm_total$HC_DepthFrom.cm > 20,]

#Model TotTPH_norm in 20to80cm profile:
hist(log(hc_dry_norm_total_20to80$TotTPH_norm))
tph_model <- lm(log(TotTPH_norm)~SiteYear, data = hc_dry_norm_total_20to80)
summary(tph_model)
tab_model(tph_model)
###################Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        5.81060    0.86409   6.725 9.71e-06 ***
#SiteYearStony1996 -0.02105    1.49664  -0.014   0.9890    
#SiteYearStony2006  3.02188    1.39330   2.169   0.0478 *  

#Same model but ANOVA-way:
aov1<-aov(log(TotTPH_norm)~SiteYear, data = hc_dry_norm_total_20to80)
summary(aov1)
TukeyHSD(aov1)

#Model TotPAH_norm in 20to80cm profile:
hist(log(hc_dry_norm_total_20to80$TotPAH_norm))
pah_model <- lm(TotPAH_norm~SiteYear , data = hc_dry_norm_total_20to80)
summary(pah_model)
tab_model(pah_model)
#####################Estimate Std. Error t value Pr(>|t|)  
#(Intercept)          5.739      3.021   1.900   0.0713 .
#SiteYearStony1996    1.952      4.272   0.457   0.6524  
#SiteYearStony2006    8.753      4.777   1.832   0.0811 .


#HC_Stock.Mgha in sediment profile:========
#Compute dry_bulk_density.gcm3 and CarbonStock.Mgha off sc data per core:
sc_core <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  mutate (Core_in.cm            = PipeLenght.cm  - CompactionIn ,#Compaction in cm to get core length in cm
          Pipe_in.cm            = PipeLenght.cm  - CompactionOut ,#Compaction in cm to get how deep we hammerred pipe in in cm
          Compaction_Correction_Value = Core_in.cm / Pipe_in.cm, 
          SliceLength.cm       = DepthTo.cm - DepthFrom.cm, #height of slice which is cylinder
          SliceVolume.cm3_corrected       = (pi*(PipeDiameter.cm/2)^2) * SliceLength.cm / Compaction_Correction_Value,  #PVC pipe of 5 cm diameter
          dry_bulk_density.gcm3_corrected = DryWeight.g /SliceVolume.cm3_corrected , #Dry bulk density corrected for compaction, C-stock normalized to g per cm3
          CarbonDensity.gcm3    = dry_bulk_density.gcm3_corrected * C_percent/100,
          CarbonStock.Mgha      = CarbonDensity.gcm3 *100 * SliceLength.cm )# gcm3 times 100 gives Mgha

#Average C_percent, Dry_Fraction, dry_bulk_density.gcm3_correcte
#As slices in sc dataset were finer that in hc dataset:
sc_core_high <- sc_core %>% 
  dplyr::group_by(SampleID_hc) %>% 
  summarise_at(vars(C_percent, Dry_Fraction, dry_bulk_density.gcm3_corrected),mean)

#Compute HC_Core_Sum.Mgha in the sediment profile:
hc_dry_core<- left_join(hc,sc_core_high, by = "SampleID_hc") #join data with mea Dry_Fraction
  
hc_dry_core_sum <- hc_dry_core %>% 
  #create function to express HC per Dry.Fraction and normalized to mean C_percent:
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x / hc_dry_core$Dry_Fraction ))%>%
  #compute HC stock based on slice parameters:
  mutate (HC_SliceLength.cm     = HC_DepthTo.cm - HC_DepthFrom.cm, #height of slice which is cylinder
          HC_SliceVolume.cm3    = square_cm2 * HC_SliceLength.cm,  #sliced square prims of 2by2cm or 3by3cm bottom
          HC_Fraction = Total / 1000000, #Fraction of total HC in each sample (mg / kg)
          HC_Density.gcm3    = dry_bulk_density.gcm3_corrected *HC_Fraction, #normalized to g per cm3
          HC_Stock.Mgha      = HC_Density.gcm3 *100 * HC_SliceLength.cm)%>%
  dplyr::group_by(Site,SiteYear, transect, Elevation) %>% #Grouping by core
  summarise(HC_Core_Sum.Mgha = sum(HC_Stock.Mgha, na.rm = T)) %>%
  dplyr::group_by (SiteYear) %>% #grouping per site
    summarise(Mean_Core_HC = mean(HC_Core_Sum.Mgha), #computing mean HC_Stock.Mgha per SiteYear
            sd_core = sd(HC_Core_Sum.Mgha, na.rm=T),
            N = n(),
            se_core = sd_core/sqrt(N))

hc_dry_core_sum

#HC_Stock.Mgha  below 20 cm========
hc_dry_core_20to80cm <- left_join(hc,sc_core_high, by = "SampleID_hc") %>%#join data with Dry Fraction
  filter(HC_DepthTo.cm >20) 

#Compute HC_Stock.Mgha inside 20to80cm profile:
hc_dry_core_sum_20to80cm  <- hc_dry_core_20to80cm %>% 
  #create function to express HC per Dry.Fraction, normalized to mean C_percent:
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x / hc_dry_core_20to80cm$Dry_Fraction ))%>%
  #compute HC stock based on slice parameters:
  mutate (HC_SliceLength.cm     = HC_DepthTo.cm - HC_DepthFrom.cm, #height of slice which is cylinder
          HC_SliceVolume.cm3    = square_cm2 * HC_SliceLength.cm,  #sliced square prims of 2by2cm or 3by3cm bottom
          HC_Fraction = Total / 1000000, #Fraction of total HC in each sample (mg / kg)
          HC_Density.gcm3    = dry_bulk_density.gcm3_corrected *HC_Fraction, #normalized to g per cm3
          HC_Stock.Mgha      = HC_Density.gcm3 *100 * HC_SliceLength.cm)

#Add-up stock of each slice:
hc_20to80cm <- hc_dry_core_sum_20to80cm %>%
  dplyr::group_by(Site,SiteYear, transect, Elevation) %>% #Grouping by core
  summarise(HC_Core_Sum.Mgha = sum(HC_Stock.Mgha, na.rm = T))%>% #summing up per core
  dplyr::group_by (SiteYear) %>% #grouping per SiteYear
  summarise(Mean_Core_HC = mean(HC_Core_Sum.Mgha), #computing mean HC_Stock.Mgha per SiteYear
            sd_core = sd(HC_Core_Sum.Mgha, na.rm=T),
            N = n(),
            se_core = sd_core/sqrt(N))

hc_20to80cm
#SiteYear  Mean_Core_HC sd_core     N se_core
# Stony1986         45.7    12.2     3    7.02
# Stony1996         33.7    26.6     3   15.3 
# Stony2006         39.8    23.8     3   13.7 

hist(log(hc_20to80cm $ HC_Core_Sum.Mgha))
tph_model <- lm(log(HC_Core_Sum.Mgha)~SiteYear, data = hc_20to80cm)
summary(tph_model)
tab_model(tph_model)

#Correlation of TPH ~ PAH (entire core):========
hc_corr_pah <- hc_dry_norm_total %>%
  select(SampleID,TotPAH_norm) %>%
  gather(Type, pah_content, -SampleID) #Prep long format for stack-plotting.

hc_corr_tph <- hc_dry_norm_total %>%
  select(SampleID,TotTPH_norm ) %>%
  gather(Type, tph_content, -SampleID) #Prep long format for stack-plotting.

hc_corr <- cbind(hc_corr_pah,hc_corr_tph)

cor.test(hc_corr$pah_content,hc_corr$tph_content)
#t = 12.431, df = 49, p-value < 2.2e-16
# cor =  0.8713543 
