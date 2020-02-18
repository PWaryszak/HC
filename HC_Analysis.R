#HC_Analysis:
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(funrar)

setwd("~/00DeakinUni/R/BCL_R/BCL/HC")

hc <- read.csv("hc.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek.csv") # read in CN data

#Compute mean dry.fraction + mean C_percent:
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars(CompactionIn, CompactionOut, C_percent, Dry_Fraction),mean)

#Join HC and CN data
#Estimate PAH and TPH per dry mass:
hc_dry <- left_join(hc,sc_high, by = "SampleID_hc") #join data with Dry Fraction

hc_dry <- hc_dry %>%
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x * hc_dry$Dry_Fraction )) #estimate mg/kg per dry mass

#Normalize PAG & TPH per Total Organic Carbon:=====
#Total Organic Carbon comes from Maria's CN-analysis of 9 replicate cores (Mean C-values per Depth reps).
hc_dry_norm <- hc_dry %>%
  mutate(Nap = Naphthalene / C_percent,
         Met = X2.Methylnaphthalene /C_percent,
         Acy = Acenaphthylene / C_percent,
         Ace = Acenaphthene / C_percent,
         Flu = Fluorene / C_percent,
         Pyr = Pyrene / C_percent,
         BaA = Benzaanthracene / C_percent,
         Chr = Chrysene / C_percent,
         BbG = Benzobfluoranthene / C_percent,
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

hc_dry_norm$TotALL_norm <- rowSums(hc_dry_norm[,42:60]) #Normalized Total Hydrocarbons
hc_dry_norm$TotPAH_norm <- rowSums(hc_dry_norm[,42:55])#Normalized Total PAH only
hc_dry_norm$TotTPH_norm <- rowSums(hc_dry_norm[,56:60]) #Normalized Total Hydrocarbons

View(hc_dry_norm)
#Example of Normalizing to %OC:
#If sediment contains 5 mg/kg of total PAHs and 0.55% OC, then
#1% OC normalised concentration = 5/0.55 = 9.1 mg/kg of total PAHs (1% OC). 

#ANOVA on total normalized TPH & PAH in relation to mangrove arrival time:======
hist(hc_dry_norm$TotTPH_norm)
tph_model <- lm(log(TotTPH_norm)~SiteYear, data = hc_dry_norm)
summary(tph_model)
tab_model(tph_model)

hist(hc_dry_norm$TotPAH_norm)
pah_model <- lm(log(TotPAH_norm)~SiteYear, data = hc_dry_norm)
summary(tph_model)
tab_model(pah_model)


#ANOVA on relatvie abundances======
hc_dry_stack_pah <- hc_dry_norm %>%
  group_by(SiteYear) %>%
  summarise_at( vars(Nap:BghiP), sum, na.rm=T) %>% #Compute sum of all PAH-s per column
  gather(PAH_Type, content, -SiteYear) #Prep long format for stack-plotting.

hist(hc_dry_stack_pah$TotTPH_norm)
tph_model <- lm(log(TotTPH_norm)~SiteYear, data = hc_dry_norm)
summary(tph_model)
tab_model(tph_model)

hist(hc_dry_norm$TotPAH_norm)
pah_model <- lm(log(TotPAH_norm)~SiteYear, data = hc_dry_norm)
summary(tph_model)
tab_model(pah_model)


#ANOVA on normalized TPH & PAH relative abundances in relation to mangrove arrival time:======
hc_dry_stack <- hc_dry_norm %>%
  group_by(SiteYear) %>%
  summarise_at( vars(Nap:BghiP), sum, na.rm=T) %>% #Compute sum of all PAH-s per column
  gather(PAH_Type, content, -SiteYear) #Prep long format for stack-plotting.


#Playing with relative abundance data======
library(funrar)
library(ade4)
data("aravo", package = "ade4")


z <-  spread(hc_dry_stack, PAH_Type, content) %>%
  remove_rownames %>%   column_to_rownames(var="SiteYear")

z_mat = as.matrix(z)
rel_mat = make_relative(z_mat)
rel_mat
data <- as.data.frame(rel_mat)

z_stack <-  rownames_to_column (.data, var = "rowname") %>%
  gather(as.data.frame(rel_mat), PAH_Type, content)
z_stack

ggplot(z_stack, aes(x = PAH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "PAH Type",y="Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "right",
        strip.text=element_text(size=18))
