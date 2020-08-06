#DRAW FIG:
#Relative abundance of PAH and TPH in relation to Total_PAH and Total_TPH, respectively.

#Load Packages=======
library("tidyverse")
library("gridExtra")
library("grid")

hc <- read.csv("hc1.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek2.csv") # read in C&N data, comes from Maria's CN-analysis of 9 replicate cores

#Compute mean dry.fraction + mean C_percent:=======
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars( C_percent, Dry_Fraction), mean)#mean values as hc samples were more coarse than sc

#Join hc and sc_high data======
hc_joined <- left_join(hc,sc_high, by = "SampleID_hc") #join data with mean C_percent and  Dry_Fraction

#Estimate PAH and TPH per dry mass by dividing content in wet sample (off hc) by mean dry fraction (off sc_high):
hc_dry <- hc_joined %>%
  #Total = sum of all TPH-s:
  mutate_at(vars(Naphthalene:Total_HC), function(x, na.rm=T)(x / hc_dry$Dry_Fraction ))%>% #estimate mg/kg per dry mass
  mutate(Total_TPH = rowSums(hc_dry[,(which(names(hc_dry)=='TPH_C_06to09'):which(names(hc_dry)=='TPH_C_37UP'))])) %>%#sum up all TPH-s to double check if it is the same as Total_HC (it is!)
  
  #PAH has to be normalized ot 0.2-10% C_percent. If outside 0.2-10% limit divide by min or max if below or above.
  #WEB: https://www.waterquality.gov.au/anz-guidelines/guideline-values/default/sediment-quality-toxicants
  mutate(OC_norm = ifelse(C_percent < 0.2, 0.2,  ifelse(C_percent > 10, 10, hc_dry$C_percent))) %>%
  
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
         TPH_37toUP = TPH_C_37UP )


#PAH relative abundance:=======
hc_dry_norm_PAH <- hc_dry %>%
  mutate_at( vars(Nap:BghiP), function(x, na.rm=T)( x / hc_dry$Total_PAH_norm *100)) %>% # % of single PAH against Total PAH
  select(Nap:BghiP)%>%
  gather(PAH_type,value) %>% 
  mutate (Site = "Stony")

pah_abund <- ggplot(data = hc_dry_norm_PAH, aes( y = value, x = Site)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", width=0.2, size = 1) +
  stat_summary(fun.y = "mean", size = 3, geom = "bar", fill = "lightblue")+
  facet_grid(~PAH_type)+
  theme_bw()+
  labs(x ="", y = "PAH Relative abundance (%)")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.title.y=element_text(size=16,colour = "black",face = "bold",),
        axis.title.x=element_text(size=16, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=11))

#ggsave(filename = "PAH_RelativeAbundance.jpeg", width = 21, height = 8,units = "cm",dpi = 600)

# Mean % PAH abundances:====
hc_dry_norm_PAH2 <- hc_dry_norm_PAH %>%
  group_by(PAH_type,Site)%>%
  summarise(AV =mean(value, na.rm=T),
            SD = sd(value),
            N = n(),
            SE= SD/sqrt(N))

hc_dry_norm_PAH2

#TPH relative abundance:====

hc_dry_TPH <- hc_dry %>%
  mutate_at( vars(TPH_06to09:TPH_37toUP), function(x, na.rm=T)( x / hc_dry$Total_TPH *100)) %>% # % of single PAH against Total PAH
  select(TPH_06to09:TPH_37toUP)%>%
  gather(TPH_type,value) %>% 
  mutate (Site = "Stony")

tph_abund <- ggplot(data = hc_dry_TPH, aes( y = value, x = Site)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", width=0.2, size = 1) +
  stat_summary(fun.y = "mean", size = 3, geom = "bar", fill = "lightblue")+
  facet_grid(~TPH_type)+
  theme_bw()+
  labs(x ="", y = "TPH Relative abundance (%)")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.title.y=element_text(size=16,colour = "black",face = "bold",),
        axis.title.x=element_text(size=16, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=11))


hc_dry_TPH_AV <- hc_dry_TPH %>%
  group_by(TPH_type)%>%
  summarise(AV =mean(value, na.rm=T),
            SD = sd(value),
            N = n(),
            SE= SD/sqrt(N))

hc_dry_TPH_AV

#Relative abundance lm:
abund_TPH <- hc_dry %>%
  mutate_at( vars(TPH_06to09:TPH_37toUP), function(x, na.rm=T)( x / hc_dry$Total_TPH *100)) %>% # % of single PAH against Total PAH
  select(TPH_06to09:TPH_37toUP, SiteYear)%>%
  gather(TPH_type,value, - SiteYear) %>% 
  mutate (Site = "Stony")

summary(lm(value~TPH_type, data = abund_TPH))#No interaction (TPH_type*SiteYear) of SiteYear with TPH_type!

abund_TPH_AV <- abund_TPH %>%
  group_by(TPH_type,SiteYear)%>%
  summarise(AV =mean(value, na.rm=T),
            SD = sd(value),
            N = n(),
            SE= SD/sqrt(N))
abund_TPH_AV

#Bind tph_stack and pah_stack together:
g <- arrangeGrob(pah_abund,tph_abund, nrow=2) #generates g

ggsave(g, filename = "FigS2_PAH_TPH_Barplot.jpeg", width = 24, height = 17,units = "cm",dpi = 600)



#Compute times above toxicity level (SQGV-Value):===========
hc_above <- hc_dry %>%
  mutate(Threshold_TPH = 280,
         Total_TPH = TPH_06to09+TPH_10to14+TPH_15to28+TPH_29to36+TPH_37toUP,
         Times_Above = Total_TPH/ Threshold_TPH,
         Sample_Above = ifelse(Times_Above >=1, 1,0))# To compute frequency of times above threshold

min(hc_above$Times_Above, na.rm = T)# 0.005847627
max(hc_above$Times_Above, na.rm = T)# 220.2039 max times above toxicity threshold

#Compute frequency of times above threshold:
good_total <- dim(subset(hc_above, !is.na(hc_above$Sample_Above))) #exclude NA-samples (43 of 51 were good records)
times_above <- dim(subset(hc_above, hc_above$Sample_Above > 0)) #34 times was above toxicity threshold
times_above*100/good_total # 86.27451  % of samples were above toxicity threshold
