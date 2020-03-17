#Load Packages=======
library("tidyverse")
library("gridExtra")
library("grid")

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

#Normalize PAH & TPH to 1% of Total Organic Carbon (C_percent):=====
#Example:If sediment contains 5 mg/kg of total PAHs and 0.55% OC, then
#1% OC normalised concentration = 5/0.55 = 9.1 mg/kg of total PAHs (1% OC). 
hc_dry_norm <- hc_dry %>%   
  mutate(Nap = Naphthalene / C_percent, #Normalized, to enable comparison with toxicity guideline values.
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

hc_dry_norm$TotALL_norm <- rowSums(hc_dry_norm[,56:77]) #Normalized Total Hydrocarbons
hc_dry_norm$TotPAH_norm <- rowSums(hc_dry_norm[,56:72]) #Normalized Total PAH only
hc_dry_norm$TotTPH_norm <- rowSums(hc_dry_norm[,73:77]) #Normalized Total Hydrocarbons

#Plot PAH relative abundance:=======
hc_dry_norm2 <- hc_dry_norm %>%
  mutate_at( vars(Nap:BghiP), function(x, na.rm=T)( x / hc_dry_norm$TotPAH_norm *100)) %>% # % of single PAH against Total PAH
  select(Nap:BghiP)%>%
  gather(PAH_type,value) %>% 
  mutate (Site = "Stony")

ggplot(data = hc_dry_norm2, aes( y = value, x = Site)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", width=0.2, size = 1) +
  stat_summary(fun.y = "mean", size = 3, geom = "bar", fill = "lightblue")+
  facet_grid(~PAH_type)+
  theme_bw()+
  labs(x ="", y = "Relative abundance (%)")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.title.y=element_text(size=16,colour = "black",face = "bold",),
        axis.title.x=element_text(size=16, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=11))

ggsave(filename = "PAH_RelativeAbundance.jpeg", 
      width = 21, 
      height = 8,
      units = "cm",
      dpi = 600)

#TPH per HC_Depth_Range :=======
#With Threshold dotted line (as per Agusti's 2019 paper)
h1 <- hc_dry_norm %>% 
  select(Year, HC_Depth_Range, TotALL_norm) %>%
  group_by(Year,HC_Depth_Range) %>%
  summarise(Mean_TPH = mean(TotALL_norm, na.rm=T),
            sd_TPH = sd(TotALL_norm, na.rm=T),
            N = n(),
            se_TPH = sd_TPH/sqrt(N))

h1_plot <- ggplot(h1,aes(x= reorder(HC_Depth_Range, desc(HC_Depth_Range)), y=Mean_TPH ,
               group = as.factor(Year)), color = reorder(as.factor(Year))) +
  labs(y = bquote('Total Petroleum Hydrocarbons   ' (mg~kg^-1)), x="Sediment depth (cm)") +
  geom_point(aes(color=as.factor(Year)),size = 2) +  geom_errorbar(aes(ymin = Mean_TPH -se_TPH, ymax = Mean_TPH+se_TPH), width = 0.2)+
  geom_line(aes(color=as.factor(Year)),size = 1.2, alpha = 0.4)+
  facet_grid( as.factor(Year) ~ .)+ 
  geom_hline(aes(yintercept = 275), linetype="dashed") + #275mg/kg is TPH toxicity threshold
  theme_bw()+
  coord_flip()+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black"),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=18))

h1_plot

ggsave(filename = "TPH_MeanContentGOOD2.jpeg", 
       width = 17, 
       height = 9,
       units = "cm",
       dpi = 600)

#PAH per Depth_Range:=======
h2 <- hc_dry_norm %>% 
  select(Year, HC_Depth_Range, TotPAH_norm) %>%
  group_by(Year,HC_Depth_Range) %>%
  summarise(Mean_PAH = mean(TotPAH_norm, na.rm=T),
            sd_PAH = sd(TotPAH_norm, na.rm=T),
            N = n(),
            se_PAH = sd_PAH/sqrt(N))

h2_plot <- ggplot(h2,aes(x= reorder(HC_Depth_Range, desc(HC_Depth_Range)), y=Mean_PAH ,
                     group = as.factor(Year)), color = reorder(as.factor(Year))) +
  labs(y = bquote('Total Polyaromatic Hydrocarbons   ' (mg~kg^-1)), x="Sediment depth (cm)") +
  #labs(x = "Depth (cm)",y="PAH (mg/kg)")+
  geom_point(aes(color=as.factor(Year)),size = 2) +  geom_errorbar(aes(ymin = Mean_PAH -se_PAH, ymax = Mean_PAH+se_PAH), width = 0.2)+
  geom_line(aes(color=as.factor(Year)),size = 1.2, alpha = 0.4)+
  facet_grid( as.factor(Year) ~ .)+ 
  geom_hline(aes(yintercept = 50), linetype="dashed") + #50 mg/kg is PAH toxicity threshold
  #or 50000/kg as per web: https://www.waterquality.gov.au/anz-guidelines/guideline-values/default/sediment-quality-toxicants
  theme_bw()+
  coord_flip()+
  #ggtitle("PAH @ Stony Creek (mg/kg)")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black"),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=18))
h2_plot

ggsave(filename = "PAH_MeanContentGOOD_DividedByDRyFraction2.jpeg", 
       width = 17, 
       height = 9,
       units = "cm",
       dpi = 600)

#Bind h1 and h2 together
grid.arrange(h1_plot, h2_plot, ncol = 1)
g <- arrangeGrob(h1_plot, h2_plot, nrow=2) #generates g

ggsave(g, filename = "PAH_TPH_MeanContentPerDepthGOOD.jpeg", 
       width = 17, 
       height = 19,
       units = "cm",
       dpi = 600)


#Stack Bar of PAH and TPH at three sites============
hc_dry_stack_pah <- hc_dry_norm %>%
  group_by(SiteYear) %>%
  summarise_at( vars(Nap:BghiP), sum, na.rm=T) %>% #Compute sum of all HC-s per column
  gather(PAH_Type, content, -SiteYear) #Prep long format for stack-plotting.

pah_stack <- ggplot(hc_dry_stack_pah, aes(x = PAH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "PAH Type",y="Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=18))

pah_stack

#Stack Bar of TPH at three sites:
hc_dry_stack_tph <- hc_dry_norm %>%
  group_by(SiteYear) %>%
  summarise_at( vars(TPH_06to09:Total_norm), sum, na.rm=T) %>% #Compute sum of all PAH-s per column
  gather(TPH_Type, content, -SiteYear) #Prep long format for stack-plotting.

tph_stack <- ggplot(hc_dry_stack_tph, aes(x = TPH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "TPH Type",y="Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        strip.text=element_text(size=18))

tph_stack
#Bind tph_stack and pah_stack together:
grid.arrange(pah_stack,tph_stack,  ncol = 1)
g <- arrangeGrob(pah_stack,tph_stack, nrow=2) #generates g

ggsave(g, filename = "PAH_TPH_Stack.jpeg", 
       width = 17, 
       height = 19,
       units = "cm",
       dpi = 600)

#Compute times above toxicity level (High-Guideline-Value):===========
hc_above <- hc_dry_norm  %>%
  mutate(Threshold_TPH = 275,
         Total_TPH = TPH_06to09+TPH_10to14+TPH_15to28+TPH_29to36+TPH_36toMore,
         Times_Above = Total_TPH/ Threshold_TPH,
         Sample_Above = ifelse(Times_Above >=1, 1,0))# To compute frequency of times above threshold

min(hc_above$Times_Above, na.rm = T)# 0.0009771783
max(hc_above$Times_Above, na.rm = T)# 47.60946:

#Compute frequency of times above threshold:
good_total <- dim(subset(hc_above, !is.na(hc_above$Sample_Above))) #exclude NA-samples (43 of 51 were good records)
times_above <- dim(subset(hc_above, hc_above$Sample_Above >0)) #34 times was above toxicity threshold
times_above*100/good_total # 79.06977 % of samples were above toxicity threshold


#HC Stock Tonnes / ha (as  per core) =====
#We need to compute dry_bulk_density.gcm3 and CarbonStock.Mgha off sc data per core:
sc_core <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  mutate (SliceLength.cm       = DepthTo.cm - DepthFrom.cm, #height of slice which is cylinder
          SliceVolume.cm3       = (pi*(PipeDiameter.cm/2)^2) * SliceLength.cm,  #PVC pipe of 5 cm diameter
          dry_bulk_density.gcm3 = DryWeight.g /SliceVolume.cm3 , #Dry bulk density
          Core_in.cm            = PipeLenght.cm  - CompactionIn ,#Compaction in cm to get core length in cm
          Pipe_in.cm            = PipeLenght.cm  - CompactionOut ,#Compaction in cm to get how deep we hammerred pipe in in cm
          Compaction_Correction_Value = Core_in.cm / Pipe_in.cm, 
          dry_bulk_density.gcm3_corrected = dry_bulk_density.gcm3 * Compaction_Correction_Value,
          CarbonDensity.gcm3    = dry_bulk_density.gcm3_corrected * C_percent/100,
          CarbonStock.Mgha      = CarbonDensity.gcm3 *100 * SliceLength.cm )# gcm3 times 100 gives Mgha

#HC cores were more coarsly sliced and we merge them with sc cores to get
#mean Dry_Fraction and C_percent per SampleID_hc:
sc_core_high <- sc_core %>% 
  dplyr::group_by(SampleID_hc) %>% 
  summarise_at(vars(C_percent, Dry_Fraction, dry_bulk_density.gcm3_corrected),mean)

View(sc_core_high)

#Join HC and CN data
hc_dry_core <- left_join(hc,sc_core_high, by = "SampleID_hc")#join data with Dry Fraction
names(hc_dry_core )

hc_dry_core_sum <- hc_dry_core %>% 
  #create function to express HC per Dry.Fraction and normalized to mean C_percent:
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x / hc_dry_core$Dry_Fraction ))%>%
  #computer HC stock based on slice parameters:
  mutate (HC_SliceLength.cm     = HC_DepthTo.cm - HC_DepthFrom.cm, #height of slice which is cylinder
          HC_SliceVolume.cm3    = square_cm2 * HC_SliceLength.cm,  #sliced square prims of 2by2cm or 3by3cm bottom
          HC_Fraction = Total / 1000000, #Fraction of total HC in each sample (mg / kg)
          HC_Density.gcm3    = dry_bulk_density.gcm3_corrected *HC_Fraction, #
          HC_Stock.Mgha      = HC_Density.gcm3 *100 * HC_SliceLength.cm )

View(hc_dry_core_sum)

#Sum up HC_Stock.Mgha per each core:
hc_end <- hc_dry_core_sum %>%
  dplyr::group_by(Site,SiteYear, transect, Elevation) %>% #Grouping by core
  summarise(HC_Core_Sum.Mgha = sum(HC_Stock.Mgha, na.rm = T))%>% #summing up per core
  dplyr::group_by (SiteYear) %>% #grouping per site
  summarise(Mean_Core_HC = mean(HC_Core_Sum.Mgha), #computing mean HC_Stock.Mgha per SiteYear
            sd_core = sd(HC_Core_Sum.Mgha, na.rm=T),
            N = n(),
            se_core = sd_core/sqrt(N))


hc_end
#####SiteYear  Mean_Core_HC sd_core   N se_core
#1 Stony1986         48.1    13.0     3    7.52
#2 Stony1996         29.0    33.9     3   19.6 
#3 Stony2006         42.8    25.4     3   14.6 