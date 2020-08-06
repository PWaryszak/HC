#DRAW FIG:
#Mean TPH and PAH content (±SE) in sediments extruded at three mangrove restoration sites (mangroves arrived in 1986,1996 and in 2006).
#The vertical dotted lined indicates SQGV toxicity level associated with adverse biological effect
#(280 mg kg-1 for of TPHs and at 50 mg kg-1 for PAHs).

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
  mutate(TotTPH = rowSums(hc_dry[,(which(names(hc_dry)=='TPH_C_06to09'):which(names(hc_dry)=='TPH_C_37UP'))])) %>%#sum up all TPH-s to double check if it is the same as Total_HC (it is!)
  
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

#Plot TPH per HC_Depth_Range :=======
#With Threshold dotted line (as per Agusti's 2019 paper)
h1 <- hc_dry %>% 
  select(Year, HC_Depth_Range, Total_HC) %>%
  group_by(Year,HC_Depth_Range) %>%
  summarise(Mean_TPH = mean(Total_HC, na.rm=T),
            sd_TPH = sd(Total_HC, na.rm=T),
            N = n(),
            se_TPH = sd_TPH/sqrt(N))

h1_plot <- ggplot(h1,aes(x= reorder(HC_Depth_Range, desc(HC_Depth_Range)), y=Mean_TPH ,
                         group = as.factor(Year)), color = reorder(as.factor(Year))) +
  labs(y = bquote('Total Petroleum Hydrocarbons   ' (mg~kg^-1)), x="Sediment depth (cm)") +
  geom_point(aes(color=as.factor(Year)),size = 2) +  geom_errorbar(aes(ymin = Mean_TPH -se_TPH, ymax = Mean_TPH+se_TPH), width = 0.2)+
  geom_line(aes(color=as.factor(Year)),size = 1.2, alpha = 0.4)+
  scale_y_continuous(limits = c(0,60000), breaks = c(0,10000,20000,30000,40000,50000,60000))+
  
  facet_grid( as.factor(Year) ~ .)+ 
  geom_hline(aes(yintercept = 275), linetype="dashed") + #275mg/kg is TPH toxicity threshold
  theme_bw()+
  coord_flip()+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black"),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.background = element_rect(fill="white"),
        strip.text=element_text(size=18))

h1_plot


#PAH per Depth_Range:=======
h2 <- hc_dry %>% 
  select(Year, HC_Depth_Range, Total_PAH_norm) %>%
  group_by(Year,HC_Depth_Range) %>%
  summarise(Mean_PAH = mean(Total_PAH_norm, na.rm=T),
            sd_PAH = sd(Total_PAH_norm, na.rm=T),
            N = n(),
            se_PAH = sd_PAH/sqrt(N))

h2_plot <- ggplot(h2,aes(x= reorder(HC_Depth_Range, desc(HC_Depth_Range)), y=Mean_PAH ,
                         group = as.factor(Year)), color = reorder(as.factor(Year))) +
  labs(y = bquote('Total Polyaromatic Hydrocarbons   ' (mg~kg^-1)), x="Sediment depth (cm)") +
  #labs(x = "Depth (cm)",y="PAH (mg/kg)")+
  geom_point(aes(color=as.factor(Year)),size = 2) +  geom_errorbar(aes(ymin = Mean_PAH -se_PAH, ymax = Mean_PAH+se_PAH), width = 0.2)+
  geom_line(aes(color=as.factor(Year)),size = 1.2, alpha = 0.4)+
  scale_y_continuous(limits = c(0,60), breaks = c(0,10,20,30,40,50,60))+
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
        strip.background = element_rect(fill="white"),
        strip.text=element_text(size=18))
h2_plot

#Arrange plots and save:
hc_plots <- arrangeGrob(h1_plot,h2_plot, nrow=2) #Bind two plots together
ggsave(hc_plots, filename = "FIG3_PAH_TPH_LINE.jpeg", width = 17, height = 19, units = "cm", dpi = 600)

