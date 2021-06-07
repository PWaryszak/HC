#DRAW FIG:
#Mean TPH and PAH content (±SE) in sediments extruded at three mangrove restoration sites (mangroves arrived in 1986,1996 and in 2006).
#The vertical dotted lined indicates SQGV toxicity level associated with adverse biological effect
#(280 mg kg-1 for of TPHs and at 50 mg kg-1 for PAHs).
#Load Packages=======

# Start the code
rm(list = ls(all = TRUE))

library("tidyverse")
library("gridExtra")
library("grid")
library(ggpubr)


#Load DATA:=======
hc <- read.csv("hc1.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek2.csv") # read in C&N data, comes from Maria's CN-analysis of 9 replicate cores

sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation, keep high only
  group_by(SampleID_hc) %>%
  summarise_at(vars( C_percent, Dry_Fraction), mean)#mean values as hc samples were more coarse than sc

#Joined hc and sc_high data:
hc_joined <- read.csv("hc_joined.csv") #joined H and OC data with C_percent and  Dry_Fraction

#Estimate PAH and TPH per dry mass by dividing content in wet sample (off hc) by mean dry fraction (off sc_high):
hc_dry <- hc_joined %>%
  #Total = sum of all TPH-s, Total_TPH is sum of TPH in dry sediments
  mutate_at(vars(Naphthalene:Total_HC), function(x, na.rm=T)(x / hc_joined$Dry_Fraction ))%>% #estimate mg/kg per dry mass
  
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
         TPH_37toUP = TPH_C_37UP ) 

#Double check:
check<- hc_dry %>% mutate(Total_TPH = rowSums(hc_dry[,(which(names(hc_dry)=='TPH_C_06to09'):which(names(hc_dry)=='TPH_C_37UP'))])) #sum up all TPH-s to double check if it is the same as Total_HC (it is!)
check$Total_TPH == hc_dry$Total_HC #All good! TRUE!

#write.csv(hc_dry, file = "HCdry.csv",row.names = F)
#TWO PANELS version Plot =========
names(hc_dry)

#Check values at 75-100:
from75up <- filter(hc_dry, HC_Depth_Range == "75_100") %>% select(Total_HC, SampleID)
from75up


#TPH plot:
p1 <- hc_dry %>% 
  select(Year,  HC_Depth_Range, Total_HC) %>%
  group_by(HC_Depth_Range) %>%
  summarise(Mean_TPH = mean(Total_HC, na.rm=T),
            sd_TPH = sd(Total_HC, na.rm=T),
            N = n(),
            se_TPH = sd_TPH/sqrt(N))

p1_plot <-  ggplot(p1,aes(x= reorder(HC_Depth_Range, desc(HC_Depth_Range)), y=Mean_TPH, group=1))+ #group = Mean_TPH, fill = reorder(as.factor(HC_Depth_Range))))+
  
  labs(y = bquote('Total Petroleum Hydrocarbons   ' (mg~kg^-1)), x="Sediment depth (cm)") +
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = Mean_TPH -se_TPH, ymax = Mean_TPH+se_TPH), width = 0.2)+
  geom_line(size = 1.2)+ #to turn line grey use , alpha = 0.3
  scale_y_continuous(limits = c(0,60000), breaks = c(0,10000,20000,30000,40000,50000,60000))+
  geom_hline(aes(yintercept = 50), linetype="dashed") + #275mg/kg is TPH toxicity threshold
  theme_bw()+
  coord_flip()+
  theme(axis.text.x=element_text(vjust=0.5,size=14, face="bold", colour = "black"),
        axis.text.y=element_text(size=14,colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        axis.title.y=element_text(size=16,colour = "black",face = "bold",),
        axis.title.x=element_text(size=16, face = "bold", hjust = 0.5),
        legend.position = "none")

p1_plot

#PAH:

p2 <- hc_dry %>% 
  select(Year, HC_Depth_Range, Total_PAH_norm) %>%
  group_by(HC_Depth_Range) %>%
  summarise(Mean_PAH = mean(Total_PAH_norm, na.rm=T),
            sd_PAH = sd(Total_PAH_norm, na.rm=T),
            N = n(),
            se_PAH = sd_PAH/sqrt(N)) %>% mutate(HC_Depth_Range2 = HC_Depth_Range) 


p2_plot <- ggplot(p2,aes(x= reorder(HC_Depth_Range2, desc(HC_Depth_Range2)), y=Mean_PAH,group = 1)) +
  labs(y = bquote('Total Polycyclic Aromatic Hydrocarbons   ' (mg~kg^-1)), x="      ") + #Sediment depth (cm) = Paste to use for vertical plot arrangment
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean_PAH -se_PAH, ymax = Mean_PAH+se_PAH), width = 0.2)+
  geom_line(size = 1.2)+ #to mke line grey use , alpha = 0.3
  scale_y_continuous(limits = c(0,51), breaks = c(0,10,20,30,40,50))+
  geom_hline(aes(yintercept = 50), linetype="dashed") + #50 mg/kg is PAH toxicity threshold
  #or 50000/kg as per web: https://www.waterquality.gov.au/anz-guidelines/guideline-values/default/sediment-quality-toxicants
  theme_bw()+
  coord_flip()+
  #scale_color_manual(values = cols)+
  
  theme(axis.text.x=element_text(vjust=0.5,size=14, face="bold", colour = "black"),
        axis.text.y=element_text(size=14,colour = "white"), #"white to use for vertical plot arrangment
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        axis.title.y=element_text(size=16,colour = "black",face = "bold",),
        axis.title.x=element_text(size=16, face = "bold", hjust = 0.5),
        legend.position = "none")
p2_plot

#Arrange plots ertically and save
p_plots  <- ggarrange(p1_plot,p2_plot, ncol=1, 
                      labels = c("a)", "b)"),
                      label.x = 0,
                      label.y = 1)
                      
p_plots
#ggsave(p_plots, filename = "FIG2_TWO_PANELS_60000.pdf", width = 17, height = 19, units = "cm", dpi = 600)


#Arrange plots horizontally as per reviewer request an save: (Turn off labels in p2_plot)
p_plots_horizon  <- ggarrange(p1_plot,p2_plot, ncol=2, 
                      labels = c("a)", "b)"),
                      label.x = 0,
                      label.y = 1)

p_plots_horizon
ggsave(p_plots_horizon, filename = "FIG2_TWO_PANELS_HORIZON.pdf", width = 32, height = 12, units = "cm", dpi = 600)


#PAH/TPH Minimum and maximum values=======
round(range(hc_dry$Total_HC),1) #1.6 61657.1
max(hc_dry$Total_HC)/280 #220 times above toxicity threshold
max_row <- filter (hc_dry, Total_HC > 61657)
max_row[1,1] #Stony.86 T1-High-30_50

range(hc_dry$Total_PAH_norm)# 0.4198813 35.5282081
