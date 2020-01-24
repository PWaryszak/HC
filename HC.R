#Load DATA=======
library(tidyverse)
hc <- read.csv("hc.csv") #read in HydroCarbon data
sc <- read.csv("StonyCreek.csv") # read in CN data

#Compute mean dry.fraction + mean C_percent:
sc_high <- filter (sc, Elevation =="High" & habitat =="mangrove") %>% #all hc cores were taken at High elevation
  group_by(SampleID_hc) %>%
  summarise_at(vars(CompactionIn, CompactionOut, C_percent, Dry_Fraction),mean)

#Join HC and CN data
#Estimate PAH and TPH per dry mass:
hc_dry <- left_join(hc,sc_high, by = "SampleID_hc") #join data with Dry Fraction
  
hc_dry <- hc_dry %>%
  mutate_at(vars(Naphthalene:Total), function(x, na.rm=T)(x * hc_dry$Dry_Fraction )) #estimate mg/kg per dry mass

#Normalize PAG & TPH per Total Organic Carbon:=====
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

hc_dry_norm$Total_norm <- rowSums(hc_dry_norm[,42:60]) #Normalized Total Hydrocarbons
hc_dry_norm$TotPAH_norm <- rowSums(hc_dry_norm[,42:55])#Normalized Total PAH only

#PAH relative abundance:=======
hc_dry_norm2 <- hc_dry_norm %>%
  mutate_at( vars(Nap:BghiP), function(x, na.rm=T)( x / hc_dry_norm$TotPAH_norm *100)) %>%
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
      width = 19, 
      height = 8,
      units = "cm",
      dpi = 600)

#Total HC per Depth_Range, emulate Agusti's 2019 paper:=======
hc2 <- hc_dry_norm %>% 
  select(Year, Depth_Range, Total_norm) %>%
  group_by(Year,Depth_Range) %>%
  summarise(Mean_TPH = mean(Total_norm, na.rm=T),
            sd_TPH = sd(Total_norm, na.rm=T),
            N = n(),
            se_TPH = sd_TPH/sqrt(N))

ggplot(hc2,aes(x= reorder(Depth_Range, desc(Depth_Range)), y=Mean_TPH ,
               group = as.factor(Year)), color = reorder(as.factor(Year))) +
  labs(x = "Depth (cm)",y="TPH (mg/kg)")+
  geom_point(aes(color=as.factor(Year)),size = 2) +  geom_errorbar(aes(ymin = Mean_TPH -se_TPH, ymax = Mean_TPH+se_TPH), width = 0.2)+
  geom_line(aes(color=as.factor(Year)),size = 1.2, alpha = 0.4)+
  facet_grid( as.factor(Year) ~ .)+ 
  geom_hline(aes(yintercept = 275), linetype="dashed") + #275mg/kg is TPH toxicity threshold
  theme_bw()+
  coord_flip()+
  #ggtitle("TPH @ Stony Creek (mg/kg)")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black"),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=18))

ggsave(filename = "TPH_MeanContent.jpeg", 
       width = 17, 
       height = 9,
       units = "cm",
       dpi = 600)


#Stack Bar of PAH at three sites============
hc_dry_stack <- hc_dry_norm %>%
  group_by(SiteYear) %>%
  summarise_at( vars(Nap:BghiP), sum, na.rm=T) %>%
  gather(PAH_Type, content, -SiteYear)

ggplot(hc_dry_stack, aes(x = PAH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "PAH Type",y="Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "right",
        strip.text=element_text(size=18))

ggsave(filename = "PAH_Stack.jpeg", 
       width = 17, 
       height = 9,
       units = "cm",
       dpi = 600)

