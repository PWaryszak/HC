#DRAW FIG 2 for HC MS:
#Relative abundance of total PAH, TPH and their components at three mangrove study sites
#(Stony Creek, Australia),
#representing three separate mangrove arrival times (Year 1986, 1996 and 2006)

#LOAD packages and DATA:========
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
#Stack Bar of PAH and TPH at three sites============
hc_dry_stack_pah <- hc_dry %>%
  group_by(SiteYear) %>%
  summarise_at( vars(Nap:Total_PAH_norm), sum, na.rm=T) %>% #Compute sum of all HC-s per column
  gather(PAH_Type, content, -SiteYear) #Prep long format for stack-plotting.

pah_stack <- ggplot(hc_dry_stack_pah, aes(x = PAH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "",y="PAH Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "none",
        strip.text=element_text(size=18))

pah_stack

#Stack Bar of TPH at three sites:
hc_dry_stack_tph <- hc_dry %>%
  group_by(SiteYear) %>%
  summarise_at( vars(TPH_06to09:TPH_37toUP, Total_TPH), sum, na.rm=T) %>% #Compute sum of all PAH-s per column
  gather(TPH_Type, content, -SiteYear) #Prep long format for stack-plotting.

hc_dry_stack_tph$TPH_Type <- factor(hc_dry_stack_tph$TPH_Type, 
                                    levels = c("TPH_06to09", "TPH_10to14" ,"TPH_15to28" ,"TPH_29to36" ,"TPH_37toUP", "Total_TPH"  ))

tph_stack <- ggplot(hc_dry_stack_tph, aes(x = TPH_Type, y = content, fill = SiteYear)) + 
  geom_bar(position = "fill",stat = "identity")+
  scale_y_continuous(labels = scales::percent_format())+
  
  labs(x = "",y="TPH Abundance (%)", fill = "Arrival Time")+
  theme(axis.text.x=element_text(vjust=0.5,size=12, face="bold", colour = "black", angle = 45),
        axis.text.y=element_text(size=8,colour = "black"),
        axis.title.y=element_text(size=14,colour = "black",face = "bold",),
        axis.title.x=element_text(size=14, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        strip.text=element_text(size=18))

tph_stack

#Bind tph_stack and pah_stack together:
g <- arrangeGrob(pah_stack,tph_stack, nrow=2) #generates g
ggsave(g, filename = "Fig2_PAH_TPH_Stack.jpeg", width = 17, height = 19,units = "cm",dpi = 600)

