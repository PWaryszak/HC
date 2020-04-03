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

#Plot PAH relative abundance:=======
hc_dry_norm_PAH <- hc_dry_norm_total %>%
  mutate_at( vars(Nap:BghiP), function(x, na.rm=T)( x / hc_dry_norm_total$TotPAH_norm *100)) %>% # % of single PAH against Total PAH
  select(Nap:BghiP)%>%
  gather(PAH_type,value) %>% 
  mutate (Site = "Stony")

ggplot(data = hc_dry_norm_PAH, aes( y = value, x = Site)) +
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

#Compute Mean values of PAH abundances:
hc_dry_norm_PAH2 <- hc_dry_norm_PAH %>%
  group_by(PAH_type)%>%
  summarise(mean(value, na.rm=T))

hc_dry_norm_PAH2


#TPH per HC_Depth_Range :=======
#With Threshold dotted line (as per Agusti's 2019 paper)
h1 <- hc_dry_norm_total %>% 
  select(Year, HC_Depth_Range, TotTPH_norm) %>%
  group_by(Year,HC_Depth_Range) %>%
  summarise(Mean_TPH = mean(TotTPH_norm, na.rm=T),
            sd_TPH = sd(TotTPH_norm, na.rm=T),
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

#PAH per Depth_Range:=======
h2 <- hc_dry_norm_total %>% 
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

hc_plots <- arrangeGrob(h1_plot,h2_plot, nrow=2) #generates g

ggsave(hc_plots, filename = "PAH_TPH_LINE.jpeg", 
       width = 17, 
       height = 19,
       units = "cm",
       dpi = 600)


#Stack Bar of PAH and TPH at three sites============
hc_dry_stack_pah <- hc_dry_norm_total %>%
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
hc_dry_stack_tph <- hc_dry_norm_total %>%
  group_by(SiteYear) %>%
  summarise_at( vars(TPH_06to09:TotTPH_norm), sum, na.rm=T) %>% #Compute sum of all PAH-s per column
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
g <- arrangeGrob(pah_stack,tph_stack, nrow=2) #generates g

ggsave(g, filename = "PAH_TPH_Stack.jpeg", 
       width = 17, 
       height = 19,
       units = "cm",
       dpi = 600)

#Compute times above toxicity level (High-Guideline-Value):===========
hc_above <- hc_dry_norm_total  %>%
  mutate(Threshold_TPH = 275,
         Total_TPH = TPH_06to09+TPH_10to14+TPH_15to28+TPH_29to36+TPH_36toMore,
         Times_Above = Total_TPH/ Threshold_TPH,
         Sample_Above = ifelse(Times_Above >=1, 1,0))# To compute frequency of times above threshold

min(hc_above$Times_Above, na.rm = T)# 0.0009771783
max(hc_above$Times_Above, na.rm = T)#55.45798

#Compute frequency of times above threshold:
good_total <- dim(subset(hc_above, !is.na(hc_above$Sample_Above))) #exclude NA-samples (43 of 51 were good records)
times_above <- dim(subset(hc_above, hc_above$Sample_Above > 0)) #34 times was above toxicity threshold
times_above*100/good_total # 72.54902 % of samples were above toxicity threshold


