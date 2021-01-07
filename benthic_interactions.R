library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats) 
library(lme4)
library(writexl) #write_xlsx(list(Sheet1=df1,Sheet2=df2),"mydata.xlsx")
library(stringr)
library(metafor)
library(car)
library(nortest)
library(GAD)
library(AICcmodavg) #load AICcmodavg package
library(DescTools) #winsorizing
library(robustHD) #winsorizing
library(MASS)
library(survival)
library(RColorBrewer)
library(patchwork)
library(rmarkdown)
library(ggpubr)
library(rstatix)
library(broom)
library(WRS2)
library(lsmeans)
library(sjPlot)
library(devtools)
library(ppcor) # this pacakge computes partial and semipartial correlations.

coral <-read_excel("Data/benthic_interactions.xlsx",sheet="Data") #to read your file
str(coral)

coral$Transect <- as.factor(coral$Transect)

master <- coral%>%
  filter(picture=="yes")%>%  #Need pics to verify information
  filter(Algae_around=="yes")%>% #only look at corals with potential to be overgrown
  filter(Life_stage=="adult")%>% #need adults only for chi square
  mutate(Chemically = case_when(Sponge_species == "Aiolochroia crassa (purple)" ~"Defended",
                                Sponge_species == "Amphimedon compressa"~"Defended",
                                Sponge_species == "Aplysina cauliformis"~"Defended",
                                Sponge_species == "Aplysina spp." ~ "Defended",
                                Sponge_species == "Aplysina fulva"~"Defended",
                                Sponge_species == "Chondrilla caribensis" ~"Variable",
                                Sponge_species == "Cliona delitrix" ~"Undefended",
                                Sponge_species == "Desmapsamma anchorata" ~"Variable",
                                Sponge_species == "Ectyoplasia ferox" ~"Defended",
                                Sponge_species == "Ircinian campana" ~"Defended",
                                Sponge_species == "Ircinia strobilina"~"Defended",
                                Sponge_species == "Mycale laevis" ~"Undefended",
                                Sponge_species == "Verongula rigida" ~"Defended"))%>%
  mutate(Sponge_genus = case_when(Sponge_species == "Aiolochroia crassa (purple)" ~"Aiolochroia",
                                  Sponge_species == "Amphimedon compressa"~"Amphimedon",
                                  Sponge_species == "Aplysina cauliformis"~"Aplysina",
                                  Sponge_species == "Aplysina spp." ~ "Aplysina",
                                  Sponge_species == "Aplysina fulva"~"Aplysina",
                                  Sponge_species == "Chondrilla caribensis" ~"Chondrilla",
                                  Sponge_species == "Cliona delitrix" ~"Cliona",
                                  Sponge_species == "Desmapsamma anchorata" ~"Desmapsamma",
                                  Sponge_species == "Ectyoplasia ferox" ~"Ectyoplasia",
                                  Sponge_species == "Ircinian campana" ~"Ircinian",
                                  Sponge_species == "Ircinia strobilina"~"Ircinia",
                                  Sponge_species == "Mycale laevis" ~"Mycale",
                                  Sponge_species == "Verongula rigida" ~"Verongula"))%>%
  filter(!(Transect=="rogue" & ID=="h"))%>% #taking out outlier, the point also adds to more than 100%
  mutate("Coral_avg_diameter" = ((Coral_max1+Coral_max2)/2))%>%
  mutate("Coral_area" = (3*3.14*(Coral_avg_diameter/2)^2))%>% #hemisphere 3*pi*r^2
  mutate("Sponge_area" = (Sponge_max1*Sponge_max2))%>% #rectangle
  mutate("Algae_area" = (Algae_max1*Algae_max2))%>% #rectangle
  mutate("Perc_sponge" = ((Sponge_area/Coral_area)*100))%>%
  mutate("Perc_algae" = ((Algae_area/Coral_area)*100))%>%
  mutate(CoralArea_category=cut(Coral_area, breaks=c(0, 10000, 40000, 120000), labels=c("small","medium","large")))%>%
  select(-picture, -Location, -Depth.ft, -Life_stage,
         -Observer, -Algae_inter, -Algae_around, -Coral_max1, 
         -Coral_max2, -Sponge_max1, -Sponge_max2, -Algae_max1,
         -Algae_max2, -Notes, -Coral_avg_diameter)%>%
  mutate(chemical_num = case_when(Chemically == "Defended" ~ 1,
                                  Chemically == "Undefended" ~ 0))
      #So that I can add up numbers to see if at least there is 
      #at least one chemically defended sponge   

#Ircinia felix??
#making data sheet where each row represents one coral
master.onecoral <- master%>%
  group_by(Date, Transect, Quadrat, ID, Type, Coral, 
           Sponge_inter,
           Inferred_algae,
           use_for_size)%>%
  summarise(Sponge_area=sum(Sponge_area, na.rm = TRUE), #na.rm is to add even if NA
            Algae_area=sum(Algae_area, na.rm = TRUE),
            Perc_sponge=sum(Perc_sponge, na.rm = TRUE),
            Perc_algae=sum(Perc_algae, na.rm = TRUE),
            Coral_area=sum(Coral_area, na.rm = TRUE),
            chemical_num=sum(chemical_num, na.rm = TRUE),
            Num_sponge=max(Num_sponge, na.rm = TRUE))%>%
  mutate(defended_sponge = case_when(chemical_num >= 1 ~"yes",
                                     chemical_num < 1 ~"no"))
                                  
################## Which species #######################################
pairs <- master%>%
  filter(Sponge_inter=="sponge_win")%>% #Only want those with sponge overgrowth
  filter(Sponge_Identified=="yes")%>% #Some IDK the ID
  group_by(Coral, Sponge_species, Sponge_genus, Chemically)%>%
  summarise(n=n())
  #mutate(perc = (n / sum(n))*100)

ggplot(data=pairs, aes(x=Chemically, y=n, fill=Coral))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.position="top")+
  xlab("Sponge species")+
  ylab("Frequency of interactions")+
  scale_fill_discrete(name = "Coral", labels = c("P. astreoides", "S. siderea"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width=15))+
 theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

################## chi square test #######################################
#exploring all data
chemical <- master.onecoral%>%
  filter(Inferred_algae %in% c("none","algae_win"))%>%
  filter(Sponge_inter %in% c("sponge_win"))

#combined
table <- table(chemical$defended_sponge, chemical$Inferred_algae)
table
chisq.test(table, correct=FALSE)

#no difference in frequency of algae interaction with defended or undefended sponges

#official chi square
chidata <- master.onecoral%>%
  filter(Date %in% c("8720", "8820", "81020", "81120"))%>% #the days prior could be used too, but right now juves are included
  filter(Inferred_algae %in% c("none","algae_win"))%>%
  filter(Sponge_inter %in% c("sponge_win", "solo"))

table <- table(chidata$Sponge_inter, chidata$Inferred_algae)
table
chisq.test(table, correct=FALSE)

#combined and number of sponges
table <- table(chidata$Num_sponge, chidata$Inferred_algae)
table
chisq.test(table, correct=FALSE)
#not enough samples

past.chi <- chidata%>%
  filter(Coral == "past")

past_table <- table(past.chi$Sponge_inter, past.chi$Inferred_algae)
past_table
chisq.test(past_table, correct=FALSE)

ssid.chi <- chidata %>%
  filter(Coral=="ssid")

ssid_table <- table(ssid.chi$Sponge_inter, ssid.chi$Inferred_algae)
ssid_table
chisq.test(ssid_table, correct=FALSE)

################## H2: Algae and Sponge Overgrowth Correlation #######################################
explore <- master.onecoral%>%
  filter(Date %in% c("81120","81220", "81320", "81420"))%>% 
  filter(use_for_size=="yes") #corals that I got correct algae measurements for 
                                #and that had algae around the coral
                                #ex. corals with inferred algae that I didn't measure
                                #is a "no" for use_for_size
 
#OFFICIAL
sizes <- master.onecoral%>%
  filter(Date %in% c("81220","81420"))%>% #didn't include "81320" b/c these were pete tile ones, #"81120" were from 10x1 transects)
  filter(use_for_size=="yes")%>%
  filter(Transect!="rogue")

#frequency of sponges?
ggplot(explore)+
  aes(Perc_sponge, Perc_algae)+ #(x,y)
  geom_point(aes(color=Coral))+
  #geom_point(aes(shape=Coral, color=CoralArea_category), size=2, stroke=1, alpha = 1)+
  labs(x="Sponge overgrowth (%)",y="Algal overgrowth (%)")+
  #scale_shape_manual(values=c(2, 1))+
  #scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"))+
  theme_classic()
  geom_smooth(aes(y=Perc_algae), method = "lm", formula = y~(exp(-0.1*x)), se=FALSE, size=1, colour="black")

#Stats

model <- lm(Perc_sponge~Perc_algae+Coral_area, data=explore)
model1 <- lm(Perc_sponge~Perc_algae*Coral_area, data=explore)
model2 <- lm(Perc_sponge~Perc_algae, data=explore)
AIC(model, model1, model2)
Anova(model, model1) #same
Anova(model,model2) #same
Anova(model1,model2) #same
#no effect of coral size

#Finding best fit for data
y <- (explore$Perc_algae)
x <- (explore$Perc_sponge)
z <- (explore$Coral_area)

plot(x,y)

fit1 <- lm(y~x)
fit2 <- lm(y~x+I(x^2)) #BEST FIT #fit2 <- lm(y~poly(x,2,raw=TRUE)) #same thing this one, raw=true helps
fit3 <- lm(y~(exp(-0.05*x)))
fit4 <- lm(y~(exp(-0.1*x)))

AIC(fit1, fit2, fit3, fit4)
#Quadratic is best fit
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)

#Testing for normality
lillie.test(x)
lillie.test(y)
#data are not normal

cor.test(x, y, method=c("pearson"), conf.level=TRUE)
cor.test(x, y, method=c("kendall"))
cor.test(x, y, method=c("spearman")) #for monotonic data

x <- log(x)
y <- log(y)
plot(x,y)

shapiro.test(x)  #normal
shapiro.test(y)  #normal
qqPlot(x)
qqPlot(y)
