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

coral <-read_excel("Data/benthic_interactions.xlsx",sheet="Data") #to read your file
str(coral)

coral$Transect <- as.factor(coral$Transect)

master <- coral%>%
  filter(picture=="yes")%>%  #Need pics to verify information
  #filter(Algae_around=="yes")%>% #only look at corals with potential to be overgrown, didn't do this in Jan trip
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
  filter(!(Transect=="rogue" & ID=="h"))%>% #taking out outlier, the point also adds to more than 100% probably incorrect
  mutate("Coral_avg_diameter" = ((Coral_max1+Coral_max2)/2))%>%
  mutate("Coral_area" = (Coral_max1*Coral_max2*3.14))%>% #ellipse
  #hemisphere 3*pi*r^2 <- only one diameter included
  mutate("Sponge_area" = (Sponge_max1*Sponge_max2))%>% #rectangle
  mutate("Algae_area" = (Algae_max1*Algae_max2))%>% #rectangle
  mutate("Perc_sponge" = ((Sponge_area/Coral_area)*100))%>%
  mutate("Perc_algae" = ((Algae_area/Coral_area)*100))%>%
  #don't need this...?? mutate("Sponge_area" = replace_na(Sponge_area,0))%>% #b/c summing these columns doesn't work with NAs
  #don't need this...?? mutate("Algae_area" = replace_na(Algae_area,0))%>% #b/c summing these columns doesn't work with NAs
  mutate("total_overgrowth" = (Sponge_area+Algae_area))%>%
  mutate(CoralArea_category=cut(Coral_area, breaks=c(0, 10000, 40000, 120000), labels=c("small","medium","large")))%>%
  select(-picture,  -Life_stage,
         -Observer, -Algae_inter, -Algae_around, -Notes, -Coral_avg_diameter)%>%
  mutate(chemical_num = case_when(Chemically == "Defended" ~ 1,
                                  Chemically == "Undefended" ~ 0))
    #So that I can add up numbers to see if at least there is 
    #at least one chemically defended sponge   

#Ircinia felix??
#making data sheet where each row represents one coral
#can not have sponge species though because there may be multiple
master.onecoral <- master%>%
  group_by(Date, Transect, Quadrat, ID, Type, Coral, 
           Sponge_inter,
           Inferred_algae,
           use_for_size, CoralArea_category, Coral_max1, Coral_max2,
           Sponge_max1, Sponge_max2, Algae_max1,
           Algae_max2)%>%
  summarise(Sponge_area=sum(Sponge_area, na.rm = TRUE), #na.rm is to add even if NA
            Algae_area=sum(Algae_area, na.rm = TRUE),
            Perc_sponge=sum(Perc_sponge, na.rm = TRUE),
            Perc_algae=sum(Perc_algae, na.rm = TRUE),
            Coral_area=sum(Coral_area, na.rm = TRUE),
            chemical_num=sum(chemical_num, na.rm = TRUE),
            Num_sponge=max(Num_sponge, na.rm = TRUE),
            total_overgrowth=sum(total_overgrowth, na.rm = TRUE))%>%
  mutate(defended_sponge = case_when(chemical_num >= 1 ~"yes",
                                     chemical_num < 1 ~"no"))
                                  
################## Which species #######################################
#some corals are replicated
pairs <- master%>%
  filter(Sponge_inter=="sponge_win")%>% #Only want those with sponge overgrowth
  filter(Sponge_Identified=="yes")%>% #Some IDK the ID
  group_by(Coral, Sponge_species, Sponge_genus, Chemically)%>%
  summarise(n=n())
  #mutate(perc = (n / sum(n))*100)

ggplot(data=pairs, aes(x=Sponge_species, y=n, fill=Coral))+
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
#with at least one defended sponge
table <- table(chemical$defended_sponge, chemical$Inferred_algae)
table
chisq.test(table, correct=FALSE)

#no difference in frequency of algae interaction with defended or undefended sponges

#official chi square
chidata <- master.onecoral%>%
  filter(Date %in% c("8720", "8820", "81020", "81120"))%>% #the days prior could be used too, but right now juves are NOT included
  filter(Inferred_algae %in% c("none","algae_win"))%>%
  filter(Sponge_inter %in% c("sponge_win", "solo"))

past.chi <- chidata%>%
  filter(Coral == "past")

ssid.chi <- chidata %>%
  filter(Coral=="ssid")

table <- table(chidata$Sponge_inter, chidata$Inferred_algae)
table
chisq.test(table, correct=FALSE)


#combined and number of sponges
table <- table(chidata$Num_sponge, chidata$Inferred_algae)
table
chisq.test(table, correct=FALSE)
#not enough samples

past_table <- table(past.chi$Sponge_inter, past.chi$Inferred_algae)
past_table
chisq.test(past_table, correct=FALSE)


ssid_table <- table(ssid.chi$Sponge_inter, ssid.chi$Inferred_algae)
ssid_table
chisq.test(ssid_table, correct=FALSE)

################## H2: Algae and Sponge Overgrowth Correlation #######################################
#taking out the corals with more than 1 sponge
#singles_only <- master%>%
#  filter(Sponge_num=="single")%>% 
#  filter(use_for_size=="yes")

#OFFICIAL
#includes corals with more than one sponge, but can't see the sponge species
explore <- master.onecoral%>%
  filter(Date %in% c("81120","81220", "81420"))%>% #didn't include "81320" b/c these were pete tile ones, #"81120" were from 10x1 transects)
  filter(use_for_size=="yes") %>%   #81320: "picked" out corals near tiles
  filter(Transect!="rogue") #were not constrained

#write_xlsx(list(data=explore),"Data/correlation.xlsx")

#Only corals that have both algae and sponge
#corals that I got correct algae measurements for 
#and that had algae around the coral
#ex. corals with inferred algae that I didn't measure
#is a "no" for use_for_size

explore1 <- explore%>%
  filter(!(Date=="81420" & Transect == "1" & Quadrat == "10" & ID == "a"))

explore2 <- explore%>%
  filter(!(Date=="81420" & Transect == "1" & Quadrat == "10" & ID == "a"))%>%
  filter(!(Date=="81420" & Transect == "4" & Quadrat == "5" & ID == "a"))

model <- lm(total_overgrowth~Coral_area, data=explore)
model1 <- lm(total_overgrowth~Coral_area, data=explore1)
model2 <- lm(total_overgrowth~Coral_area, data=explore2)

AIC(model, model1, model2)
Anova(model, model1, type=3)
Anova(model1, model2, type=3)

summary(model2, type=3)

explore_graph <- explore2%>%
  group_by(Date, Transect, Quadrat, ID, Coral_area)%>%
  gather("Type", "overgrowth", c(Sponge_area, Algae_area)) #to seperate sponge and algal growth

out <- as.data.frame(augment(model))
cooks <- cooks.distance(model)
plot(cooks)
#4/n(16) = 0.25
#y=66680, x=171044.239, Date TQID: 81420 1,10,a
#y=6946, x=104350.639, Date TQID: 81420, 4,5,a

# Step 1: Call the pdf command to start the plot
#pdf(file = "Figs/correlation.pdf",   # The directory you want to save the file in
#    width = 7, # The width of the plot in inches
#    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code

#CORRELATION
ggplot(explore)+
  aes(Perc_sponge, Perc_algae)+ #(x,y)
  geom_point(aes(shape=Coral, color=CoralArea_category), size=2, stroke=1, alpha = 1)+
  labs(x="Sponge Overgrowth (%)",y="Algal Overgrowth (%)")+
  scale_shape_manual(values=c(2, 1),
                     name = "Coral Species",
                     labels = c((expression(paste(italic("Porites astreoides")))),
                               expression(paste(italic("Sidastrea siderea")))))+
  scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"),
                     name="Coral Size", 
                     labels = c("small", "medium", "large"))+
  theme_classic()

# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()


# Step 1: Call the pdf command to start the plot
pdf(file = "Figs/coral_size.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code

#frequency of sponges?
graph1 <- ggplot(explore)+
  aes(Coral_area, total_overgrowth, color=Coral)+ #(x,y)
  geom_point(aes(color=Coral))+
  #geom_point(aes(shape=Coral, color=CoralArea_category), size=2, stroke=1, alpha = 1)+
  labs(x="Coral Area",y="Total Area of Overgrowth \n from Algae and Sponges")+
  #scale_shape_manual(values=c(2, 1))+
  #scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"))+
  theme_classic()+
  geom_smooth(method = "lm", se=FALSE, size=1, colour="black")+
  scale_color_manual(values=c("#F8A42F", "#FF4605"),name="Coral Species", 
                     labels = c((expression(paste(italic("Porites astreoides")))),
                                expression(paste(italic("Sidastrea siderea")))))+
  theme(legend.position="top")
  #theme(legend.text.align = 0,
  #      legend.position = c(0.8, 0.2))

graph1

# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()



# Step 1: Call the pdf command to start the plot
pdf(file = "Figs/coral_type.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code

graph2 <- ggplot(explore_graph)+
  aes(Coral_area, overgrowth, group=Type)+ #(x,y)
  geom_point(aes(color=Type))+
  geom_smooth(aes(group=Type, color = factor(Type)), method="lm", se=FALSE, size=1)+
  #geom_point(aes(shape=Coral, color=CoralArea_category), size=2, stroke=1, alpha = 1)+
  labs(x="Coral Area",y="Area of Overgrowth \n from Algae and Sponges")+
  scale_color_manual(values=c("#006400", "#964B00"),name="Group", 
                     labels = c("Algae", "Sponge"))+
  theme_classic()+
  theme(legend.position="top")

graph2
# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()

# Step 1: Call the pdf command to start the plot
pdf(file = "Figs/combined.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code

(graph1)|(graph2)

# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()


ggscatter(
  data=explore_graph, x = "Coral_area", y = "overgrowth",
  color = "Type", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Type)
  )

#algae: y=0.39x-3600, R2: 0.83
#sponge: y=0.038x+400, R2: 0.055


model <- lm(overgrowth~Coral_area+Type, data=explore_graph)
model1 <- lm(overgrowth~Coral_area*Type, data=explore_graph) #there is an interaction
model2 <- lm(overgrowth~Coral_area, data=explore_graph)      #means there is more algal overgrowth, but also has to do with sponge overgrowth

AIC(model, model1, model2)
Anova(model, model1) 
Anova(model,model2)
Anova(model1,model2) #the interaction is significant

summary(model1)
Anova(model1, type=3)
#algae: y=0.15x-880, R2: 0.34
#sponge: y=0.057x+59, R2: 0.23

#Stats
model <- lm(Perc_sponge~Perc_algae+Coral_area, data=explore)
model1 <- lm(Perc_sponge~Perc_algae*Coral_area, data=explore)
model2 <- lm(Perc_sponge~Perc_algae, data=explore) #no interaction with coral_area
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
