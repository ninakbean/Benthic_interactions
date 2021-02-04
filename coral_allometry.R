library(popbio)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpmisc) #part of ggplot2? Formula on graoh
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
library(growthrates)

field <-read_excel("Data/coral_allometry.xlsx",sheet="Field") #to read your file
photo<-read_excel("Data/coral_allometry.xlsx",sheet="Photo") #to read your file
error <- read_excel("Data/coral_allometry.xlsx", sheet = "error")

#### LOOKING AT MEASUREMENT ERROR OF THE SAME CORAL AT THE SAME TIME ##########
error_coral <- error%>%
  filter(other.SA.current!="inferred")%>% #don't want inferred algal growth atm
  group_by(coral.ID, type)%>% #not including error with each algal measurement
  summarise(coral.perm.mm=sum(coral.perm.mm, na.rm=TRUE),
            other.perm.mm=sum(other.perm.mm, na.rm=TRUE),
            coral.max.mm=sum(coral.max.mm, na.rm=TRUE),
            coral.min.mm = sum(coral.min.mm, na.rm=TRUE),
            coral.height.mm = sum(coral.height.mm, na.rm=TRUE),
            other.max.mm = sum(other.max.mm, na.rm=TRUE),
            other.min.mm = sum(other.min.mm, na.rm=TRUE))%>%
  mutate(coral.SA = (2*3.14*coral.height.mm*(coral.max.mm/2)*(coral.min.mm/2)),
         other.SA = other.max.mm*other.min.mm)
         
analyzing_error <- error_coral %>%
  group_by(coral.ID)%>% #looking at error with coral perm, coral SA, other perm, other SA
  summarise(coral.SA.avg=mean(coral.SA),
            coral.SA.sd=sd(coral.SA),
            coral.SA.se=(sd(coral.SA)/sqrt(length(coral.SA))),
            other.SA.avg=mean(other.SA),
            other.SA.sd=sd(other.SA),
            other.SA.se=(sd(other.SA)/sqrt(length(other.SA))),
            coral.perm.avg=mean(coral.perm.mm),
            coral.perm.sd=sd(coral.perm.mm),
            coral.perm.se=(sd(coral.perm.mm)/sqrt(length(coral.perm.mm))),
            other.perm.avg=mean(other.perm.mm),
            other.perm.sd=sd(other.perm.mm),
            other.perm.se=(sd(other.perm.mm)/sqrt(length(other.perm.mm))))

################### SETTING UP ##############################################

#making data sheet with just coral specs
coral <- field%>%
  select(coral.ID, coral.perm.mm,	coral.max.mm,	coral.min.mm,	coral.height.mm, algae.SA.type, sponge.SA.type)%>%
  group_by(coral.ID, algae.SA.type, sponge.SA.type)%>%
  summarise(coral.perm.mm=sum(coral.perm.mm, na.rm=TRUE),
            coral.max.mm=sum(coral.max.mm, na.rm=TRUE),
            coral.min.mm = sum(coral.min.mm, na.rm=TRUE),
            coral.height.mm = sum(coral.height.mm, na.rm=TRUE))%>%
  mutate(coral.SA = (2*3.14*coral.height.mm*(coral.max.mm/2)*(coral.min.mm/2)),
         P_SA = (coral.perm.mm/coral.SA))%>%
  filter(algae.SA.type!="NA") #to get rid of blank columns

################### CHI SQUARE ##############################################

#Chi square table:
table <- table(coral$sponge.SA.type, coral$algae.SA.type)
table

#2021 data set
#             algae_win none
#none              135   23
#sponge_win         8    1

#2020 data set
#             algae_win none
#none              19   22
#sponge_win         5    2

chi <- matrix(c((135+19), (23+22), (8+5), (1+2)), ncol=2,byrow=TRUE)
colnames(chi) <- c("none","sponge win")
rownames(chi) <- c("algae win","none")
chi <- as.table(chi)
chi

chisq.test(chi, correct=FALSE)$expected
#one expected cell is less than 5
chisq.test(chi, correct=FALSE)
#no interaction between coral and algae
fisher.test(chi)
#no significant relationship between variables
#knowing the value of one DOES NOT helps to predict the value of the other

#percent algal overgrowth
154/(154+45+13+3)
#0.7162791

#percent sponge overgrowth
3/(154+45+13+3)
#0.01395349

#percent both overgrowth
45/(154+45+13+3)
#0.2093023

#percent no overgrowth
13/(154+45+13+3)
#0.06046512

#total number
(154+45+13+3)
#215

################### SETTING UP ##############################################

#delete these columns
field <- field%>%
  select(-coral.perm.mm, -perm.loc, -coral.max.mm, -coral.min.mm, -coral.height.mm)

#just so I can read the numbers but with them duplicated down
data <- merge(field, coral, by="coral.ID")

#taking out data I don't want
combined <- data%>%
  filter(data!= "photo")%>% #only include field measurements
  filter(coral.species!="SSID")%>% #only sampling PAST
  filter(other.SA.current!="inferred")%>% #don't want inferred algal growth atm
  filter(!coral.ID %in% c(25,28, 65,67))%>%  #need to fix these samples
  filter(coral.max.mm>=40)

#making it so that each row is a coral, adding more specs, includes sponges
master <- combined%>%
  group_by(date, coral.ID, other.SA.type, other.perm.type)%>%
  summarise(coral.max.mm=mean(coral.max.mm, na.rm=TRUE),
            coral.min.mm=mean(coral.min.mm, na.rm=TRUE),
            coral.height.mm=mean(coral.height.mm, na.rm=TRUE),
            coral.SA=mean(coral.SA, na.rm=TRUE),
            coral.perm.mm=sum(coral.perm.mm, na.rm=TRUE),
            other.max.mm=sum(other.max.mm, na.rm=TRUE),
            other.min.mm=sum(other.min.mm, na.rm=TRUE),
            other.perm.mm=sum(other.perm.mm, na.rm = TRUE))%>%
  mutate(other.SA = (other.max.mm*other.min.mm),
         other.SA.perc = (other.SA/coral.SA),
         other.perm.perc=(other.perm.mm/coral.perm.mm),
         P_SA=(coral.perm.mm/coral.SA),
         log.other.SA.perc = log(other.SA.perc+0.001),
         log.coral.SA = log(coral.SA+0.001))

################### AREA ##############################################

#Hypothesis 1: SA
# Step 1: Call the pdf command to start the plot
pdf(file = "Figs/other_SA.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code

ggplot(master)+
  aes(y=log(other.SA.perc+0.001), x=log(coral.SA+0.001))+ #coral.max.mm, (x,y) #all 6cm lower have more algal overgrowth
  geom_point(aes())+
  labs(x=(expression(paste("Log (coral area"," ", (mm^2), " ", "+ 0.001)"))),y=(expression(paste("Log (algal area"," ", (mm^2), " ", "+ 0.001)"))))+
  theme_classic(base_size=12)+
  geom_smooth(method="lm", color="black", size=1)

# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()

model <- lm(log(other.SA.perc+0.001)~log(coral.SA+0.001), data=master) #can't log 0s

out <- as.data.frame(augment(model))
cooks <- cooks.distance(model)
plot(cooks)
#after taking out 6 outliars, the slope looks the same, and p value
#Only the intercept looks like it changed a lot
#4/n(129) = 0.03100775
#y=-6.907755, x= 9.109768, ID= 36 
#y=-6.907755, x=9.466860, ID= 5464 
#y=-1.928448, x=8.234300, ID= 93 
#y=-2.337783, x=8.745125, ID= 111 
#y=-6.907755, x=10.472346, ID= 126 
#y=-2.639150, x=10.131420, ID= 37 NOT THIS ONE when taking sponges out

h1.1 <- SA_perc%>%
  filter(!coral.ID %in% c(36, 5464, 93, 111, 126, 37))

model1 <- lm(log(other.SA.perc+0.001)~log(coral.SA+0.001), data=h1.1) #can't log 0s

AIC(model, model1)
Anova(model, model1, type=3)

#(y=-0.44x-0.04, R2=0.27, p<0.001)

summary(model, type=II)
summary(model1, type=II)

################### PERIMETER ##############################################
allometry <- master%>%
  filter(date %in% c(1182021, 1192021, 1202021, 1262021, 1272021))%>%  #need to fix these samples
  #filter(other.perm.perc!=0)%>% #just for now!!!
  filter(coral.max.mm>40)%>%
  mutate(CoralArea_category=cut(coral.max.mm, breaks=c(0, 70, 110, 300), labels=c("small","medium","large")))

#to create size classes
hist(allometry$coral.max.mm) 

ggplot(allometry)+
  aes(y=log(other.perm.perc+0.001), x=log(coral.perm.mm+0.001))+ #(x,y)
  geom_point(aes(color= coral.SA))+
  theme_classic()+
  geom_smooth(method="lm")

model1 <- lm(log(other.perm.perc+0.001)~log(coral.perm.mm+0.001), data=allometry) #can't log 0s
summary(model1, type=II)

#(y=-0.13x-2, R2: 0.007, p=0.59)

#Hypothesis 1: SA
# Step 1: Call the pdf command to start the plot
pdf(file = "Figs/P_SA.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code

#Hypothesis 2: perm
ggplot(allometry.no.outlier3)+
  aes(x=(P_SA), y=(other.perm.perc))+ #(x,y)
  geom_point(aes(color=CoralArea_category))+
  labs(x="Coral perimeter:surface area ratio",y="Coral perimeter overgrowth \n by algae and sponges (%)")+
  scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"),
                     name="Coral Size", 
                     labels = c("small", "medium", "large"))+
  theme_classic()+
  geom_smooth(method="lm", color="black")

# Step 3: Run dev.off() to create the file!
#“null device”<- saying now create plots in the main R plotting window again.
dev.off()

#took out 3 outliers

model <- lm((other.perm.perc)~(P_SA), data=allometry)

out <- as.data.frame(augment(model))
cooks <- cooks.distance(model)
plot(cooks)
#4/n(43) = 0.09302326
#y=0.309178744, x=0.0320328313, ID=146
#y=0.090000000, x=0.0215766489, ID=41
#y=0.252525253, x=0.0004015495, ID=53

allometry.no.outlier1 <- allometry%>%
  filter(!coral.ID %in% c(146))

allometry.no.outlier2 <- allometry%>%
  filter(!coral.ID %in% c(146, 41))

allometry.no.outlier3 <- allometry%>%
  filter(!coral.ID %in% c(146, 41, 53))

model.no.outlier1 <-  lm(other.perm.perc~P_SA, data=allometry.no.outlier1)
model.no.outlier2 <-  lm(other.perm.perc~P_SA, data=allometry.no.outlier2)
model.no.outlier3 <-  lm(other.perm.perc~P_SA, data=allometry.no.outlier3)

Anova(model, model.no.outlier1, type=3)
Anova(model.no.outlier1, model.no.outlier2, type=3)
Anova(model.no.outlier2, model.no.outlier3, type=3)

AIC(model, model.no.outlier1, model.no.outlier2, model.no.outlier3)

summary(model, type=II)
summary(model.no.outlier3, type=II)

#(y=7.66x+0.06, R2=0.37, p<0.001)


ggplot(allometry.no.outlier3)+
  aes(x=(P_SA), y=(other.perm.perc))+ #(x,y)
  geom_point(aes(color=coral.max.mm))+
  #labs(x="Sponge Overgrowth (%)",y="Algal Overgrowth (%)")+
  theme_classic()+
  geom_smooth(method="lm")
















#NOTES
#Dictyonella funicularis, 
#erksa, Klaus rus: 1981
#small brown porites is usually porites porites
#Porites asteroides is yellow morph
