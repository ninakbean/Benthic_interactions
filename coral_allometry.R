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

#Chi square table:
table <- table(coral$sponge.SA.type, coral$algae.SA.type)
table

#2021 data set
#             algae_win none
#none              117   20
#sponge_win         8    1

#2020 data set
#             algae_win none
#none              19   22
#sponge_win         5    2

chi <- matrix(c((117+19), (20+22), (8+5), (1+2)), ncol=2,byrow=TRUE)
colnames(chi) <- c("none","sponge win")
rownames(chi) <- c("algae win","none")
chi <- as.table(chi)
chi

chisq.test(chi, correct=FALSE)
#no interaction between coral and algae

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
  filter(coral.max.mm>40)

#making it so that each row is a coral, adding more specs
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
         P_SA=coral.perm.mm/coral.SA)

#Hypothesis 1
hypothesis1 <- master%>%
  filter(date %in% c(1222021, 1232021))

ggplot(hypothesis1)+
  aes(y=log(other.SA.perc), x=log(coral.SA))+ #coral.max.mm, (x,y) #all 6cm lower have more algal overgrowth
  geom_point(aes(color=other.SA.type))+
  #labs(x="Sponge Overgrowth (%)",y="Algal Overgrowth (%)")+
  theme_classic()+
  geom_smooth(method="lm")

plot <- lm(log(other.SA.perc)~log(coral.SA), data=frequency) #can't log 0s
summary(plot, type=III)

#Hypothesis 2
allometry <- master%>%
  filter(date %in% c(1182021, 1192021, 1202021))%>%  #need to fix these samples
  filter(other.perm.perc!=0)%>% #just for now!!!
  filter(coral.max.mm>40)

ggplot(allometry)+
  aes(x=log(P_SA), y=log(other.perm.perc))+ #(x,y)
  geom_point(aes(color=coral.max.mm))+
  #labs(x="Sponge Overgrowth (%)",y="Algal Overgrowth (%)")+
  theme_classic()+
  geom_smooth(method="lm")

model <- lm(log(other.perm.perc)~log(P_SA), data=allometry)
summary(model, type=III)

#exploring
ggplot(master)+
  aes(y=log(other.perm.perc), x=log(coral.perm.mm))+ #(x,y)
  geom_point(aes(color= coral.SA))+
  theme_classic()+
  geom_smooth(method="lm")


#Chi square table










#NOTES
#Dictyonella funicularis, 
#erksa, Klaus rus: 1981
#small brown porites is usually porites porites
#Porites asteroides is yellow morph
