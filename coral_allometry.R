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

################### SETTING UP ##############################################
#making data sheet with just coral specs
coral <- field%>%
  select(coral.ID, coral.perm.mm,	coral.max.mm,	coral.min.mm,	coral.height.mm, 
         algae.SA.type, sponge.SA.type, coral.species)%>%
  group_by(coral.ID, algae.SA.type, sponge.SA.type, coral.species)%>%
  summarise(coral.perm.mm=sum(coral.perm.mm, na.rm=TRUE),
            coral.max.mm=sum(coral.max.mm, na.rm=TRUE),
            coral.min.mm = sum(coral.min.mm, na.rm=TRUE),
            coral.height.mm = sum(coral.height.mm, na.rm=TRUE))%>%
  mutate(coral.SA.ellipsoid = (2*3.14*coral.height.mm*(coral.max.mm/2)*(coral.min.mm/2)),
         coral.SA.ellipse = (3.14*coral.max.mm*coral.min.mm),
         P_SA = (coral.perm.mm/coral.SA.ellipsoid))%>%
  filter(algae.SA.type!="NA")%>% #to get rid of blank columns with no info on it
                              #because I didn't copy SA.type info for all cells of a coral
  ungroup(coral.ID, algae.SA.type, sponge.SA.type, coral.species)%>%
  mutate(CoralArea_category=cut_number(coral.SA.ellipsoid, 3, labels=c("small","medium","large")))

#cut: when not given explicit break points divides values into bins of same width,
#usually won't contain an equal number of items
#cut_number: cuts it so same number in each group
#to optimize bins for the least variance...devtools::install_github("moodymudskipper/cutr")...table(cutr::smart_cut(x, list(4, "balanced"), "g"))

str(coral)
#delete these columns

field_edit <- field%>%
  select(-coral.perm.mm, -perm.loc, -coral.max.mm, -coral.min.mm, -coral.height.mm,
         -algae.SA.type, -sponge.SA.type, -coral.species)

#just so I can read the numbers but with them duplicated down
data <- merge(field_edit, coral, by="coral.ID")

#taking out data I don't want
combined <- data%>%
  filter(data!= "photo")%>% #only include field measurements, doesn't do anything now
 # filter(coral.species!="SSID")%>%  #only sampling PAST
  filter(!coral.ID %in% c(65,67))%>%  #65, no height measure, look at pics?
                                      #67, really big fragmented colony, only inferred colony size? reconfirm with pics
  filter(coral.max.mm>=40) #only looking at adults

#making it so that each row is a coral, adding more specs, includes sponges
master <- combined%>%
  group_by(date, coral.ID, other.SA.type, CoralArea_category, other.perm.type, other.SA.current, algae.SA.type, sponge.SA.type,
           coral.species)%>%
  summarise(coral.max.mm=mean(coral.max.mm, na.rm=TRUE),
            coral.min.mm=mean(coral.min.mm, na.rm=TRUE),
            coral.height.mm=mean(coral.height.mm, na.rm=TRUE),
            coral.SA.ellipsoid=mean(coral.SA.ellipsoid, na.rm=TRUE),
            coral.SA.ellipse=mean(coral.SA.ellipse, na.rm=TRUE),
            coral.perm.mm=mean(coral.perm.mm, na.rm=TRUE),
            other.max.mm=sum(other.max.mm, na.rm=TRUE),
            other.min.mm=sum(other.min.mm, na.rm=TRUE),
            other.perm.mm=sum(other.perm.mm, na.rm = TRUE))%>%
mutate(other.SA = (other.max.mm*other.min.mm),
         other.SA.perc.ellipsoid = ((other.SA/coral.SA.ellipsoid)*100),
         other.SA.perc.ellipse = ((other.SA/coral.SA.ellipse)*100),
         other.perm.perc=((other.perm.mm/coral.perm.mm)*100),
         P_SA=(coral.perm.mm/coral.SA.ellipsoid))%>%
  ungroup(date, coral.ID, other.SA.type, CoralArea_category, other.perm.type, other.SA.current, algae.SA.type, sponge.SA.type,
          coral.species)

str(master)

#splitting up data sheet so I can parse out the sponges
SA <- master%>%
  select(date, coral.ID, coral.SA.ellipsoid, coral.SA.ellipse, coral.max.mm,
         coral.perm.mm, P_SA, 
         other.SA.type, other.SA.current, other.SA, other.SA.perc.ellipsoid,
         other.SA.perc.ellipse,
         algae.SA.type, sponge.SA.type, CoralArea_category,
         coral.species)%>%
  filter(other.SA.current!="inferred")%>% #don't want inferred algal growth atm
  filter(other.SA.type!="none")%>% #taking out corals without any SA overgrowth
  group_by(date, coral.ID, coral.max.mm,coral.SA.ellipsoid, coral.SA.ellipse, coral.perm.mm, P_SA, other.SA.type,
           algae.SA.type, sponge.SA.type, CoralArea_category,
           coral.species)%>%
  summarise(other.SA=sum(other.SA, na.rm=TRUE),
          other.SA.perc.ellipsoid=sum(other.SA.perc.ellipsoid, na.rm=TRUE),
          other.SA.perc.ellipse=sum(other.SA.perc.ellipse,na.rm=TRUE))%>%
  ungroup(date, coral.ID,coral.max.mm, CoralArea_category, coral.SA.ellipsoid, coral.SA.ellipse, coral.perm.mm, P_SA, other.SA.type,
          algae.SA.type, sponge.SA.type,
          coral.species)

#I just want a data sheet where both algae and sponge are overgrowing
#If I want to combine data from last trip, I need to assume the coral is an ellipse
SA.correlation <- SA%>%
  filter(algae.SA.type=="algae_win" & sponge.SA.type=="sponge_win")%>%
 # select(other.SA)
  select(date, coral.ID, CoralArea_category, other.SA.type, coral.SA.ellipse, other.SA,
         coral.species)%>% #took out other.SA.perc.ellipse,,
  filter(!coral.ID %in% c(91, 77)) #91 & 77 is inferred sponge overgrowth

#Making data from wide to long
#algae=SA of algae on coral
#sponge=SA of sponge on coral
SA.correlation.wide <- spread(SA.correlation, other.SA.type, other.SA)

SA.2021 <- SA%>%
  filter(!date %in% c("81120","81220", "81420"))%>%
  filter(coral.species=="PAST")
  #filter(other.SA.type=="algae") #taking sponges out because not enough (5)

perm <- master%>%
  select(date, coral.ID, CoralArea_category, coral.SA.ellipsoid, coral.perm.mm, P_SA, 
         other.perm.type, other.perm.mm,other.perm.perc,
         coral.species)%>%
  filter(other.perm.type!="NA")%>% #taking out corals I didn't get perimeters for
  group_by(date, coral.ID,CoralArea_category, coral.SA.ellipsoid, coral.perm.mm, P_SA, other.perm.type,
           coral.species)%>%
  summarise(other.perm.mm=sum(other.perm.mm, na.rm=TRUE),
            other.perm.perc=sum(other.perm.perc, na.rm=TRUE))
  #filter(other.perm.type=="algae") #taking sponges out because not enough (8)

################### CHI SQUARE ##############################################
#2020 data set of PAST
#             algae_win none
#none              19   22
#sponge_win         5    2

#2020 data set (combined)
#             algae_win none
#solo              45   25
#sponge_win        22    3

#Chi square table:
table <- table(coral$sponge.SA.type, coral$algae.SA.type)
table

#COMBINED
chi <- matrix(c((135+45), (23+25), (8+22), (1+3)), ncol=2,byrow=TRUE)
colnames(chi) <- c("none","sponge win")
rownames(chi) <- c("algae win","none")
chi <- as.table(chi)
chi
chisq.test(chi, correct=FALSE)$expected
#one expected cell is less than 5
chisq.test(chi, correct=FALSE)
#no interaction between coral and algae when combined

#JUST PAST
chi <- matrix(c((135+19), (23+22), (8+5), (1+2)), ncol=2,byrow=TRUE)
colnames(chi) <- c("none","sponge win")
rownames(chi) <- c("algae win","none")
chi <- as.table(chi)
chi
chisq.test(chi, correct=FALSE)$expected
#one expected cell is less than 5
fisher.test(chi) #for small expected freq
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

################### CORRELATION ##############################################
#pdf(file = "Figs/correlation.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

#CORRELATION
  ggplot(SA.correlation.wide)+
    aes((algae),(sponge))+ #(x,y) color=CoralArea_category
    geom_point(aes(color=CoralArea_category, shape=coral.species))+
    labs(x=(expression(paste("Sponge Overgrowth"," ", (mm^2)))),y=(expression(paste("Algal Overgrowth"," ", (mm^2)))))+
  theme_classic()+
  geom_smooth(method="lm", color="black")
  
dev.off()

#percent of an ellipse coral
#low number of points because last trip, lots were also SSID, and some were off transect
#This trip had few, such a low percentage

model <- lm((sponge)~(algae),data=SA.correlation.wide)
summary(model, type=II)

#check assumptions

model.diag.metrics <- augment(model)
sresid <- studres(model) #studentized residuals
shapiro.test(sresid) # non normal

#Residuals vs Fitted
#linearity assumption
#show no fitted pattern, horizontal at zero
plot(model, 1)

#Normal Q-Q. Used to examine whether the residuals are normally distributed. 
#It’s good if residuals points follow the straight dashed line.
plot(model,2)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.diag.metrics$.std.resid)

#Homogeneity of variance of the residuals (homoscedasticity)
#scale-location plot, also known as the spread-location plot.
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points. 
plot(model, 3)

#Residuals vs Leverage. Used to identify influential cases
#Observations whose standardized residuals are greater than 3 in absolute value are possible outliers (James et al. 2014)
#Standardized residuals interpreted as the # of SE away from the regression line.
#A value of this statistic above 2(p + 1)/n indicates an observation with high leverage (P. Bruce and Bruce 2017); 
#where, p is the number of predictors and n is the number of observations.
#outlying values are generally located at the upper right corner or at the lower right corner. 

# Cook's distance
plot(model, 4, id.n = 3) #id.n to label top 5
# Residuals vs LeverageL influential obs
plot(model, 5, id.n = 5)

#high influence if Cook’s distance exceeds 4/(n - p - 1)(P. Bruce and Bruce 2017)
#where n is the number of observations and p the number of predictor variables.
#4/9-1-1
#4/7=0.5714

#y=-1.527215, x= 2.2951958, ID= 43

SA.correlation.wide.1 <- SA.correlation.wide%>%
  filter(!coral.ID %in% c(43))

model1 <- lm(log(sponge)~log(algae),data=SA.correlation.wide.1)

AIC(model, model1)
Anova(model, model1, type=3)

summary(model, type=II)
summary(model1, type=II)
#(y=-0.44x-0.04, R2=0.27, p<0.001)
#outlier changed overall result, but not significantly?


################### AREA PERCENT ##############################################
#Hypothesis 1: SA

# pdf(file = "Figs/other_SA_perc_log.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches
# 
ggplot(SA.2021)+
  aes(y=(other.SA.perc.ellipsoid), x=(coral.SA.ellipsoid))+ #coral.max.mm, (x,y) #all 6cm lower have more algal overgrowth
  geom_point(aes(color=other.SA.type))+ #color=CoralArea_category,
  labs(x=(expression(paste("Log coral (ellipsoid) area"," ", (mm^2)))),y="Log Algal area (%)")+
  theme_classic(base_size=12)+
  geom_smooth(method="lm", color="black",size=0.5, aes(color="Exp Model"), formula= (y ~ (exp(-0.005*x))))+
  annotate("text", x = 4e+06, y = 25, label = "y=e^-0.005x")

dev.off()

model <- lm(log(other.SA.perc.ellipsoid)~log(coral.SA.ellipsoid), data=SA.2021) #can't log 0s
summary(model, type=II)

y <- SA.2021$other.SA.perc.ellipsoid
x <- SA.2021$coral.SA.ellipsoid

fit3 <- lm(y~(exp(-0.005*x)))
summary(fit3)

#y=-0.74x+7.93

################### AREA OG VALUES ##############################################
#OG VALUES
#pdf(file = "Figs/other_SA_raw.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

ggplot(SA.2021)+
  aes(y=log10(other.SA), x=log10(coral.SA.ellipsoid))+ #coral.max.mm, (x,y) #all 6cm lower have more algal overgrowth
  geom_point(aes(color=other.SA.type))+
  labs(x=(expression(paste("Coral area"," ", (mm^2)))),y=(expression(paste("Algal area"," ", (mm^2)))))+
  theme_classic(base_size=12)+
  geom_smooth(method="lm", color="black")

model <- lm(log10(other.SA)~log10(coral.SA.ellipsoid), data=SA.2021) #can't log 0s
summary(ancova)

#box cox
boxcox(model) #told me to use logs on both axes
boxcox(model, lambda=seq(-0.5,0.5,0.01)) #told me to use logs on both axes

#histogram visualization
hist(SA.2021$other.SA)
y<- log10(SA.2021$other.SA)
hist(y)

hist(SA.2021$coral.SA.ellipsoid)
y<- log10(SA.2021$coral.SA.ellipsoid)
hist(y)

SA.2021.wins <- SA.2021%>%
  mutate(log_other.SA.wins = Winsorize(log10(other.SA), na.rm=TRUE))

#other models
model.wins <- lm(log_other.SA.wins~log10(coral.SA.ellipsoid), data=SA.2021.wins) #can't log 0s

AIC(model.wins, model) #winsorized better fit
anova(model.wins, model) #winsorized models are sig different
summary(model, type=II)
summary(model.wins, type=II) #winsorized model interpretation is the same

#linear regression assumptions
model.diag.metrics <- as.data.frame(augment(model))
sresid <- studres(model) #studentized residuals
shapiro.test(sresid) # non normal


#linearity assumption
plot(model, 1)

#Normal Q-Q 
plot(model,2)
shapiro_test(model.diag.metrics$.std.resid)

#Homogeneity of variance of the residuals (homoscedasticity)
plot(model, 3)

#standardized residuals > 3 in absolute value are possible outliers (James et al. 2014)
#Standardized residuals interpreted as the # of SE away from the regression line.

#high influence if Cook’s distance exceeds 4/(n - p - 1)(P. Bruce and Bruce 2017)
#where n is the number of observations and p the number of predictor variables.
#4/111-1-1
#4/109=0.037

# Cook's distance
plot(model, 4, id.n = 3) #id.n to label top 5
# Residuals vs LeverageL influential obs
plot(model, 5, id.n = 5)

#other.SA, coral.SA.ellipsoid
#ylog=1.591065,  xlog=6.112369, ID= 5468
#y=39.00004, x=1295296

#ylog=4.015360, xlog=6.314012
#y=10360.01, x=2060687, ID= 82

SA.2021.outliers <- SA.2021%>%
  filter(!coral.ID %in% c(5468, 82))

model1 <- lm(log10(other.SA)~log10(coral.SA.ellipsoid), data=SA.2021.outliers) #can't log 0s

AIC(model, model1, model.wins)
Anova(model, model1, type=3) #taking outlier has significant effect

summary(model, type=II)
summary(model1, type=II)
#no difference in interpretation
#Keep whole data set

#ANCOVA for algae and sponges
ancova <- lm(log10(other.SA)~other.SA.type*log10(coral.SA.ellipsoid), data=SA.2021) #can't log 0s
summary(ancova)
#Coral size has a significant effect
#Slopes of algae and sponge are the same because no interaction
anova(ancova)

ancova1 <- lm(log10(other.SA)~other.SA.type+log10(coral.SA.ellipsoid), data=SA.2021) #can't log 0s
anova(ancova, ancova1)
#model simplification was justified because it caused a negligible reduction in the explanatory power of the model
#taking out sponge as a factor
ancova2 <- lm(log10(other.SA)~log10(coral.SA.ellipsoid), data=SA.2021)
anova(ancova1, ancova2) #there is a significant effect of sponge with the intercept
summary(ancova1)
anova(ancova1) #overall no effect?
step(ancova)

#data sheet with no sponges because only 6 sponges and because 
#sponges have same intercept and slope as algae
SA.2021.algae <- SA.2021 %>%
  filter(other.SA.type=="algae")

# pdf(file = "Figs/other_SA_raw_algae.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches
  
ggplot(SA.2021.algae)+
  aes(y=log10(other.SA), x=log10(coral.SA.ellipsoid))+ #coral.max.mm, (x,y) #all 6cm lower have more algal overgrowth
  geom_point()+
  labs(x=(expression(paste("Log coral area"," ", (mm^2)))),y=(expression(paste("Log algal area"," ", (mm^2)))))+
  theme_classic(base_size=12)+
  geom_smooth(method="lm", color="red",size=0.5, formula= (y~x))+
  annotate("text", x = 5, y = 4.5, label = "y=0.26x+1.44, R2=0.06, p<0.01")

#geom_smooth(method="lm", color="red",size=0.5, formula= (log10(y) ~ log10(x)))

dev.off()

model <- lm(log10(other.SA)~log10(coral.SA.ellipsoid), data=SA.2021.algae) #can't log 0s

#box cox
boxcox(model) #told me to use logs on both axes
boxcox(model, lambda=seq(-0.5,0.5,0.01)) #told me to use logs on both axes

#histogram visualization
hist(SA.2021.algae$other.SA)
y<- log10(SA.2021.algae$other.SA)
hist(y)

hist(SA.2021.algae$coral.SA.ellipsoid)
y<- log10(SA.2021.algae$coral.SA.ellipsoid)
hist(y)

SA.2021.algae.wins <- SA.2021.algae%>%
  mutate(log_other.SA.wins = Winsorize(log10(other.SA), na.rm=TRUE))

#other models
model.wins <- lm(log_other.SA.wins~log10(coral.SA.ellipsoid), data=SA.2021.algae.wins) #can't log 0s

AIC(model.wins, model) #winsorized better fit
anova(model.wins, model) #winsorized models are sig different
summary(model, type=II)
summary(model.wins, type=II) #winsorized model interpretation is the same

#linear regression assumptions
model.diag.metrics <- as.data.frame(augment(model))
sresid <- studres(model) #studentized residuals
shapiro.test(sresid) # non normal

#linearity assumption
plot(model, 1)

#Normal Q-Q 
plot(model,2)
shapiro_test(model.diag.metrics$.std.resid)

#Homogeneity of variance of the residuals (homoscedasticity)
plot(model, 3)

#standardized residuals > 3 in absolute value are possible outliers (James et al. 2014)
#Standardized residuals interpreted as the # of SE away from the regression line.

#high influence if Cook’s distance exceeds 4/(n - p - 1)(P. Bruce and Bruce 2017)
#where n is the number of observations and p the number of predictor variables.
#4/111-1-1
#4/109=0.037

# Cook's distance
plot(model, 4, id.n = 3) #id.n to label top 5
# Residuals vs LeverageL influential obs
plot(model, 5, id.n = 5)

#other.SA, coral.SA.ellipsoid
#ylog=1.591065,  xlog=6.112369, ID= 5468
#y=39.00004, x=1295296

#ylog=4.015360, xlog=6.314012
#y=10360.01, x=2060687, ID= 82

#ylog=1.447158, xlog=5.686616
#y=28, x=485977.3, ID= 86

SA.2021.algae.outliers <- SA.2021.algae%>%
  filter(!coral.ID %in% c(5468, 82, 86))

model1 <- lm(log10(other.SA)~log10(coral.SA.ellipsoid), data=SA.2021.algae.outliers) #can't log 0s

AIC(model, model1, model.wins)
Anova(model, model1, type=3) #taking outlier has significant effect

summary(model, type=II) #USE THIS
summary(model1, type=II)
#no difference in interpretation
#Keep whole data set

################### AREA BAR PLOT ##############################################
SA.2021.bar <- SA.2021%>%
  group_by(CoralArea_category, other.SA.type)%>%
  summarise(other.SA.perc.ellipsoid.mean=mean(other.SA.perc.ellipsoid, na.rm=TRUE),
            other.SA.perc.ellipsoid.std=sd(other.SA.perc.ellipsoid, na.rm=TRUE),
            other.SA.perc.ellipsoid.se=(sd(other.SA.perc.ellipsoid)/sqrt(length(other.SA.perc.ellipsoid))))

#BAR PLOT
ggplot(SA.2021.bar)+
  aes(CoralArea_category, other.SA.perc.ellipsoid.mean, fill=other.SA.type)+ #(x,y)
  geom_bar(width = 0.8, position = position_dodge(width=1), stat="identity")+
  geom_errorbar(aes(ymin=other.SA.perc.ellipsoid.mean-other.SA.perc.ellipsoid.se, 
                    ymax=other.SA.perc.ellipsoid.mean+other.SA.perc.ellipsoid.se, 
                    width=0.1),position=position_dodge(1))+
  theme_classic(base_size=12)+
  labs(x="Coral size",y="Surface area overgrown (%)")

coral_ranges <- SA.2021 %>% 
  group_by(CoralArea_category) %>% 
  summarise(max.SA = max(coral.SA.ellipsoid),
            min.SA = min(coral.SA.ellipsoid),
            max.diam = max(coral.max.mm),
            min.diam = min(coral.max.mm))

################### PERIMETER ##############################################
#pdf(file = "Figs/perm_raw.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

ggplot(perm)+
  aes(y=(other.perm.mm), x=(coral.perm.mm))+ #(x,y)
  geom_point(aes(color=other.perm.type))+
  theme_classic()+
  geom_smooth(method="lm", color="black")+
  labs(x=(expression(paste("Coral Perimeter"," ", (mm)))),y=(expression(paste("Algae Overgrowth on Perimeter"," ", (mm)))))

dev.off()

#pdf(file = "Figs/perm_perc.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

ggplot(perm)+
  aes(y=(other.perm.perc), x=(coral.perm.mm))+ #(x,y)
  geom_point(aes(color=other.perm.type))+
  theme_classic()+
  geom_smooth(method="lm", color="black")+
  labs(x=(expression(paste("Coral Perimeter"," ", (mm)))),y="Algae Overgrowth on Perimeter (%)")

dev.off()

model <- lm((other.perm.perc)~(coral.perm.mm), data=perm) #can't log 0s
summary(model, type=II)

#y=-0.08+63.08

model.interaction <- lm((other.perm.mm)~(coral.perm.mm)+coral.SA.ellipsoid, data=perm) #can't log 0s
model <- lm((other.perm.mm)~(coral.perm.mm), data=perm) #can't log 0s
summary(model, type=II)
summary(model.interaction, type=II)

#(y=-0.13x-2, R2: 0.007, p=0.59)

#check assumptions

model.diag.metrics <- augment(model)
sresid <- studres(model) #studentized residuals
shapiro.test(sresid) # non normal

#Residuals vs Fitted
#linearity assumption
#show no fitted pattern, horizontal at zero
plot(model, 1)

#Normal Q-Q. Used to examine whether the residuals are normally distributed. 
#It’s good if residuals points follow the straight dashed line.
plot(model,2)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.diag.metrics$.std.resid)

#Homogeneity of variance of the residuals (homoscedasticity)
#scale-location plot, also known as the spread-location plot.
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points. 
plot(model, 3)

#Residuals vs Leverage. Used to identify influential cases
#Observations whose standardized residuals are greater than 3 in absolute value are possible outliers (James et al. 2014)
#Standardized residuals interpreted as the # of SE away from the regression line.
#A value of this statistic above 2(p + 1)/n indicates an observation with high leverage (P. Bruce and Bruce 2017); 
#where, p is the number of predictors and n is the number of observations.
#outlying values are generally located at the upper right corner or at the lower right corner. 

# Cook's distance
plot(model, 4, id.n = 3) #id.n to label top 5
# Residuals vs LeverageL influential obs
plot(model, 5, id.n = 5)

#high influence if Cook’s distance exceeds 4/(n - p - 1)(P. Bruce and Bruce 2017)
#where n is the number of observations and p the number of predictor variables.
#4/49-1-1
#4/47=0.085

#y=453, x= 640, ID= 5470

perm.1 <- perm%>%
  filter(!coral.ID %in% c(5470))

model1 <- lm((other.perm.mm)~(coral.perm.mm), data=perm.1) #can't log 0s

AIC(model, model1)
Anova(model, model1, type=3)

summary(model, type=II)
summary(model1, type=II)
#(y=0.35x+4.21, R2=0.32, p<0.001)
#outlier changed overall result, but mostly just the intercept, but significantly different

################### P:SA ##############################################
#Hypothesis 1: SA
# Step 1: Call the pdf command to start the plot
#pdf(file = "Figs/P_SA.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code

perm.edit <- perm%>%
mutate(P_SA_cat=cut(P_SA, breaks=c(0, 0.0025,0.005,0.02), labels=c("small","medium","large")))%>%
  group_by(P_SA_cat)%>%
  summarise(other.perm.perc.mean=mean(other.perm.perc, na.rm=TRUE),
            other.perm.perc.std=sd(other.perm.perc, na.rm=TRUE),
            other.perm.perc.se=(sd(other.perm.perc)/sqrt(length(other.perm.perc))))


#Hypothesis 2: perm
ggplot(perm.edit)+
  aes(P_SA_cat, other.perm.perc.mean)+ #(x,y)
  geom_bar(width = 0.7, position = position_dodge(width=0.7), stat="identity")+
  geom_errorbar(aes(ymin=other.perm.perc.mean-other.perm.perc.se, ymax=other.perm.perc.mean+other.perm.perc.se, width=0.1),position=position_dodge(1))+
  theme_classic(base_size=12)+
  labs(x="Coral perimeter:surface area ratio",y="Coral perimeter overgrowth \n by algae")
  

# Step 1: Call the pdf command to start the plot
#pdf(file = "Figs/P_SA_perc.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code

ggplot(perm)+
  aes((P_SA), (other.perm.perc))+ #(x,y)
  geom_point(aes(color=CoralArea_category, shape=other.perm.type))+
  labs(x="Coral perimeter:surface area ratio",y="Coral perimeter overgrowth \n by algae (%)")+
  scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"),
                     name="Coral Size", 
                     labels = c("small", "medium", "large"))+
  theme_classic()+
  geom_smooth(method="lm", color="black")
  
dev.off()


#pdf(file = "Figs/P_SA_raw.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 4) # The height of the plot in inches

ggplot(perm)+
  aes((P_SA), (other.perm.mm))+ #(x,y)
  geom_point(aes(color=CoralArea_category, shape=other.perm.type))+
  labs(x="Coral perimeter:surface area ratio",y="Coral perimeter overgrowth \n by algae (mm)")+
  scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"),
                     name="Coral Size", 
                     labels = c("small", "medium", "large"))+
  theme_classic()+
  geom_smooth(method="lm", color="black")

dev.off()



#took out 3 outliers

model.OG <- lm((other.perm.mm)~(P_SA), data=perm)
bc <- boxcox((other.perm.mm)~(P_SA), data=perm)
(lambda <- bc$x[which.max(bc$y)])

model <- lm((other.perm.mm^0.3)~(P_SA), data=perm)


#check assumptions

model.diag.metrics <- augment(model)
sresid <- studres(model) #studentized residuals
shapiro.test(sresid) # non normal

#Residuals vs Fitted
#linearity assumption
#show no fitted pattern, horizontal at zero
plot(model, 1)

#Normal Q-Q. Used to examine whether the residuals are normally distributed. 
#It’s good if residuals points follow the straight dashed line.
plot(model,2)

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.diag.metrics$.std.resid)

#Homogeneity of variance of the residuals (homoscedasticity)
#scale-location plot, also known as the spread-location plot.
#This plot shows if residuals are spread equally along the ranges of predictors. 
#It’s good if you see a horizontal line with equally spread points. 
plot(model, 3)

#Residuals vs Leverage. Used to identify influential cases
#Observations whose standardized residuals are greater than 3 in absolute value are possible outliers (James et al. 2014)
#Standardized residuals interpreted as the # of SE away from the regression line.
#A value of this statistic above 2(p + 1)/n indicates an observation with high leverage (P. Bruce and Bruce 2017); 
#where, p is the number of predictors and n is the number of observations.
#outlying values are generally located at the upper right corner or at the lower right corner. 

# Cook's distance
plot(model, 4, id.n = 3) #id.n to label top 5
# Residuals vs LeverageL influential obs
plot(model, 5, id.n = 5)

#high influence if Cook’s distance exceeds 4/(n - p - 1)(P. Bruce and Bruce 2017)
#where n is the number of observations and p the number of predictor variables.
#4/41-1-1
#4/39=0.10256

#(x=log(P_SA), y=log(other.perm.mm))

#y=6.115892, x= -8.925084, ID= 5470
#y=1.609438, x= -6.529310, ID= 5
#y=4.941642, x= -4.168126, ID= 36

#off of transformed data
#P_SA=0.0154812456, ID= 36
#P_SA=0.0001330103, ID= 5470

P_SA.1 <- perm%>%
  filter(!coral.ID %in% c(36))

P_SA.2 <- perm%>%
  filter(!coral.ID %in% c(36, 5470))


model1 <- lm((other.perm.mm^0.3)~(P_SA), data=P_SA.1) #can't log 0s
model2 <- lm((other.perm.mm^0.3)~(P_SA), data=P_SA.2) #can't log 0s

AIC(model, model1, model2)
Anova(model, model1, type=3)
Anova(model1, model2, type=3)
Anova(model2, model3, type=3)

summary(model, type=II)
summary(model2, type=II)
#(y=12.88x+3.73, R2=0.002, p=0.79)
#outlier did change overall result and significance

#### LOOKING AT MEASUREMENT ERROR OF THE SAME CORAL AT THE SAME TIME ##########
coral.error <- error%>%
  select(coral.ID, type, coral.perm.mm,	coral.max.mm,	coral.min.mm,	coral.height.mm)%>%
  group_by(coral.ID, type)%>%
  summarise(coral.perm.mm=sum(coral.perm.mm, na.rm=TRUE),
            coral.max.mm=sum(coral.max.mm, na.rm=TRUE),
            coral.min.mm = sum(coral.min.mm, na.rm=TRUE),
            coral.height.mm = sum(coral.height.mm, na.rm=TRUE))%>%
  mutate(coral.SA.ellipsoid = (2*3.14*coral.height.mm*(coral.max.mm/2)*(coral.min.mm/2)),
         coral.SA.ellipse = (3.14*coral.max.mm*coral.min.mm),
         P_SA = (coral.perm.mm/coral.SA.ellipsoid))

coral.error$coral.ID <- as.factor(coral.error$coral.ID)

coral.error.graph <- coral.error%>%
  select(coral.ID, coral.SA.ellipsoid, coral.perm.mm, P_SA)%>%
  group_by(coral.ID)%>%
  summarise(coral.SA.ellipsoid.mean=mean(coral.SA.ellipsoid, na.rm=TRUE),
            coral.SA.ellipsoid.std=sd(coral.SA.ellipsoid, na.rm=TRUE),
            coral.SA.ellipsoid.se=(sd(coral.SA.ellipsoid)/sqrt(length(coral.SA.ellipsoid))),
            coral.perm.mm.mean=mean(coral.perm.mm, na.rm=TRUE),
            coral.perm.mm.std=sd(coral.perm.mm, na.rm=TRUE),
            coral.perm.mm.se=(sd(coral.perm.mm)/sqrt(length(coral.perm.mm))),
            P_SA.mean=mean(P_SA, na.rm=TRUE),
            P_SA.std=sd(P_SA, na.rm=TRUE),
            P_SA.se=(sd(P_SA)/sqrt(length(P_SA))))

#The error would be of the coral area
ggplot(coral.error.graph)+
  aes(y=coral.SA.ellipsoid.mean, x=coral.ID)+
  geom_bar(width = 0.7, position = position_dodge(width=0.7), stat="identity")+
  geom_errorbar(aes(ymin=coral.SA.ellipsoid.mean-coral.SA.ellipsoid.se, ymax=coral.SA.ellipsoid.mean+coral.SA.ellipsoid.se, width=0.1),position=position_dodge(1))+
  theme_classic(base_size=12)+
  labs(x="Coral ID",y="Coral surface area")

model <- aov(coral.SA.ellipsoid~coral.ID, data=coral.error)
summary(model)

#Variation between
(2.350e+11 - 2.478e+10 )/2
#1.0511e+11

#total
1.0511e+11 + 2.478e+10 
#1.2989e+11

#Variation from ID
1.0511e+11/1.2989e+11
#0.81% of variation is explained by between group surface area
#0.81*0.81 = 0.6561 error

ggplot(coral.error.graph)+
  aes(y=P_SA.mean, x=coral.ID)+
  geom_bar(width = 0.7, position = position_dodge(width=0.7), stat="identity")+
  geom_errorbar(aes(ymin=P_SA.mean-P_SA.se, ymax=P_SA.mean+P_SA.se, width=0.1),position=position_dodge(1))+
  theme_classic(base_size=12)+
  labs(x="Coral ID",y="Coral surface area")


model <- aov(P_SA~coral.ID, data=coral.error)
summary(model)

#Variation between
(1.412e-07  - 1.925e-08)/2
#6.0975e-08

#total
6.0975e-08 + 1.925e-08 
#8.0225e-08

#Variation from ID
6.0975e-08/8.0225e-08
#0.76 of variation is explained by between group P_SA


ggplot(coral.error.graph)+
  aes(y=coral.perm.mm.mean, x=coral.ID)+
  geom_bar(width = 0.7, position = position_dodge(width=0.7), stat="identity")+
  geom_errorbar(aes(ymin=coral.perm.mm.mean-coral.perm.mm.se, ymax=coral.perm.mm.mean+coral.perm.mm.se, width=0.1),position=position_dodge(1))+
  theme_classic(base_size=12)+
  labs(x="Coral ID",y="Coral surface area")

model <- aov(coral.perm.mm~coral.ID, data=coral.error)
summary(model)

#Variation between
(39740 - 83  )/2
#19828.5

#total
19828.5 + 83
#19911.5

#Variation from ID
19828.5/19911.5
#0.995 of variation is explained by between group perimeter











#NOTES
#Dictyonella funicularis, 
#erksa, Klaus rus: 1981
#small brown porites is usually porites porites
#Porites asteroides is yellow morph


## ERROR IN OVERGROWTH MEASUREMENTS
#delete these columns
field_edit.error <- error%>%
  select(-coral.perm.mm, -perm.loc, -coral.max.mm, -coral.min.mm, -coral.height.mm)

#just so I can read the numbers but with them duplicated down
data.error <- merge(field_edit.error, coral.error, by=c("coral.ID","type"))

#making it so that each row is a coral, adding more specs, includes sponges
master.error <- data.error%>%
  group_by(coral.ID, type, other.SA.type, other.perm.type, other.SA.current)%>%
  summarise(coral.max.mm=mean(coral.max.mm, na.rm=TRUE),
            coral.min.mm=mean(coral.min.mm, na.rm=TRUE),
            coral.height.mm=mean(coral.height.mm, na.rm=TRUE),
            coral.SA.ellipsoid=mean(coral.SA.ellipsoid, na.rm=TRUE),
            coral.SA.ellipse=mean(coral.SA.ellipse, na.rm=TRUE),
            coral.perm.mm=mean(coral.perm.mm, na.rm=TRUE),
            other.max.mm=sum(other.max.mm, na.rm=TRUE),
            other.min.mm=sum(other.min.mm, na.rm=TRUE),
            other.perm.mm=sum(other.perm.mm, na.rm = TRUE))%>%
  mutate(other.SA = (other.max.mm*other.min.mm),
         other.SA.perc.ellipsoid = ((other.SA/coral.SA.ellipsoid)*100),
         other.SA.perc.ellipse = ((other.SA/coral.SA.ellipse)*100),
         other.perm.perc=((other.perm.mm/coral.perm.mm)*100),
         P_SA=(coral.perm.mm/coral.SA.ellipsoid))%>%
  ungroup(coral.ID, type, other.SA.type, other.perm.type, other.SA.current)

str(master.error)

#splitting up data sheet so I can parse out the sponges
SA.error <- master.error%>%
  select(coral.ID, type, other.SA.type, other.SA.current, coral.max.mm, coral.min.mm, coral.height.mm, 
         coral.SA.ellipsoid, coral.SA.ellipse,
         other.max.mm, other.min.mm, other.SA, other.SA.perc.ellipsoid, other.SA.perc.ellipse)%>%
  filter(other.SA.current!="inferred")%>%
  filter(other.SA.type!="none")%>% #taking out corals without any SA overgrowth
  group_by(coral.ID, type, other.SA.type, coral.max.mm, coral.min.mm, coral.height.mm, 
           coral.SA.ellipsoid, coral.SA.ellipse)%>%
  summarise(other.SA=sum(other.SA, na.rm=TRUE),
            other.SA.mean=mean(other.SA, na.rm=TRUE),
            other.SA.std=sd(other.SA, na.rm=TRUE),
            other.SA.se=(sd(other.SA)/sqrt(length(other.SA))),
            other.SA.perc.ellipsoid=sum(other.SA.perc.ellipsoid, na.rm=TRUE),
            other.SA.perc.ellipse=sum(other.SA.perc.ellipse,na.rm=TRUE),
            other.max.mm=sum(other.max.mm, na.rm=TRUE),
            other.min.mm=sum(other.min.mm, na.rm=TRUE))%>%
  ungroup(coral.ID, type, other.SA.type, coral.max.mm, coral.min.mm, coral.height.mm, 
          coral.SA.ellipsoid, coral.SA.ellipse)

SA.error$coral.ID <- as.factor(SA.error$coral.ID)

perm.error <- master.error%>%
  select(coral.ID, type, other.perm.type, coral.perm.mm, other.perm.mm, other.perm.perc,
         P_SA)%>%
  group_by(coral.ID, type, other.perm.type, coral.perm.mm, 
           P_SA)%>%
  summarise(other.perm.mm=sum(other.perm.mm, na.rm=TRUE),
            other.perm.perc=sum(other.perm.perc, na.rm=TRUE),
            other.perm.mean=mean(other.perm.mm, na.rm=TRUE),
            other.perm.std=sd(other.perm.mm, na.rm=TRUE),
            other.perm.se=(sd(other.perm.mm)/sqrt(length(other.perm.mm))))%>%
  filter(other.perm.type %in% c("algae")) #taking out corals without any SA overgrowth




