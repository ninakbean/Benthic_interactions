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
library(nlcor)
library(ppcor) # this pacakge computes partial and semipartial correlations.

x <- seq(0,3*pi,length.out=100)
y <- sin(x)
plot(x,y,type="l")
# linear correlation is small
cor(x,y)
# [1] 6.488616e-17
# nonlinear correlation is more representative
nlcor(x,y, plt = T)
# $cor.estimate
# [1] 0.9774
# $adjusted.p.value
# [1] 1.586302e-09
# $cor.plot

#8/6: all corals accounted for in quads 1-2 of KW. After, did solo & interactions of all specoes. 
#Did not record algae

#8/7: Recorded past and ssid solo and interaction of ALL sizes
#5 overgrowths: three 10x1

#8/8: all ssid and past solo and inter that were not juve. some juve pics taken...
#10 over growths: five 10x1

#8/10: all solo and sponge interactions of adults (use for contingency tabe)
#4 over growths:  three 10x1

#8/11: Just looked at overgrowth scentarios and solos. made measurements (use for contingency table)
#9 overgrowths: three 10x1

#8/12: targeted overgrown corals to fill in correlation. made measurements 
#Did not fill in more of contingency table
#25 overgrowths: six 25X1

setwd("/Users/nina/Projects/St John 2020/Research_Data")
coral <-read_excel("St_john_2020.xlsx",sheet="Data") #to read your file
str(coral)

coral$Transect <- as.factor(coral$Transect)

master <- coral%>%
  filter(Coral!="NA")%>% #take out no corals
  filter(Coral %in% c("ssid","past"))%>%
  filter(picture=="yes")  #Need pics to verify information
                          #All of them are taken out from other filters, but taking them out now just in case

################## Frequency of interactions #######################################
pairs <- master%>%
  filter(Sponge_inter=="sponge_win")%>% #Only want those with sponge overgrowth
  filter(Num_sponge=="1")%>% #Don't know what to do when corals have more than 1 sponge
  filter(Sponge_Identified=="yes")%>%
  group_by(Coral, Sponge_1)%>%
  summarise(n=n())
  #mutate(perc = (n / sum(n))*100)

ssid_pairs <- pairs%>%
  filter(Coral=="ssid")%>%
  unite(pair,Coral,Sponge_1,sep="_")

past_pairs <- pairs%>%
  filter(Coral=="past")%>%
  unite(pair,Coral,Sponge_1,sep="_")


ggplot(data=ssid_pairs, aes(x=pair, y=n)) +
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 8, angle = 50, hjust = 1))

ggplot(data=past_pairs, aes(x=pair, y=n)) +
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 8, angle = 50, hjust = 1))

ggplot(data=pairs, aes(x=Sponge_1, y=n, fill=Coral))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(legend.position="top")+
  xlab("Sponge species")+
  ylab("Frequency of interactions")+
  scale_fill_discrete(name = "Coral", labels = c("P. astreoides", "S. siderea"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width=15))+
 theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

################## chi square test #######################################
chidata <- master%>%
  filter(Life_stage=="adult")%>% #need adults only for chi square
  filter(Date %in% c("8720", "8820", "81020", "81120"))%>% #the days prior could be used too, but right now juves are included
  #filter(Algae_around=="yes")%>% #only look at corals with potential to be overgrown
  filter(Inferred_algae %in% c("none","algae_win"))%>%
  filter(Sponge_inter %in% c("sponge_win", "solo"))

past <- chidata%>%
  filter(Coral == "past")

ssid <- chidata %>%
  filter(Coral=="ssid")

#filter(Inferred_algae %in% c("none","algae_win"))%>% look at inffered algae

#all corals
chi <- chidata%>%
  dplyr::count(Sponge_inter, Inferred_algae, sort = TRUE)
chi

#split into ssid and past
chisep <- chidata%>%
  dplyr::count(Coral, Sponge_inter, Inferred_algae, sort = TRUE)
chisep

#Use Yates' continuity correction for:
#This formula is chiefly used when at least one cell of the table has an expected count smaller than 5. 

table <- table(chidata$Sponge_inter, chidata$Inferred_algae)
table
chisq.test(table, correct=FALSE)

#run tables seperately
past_table <- table(past$Sponge_inter, past$Inferred_algae)
past_table
chisq.test(past_table, correct=FALSE)

ssid_table <- table(ssid$Sponge_inter, ssid$Inferred_algae)
ssid_table
chisq.test(ssid_table, correct=FALSE)

################## H2: Algae and Sponge Overgrowth Correlation #######################################
#good info here!: https://datascience.stackexchange.com/questions/64260/pearson-vs-spearman-vs-kendall
#pearson = for normal data
#kendall = not normal, rank based
#spearman = not normal, rank based

#spearman: There is a error message because R cannot compute exact p values 
#(the test is based on ranks, we have few cars with the same hp or wt).

sizes <- master%>%
  filter(Date %in% c("81120", "81220", "81320", "81420"))%>%
  filter(use_for_size=="yes")%>%
  mutate("Coral_avg" = ((Coral_max1+Coral_max2)/2))%>%
  mutate("Coral_area" = (2*3.14*(Coral_avg/2)^2))%>%
  mutate("Sponge_area" = (Sponge_max1*Sponge_max2))%>%
  mutate("Algae_area" = (Algae_max1*Algae_max2))%>%
  mutate("Perc_sponge" = ((Sponge_area/Coral_area)*100))%>%
  mutate("Perc_algae" = ((Algae_area/Coral_area)*100))%>%
  filter(Sponge_area!="0")%>%
  filter(Algae_area!="0")%>%
  mutate("Total overgrowth" = (Perc_sponge+Perc_algae))%>%
  filter(!(Transect=="rogue" & ID=="h"))%>% #taking out outlier, the point also adds to more than 100%
  mutate(CoralArea_category=cut(Coral_area, breaks=c(0, 10000, 40000, 120000), labels=c("small","medium","large")))
  
hist(sizes$Coral_area, breaks=10)


ggplot(sizes)+
  aes(Perc_sponge, Perc_algae, group=Coral)+ #(x,y)
  geom_point(aes(shape=Coral, color=CoralArea_category), size=2, stroke=1, alpha = 1)+
  labs(x="Sponge overgrowth (%)",y="Algal overgrowth (%)")+
  scale_shape_manual(values=c(2, 1))+
  scale_color_manual(values=c("#F8A42F", "#FF4605", "#860111"))+
  theme_classic()
  geom_smooth(aes(y=Perc_algae), method = "lm", formula = y~(exp(-0.1*x)), se=FALSE, size=1, colour="black")

#size bin
#lots of leverage on the right
#plot CI around the quadratic fit. 

  scale_color_gradient(low="blue", high="light blue")

  scale_color_gradientn(colours = rainbow(5))
 
hist(sizes$Coral_area)

#geom_smooth(aes(y=Perc_algae), method = "lm", formula = y ~ poly(x, 2), se=FALSE, size=1, colour="black")

#test for finding out a,b, and c from a quadratic equation
x <- seq(from=1, to=5, length.out=100)
y <- 4 - 0.6*x + 0.1*x^2 
ll <-(lm(y~x+I(x^2)))  
summary(ll)
## done 

#Finding best fit for data
y <- (sizes$Perc_algae)
x <- (sizes$Perc_sponge)
z <- sizes$Coral_area
x <-((0.06213*x^2)-1.96196*x+19.57629) #linearizing (after taking out outlier)

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
#a= 0.06213
#b= -1.96196
#c= 19.57629

#looking for outlier from quadratic
plot(cooks.distance(fit2))
cooks.distance(fit2)
any(cooks.distance(fit2) > 1) #proposed thresold
any(cooks.distance(fit2) > (4 * nrow(sizes)))  #proposed thresold

shapiro.test(x) 
shapiro.test(y)
lillie.test(x)
lillie.test(y)
plot(x,y)

pcor.test(x, y, z, method = c("pearson"))
pcor.test(x, y, z, method = c("kendall"))
pcor.test(x, y, z, method = c("spearman"))

cor.test(x, y, method=c("pearson"), conf.level=TRUE)
cor.test(x, y, method=c("kendall"))
cor.test(x, y, method=c("spearman")) #for monotonic data
nlcor(x,y, refine=0.5, plt = T) #same result as pearsons..

#Keeping outlier and linearizing
y <- sizes$Perc_algae
x <- sizes$Perc_sponge
x <-((0.06213*x^2)-1.96196*x+19.57629) 
plot(x,y)

fit1 <- lm(y~x+z)
fit1.2 <- lm(y~x)

fit2 <- lm(y~x+I(x^2)+z) 
fit2.2 <- lm(y~x+I(x^2)) #BEST FIT #fit2 <- lm(y~poly(x,2,raw=TRUE)) #same thing this one, raw=true helps

fit3 <- lm(y~(exp(-0.05*x))+z)
fit3.2 <- lm(y~(exp(-0.05*x)))

fit4 <- lm(y~(exp(-0.02*x))+z)
fit4.2 <- lm(y~(exp(-0.02*x)))

AIC(fit1, fit1.2, fit2, fit2.2, fit3, fit3.2, fit4, fit4.2) #quad still better fit...

# Transforming data
PTx <- powerTransform(x, family="bcPower")
summary(PTx)
PTy <- powerTransform(y, family="bcPower")
summary(PTy)
#isixsigma.com/tools-templates/normality/making-data-normal-using-box-cox-power-transformation/

#Data were not normal so I transformed
y <- sizes$Perc_algae
x <- sizes$Perc_sponge
x <-((0.06213*x^2)-1.96196*x+19.57629) 
z <- sizes$Coral_area
plot(x,y)

x <- log(x)
y <- log(y)

shapiro.test(x) 
shapiro.test(y)
lillie.test(x)
lillie.test(y)
qqPlot(x)
qqPlot(y)
plot(x,y)

#look at residuals
#shapiro.test(residuals(mod2))
#lillie.test(residuals(mod2))
#leveneTest(master_biomass$log_resp, master_biomass$Reproduction_mode) #p-value = 4.376e-05, #homogeneity of variances

pcor.test(x, y, z, method = c("pearson"))
pcor.test(x, y, z, method = c("kendall"))
pcor.test(x, y, z, method = c("spearman"))

cor.test(x, y, method=c("pearson"), conf.level=TRUE)
cor.test(x, y, method=c("kendall"))
cor.test(x, y, method=c("spearman")) #for monotonic data
nlcor(x, y, refine=0.5, plt = T) #same result as pearsons


#### PRELIM TO GET SPECIES ##################
#unite("Interaction", Coral, Interaction, Transect, sep="_", na.rm=TRUE)
  
  
  ggplot(graph, aes(x=Interaction, y=n))+
    geom_bar(stat="identity")+
    theme_classic()+
    scale_x_discrete(labels = function(labels) {
      fixedLabels <- c()
      for (l in 1:length(labels)) {
        fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
      }
      return(fixedLabels)})





require(MASS) ## for mvrnorm()
set.seed(1)
dat <- mvrnorm(1000, mu = c(4,5), Sigma = matrix(c(1,0.8,1,0.8), ncol = 2))
dat <- data.frame(dat)
names(dat) <- c("X","Y")
plot(dat)
res <- resid(mod <- lm(Y ~ X, data = dat))
res.qt <- quantile(res, probs = c(0.05,0.95))
want <- which(res >= res.qt[1] & res <= res.qt[2])
plot(dat, type = "n")

points(dat[-want,], col = "black", pch = 21, bg = "black", cex = 0.8)
points(dat[want,], col = "red", pch = 21, bg = "red", cex = 0.8)
abline(mod, col = "blue", lwd = 2)

ares <- abs(res)
absres.qt <- quantile(ares, prob = c(.9))
abswant <- which(ares <= absres.qt)
## plot - virtually the same, but not quite
plot(dat, type = "n")
points(dat[-abswant,], col = "black", pch = 21, bg = "black", cex = 0.8)
points(dat[abswant,], col = "red", pch = 21, bg = "red", cex = 0.8)
abline(mod, col = "blue", lwd = 2)
head(cooks.distance(mod))
any(cooks.distance(mod) > 1)
any(cooks.distance(mod) > (4 * nrow(dat)))
