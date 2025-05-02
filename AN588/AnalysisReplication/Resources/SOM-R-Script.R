# SOM File 1
# R script: Data analysis for the manuscript entitled, "The evolution of the adolescent 
#   growth spurt: urinary biomarkers of bone turnover in wild chimpanzees".


#Open packages

library("dplyr")
library("tidyr")
library("car")
library("pbkrtest")
library("lubridate")
library("XML")
library("arm")
library("MuMIn")
library("MASS")
library("date")
library("chron")
library("EMAtools")
library("mgcv")
library("gratia")
library("quantmod")
library("segmented")
library("ggplot2")
library("JWileymisc")
library("multilevelTools")



#Import datasets
collag <- read.csv("Collagen_Ngogo")
osteoc <- read.csv("Osteocalcin_Ngogo")

#Reformat as tables
collag <- tibble::as_tibble(collag)
osteoc <- tibble::as_tibble(osteoc)


#Filter questionable samples

collag<-filter(collag,!individualID=="Morton(x)Morton??(y)")
osteoc<-filter(osteoc,!individualID=="Morton(x)Morton??(y)")

#Format time variable
collag$time<-chron(times = collag$collTime,format="h:m:s")
osteoc$time<-chron(times = osteoc$collTime,format="h:m:s")

#Z-transform continuous variables
collag$timeSt<-(collag$time - mean(collag$time, na.rm=TRUE))/sd(collag$time, na.rm=TRUE)
collag$ageSt<-(collag$sampleAge - mean(collag$sampleAge, na.rm=TRUE))/sd(collag$sampleAge, na.rm=TRUE)

osteoc$timeSt<-(osteoc$time - mean(osteoc$time, na.rm=TRUE))/sd(osteoc$time, na.rm=TRUE)
osteoc$ageSt<-(osteoc$sampleAge - mean(osteoc$sampleAge, na.rm=TRUE))/sd(osteoc$sampleAge, na.rm=TRUE)


#Select variables
collag1<-dplyr::select(collag,individualID,Sex,NTXnM_SG,timeSt,ageSt,sampleAge)
osteoc1<-dplyr::select(osteoc,individualID,Sex,uMidOC_SG,timeSt,ageSt,sampleAge)

collag1<-na.omit(collag1)
osteoc1<-na.omit(osteoc1)








#COLLAGEN

#Calculating individual means

collagAvg<-collag1 %>% group_by(individualID,Sex) %>% summarise_all(funs(mean))

collagAvg$Sex=as.factor(collagAvg$Sex)


powerTransform(collagAvg$NTXnM_SG,family="bcPower")
#Result - 0.2448516



#ALL INDIVIDUALS

gam_collag_all<-gam((NTXnM_SG^0.2448516)/0.2448516 ~ Sex + s(sampleAge,bs="cr"), method="REML", data=collagAvg)
summary(gam_collag_all)

gam_collag_all<-gam((NTXnM_SG^0.2448516)/0.2448516 ~ s(sampleAge,bs="cr"), method="REML", data=collagAvg)
summary(gam_collag_all)


plot(gam_collag_all, residuals=TRUE)

gam.check(gam_collag_all, k.rep=1000)


#FINDING PEAKS IN GAM FIT

fits1 <- gam_collag_all$fitted.values #getting fits; you can take residuals

df1 <- cbind(collagAvg$sampleAge, collagAvg$NTXnM_SG, fits1)
df1 <- as.data.frame(df1)

index <- which.max(df1$fits) #very useful command
row.in.df <- df1[index,]
AgeOfMaxCollagen <- collagAvg$sampleAge[which.max(fits1)]

AgeOfMaxCollagen


##Using derivatives to identify points of decrease

collag_all_d1 <- derivatives(gam_collag_all,  term = "s(sampleAge)", type = "central")

collag_all_d1 <- collag_all_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

collag_all_d1_plot <- draw(collag_all_d1) +
  geom_line(data = collag_all_d1, aes(x = data, y = change), lwd = 2)

with(collag_all_d1, diff(is.na(change) * 1))



#Using segmented regression to identify breakpoints

lm_collag <- lm((NTXnM_SG^0.2448516)/0.2448516 ~ sampleAge,data=collagAvg)
summary(lm_collag)

lm_collag_segmented_all<-segmented(lm_collag, seg.Z=~sampleAge, npsi=3)
summary(lm_collag_segmented_all)

print(lm_collag_segmented_all)





#MALES ONLY

collagAvgMales <- filter(collagAvg, Sex=="Male")

powerTransform(collagAvgMales$NTXnM_SG,family="bcPower")
#Result - 0.1585482

gam_collag_males<-gam((NTXnM_SG^0.1585482)/0.1585482 ~ s(sampleAge,bs="cr"), method="REML", data=collagAvgMales)
summary(gam_collag_males)


plot(gam_collag_males, residuals=TRUE)

gam.check(gam_collag_males, k.rep=1000)


##FINDING PEAKS IN GAM FIT

fits2 <- gam_collag_males$fitted.values #getting fits, you can take residuals,

df2 <- cbind(collagAvgMales$sampleAge, collagAvgMales$NTXnM_SG, fits2)
df2 <- as.data.frame(df2)

index <- which.max(df2$fits) #very usefull command
row.in.df <- df2[index,]
AgeOfMaxCollagen_Males <- collagAvgMales$sampleAge[which.max(fits2)]

AgeOfMaxCollagen_Males



##Using derivatives to identify periods of change

collag_males_d1 <- derivatives(gam_collag_males,  term = "s(sampleAge)", type = "central")

collag_males_d1 <- collag_males_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

draw(collag_males_d1) +
  geom_line(data = collag_males_d1, aes(x = data, y = change), lwd = 2)

with(collag_males_d1, diff(is.na(change) * 1))



#Using segmented regression to identify breakpoints

lm_collag_males<-lm((NTXnM_SG^0.1585482)/0.1585482 ~ sampleAge, data=collagAvgMales)
summary(lm_collag_males)

lm_collag_segmented_males<-segmented(lm_collag_males, seg.Z=~sampleAge, npsi=3)
print(lm_collag_segmented_males)





#FEMALES

collagAvgFemales <- filter(collagAvg, Sex=="Female")

powerTransform(collagAvgFemales$NTXnM_SG,family="bcPower")
#Result - 0.3197467


gam_collag_females<-gam((NTXnM_SG^0.3197467)/0.3197467 ~ s(sampleAge,bs="cr"), method="REML", data=collagAvgFemales)
summary(gam_collag_females)


plot(gam_collag_females, residuals=TRUE)

gam.check(gam_collag_females, k.rep=1000)


#FINDING PEAKS IN GAM FIT

fits3 <- gam_collag_females$fitted.values #getting fits, you can take residuals,

df3 <- cbind(collagAvgFemales$sampleAge, collagAvgFemales$NTXnM_SG, fits3)
df3 <- as.data.frame(df3)

index <- which.max(df3$fits) #very usefull command
row.in.df <- df3[index,]
AgeOfMaxCollagen_Females <- collagAvgFemales$sampleAge[which.max(fits3)]

AgeOfMaxCollagen_Females


#Using derivatives to identify periods of change

collag_females_d1 <- derivatives(gam_collag_females,  term = "s(sampleAge)", type = "central")

collag_females_d1 <- collag_females_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

draw(collag_females_d1) +
  geom_line(data = collag_females_d1, aes(x = data, y = change), lwd = 2)

with(collag_females_d1, diff(is.na(change) * 1))


#Using segmented regression to identify breakpoints

lm_collag_females<-lm((NTXnM_SG^0.3197467)/0.3197467 ~ sampleAge, data=collagAvgFemales)
summary(lm_collag_females)

lm_collag_segmented_females <- segmented(lm_collag_females, seg.Z=~sampleAge, npsi=2)
print(lm_collag_segmented_females)









#OSTEOCALCIN

#Generating individual means

osteocAvg<-osteoc1 %>% group_by(individualID,Sex) %>% summarise_all(funs(mean))

osteocAvg$Sex=as.factor(osteocAvg$Sex)


powerTransform(osteocAvg$uMidOC_SG,family="bcPower")
#Result - 0.3974672


#GAM - All individuals

gam_OC_all_1<-gam((uMidOC_SG^0.3974672)/0.3974672 ~ Sex + s(sampleAge,bs="cr"), method="REML", data=osteocAvg)
summary(gam_OC_all_1)


gam_OC_all_2<-gam((uMidOC_SG^0.3974672)/0.3974672 ~ s(sampleAge,bs="cr"), method="REML", data=osteocAvg)
summary(gam_OC_all_2)


plot(gam_OC_all_2, residuals=TRUE)

gam.check(gam_OC_all_2, k.rep=1000)


#FINDING PEAKS IN GAM FIT

fits4 <- gam_OC_all_2$fitted.values #getting fits, you can take residuals,

df4 <- cbind(osteocAvg$sampleAge, osteocAvg$uMidOC_SG, fits4)
df4 <- as.data.frame(df4)

index <- which.max(df4$fits) #very usefull command
row.in.df <- df4[index,]
AgeOfMaxOsteocalcin <- osteocAvg$sampleAge[which.max(fits4)]

AgeOfMaxOsteocalcin


#Using derivatives to identify periods of change

osteoc_all_d1 <- derivatives(gam_OC_all_2,  term = "s(sampleAge)", type = "central")

osteoc_all_d1 <- osteoc_all_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

draw(osteoc_all_d1) +
  geom_line(data = osteoc_all_d1, aes(x = data, y = change), lwd = 2)

with(osteoc_all_d1, diff(is.na(change) * 1))


#Using segmented regression to identify breakpoints

lm_OC_all<-lm((uMidOC_SG^0.3974672)/0.3974672 ~ sampleAge, data=osteocAvg)
summary(lm_OC_all)

lm_OC_segmented_all <- segmented(lm_OC_all, seg.Z=~sampleAge, npsi=3)
print(lm_OC_segmented_all)



#MALES

osteocAvgMales <- filter(osteocAvg, Sex=="Male")

powerTransform(osteocAvgMales$uMidOC_SG,family="bcPower")
#Result - 0.4060278  


gam_OC_males<-gam((uMidOC_SG^0.4060278)/0.4060278 ~ s(sampleAge,bs="cr"), method="REML", data=osteocAvgMales)
summary(gam_OC_males)


plot(gam_OC_males, residuals=TRUE)

gam.check(gam_OC_males, k.rep=1000)



##FINDING PEAKS IN GAM FIT

fits5 <- gam_OC_males$fitted.values #getting fits, you can take residuals,

df5 <- cbind(osteocAvgMales$sampleAge, osteocAvgMales$uMidOC_SG, fits5)
df5 <- as.data.frame(df5)

index <- which.max(df5$fits) #very usefull command
row.in.df <- df5[index,]
AgeOfMaxOsteocalcin_Males <- osteocAvgMales$sampleAge[which.max(fits5)]

AgeOfMaxOsteocalcin_Males



#Using derivatives to identify periods of change

osteoc_males_d1 <- derivatives(gam_OC_males,  term = "s(sampleAge)", type = "central")

osteoc_males_d1 <- osteoc_males_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

draw(osteoc_males_d1) +
  geom_line(data = osteoc_males_d1, aes(x = data, y = change), lwd = 2)

with(osteoc_males_d1, diff(is.na(change) * 1))


#Using segmented regression to identify breakpoints

lm_OC_males<-lm((uMidOC_SG^0.4060278)/0.4060278 ~ sampleAge, data=osteocAvgMales)
summary(lm_OC_males)

lm_OC_segmented_males <- segmented(lm_OC_males, seg.Z=~sampleAge, npsi=3)
print(lm_OC_segmented_males)




#FEMALES

osteocAvgFemales <- filter(osteocAvg, Sex=="Female")

powerTransform(osteocAvgFemales$uMidOC_SG,family="bcPower")
#Result - 0.316367  


gam_OC_females<-gam((uMidOC_SG^0.316367)/0.316367 ~ s(sampleAge,bs="cr"), method="REML", data=osteocAvgFemales)
summary(gam_OC_females)


plot(gam_OC_females, residuals=TRUE)

gam.check(gam_OC_females, k.rep=1000)


#FINDING PEAKS IN GAM FIT

fits6 <- gam_OC_females$fitted.values #getting fits, you can take residuals,

df6 <- cbind(osteocAvgFemales$sampleAge, osteocAvgFemales$uMidOC_SG, fits6)
df6 <- as.data.frame(df6)

index <- which.max(df6$fits) #very usefull command
row.in.df <- df6[index,]
AgeOfMaxOsteocalcin_Females <- osteocAvgFemales$sampleAge[which.max(fits6)]

AgeOfMaxOsteocalcin_Females



#Using derivatives to identify periods of change

osteoc_females_d1 <- derivatives(gam_OC_females,  term = "s(sampleAge)", type = "central")

osteoc_females_d1 <- osteoc_females_d1 %>%
  mutate(change = if_else(lower < 0 & upper > 0, NA_real_, derivative))

draw(osteoc_females_d1) +
  geom_line(data = osteoc_females_d1, aes(x = data, y = change), lwd = 2)

with(osteoc_females_d1, diff(is.na(change) * 1))


#Using segmented regression to identify breakpoints

lm_OC_females<-lm((uMidOC_SG^0.316367)/0.316367 ~ sampleAge, data=osteocAvgFemales)
summary(lm_OC_females)

lm_OC_segmented_females <- segmented(lm_OC_females, seg.Z=~sampleAge, npsi=3)
print(lm_OC_segmented_females)







#FIGURES

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)

#Color-Blind Palette for Steroids
cbPalette <- c("#009E73", "#CC79A7", "#de2d26")


# retrieve the fitted data for segmented regressions
lm_collag_segmented_all_fitted <- fitted(lm_collag_segmented_all)
collag_segmented_all_model <- data.frame(Distance = collagAvg$sampleAge, Elevation = lm_collag_segmented_all_fitted)
collag_segmented_all_model_2 <- data.frame(Distance = c(2.46, 4.4, 11.6, 18.5, 66), Elevation = c(8.00, 9.14, 10.03, 7.67, 7.42))

lm_OC_segmented_all_fitted <- fitted(lm_OC_segmented_all)
OC_segmented_all_model <- data.frame(Distance = osteocAvg$sampleAge, Elevation = lm_OC_segmented_all_fitted)
OC_segmented_all_model_2 <- data.frame(Distance = c(2.46, 4.7, 12, 13.4, 66), Elevation = c(1.096, 1.55, 1.58, 1.14, 0.74))


lm_collag_segmented_males_fitted <- fitted(lm_collag_segmented_males)
collag_segmented_males_model <- data.frame(Distance = collagAvgMales$sampleAge, Elevation = lm_collag_segmented_males_fitted)
collag_segmented_males_model_2 <- data.frame(Distance = c(2.46, 12.0, 14.1, 40.6, 52.062), Elevation = c(9.83, 11.68, 9.83, 9.12, 9.97))


lm_collag_segmented_females_fitted <- fitted(lm_collag_segmented_females)
collag_segmented_females_model <- data.frame(Distance = collagAvgFemales$sampleAge, Elevation = lm_collag_segmented_females_fitted)
collag_segmented_females_model_2 <- data.frame(Distance = c(2.93, 13.0, 18.7, 66), Elevation = c(9.19, 9.43, 7.10, 6.85))


lm_OC_segmented_males_fitted <- fitted(lm_OC_segmented_males)
OC_segmented_males_model <- data.frame(Distance = osteocAvgMales$sampleAge, Elevation = lm_OC_segmented_males_fitted)
OC_segmented_males_model_2 <- data.frame(Distance = c(2.46, 11.6, 13.2, 40.7, 52.062), Elevation = c(1.17, 1.65, 1.087, 1.10, 0.899))


lm_OC_segmented_females_fitted <- fitted(lm_OC_segmented_females)
OC_segmented_females_model <- data.frame(Distance = osteocAvgFemales$sampleAge, Elevation = lm_OC_segmented_females_fitted)
OC_segmented_females_model_2 <- data.frame(Distance = c(2.93, 6.9, 13.9, 47.1, 66), Elevation = c(1.81, 2.35, 1.68, 0.96, 1.49))






#Figure 1 - All points

#COLLAGEN - GAM

predict_collagAvg <- predict(gam_collag_all, collagAvg, se.fit=TRUE, interval="confidence")

collagAll_plot_data <- cbind(collagAvg, predict_collagAvg)

a1 <- ggplot(collagAll_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.2448516)/0.2448516)) +
  theme_bw()

a2<-a1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                     size=1.3, linetype="dashed", color="cornflowerblue") + 
                     geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
                    size=1.3, linetype="dashed", color="cornflowerblue") + 
                    geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
a2

a3<-a2+labs(title="Collagen", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(a3)



#COLLAGEN - SEGMENTED REGRESSION

b1 <- ggplot(collagAll_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.2448516)/0.2448516)) +
  theme_bw()

b2<-b1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
b2

b3<-b2+labs(title="", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=collag_segmented_all_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(b3)






#OSTEOCALCIN - GAM

predict_osteocAvg <- predict(gam_OC_all_2, osteocAvg, se.fit=TRUE, interval="confidence")

osteocAll_plot_data <- cbind(osteocAvg, predict_osteocAvg)

c1 <- ggplot(osteocAll_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.3974672)/0.3974672)) +
  theme_bw()

c2<-c1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                     size=1.3, linetype="dashed", color="cornflowerblue") + 
                     geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
                     size=1.3, linetype="dashed", color="cornflowerblue") + 
                     geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
c2

c3<-c2+labs(title="Osteocalcin", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(c3)



#OSTEOCALCIN - SEGMENTED REGRESSION

d1 <- ggplot(osteocAll_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.3974672)/0.3974672)) +
  theme_bw()

d2<-d1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
d2

d3<-d2+labs(title="", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=OC_segmented_all_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(d3)



compiled_figure_1 <- ggarrange(a3,c3,b3,d3, nrow=2, ncol = 2)

compiled_figure_1

compiled_figure_1 <- ggarrange(a3,c3,b3,d3, nrow=2, ncol = 2,
                               labels = c('A','B','C','D'), font.label= list(size=15,face="bold"), 
                               hjust=-20, vjust=4.5 )






#FIGURE 2 - MALES ONLY

##COLLAGEN - GAM

predict_collagAvg_males <- predict(gam_collag_males, collagAvgMales, se.fit=TRUE, interval="confidence")

collagMales_plot_data <- cbind(collagAvgMales, predict_collagAvg_males)

e1 <- ggplot(collagMales_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.1585482)/0.1585482)) +
  theme_bw()

e2<-e1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                                                   size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
            size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
e2

e3<-e2+labs(title="Collagen", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(e3)



#COLLAGEN - SEGMENTED REGRESSION

f1 <- ggplot(collagMales_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.1585482)/0.1585482)) +
  theme_bw()

f2<-f1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
f2

f3<-f2+labs(title="", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=collag_segmented_males_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(f3)



#OSTEOCALCIN - GAM

predict_osteocAvg_males <- predict(gam_OC_males, osteocAvgMales, se.fit=TRUE, interval="confidence")

osteocMales_plot_data <- cbind(osteocAvgMales, predict_osteocAvg_males)

g1 <- ggplot(osteocMales_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.4060278)/0.4060278)) +
  theme_bw()

g2<-g1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                                                   size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
            size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
g2

g3<-g2+labs(title="Osteocalcin", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(g3)



#OSTEOCALCIN - SEGMENTED REGRESSION

h1 <- ggplot(osteocMales_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.4060278)/0.4060278)) +
  theme_bw()

h2<-h1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
h2

h3<-h2+labs(title="", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=OC_segmented_males_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(h3)



compiled_figure_2 <- ggarrange(e3,g3,f3,h3, nrow=2, ncol = 2)

compiled_figure_2










#FIGURE 3 - FEMALES ONLY

##COLLAGEN

predict_collagAvg_females <- predict(gam_collag_females, collagAvgFemales, se.fit=TRUE, interval="confidence")

collagFemales_plot_data <- cbind(collagAvgFemales, predict_collagAvg_females)

i1 <- ggplot(collagFemales_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.3197467)/0.3197467)) +
  theme_bw()

i2<-i1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                                                   size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
            size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
i2

i3<-i2+labs(title="Collagen", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(i3)



#COLLAGEN - SEGMENTED REGRESSION

j1 <- ggplot(collagFemales_plot_data, aes(x = sampleAge, y = (NTXnM_SG^0.3197467)/0.3197467)) +
  theme_bw()

j2<-j1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
j2

j3<-j2+labs(title="", y="Collagen (nM corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=collag_segmented_females_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(j3)




##OSTEOCALCIN

#OSTEOCALCIN - GAM

predict_osteocAvg_females <- predict(gam_OC_females, osteocAvgFemales, se.fit=TRUE, interval="confidence")

osteocFemales_plot_data <- cbind(osteocAvgFemales, predict_osteocAvg_females)

k1 <- ggplot(osteocFemales_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.316367)/0.316367)) +
  theme_bw()

k2<-k1 + geom_smooth(aes(x = sampleAge, y = fit),method="gam",se=TRUE,size=1.3, 
                     color="darkblue") + geom_line(aes(x = sampleAge, y = fit+(1.96*se.fit)),
                                                   size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_line(aes(x = sampleAge, y = fit-(1.96*se.fit)),
            size=1.3, linetype="dashed", color="cornflowerblue") + 
  geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
k2

k3<-k2+labs(title="Osteocalcin", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans"))
print(k3)



#OSTEOCALCIN - SEGMENTED REGRESSION

l1 <- ggplot(osteocFemales_plot_data, aes(x = sampleAge, y = (uMidOC_SG^0.316367)/0.316367)) +
  theme_bw()

l2<-l1 + geom_vline(xintercept = c(5,10), linetype="dotted") + geom_point()
l2

l3<-l2+labs(title="", y="Osteocalcin (ng/ml corr. SG) [Transformed]", x="Age (Years)") +
  theme(plot.title=element_text(size=25,face="bold",hjust=0.5)) +
  theme(axis.title.x=element_text(size=16,vjust=0,face="bold")) +
  theme(axis.title.y=element_text(size=16,face="bold",vjust=2.5)) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  scale_x_continuous(breaks=c(5,10,20,30,40,50,60), labels=c("5","10","20","30","40","50","60")) +
  theme(axis.text.x=element_text(colour="black",size=16,)) +
  theme(axis.text.y=element_text(colour="black",angle = 90,size=16,hjust=0.5,vjust=2)) +
  theme(text=element_text(family="sans")) +
  geom_line(data=OC_segmented_females_model_2, aes(x = Distance, y = Elevation), color="red", size=1.5, linetype="solid")
print(l3)



compiled_figure_3 <- ggarrange(i3,k3,j3,l3, nrow=2, ncol = 2)

compiled_figure_3
