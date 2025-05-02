library("dplyr")
library("pbkrtest")
library("lubridate")
library("XML")
library("arm")
library("MuMIn")
library("MASS")
library("date")
library("chron")
# EMA tools has been deprecated for R 4.4.3;
# lme4 apparently can peform the same mixed
# effects modeling so I'm substituting it in
# library("EMAtools")
library("lme4")
library("mgcv")
library("gratia")
library("quantmod")
library("segmented")
library("ggplot2")
library("JWileymisc")
library("multilevelTools")
library("tidyverse")

# adding ggpubr because idk how they made the plot grids without it
library("ggpubr")

#Import datasets
collag <- read.csv("https://github.com/bygentry/AN588_DataAnalysisRecreation_gentrym1/raw/refs/heads/main/Resources/Ngogo-collagen.csv")
osteoc <- read.csv("https://github.com/bygentry/AN588_DataAnalysisRecreation_gentrym1/raw/refs/heads/main/Resources/Ngogo-osteocalcin.csv")

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


