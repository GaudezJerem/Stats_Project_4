####################################################################################################################
# Author: Jeremie Gaudez
#Date: 06/12/2019
#This code takes a specific processed mass spectrometry dataset with heavy/light phosphorylation ratios
#And tries to fit a sigmoidal model to one of its substrates
####################################################################################################################


library(reshape2) #Library to use melt
library(tidyverse)

# Don't forget to set the working directory to the dataset folder!
Dataset_Clean=read.csv("Table-S2-fitted-dataset.csv")
#Get only PRC1 data for validation
PRC1_sub=subset(Dataset_Clean, Gene_names %in% "PRC1"&Position %in% "481")

#Get PRC1 data which is interesting (by Phosphorylation position
#for control data
PRC1_Melt=melt(PRC1_sub, id="Gene_names", 
measure=grep("Control_", names(PRC1_sub)))

#Create a new column with error data and subsequently removing the error rows
for (var in PRC1_Melt$variable){
  if(grepl("error", var)){
    PRC1_Melt$ERROR=PRC1_Melt[(1+length(PRC1_Melt[,1])/2):length(PRC1_Melt[,1]), 3]
    PRC1_Melt=subset(PRC1_Melt[1:(length(PRC1_Melt[,1])/2),])
    break
  }
}


#Creating a time column
index = 1:length(PRC1_Melt[,1])
for(i in index){
  if (grepl("_0", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=0
  }
  if (grepl("2_5", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=2.5
  }
  if (grepl("l_5", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=5
  }
  if (grepl("7_5", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=7.5
  }
  if (grepl("0", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=0
  }
  if (grepl("10", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=10
  }
  if (grepl("20", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=20
  }
  if (grepl("30", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=30
  }
  if (grepl("45", PRC1_Melt$variable[i])){
    PRC1_Melt$TIME[i]=45
  }
}
PRC1_Melt=select(PRC1_Melt,-Gene_names, -variable)
#Plotting
ggplot(data=PRC1_Melt, aes(x=TIME, y=value))+
  geom_point(na.rm=TRUE)

#Package for sigmoidial modelling
library(sicegar)

#Tranformating data for sigmoidal modellimg
dataInput=data.frame(PRC1_Melt$TIME, 
                     PRC1_Melt$value)
names(dataInput)=c("time", "intensity")

#Model fitting
sig_mod=multipleFitFunction(dataInput, model="sigmoidal")
t(sig_mod)
fitAndCategorize(dataInput)

#With inverting data
inverseData=mutate(PRC1_Melt, 1/value, #(1/value)*(ERROR/value), 
                   value=NULL)#, ERROR=NULL)
names(inverseData)=c("time", "intensity")#, "error")
ggplot(data=inverseData, aes(x=time, y=intensity))+
  geom_point()#+
  #geom_errorbar(aes(ymin=intensity-error, ymax=intensity+error))

#inverseDataMod=select(inverseData, -error)
sig_mod=multipleFitFunction(dataInput=inverseData, model="sigmoidal")
t(sig_mod)

FACMod=fitAndCategorize(inverseData)
FACMod[[2]][8]
################################################################
