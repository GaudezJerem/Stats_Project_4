####################################################################################################################
#Author: Jeremie Gaudez
#Date: 06/12/2019
#This code tries to fit a sigmoidal dephosphorylation kinetcics model to phosphorylation data over time 
#obtained by mass spectrometry. In this case, it returns a dataframe with candidate B55 substrates, B55 being a 
#phosphatase that follows sigmoidal kinetics. 
#It then tests the accuracy of the model by checking how many of the B55 substrates were correctly identified 
#and how many identified substrates are false positives. 
##################################################################################################################


library(reshape2) #Library to use melt
library(tidyverse)

#Package for sigmoidal fit
library(sicegar)

#Import the normalised dataset - set the correct working directory!
Dataset_Clean=read.csv("Table-S2-fitted-dataset.csv")

#Create output dataframe
Prots_Model_fit_DF=data.frame(Proteins=character(),
                              ID=character(),
                              Position=integer(), 
                              P_value=numeric())

#Loop that goes through every protein (takes a long time)
for(p in 1:length(Dataset_Clean[,1])){
  #Assigning Output variabbles for certein protein characteristics
  Name=Dataset_Clean[p,2]
  Pos=Dataset_Clean[p,4]
  id=Dataset_Clean[p,1]
  
  #Reshape the dataset to isolate contol data in a column
  Melted_Control=melt(Dataset_Clean[p,], id="Gene_names", 
                      measure=grep("Control_", names(Dataset_Clean)))
  
  #Looping through control variables to creata a new "TIME" column in minutes
  index = 1:length(Melted_Control[,1])
  for(i in index){
    if (grepl("_0", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=0
    }
    if (grepl("2_5", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=2.5
    }
    if (grepl("l_5", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=5
    }
    if (grepl("7_5", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=7.5
    }
    if (grepl("10", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=10
    }
    if (grepl("20", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=20
    }
    if (grepl("30", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=30
    }
    if (grepl("45", Melted_Control$variable[i])){
      Melted_Control$TIME[i]=45
    }
  }
  #Since the data represents an inverted sigmoidal curve and the sicegar
  #Package doesn't like it, inverting the data
  inverseData=mutate(Melted_Control, 1/value, 
                     value=NULL)
  
  #Selecting only the TIME and value columns
  inverseData=select(inverseData, -Gene_names, -variable)
  #Renaming for compatibility with sicegar package
  names(inverseData)=c("time", "intensity")
  inverseData=na.omit(inverseData)
  
  #Modelling with various autodetermined parameters
  #Sometimes errors arise, in that case skip to the next protein
  skip_to_next=FALSE
  FACMod=tryCatch({fitAndCategorize(inverseData)
  }, error=function(e){
    message(e)
    skip_to_next <<- TRUE
  })
  if(skip_to_next){next}
  
  #Sometimes NA is returned for the model fitting, in that case skip
  if(anyNA(FACMod[2])){
    message("NA for prot ", id, " ", Name)
    next
  }
  
  #Model returns mutiple p values, including one for slope
  #If that p value is under 0.01 (more stringent because of 7000 proteins)
  #Then add to the output dataframe
  if(FACMod[[2]][8]<0.01){
    binding_DF=data.frame(Proteins=Name, ID=id, Position=Pos, 
                          P_value = FACMod[[2]][8])
    Prots_Model_fit_DF=rbind(Prots_Model_fit_DF, binding_DF)
  }
}

#Save the dataframe
write.csv(Prots_Model_fit_DF, file="Prots_model_P_val.csv")

#New dataset from table
Dataset_Analysed=read.csv("Prots_model_P_val.csv")

#P-value correction if needed
Data_Adj=p.adjust(Dataset_Analysed$slopeParam_Pr_t, method = "BH")
Dataset_Analysed$ADJ_P_val=Data_Adj
Filtered_Data=filter(Dataset_Analysed, ADJ_P_val<0.01)

#Isolate protein hits from whole dataset
Com_Data_Sub=subset(Dataset_Clean, 
                    RecordsIndex %in% Filtered_Data$ID)

#Removing Cross-Mixing columns (superfluous) and replacing bad values
Com_Data_Sub=select(Com_Data_Sub, -grep("CM", names(Com_Data_Sub)))
Com_Data_Sub[Com_Data_Sub=="#VALUE!"]=NaN

#Making all value columns the same type
Transformed=apply(Com_Data_Sub[5:28], 2, as.double)

#Combining back to a single dataframe with all information
Com_Data_Sub=cbind(Com_Data_Sub[1:4], Transformed)

# Making one column with all information for each protein
Isolated_Prots=pivot_longer(Com_Data_Sub, Control_0:GWL_45)
Isolated_Prots=as.data.frame(Isolated_Prots)

#Adding two new columns, one for time and one for Condition
index = 1:length(Isolated_Prots[,1])
for(i in index){
  if (grepl("_0", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=0
  }
  if (grepl("_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=5
  }
  if (grepl("2_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=2.5
  }
  if (grepl("7_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=7.5
  }
  if (grepl("10", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=10
  }
  if (grepl("20", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=20
  }
  if (grepl("30", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=30
  }
  if (grepl("45", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=45
  }
  if (grepl("Control_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="Control"
  }
  if (grepl("B55_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="B55_KD"
  }
  if (grepl("GWL_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="GWL_KD"
  }
}

#Saving file
write.csv(Isolated_Prots, file="Proteins_Time_Condition")

ggplot(filter(Isolated_Prots, Condition=="Control"), 
       aes(x=TIME, y=value, color=Gene_names))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("Time(minutes)")+
  ggtitle("Dephosphorylation over time")+
  geom_point()+
  geom_line()+
  geom_text(data=subset(Isolated_Prots, TIME==45&Condition=="Control"), 
            aes(label=Gene_names), check_overlap = TRUE, nudge_y = 0.01)

#Evaluating model accuracy
#B_55 substrates
B55_substrates=read.csv("JCB_201606033_TableS3.csv")
Accu=0
for(site in B55_substrates$RecordsIndex){
  if (site %in% Isolated_Prots$RecordsIndex){
    Accu = Accu + 1
  }
}
Accu_p_Confirmed = Accu / length(B55_substrates$RecordsIndex)

#Total accuracy
Accu_p_Model = Accu / length(Com_Data_Sub$RecordsIndex)

#How many false positives are Independent Substrates
Ind_substrates=read.csv("JCB_201606033_TableS4.csv")
Accu_ind = 0
for(site in Ind_substrates$ï..RecordsIndex){
  if (site %in% Isolated_Prots$RecordsIndex){
    Accu_ind = Accu_ind + 1
  }
}

#New plot with all bad hits

#Isolating false positives
B55_substrates=read.csv("JCB_201606033_TableS3.csv")
Com_Data_Sub=subset(Filtered_Data, 
                    !ID %in% B55_substrates$RecordsIndex)

#Create dataest with all false positives values
Com_Data_Sub=subset(Dataset_Clean, 
                    RecordsIndex %in% Com_Data_Sub$ID)
Com_Data_Sub=select(Com_Data_Sub, -grep("CM", names(Com_Data_Sub)))
Com_Data_Sub[Com_Data_Sub=="#VALUE!"]=NaN

#Making all value columns the same type
Transformed=apply(Com_Data_Sub[5:28], 2, as.double)

#Combining back to a single dataframe with all information
Com_Data_Sub=cbind(Com_Data_Sub[1:4], Transformed)

# Making one column with all information for each protein
Isolated_Prots=pivot_longer(Com_Data_Sub, Control_0:GWL_45)
Isolated_Prots=as.data.frame(Isolated_Prots)

#Adding two new columns, one for time and one for Condition
index = 1:length(Isolated_Prots[,1])
for(i in index){
  if (grepl("_0", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=0
  }
  if (grepl("_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=5
  }
  if (grepl("2_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=2.5
  }
  if (grepl("7_5", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=7.5
  }
  if (grepl("10", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=10
  }
  if (grepl("20", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=20
  }
  if (grepl("30", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=30
  }
  if (grepl("45", Isolated_Prots$name[i])){
    Isolated_Prots$TIME[i]=45
  }
  if (grepl("Control_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="Control"
  }
  if (grepl("B55_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="B55_KD"
  }
  if (grepl("GWL_", Isolated_Prots$name[i])){
    Isolated_Prots$Condition[i]="GWL_KD"
  }
}

#Saving file
write.csv(Isolated_Prots, file="Proteins_Time_Condition")

ggplot(filter(Isolated_Prots, Condition=="Control"), 
       aes(x=TIME, y=value, color=Gene_names))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlab("Time(minutes)")+
  ggtitle("Dephosphorylation over time - False positives")+
  geom_point()+
  geom_line()+
  geom_text(data=subset(Isolated_Prots, TIME==45&Condition=="Control"), 
            aes(label=Gene_names), check_overlap = TRUE, nudge_y = 0.01)
