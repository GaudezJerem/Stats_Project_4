require(tidyverse)
require(cowplot)
  
# import control only data
RAWlongdata <- as_tibble(read.csv("S2_fitted_dataset_short.csv"))
# import 'high confidence' hits
RAWhighconf <- as_tibble(read.csv("S3_high_confidence_B55_sub.csv"))

# filter out non and low confidence substrates
Substrates <- semi_join(RAWlongdata, RAWhighconf, by="RecordsIndex")
# filter out pellet fraction substrates
SupernatantOnly <- filter(Substrates, str_detect(RecordsIndex, "_S"))

# basic residue count function and add to table
CountBasic <- function(sequence){
    STRseq <- toString(sequence)
    Lysine <- str_count(STRseq, pattern="K")
    Arginine <- str_count(STRseq, pattern="R")
    BasicRes <- Lysine + Arginine
    return(BasicRes)
  }
CBVector <- Vectorize(CountBasic, vectorize.args="sequence")

SupernatantOnly <- add_column(SupernatantOnly,
                              "Basic_Residues"=CBVector(SupernatantOnly$Sequence_window),
                              .after="Sequence_window")

# wide format to long, keeping errors
SupernatantLong <- gather(SupernatantOnly, temporarycondition, temporaryvalues,
                          Control_0:Sim_Control_40, factor_key=TRUE)

# generate condition and timepoint columns  
ConditionGen <- function(condition){
    if (startsWith(toString(condition), "Sim_")){
      actual <- "Sim_Control"}
    else{
      if (startsWith(toString(condition), "error_"))
        actual <- "stderr"
      else{
        actual <- "Mean_Control"}}
    return(actual)
  }
CGVector <- Vectorize(ConditionGen, vectorize.args="condition")
TimeGen <- function(condition){
    if (endsWith(toString(condition), "l_0")){
      time <- 0.0}
    else{
      if (endsWith(toString(condition), "l_2_5"))
        time <- 2.5
      else{
        if (endsWith(toString(condition), "l_5"))
          time <- 5.0
        else{
          if (endsWith(toString(condition), "l_7_5"))
            time <- 7.5
          else{
            if (endsWith(toString(condition), "l_10"))
              time <- 10.0
            else{
              if (endsWith(toString(condition), "l_20"))
                time <- 20.0
              else{
                if (endsWith(toString(condition), "l_30"))
                  time <- 30.0
                else{
                  if (endsWith(toString(condition), "l_40"))
                  time <- 40.0
                  else
                    time <- 45.0}}}}}}}
    return(time)
  }
TGVector <- Vectorize(TimeGen, vectorize.args="condition")

# add columns to table
SupernatantLong <- add_column(SupernatantLong,
                              "Condition"=CGVector(SupernatantLong$temporarycondition),
                              .before="temporaryvalues")
SupernatantLong <- add_column(SupernatantLong,
                              "Time_mins"=TGVector(SupernatantLong$temporarycondition),
                              .before="temporaryvalues")

# errors as own column
SupernatantNOERR <- filter(SupernatantLong, Condition!="stderr") 
  # table of errors
Stderr <- filter(SupernatantLong, Condition=="stderr")
Stderr <- dplyr::select(Stderr, 1, 11, 14) 
Stderr <- rename(Stderr, StdErr=temporaryvalues, temporaryconditionedit=temporarycondition)
  # function creating common variable for table joining
TempEditor <- function(condition){
    if (endsWith(toString(condition), "Control_0")){
      FinalCon <- "Control_0"}
    else{
      if (endsWith(toString(condition), "Control_2_5"))
        FinalCon <- "Control_2_5"
      else{
        if (endsWith(toString(condition), "Control_5"))
          FinalCon <- "Control_5"
        else{
          if (endsWith(toString(condition), "Control_7_5"))
            FinalCon <- "Control_7_5"
          else{
            if (endsWith(toString(condition), "Control_10"))
              FinalCon <- "Control_10"
            else{
              if (endsWith(toString(condition), "Control_20"))
                FinalCon <- "	Control_20"
              else{
                if (endsWith(toString(condition), "Control_30"))
                  FinalCon <- "Control_30"
                else{
                  FinalCon <- "Control_45"}}}}}}}
    return(FinalCon)
  }
  TEVector <- Vectorize(TempEditor, vectorize.args="condition")
  # adding column 
  Stderr <- add_column(Stderr, "temporarycondition2"=TEVector(Stderr$temporaryconditionedit),
           .after="temporaryconditionedit")
Stderr <- dplyr::select(Stderr, -2)
Stderr <- rename(Stderr, RecordsIndex2=RecordsIndex)
# removing temp columns
FinalData <- left_join(SupernatantNOERR, Stderr, by=c("RecordsIndex"="RecordsIndex2",
                                                      "temporarycondition"="temporarycondition2"))
FinalData <- rename(FinalData, Experimental_values=temporaryvalues)
FinalData <- dplyr::select(FinalData, -11)

FinalData <- FinalData[c(1,2,3,4,5,7,8,9,11,6,10,12,13,14)]

# rate of change 
RateofChange <- function(recordindex){
  SinglePosition <- filter(FinalData, FinalData[1]==toString(recordindex) & 
                             FinalData[9]=="Mean_Control")
  LinReg <- lm(Experimental_values ~ Time_mins, data=SinglePosition)
  Slope <- LinReg$coefficients[2]
  return(Slope)
}
RoCVector <- Vectorize(RateofChange, vectorize.args="recordindex")

FinalData <- add_column(FinalData,
                        "rate_of_change"=RoCVector(FinalData$RecordsIndex),
                        .after="k_b55")

# plot of RoC vs kB55 coloured by Basic Residues
RoCkB55Reg <- lm(k_b55 ~ rate_of_change, data=FinalData)
RoCvsk_b55 <- ggplot(FinalData, aes(x=rate_of_change, y=k_b55, colour=Basic_Residues)) +
  geom_point() +
  geom_abline(intercept=RoCkB55Reg$coefficients[1],slope=RoCkB55Reg$coefficients[2]) +
  theme_bw() +
  labs(x="Rate of Change (H/Lmins)", y="k_b55", title="kB55 vs Rate of Change")

RoCkBasicReg <- lm(rate_of_change ~ Basic_Residues, data=FinalData)
RoCvsBasic <- ggplot(FinalData, aes(x=Basic_Residues, y=rate_of_change,
                                    colour=Basic_Residues)) +
  geom_point() +
  geom_abline(intercept=RoCkBasicReg$coefficients[1],slope=RoCkBasicReg$coefficients[2]) +
  theme_bw() +
  labs(x="Basic Residues", y="Rate of Change (H/Lmins)",
       title="Number of Basic Residues vs Rate of Change")

kB55BasicReg <- lm(k_b55 ~ Basic_Residues, data=FinalData)
kB55vsBasic <- ggplot(FinalData, aes(x=Basic_Residues, y=rate_of_change, 
                                     colour=Basic_Residues)) +
  geom_abline(intercept=kB55BasicReg$coefficients[1],slope=RoCkBasicReg$coefficients[2]) +
  geom_point() +
  theme_bw() +
  labs(x="Basic Residues", y="kB55",
       title="Number of Basic Residues vs Rate of Change")

plot_grid(RoCvsBasic, kB55vsBasic)