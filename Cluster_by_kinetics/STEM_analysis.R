#Analysis of STEM results

###########################################
#Compute significance of overlap between the lists from STEM and 
#classification from the paper
#for control condition

Control_summary <- read.csv("C:/Users/Anya/Dropbox/PhD Oxford 2/Stats/Data/STEM analysis/Can we differentiate between B55 dep and indep within condition/Results/Control_summary.csv", header=TRUE)

#Generate vectors of characters of protein names. 
Model_7 <- as.character(Control_summary[,1])
Model_13 <- as.character(Control_summary[,2])
Model_41 <- as.character(Control_summary[,3])
Model_6 <- as.character(Control_summary[,4])
Model_2 <- as.character(Control_summary[,5])

B55_INDEP <- as.character(Control_summary[,7])
B55_DEP <- as.character(Control_summary[,6])

#find:
# 1. genes in common
# 2. genes in each set
#3. total number of genes
#Then compute the hypergeometric test

p_model7_B55DEP <- phyper(length(intersect(Model_7, B55_DEP)), 257, 7391-257, 203, lower.tail = FALSE)
p_model7_B55INDEP <- phyper(length(intersect(Model_7, B55_INDEP)), 257, 7391-257, 153, lower.tail = FALSE)

p_model13_B55DEP <- phyper(length(intersect(Model_13, B55_DEP)), 103, 7391-103, 203, lower.tail = FALSE)
p_model13_B55INDEP <- phyper(length(intersect(Model_13, B55_INDEP)), 103, 7391-103, 153, lower.tail = FALSE)

p_model41_B55DEP <- phyper(length(intersect(Model_41, B55_DEP)), 69, 7391-69, 203, lower.tail = FALSE)
p_model41_B55INDEP <- phyper(length(intersect(Model_41, B55_INDEP)), 69, 7391-69, 153, lower.tail = FALSE)

p_model6_B55DEP <- phyper(length(intersect(Model_6, B55_DEP)), 89, 7391-89, 203, lower.tail = FALSE)
p_model6_B55INDEP <- phyper(length(intersect(Model_6, B55_INDEP)), 89, 7391-89, 153, lower.tail = FALSE)

p_model2_B55DEP <- phyper(length(intersect(Model_2, B55_DEP)), 76, 7391-76, 203, lower.tail = FALSE)
p_model2_B55INDEP <- phyper(length(intersect(Model_2, B55_INDEP)), 76, 7391-76, 153, lower.tail = FALSE)

#Do the same as above for Greatwall

GWL_summary <- read.csv("C:/Users/Anya/Dropbox/PhD Oxford 2/Stats/Data/STEM analysis/Can we differentiate between B55 dep and indep within condition/Results/GWL_summary.csv")

Model_7GWL <- as.character(GWL_summary[,1])
Model_37GWL <- as.character(GWL_summary[,2])
Model_33GWL <- as.character(GWL_summary[,3])
Model_10GWL <- as.character(GWL_summary[,4])

B55_INDEP <- as.character(GWL_summary[,6])
B55_DEP <- as.character(GWL_summary[,5])


#find:
# 1. genes in common
# 2. genes in each set
#3. total number of genes

p_model7GWL_B55DEP <- phyper(length(intersect(Model_7GWL, B55_DEP)), 247, 7391-247, 203, lower.tail = FALSE)
p_model7GWL_B55INDEP <- phyper(length(intersect(Model_7GWL, B55_INDEP)), 247, 7391-247, 153, lower.tail = FALSE)

p_model37GWL_B55DEP <- phyper(length(intersect(Model_37GWL, B55_DEP)), 81, 7391-81, 203, lower.tail = FALSE)
p_model37GWL_B55INDEP <- phyper(length(intersect(Model_37GWL, B55_INDEP)), 81, 7391-81, 153, lower.tail = FALSE)

p_model33GWL_B55DEP <- phyper(length(intersect(Model_33GWL, B55_DEP)), 72, 7391-69, 203, lower.tail = FALSE)
p_model33GWL_B55INDEP <- phyper(length(intersect(Model_33GWL, B55_INDEP)), 72, 7391-69, 153, lower.tail = FALSE)

p_model10GWL_B55DEP <- phyper(length(intersect(Model_10GWL, B55_DEP)), 60, 7391-89, 203, lower.tail = FALSE)
p_model10GWL_B55INDEP <- phyper(length(intersect(Model_10GWL, B55_INDEP)), 60, 7391-89, 153, lower.tail = FALSE)

#create colour coded image
#control

my.data <- data.matrix(Summary_phyper_control[,2:3])
rownames(my.data) <- c("Model 7", "Model 13", "Model 41", "Model 6", "Model 2")
colnames(my.data) <- c("B55-dependent", "B55-independent")

x = 1:ncol(my.data)
y = 1:nrow(my.data)
centers <- expand.grid(y,x)

#make the plot margins a little bigger
par(mar = c(2,7,4,2))

image(x, y, t(my.data),
      col = c(rgb(0,0.6, 0, 0.3), rgb(0,0.8,0,0.3),rgb(0.4,0.8,0,0.3), rgb(0.8,0.8,0,0.3), rgb(1,0.2,0,0.3)),
      breaks = c(min(my.data), 0.001, 0.01, 0.05, 0.1, 1),
      xaxt = 'n', 
      yaxt = 'n', 
      xlab = '', 
      ylab = '',
      ylim = c(max(y) + 0.5, min(y) - 0.5)
)
text(centers[,2], centers[,1], c(my.data), col= "black")
#add margin text
mtext(paste(attributes(my.data)$dimnames[[2]]), at=1:ncol(my.data), padj = -1)
mtext(attributes(my.data)$dimnames[[1]], at=1:nrow(my.data), side = 2, las = 1, adj = 1.2)

#add black lines
abline(h=y + 0.5)
abline(v=x + 0.5)

#Greatwall

my.data.GWL <- data.matrix(Summary_phyper_GWL[,2:3])
rownames(my.data.GWL) <- c("Model 7", "Model 37", "Model 33", "Model 10")
colnames(my.data.GWL) <- c("B55-dependent", "B55-independent")

x = 1:ncol(my.data.GWL)
y = 1:nrow(my.data.GWL)
centers <- expand.grid(y,x)

#make the plot margins a little bigger
par(mar = c(2,7,4,2))

image(x, y, t(my.data.GWL),
      col = c(rgb(0,0.6, 0, 0.3), rgb(0.2,0.8,0,0.3),rgb(0.6,0.8,0,0.3), rgb(1,0.2,0,0.3)),
      breaks = c(min(my.data.GWL), 0.001, 0.01, 0.05, 0.5),
      xaxt = 'n', 
      yaxt = 'n', 
      xlab = '', 
      ylab = '',
      ylim = c(max(y) + 0.5, min(y) - 0.5)
)
text(centers[,2], centers[,1], c(my.data.GWL), col= "black")
#add margin text
mtext(paste(attributes(my.data.GWL)$dimnames[[2]]), at=1:ncol(my.data.GWL), padj = -1)
mtext(attributes(my.data.GWL)$dimnames[[1]], at=1:nrow(my.data.GWL), side = 2, las = 1, adj = 1.2)

#add black lines
abline(h=y + 0.5)
abline(v=x + 0.5)

#################################################
#GO Plotting
library(stringr)
library(tidyverse)
library(ggrepel)

#Greatwall
GO_GWL <- read.csv("C:/Users/Anya/Dropbox/PhD Oxford 2/Stats/Data/k means/Can we differentiate between B55 dep and indep within condition/Results/GO_GWL.csv")
GO_GWL[,9]<-as.numeric(paste(GO_GWL[,9]))

ggplot(GO_GWL, aes(GO_GWL[,1], -log(GO_GWL$Corrected.p.value, 2)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05, 2), linetype = "dashed", colour = "red")+
  geom_text_repel(data = subset(GO_GWL, -log(GO_GWL$Corrected.p.value, 2) > -log(0.05, 2)), aes(data[,1], -log(data[,9],2), label = str_wrap(Category.Name, 20)), size = 3)+
  ggtitle("siRNA Greatwall")+
  xlab("Model ID")+
  ylab("-log(corrected p-value)")

#B55
GO_B55 <- read.csv("C:/Users/Anya/Dropbox/PhD Oxford 2/Stats/Data/STEM analysis/Can we differentiate between B55 dep and indep within condition/Results/GO_B55.csv")
GO_B55[,9]<-as.numeric(paste(GO_B55[,9]))

ggplot(GO_B55, aes(GO_B55[,1], -log(GO_B55$Corrected.p.value, 2)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05, 2), linetype = "dashed", colour = "red")+
  geom_text_repel(data = subset(GO_B55, -log(GO_B55$Corrected.p.value, 2) > -log(0.05, 2)), aes(data[,1], -log(data[,9],2), label = str_wrap(Category.Name, 20)), size = 3)+
  ggtitle("siRNA B55")+
  xlab("Model ID")+
  ylab("-log(corrected p-value)")

#Control
GO_Control <- read.csv("C:/Users/Anya/Dropbox/PhD Oxford 2/Stats/Data/STEM analysis/Can we differentiate between B55 dep and indep within condition/Results/GO_Control.csv")
GO_Control[,9]<-as.numeric(paste(GO_Control[,9]))

ggplot(GO_Control, aes(GO_Control[,1], -log(GO_Control$Corrected.p.value, 2)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05, 2), linetype = "dashed", colour = "red")+
  geom_text_repel(data = subset(GO_Control, -log(GO_Control$Corrected.p.value, 2) > -log(0.05, 2)), aes(data[,1], -log(data[,9],2), label = str_wrap(Category.Name, 20)), size = 3)+
  ggtitle("siRNA Control")+
  xlab("Model ID")+
  ylab("-log(corrected p-value)")