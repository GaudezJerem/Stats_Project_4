# Setting working directory

path <- rstudioapi::getActiveDocumentContext()$path %>%
  strsplit("/") %>%
  unlist() %>%
  .[-length(.)] %>%
  paste(., collapse = "/")
setwd(path)
getwd()

#####################################################
#                   Packages                        #
#####################################################

LoadLibraries <- function(){
  library(magrittr)
  library(ggplot2)
  library(wesanderson)
  library(dplyr)
  library(car)
  library(dotwhisker)
  library(broom)
  library(limma)
  library(readxl)
  library(edgeR)
  library(MASS)
  library(reshape2)
  library(reshape)
}

LoadLibraries()

#####################################################
#                     Import                        #
#####################################################

file <- "../2_data_input/Other data/Table S1 - supernatant.xlsx"

data <- read_excel(file)

#####################################################
#                     Data-Tidy                     #
#####################################################

# Removing rows that have NA's by subsetting to tidied data

## Tidied data that is used to subset messy data
tidied_data <- read_excel("../2_data_input/Other data/Table S2 - fitted dataset.xlsx")

## Need to remove the rows from pellet in tidied as only using supernatent data

source("../3_functions/var_names.R")

pel_vs_sup <- var_names(colData = tidied_data$RecordsIndex,
                       pattern = "_", 
                       begin = 1, 
                       end = 1)

tidied_data$pel_sup <- pel_vs_sup

tidied_data_subset <- subset(tidied_data, tidied_data$pel_sup == "S")

tidied_data <- tidied_data_subset

## Extracting peptide positions for subsetting messy with tidy

tidied_pos <- tidied_data$Position

data_pos <- rep("NA", length(data$`Positions within proteins`))

position <- data$`Positions within proteins`

##  Selecting the first protein from potential alternatice splices
for (i in 1:length(data_pos)){
  
  split <- strsplit(position[i], ";")
  vector <- split[[1]]
  first <- vector[1]
  data_pos[i] <- first

}

## Subsetting messy data

tidy_IDs <- paste(tidied_data$Protein, 
                  tidied_pos,
                  sep = "_")

row_IDs <- paste(data$Protein,
                 data_pos,
                 sep = "_")

data$ID <- row_IDs
tidied_data$ID <- tidy_IDs

# The two data-sets are now matching and can be compared 
data_subset <- subset(data, data$ID %in% tidied_data$ID)
tidied_subset <- subset(tidied_data, tidied_data$ID %in% data$ID)

ID_1 <- data_subset$ID
ID_2 <- tidied_subset$ID

# One of the rows has a duplicated ID in data and cannot be used
ID_1 <- sort(ID_1)
ID_2 <- sort(ID_2)

# duplicate found
duplicates <- ID_1[duplicated(ID_1)]

data_subset <- data_subset[data_subset$ID != duplicates, ]
tidied_subset <- tidied_subset[tidied_subset$ID != duplicates, ]

# Ordering two data-sets to match
data_sorted <- data_subset[order(data_subset$ID), ]
tidied_sorted <- tidied_subset[order(tidied_subset$ID), ]

#####################################################
#                   Cross Mixing                    #
#####################################################

# Need to transform the data by multipling by CM scaler in tidied_data

# Only want columns of phospho info for b55 KD and control

phos_data_control <- data_sorted[, 19:66]
phos_data_B55 <- data_sorted[, 99:130]

# Need to multiply each value in the phos_data_B55 by the mean CM value in the tidied data
## The dataframes have been organised to have matching rows

CM_B55 <- tidied_sorted$CM_Con_L_B55_H_

B_55_CM <- phos_data_B55

# Sweep function multiplies everything in each row by the corresponding value in vector
B_55_CM <- sweep(B_55_CM, MARGIN = 1, CM_B55, '*')

phos_data_B55 <- B_55_CM

#####################################################
#                  Data-Labelling                   #
#####################################################

# Labelling colnames
## Desired result:
## Control_tp0_rep1, Control_tp2.5_rep1

## Number of time points
tps <- c("tp0", "tp2.5", "tp5", "tp7.5", "tp10", "tp20", "tp30", "tp45")

## Sets of times
control_sets <- length(phos_data_control)/length(tps)
B55_sets <- length(phos_data_B55)/length(tps)

## Control reps

control_reps <- list()
for (i in 1:control_sets){
  
  rep_num <- paste("rep_", i,
                   sep = "")
  
  rep_tps <- rep(rep_num, length(tps))
  
  control_reps[[i]] <- rep_tps
  
  }

## Flatten list
cnt_reps <- unlist(control_reps)

## Reps for control
cnt_tps <- rep(tps, control_sets)

## Factor label for control
cnt_label <- rep("Control", length(colnames(phos_data_control)))

## Making data label
control_colnames <- rep("NA", length(colnames(phos_data_control)))

for (i in 1:length(colnames(phos_data_control))) {
  
  new_label <- paste(cnt_label[i], "_", cnt_tps[i], "_", cnt_reps[i],
                     sep = "")
  
  control_colnames[i] <- new_label
  
  }

# Updating colnames
## Replacing colnames with new labels

colnames(phos_data_control) <- control_colnames

# Repeating above steps for B55

## Control reps

B55_reps <- list()
for (i in 1:B55_sets){
  
  rep_num <- paste("rep_", i,
                   sep = "")
  
  rep_tps <- rep(rep_num, length(tps))
  
  B55_reps[[i]] <- rep_tps
  
}

## Flatten list
b55_reps <- unlist(B55_reps)

## Reps for B55
b55_tps <- rep(tps, B55_sets)

## Factor label for B55
b55_label <- rep("B55", length(colnames(phos_data_B55)))

## Making data label
B55_colnames <- rep("NA", length(colnames(phos_data_B55)))

for (i in 1:length(colnames(phos_data_B55))) {
  
  new_label <- paste(b55_label[i], "_", b55_tps[i], "_", b55_reps[i],
                     sep = "")
  
  B55_colnames[i] <- new_label
  
}

# Updating colnames
## Replacing colnames with new labels

colnames(phos_data_B55) <- B55_colnames

#####################################################
#                    Data-Frame                     #
#####################################################

# Combining two data-sets
phos_data <- cbind(phos_data_control, phos_data_B55)

# Renaming variables for consistency

names(phos_data) <- names(phos_data) %>%
  gsub(pattern = "tp2.5",
       replacement = "tp02.5") %>%
  gsub(pattern = "tp5",
       replacement = "tp05.0") %>%
  gsub(pattern = "tp7.5",
       replacement = "tp07.5") %>%
  gsub(pattern = "tp10",
       replacement = "tp10.0") %>%
  gsub(pattern = "tp20",
       replacement = "tp20.0") %>%
  gsub(pattern = "tp30",
       replacement = "tp30.0") %>%
  gsub(pattern = "tp45",
       replacement = "tp45.0") %>%
  gsub(pattern = "tp0_",
       replacement = "tp00.0_")

# Re-ordering columns for readibility
## 0 - 45 and from Control - B55

## Ordered by numbers, decrease = FALSE
## Change boolean depending on desired time order
phos_data <- colnames(phos_data) %>%
  {phos_data[ , order(.)]}

## Ordered by first characters, decrease = TRUE
## Substr used to select first two charactors of 'Control' and 'B55' for re-ordering
## Change substr start and stop argument values if conditions are positioned differently
## Change boolean depending on desired condition ordering
phos_data <- colnames(phos_data) %>%
  {substr(., 1, 2)} %>%
  {phos_data[ , order(., decreasing = TRUE)]}

#####################################################
#                     DESEQ DATA                    #
#####################################################

# Need to extract variable names for deseq data set

# Function used to extract names (see README for info)
source("../3_functions//var_names.R")

#####################################################
#                     Condition                     #
#####################################################

phos_data_col <- colnames(phos_data)

condition <- var_names(colData = phos_data_col,
                       pattern = "_t", 
                       begin = FALSE, 
                       end = -1)


# Transforming from chr to factor type
condition <- as.factor(condition)
levels(condition)
## Defining Control as base and B55 as treatment
condition <- factor(condition, 
                    levels(condition)[2:1])

#####################################################
#                       Hours                       #
#####################################################

hours <- var_names(colData = phos_data_col, 
                   pattern = "_r", 
                   begin = -4, 
                   end = -1)

##Transforming from chr to factor
str(hours)
hours <- as.factor(hours)
levels(hours)

#####################################################
#                        Reps                       #
#####################################################

# Reps:

rep <- var_names(colData = phos_data_col, 
                 pattern = "rep", 
                 begin = 4, 
                 end = 4)

##Transforming from chr to factor
str(rep)
rep <- as.factor(rep)


#####################################################
#                     colData                       #
#####################################################

# The colData now describes the samples correctly

## Variables should correspond to the correct sample
phos_colData <- data.frame(row.names = phos_data_col, 
                              condition, 
                              hours, 
                              rep)

#####################################################
#                     DGEList                       #
#####################################################


# Matrix is needed for construction of DESeq data class
## Transforming to a matrix
phos_counts <- as.matrix(phos_data)
## Setting NAs as zero
phos_counts[is.na(phos_counts)] <- 0

#Creating DGEList for analysis

#Samples
#Can use previously created phos_colData with the sample information.

#Counts
#Can use previously tidied up counts matrix.

#Genes
#Can use the protein+pos ID in the rows of count data
genetable <- data_subset$ID

# DGElist
phos_dge <- DGEList(counts = phos_counts,
                       samples = phos_colData,
                       genes = genetable)

names(phos_dge)

phos_dge$samples
head(phos_dge$counts)
head(phos_dge$genes)

#####################################################
#                Differential Expression            #
#####################################################

z <- factor(paste(condition, hours, sep = "_"))
design1 <- model.matrix(~0 + z)
colnames(design1) <- gsub("z", "", colnames(design1))

# Contrast Matrix

contr.matrix <- makeContrasts(
  tp00.0 = B55_00.0 - Control_00.0,
  tp02.5 = B55_02.5 - Control_02.5,
  tp05.0 = B55_05.0 - Control_05.0,
  tp07.5 = B55_07.5 - Control_07.5,
  tp10.0 = B55_10.0 - Control_10.0,
  tp20.0 = B55_20.0 - Control_20.0,
  tp30.0 = B55_30.0 - Control_30.0,
  tp45.0 = B55_45.0 - Control_45.0,
  levels = design1
)

# Fitting Linear models

## Trying to fit phos_data that still has NAs
vfit <- lmFit(phos_data, design1)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

## Fitting Ebayes

vfit2 <- eBayes(vfit)

#####################################################
#             Exporting Limma Results               #
#####################################################

# decideTests performs multiple testing with FDR adjustment
## produces a matrix of 0's, 1's, -1s according to no DE, up DE, and down DE, respectively.
dt <-decideTests(vfit2, 
                 method = "global", 
                 adjust.method = "BH",
                 p.value = 0.05)
summary(dt)

# writing a table of treat and decidetests outputs
write.fit(vfit2, dt, file="../4_data_output/decide_tests.csv")
limma_results <- read.csv(file="../4_data_output/decide_tests.csv",
                          sep = "")

up_down <- limma_results[, 28:35]
up_down_sums <- as.vector(rowSums(up_down) > 0)
limma_results$up_down <- ifelse(up_down_sums, "up", "down")

#####################################################
#                 Extracting P-values               #
#####################################################

# Library for combining p-values
library(scran)

# Creating matrix of p.values generated from treat 
treat_pvalues <- vfit2$p.value

# Applying simes method to combine the P values across the contrasts for each gene
treat_pval.list <- lapply(1:ncol(treat_pvalues), FUN=function(i) { treat_pvalues[,i] })
treat_per.gene <- do.call(scran::combinePValues, c(treat_pval.list, list(method="simes")))
# Applying adjustment method 'BH' to these Adj-pvalues for each gene
treat_adj.per.gene <- p.adjust(treat_per.gene, method="BH")
# Extracting adj-pvalues <= 0.01
treat_sig.per.gene <- treat_adj.per.gene <= 0.05

# Adding data to master data

data_sorted$sig <- treat_sig.per.gene
data_sorted$adj_p_value <- treat_adj.per.gene
data_sorted$sig[is.na(data_sorted$sig)] <- FALSE
data_sorted$up_down <- limma_results$up_down

# Extracting proteins that have significant difference 

phos_deg <- subset(data_sorted, data_sorted$sig == TRUE)

# Extracting proteins sig high phos when B55 is knocked down
## This would meet expected profile change if substrate of B55

phos_deg <- subset(phos_deg, phos_deg$up_down == "up")

phos_deg <- phos_deg[order(phos_deg$adj_p_value, decreasing = FALSE),]

write.csv(phos_deg, "../4_data_output/phos_degs.csv")

#####################################################
#            Comparing with model data              #
#####################################################

# Comparing list of DE proteins with high confidence proteins
## Will produce set of even more high confident substrates

file2 <- "../2_data_input/Other data/Table S3 - high confidence B55 substrates.xlsx"
high_conf <- read_excel(file2)

# Subsetting for proteins shared between models
shared_phos_de <- subset(phos_deg, phos_deg$ID %in% high_conf$ID)
## 161 proteins found in both data sets

write.csv(shared_phos_de, "../4_data_output/shared_phos_degs.csv")

#####################################################
#                    GO-analysis                    #
#####################################################

# Extracting the gene names for GO-analysis

de_genes <- shared_phos_de$Protein
background_genes <- data_sorted$Protein

write.csv(de_genes, "../4_data_output/enriched_proteins.csv")

write.csv(background_genes, "../4_data_output/background_proteins.csv")


#####################################################
#                 Example Profiles                  #
#####################################################

# Plotting profiles of examples proteins enriched in GO-terms for mitosis/chromosome sep

gene_count_data <- phos_data
gene_count_data$protein <- data_sorted$Protein
gene_count_data$ID <- data_sorted$ID

# Subsetting count data for only ones sig and in model data = 161 high high conf genes
gene_count_data <- subset(gene_count_data, gene_count_data$ID %in% shared_phos_de$ID)

# Plot functions
source("../3_functions/plot_functions.R")

plot.median <- function(x) {
  m <- mean(x)
  c(y = m, ymin = m, ymax = m)
}

# Subsetting for gene of interest
## Using ID and refering to table to select correct gene of interest
sig_CDCA5 <- subset(gene_count_data, gene_count_data$ID == "Q96FF9_159")
sig_NUMA1 <- subset(gene_count_data, gene_count_data$ID == "Q14980_2000")
sig_SFPQ <- subset(gene_count_data, gene_count_data$ID == "P23246_687")
sig_CDC20 <- subset(gene_count_data, gene_count_data$ID == "Q12834_41")

# Removing un-needed data
sig_CDCA5$protein <- NULL
sig_CDCA5$ID <- NULL

sig_NUMA1$protein <- NULL
sig_NUMA1$ID <- NULL

sig_SFPQ$protein <- NULL
sig_SFPQ$ID <- NULL

sig_CDC20$protein <- NULL
sig_CDC20$ID <- NULL


# Melting data to suitable structure for ggplot
melt_CDCA5 <- melt(sig_CDCA5)
melt_NUMA1 <- melt(sig_NUMA1)
melt_SFPQ <- melt(sig_SFPQ)
melt_CDC20 <- melt(sig_CDC20)

## Extracting descriptions

## Condition:
condition1 <- var_names(colData = melt_CDCA5$variable,
                       pattern = "_t", 
                       begin = FALSE, 
                       end = -1)

condition1 <- as.factor(condition1)

levels(condition1)

condition1 <- factor(condition1, 
                    levels(condition1)[2:1])

# Time:

hours1 <- var_names(colData = melt_CDCA5$variable, 
                   pattern = "_r", 
                   begin = -4, 
                   end = -1)

# Adding variables to data

melt_CDCA5$hours <- hours1
melt_CDCA5$condition <- condition1
melt_CDCA5$gene <- rep("CDCA5 pT159", length(condition1))

melt_NUMA1$hours <- hours1
melt_NUMA1$condition <- condition1
melt_NUMA1$gene <- rep("NUMA1 pT2000", length(condition1))

melt_SFPQ$hours <- hours1
melt_SFPQ$condition <- condition1
melt_SFPQ$gene <- rep("SFPQ pT687", length(condition1))

melt_CDC20$hours <- hours1
melt_CDC20$condition <- condition1
melt_CDC20$gene <- rep("CDC20 pS41", length(condition1))

# Combining data-sets

profiles <- rbind(melt_CDCA5,
                  melt_NUMA1,
                  melt_SFPQ,
                  melt_CDC20)

# Plotting

p <- ggplot(profiles,
            aes(x = hours,
                y = value,
                colour = condition))

p <- p + stat_summary(fun.data="plot.median", 
                                   geom = "line",
                                   alpha = 0.8,
                                   size=0.9,
                                   aes(group = condition))

p <- p + geom_point(size = 2,
                  alpha = 0.8)

p <- p + theme_Publication()

p <- p + scale_color_manual(values=c("#118AB2", "#FFD166"))

p <- p + ylab("(H/L)")

p <- p + xlab('Time (min)')

p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


p <- p + facet_wrap(~ gene, 
                    scales = "free_y")

p <- p + theme(legend.key = element_blank(), 
               strip.background = element_rect(colour="#FFFFFF", fill="#FFFFFF") ) 

p


png("profiles.png", units="cm", width=20, height=15, res=300)
p
dev.off()

# MDS

mds <- plotMDS(phos_dge$counts, plot = FALSE)

## Converting to dataframa
mds_x_data <- as.data.frame(mds$cmdscale.out)
## Adding column names
colnames(mds_x_data) <- c("MD1", "MD2")

## Converting to long format with mds1 and mds2 under same column
## Allows plotting with facets
mds_x_data <- gather(data = mds_x_data,
                     mds,
                     value)

## Sample information added for plot
mds_x_data$sample <- colnames(phos_counts)
mds_x_data$condition <- condition
mds_x_data$hours <- hours

# Plotting MDS over time

mds_plot <- ggplot(data = mds_x_data, 
                   aes(x = hours, 
                       y = value, 
                       group = mds, 
                       colour = condition))

mds_plot <- mds_plot + xlab("Time (min)")

mds_plot <- mds_plot + ylab("Leading LogFC")
mds_plot <- mds_plot + facet_grid(rows = vars(mds), 
                                  scales = "free")

mds_plot <- mds_plot + theme_Publication() #+ scale_colour_Publication()

mds_plot <- mds_plot + scale_color_manual(values=c("#118AB2", "#FFD166"))

mds_plot <- mds_plot + stat_summary(fun.data="plot.median", 
                                    geom = "line",
                                    alpha = 0.8,
                                    width=0.1, size=0.9,
                                    aes(group = condition))

mds_plot <- mds_plot + geom_point(size = 2.5,
                                  alpha = 0.8)

mds_plot <- mds_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

mds_plot <- mds_plot + theme(legend.key = element_blank(), 
                             strip.background = element_rect(colour="#FFFFFF", fill="#FFFFFF") )                   


mds_plot

png("mds.png", units="cm", width=20, height=15, res=300)
mds_plot
dev.off()

# Shares genes plot 

library(UpSetR)

# genes in datas

high_id <- paste(high_conf$Protein, high_conf$Position, sep = "_")
high_conf$ID <- high_id

my_genes <- as.vector(phos_deg$ID)
their_genes <- as.vector(high_conf$ID)

identifier <- c(my_genes, their_genes)
identifier <- unique(identifier)

# ifelse
match <-  identifier %in% as.vector(phos_deg$ID)
match <-  ifelse(match, 1, 0)

match2 <- identifier %in% as.vector(high_conf$ID)
match2 <- ifelse(match2, 1, 0)

# Data for input
upset_input <- as.data.frame(identifier)
upset_input$limma <- match
upset_input$model <- match2

# Upset plot
library(wesanderson)

png("upset.png", units="cm", width=15, height=15, res=300)
upset(upset_input, 
      nsets = 2,
      sets.bar.color = c("#118AB2", "#FFD166"),
      text.scale = 2,
      set_size.numbers_size = 0.5)
dev.off()

