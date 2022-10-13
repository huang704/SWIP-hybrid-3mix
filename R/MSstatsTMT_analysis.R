# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("MSstatsTMT")

library(MSstatsTMT) #v2.2.0
?MSstatsTMT

library(tidyr)
library(dplyr)

source('modelEvaluation.R')
source("groupComparisonTMT.R")
source("utils_group_comparison.R")

######################################################################################
############################### run MSstatsTMT analysis ###########################
######################################################################################
load("Data/pd_psm.rda") 
load("Data/pd_annotation.rda")

processed.input <- PDtoMSstatsTMTFormat(pd_psm,
                                    pd_annotation, 
                                    which.proteinid = "Master.Protein.Accessions")

save(processed.input, file="TMT-Swip-PD-processed-input.rda")

quant.pd <- proteinSummarization(data = processed.input,
                                 method = "msstats",
                                 global_norm = TRUE, # perform global normalization using 'norm'
                                 reference_norm = TRUE,
                                 MBimpute = FALSE # not impute missing values
                                 )

save(quant.pd, file = "TMT-Swip-PD-quantResult.rda")

######################################################################################
############################### run statistical analysis ###########################
######################################################################################
# run statistical testing for swip dataset
load("TMT-Swip-PD-quantResult.rda")
data <- quant.pd$ProteinLevelData
data <- data %>% filter(!Condition %in% c("Norm", "Empty"))

#make sure the biological replicate is the combination of mixture and group
data$BioReplicate <- paste(data$Mixture,
                           matrix(unlist(strsplit(as.character(data$Condition) , "\\.")), byrow = TRUE, ncol = 2)[,1],
                           sep="_")
data$Condition <- as.factor(as.character(data$Condition))
data$Condition <- gsub(" ", "", data$Condition)

# perform contrast comparison
levels(as.factor(data$Condition))
# 'Norm' will be removed during tesing and should be not considered in the contrast
comparison1 <-matrix(c(-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,
                       1/7,1/7,1/7,1/7,1/7,1/7,1/7),nrow=1)

comparison2 <-matrix(c(0.5,0,0,0,0,0,-0.5,
                       0.5,0,0,0,0,0,-0.5),nrow=1)

comparison3 <-matrix(c(0.5,0,0,0,0,-0.5,0,
                       0.5,0,0,0,0,-0.5,0),nrow=1)

comparison4 <-matrix(c(0.5,0,0,0,-0.5,0,0,
                       0.5,0,0,0,-0.5,0,0),nrow=1)

comparison5 <-matrix(c(0.5,0,0,-0.5,0,0,0,
                       0.5,0,0,-0.5,0,0,0),nrow=1)

comparison6 <-matrix(c(0.5,0,-0.5,0,0,0,0,
                       0.5,0,-0.5,0,0,0,0),nrow=1)

comparison7 <-matrix(c(0.5,-0.5,0,0,0,0,0,
                       0.5,-0.5,0,0,0,0,0),nrow=1)

comparison8 <-matrix(c(0,0,0,0,0,-0.5,0.5,
                       0,0,0,0,0,-0.5,0.5),nrow=1)

comparison9 <-matrix(c(0,0,0,0,-0.5,0,0.5,
                       0,0,0,0,-0.5,0,0.5),nrow=1)

comparison10 <-matrix(c(0,0,0,-0.5,0,0,0.5,
                       0,0,0,-0.5,0,0,0.5),nrow=1)

comparison11 <-matrix(c(0,0,-0.5,0,0,0,0.5,
                       0,0,-0.5,0,0,0,0.5),nrow=1)

comparison12 <-matrix(c(0,-0.5,0,0,0,0,0.5,
                       0,-0.5,0,0,0,0,0.5),nrow=1)

comparison13 <-matrix(c(0,0,0,0,-0.5,0.5,0,
                        0,0,0,0,-0.5,0.5,0),nrow=1)

comparison14 <-matrix(c(0,0,0,-0.5,0,0.5,0,
                        0,0,0,-0.5,0,0.5,0),nrow=1)

comparison15 <-matrix(c(0,0,-0.5,0,0,0.5,0,
                        0,0,-0.5,0,0,0.5,0),nrow=1)

comparison16 <-matrix(c(0,-0.5,0,0,0,0.5,0,
                        0,-0.5,0,0,0,0.5,0),nrow=1)

comparison17 <-matrix(c(0,0,0,-0.5,0.5,0,0,
                        0,0,0,-0.5,0.5,0,0),nrow=1)

comparison18 <-matrix(c(0,0,-0.5,0,0.5,0,0,
                        0,0,-0.5,0,0.5,0,0),nrow=1)

comparison19 <-matrix(c(0,-0.5,0,0,0.5,0,0,
                        0,-0.5,0,0,0.5,0,0),nrow=1)

comparison20 <-matrix(c(0,0,-0.5,0.5,0,0,0,
                        0,0,-0.5,0.5,0,0,0),nrow=1)

comparison21 <-matrix(c(0,-0.5,0,0.5,0,0,0,
                        0,-0.5,0,0.5,0,0,0),nrow=1)

comparison22 <-matrix(c(0,-0.5,0.5,0,0,0,0,
                        0,-0.5,0.5,0,0,0,0),nrow=1)

comparison23 <-matrix(c(0,0,0,0,0,0,0,
                       1,0,0,0,0,0,-1),nrow=1)

comparison24 <-matrix(c(0,0,0,0,0,0,0,
                       1,0,0,0,0,-1,0),nrow=1)

comparison25 <-matrix(c(0,0,0,0,0,0,0,
                       1,0,0,0,-1,0,0),nrow=1)

comparison26 <-matrix(c(0,0,0,0,0,0,0,
                       1,0,0,-1,0,0,0),nrow=1)

comparison27 <-matrix(c(0,0,0,0,0,0,0,
                       1,0,-1,0,0,0,0),nrow=1)

comparison28 <-matrix(c(0,0,0,0,0,0,0,
                       1,-1,0,0,0,0,0),nrow=1)

comparison29 <-matrix(c(0,0,0,0,0,0,0,
                       0,0,0,0,0,-1,1),nrow=1)

comparison30 <-matrix(c(0,0,0,0,0,0,0,
                       0,0,0,0,-1,0,1),nrow=1)

comparison31 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,0,-1,0,0,1),nrow=1)

comparison32 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,-1,0,0,0,1),nrow=1)

comparison33 <-matrix(c(0,0,0,0,0,0,0,
                        0,-1,0,0,0,0,1),nrow=1)

comparison34 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,0,0,-1,1,0),nrow=1)

comparison35 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,0,-1,0,1,0),nrow=1)

comparison36 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,-1,0,0,1,0),nrow=1)

comparison37 <-matrix(c(0,0,0,0,0,0,0,
                        0,-1,0,0,0,1,0),nrow=1)

comparison38 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,0,-1,1,0,0),nrow=1)

comparison39 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,-1,0,1,0,0),nrow=1)

comparison40 <-matrix(c(0,0,0,0,0,0,0,
                        0,-1,0,0,1,0,0),nrow=1)

comparison41 <-matrix(c(0,0,0,0,0,0,0,
                        0,0,-1,1,0,0,0),nrow=1)

comparison42 <-matrix(c(0,0,0,0,0,0,0,
                        0,-1,0,1,0,0,0),nrow=1)

comparison43 <-matrix(c(0,0,0,0,0,0,0,
                        0,-1,1,0,0,0,0),nrow=1)

comparison44 <-matrix(c(1,0,0,0,0,0,-1,
                       0,0,0,0,0,0,0),nrow=1)

comparison45 <-matrix(c(1,0,0,0,0,-1,0,
                       0,0,0,0,0,0,0),nrow=1)

comparison46 <-matrix(c(1,0,0,0,-1,0,0,
                       0,0,0,0,0,0,0),nrow=1)

comparison47 <-matrix(c(1,0,0,-1,0,0,0,
                       0,0,0,0,0,0,0),nrow=1)

comparison48 <-matrix(c(1,0,-1,0,0,0,0,
                       0,0,0,0,0,0,0),nrow=1)

comparison49 <-matrix(c(1,-1,0,0,0,0,0,
                       0,0,0,0,0,0,0),nrow=1)

comparison50 <-matrix(c(0,0,0,0,0,-1,1,
                       0,0,0,0,0,0,0),nrow=1)

comparison51 <-matrix(c(0,0,0,0,-1,0,1,
                       0,0,0,0,0,0,0),nrow=1)

comparison52 <-matrix(c(0,0,0,-1,0,0,1,
                        0,0,0,0,0,0,0),nrow=1)

comparison53 <-matrix(c(0,0,-1,0,0,0,1,
                        0,0,0,0,0,0,0),nrow=1)

comparison54 <-matrix(c(0,-1,0,0,0,0,1,
                        0,0,0,0,0,0,0),nrow=1)

comparison55 <-matrix(c(0,0,0,0,-1,1,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison56 <-matrix(c(0,0,0,-1,0,1,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison57 <-matrix(c(0,0,-1,0,0,1,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison58 <-matrix(c(0,-1,0,0,0,1,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison59 <-matrix(c(0,0,0,-1,1,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison60 <-matrix(c(0,0,-1,0,1,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison61 <-matrix(c(0,-1,0,0,1,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison62 <-matrix(c(0,0,-1,1,0,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison63 <-matrix(c(0,-1,0,1,0,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison64 <-matrix(c(0,-1,1,0,0,0,0,
                        0,0,0,0,0,0,0),nrow=1)

comparison65 <-matrix(c(-1,0,0,0,0,0,0,
                       1,0,0,0,0,0,0),nrow=1)

comparison66 <-matrix(c(0,-1,0,0,0,0,0,
                       0,1,0,0,0,0,0),nrow=1)

comparison67 <-matrix(c(0,0,-1,0,0,0,0,
                       0,0,1,0,0,0,0),nrow=1)

comparison68 <-matrix(c(0,0,0,-1,0,0,0,
                       0,0,0,1,0,0,0),nrow=1)

comparison69 <-matrix(c(0,0,0,0,-1,0,0,
                       0,0,0,0,1,0,0),nrow=1)

comparison70 <-matrix(c(0,0,0,0,0,-1,0,
                       0,0,0,0,0,1,0),nrow=1)

comparison71 <-matrix(c(0,0,0,0,0,0,-1,
                       0,0,0,0,0,0,1),nrow=1)

comparison72 <-matrix(c(-0.5,0,0,0,0,0,0.5,
                       0.5,0,0,0,0,0,-0.5),nrow=1)

comparison73 <-matrix(c(-0.5,0,0,0,0,0.5,0,
                       0.5,0,0,0,0,-0.5,0),nrow=1)

comparison74 <-matrix(c(-0.5,0,0,0,0.5,0,0,
                       0.5,0,0,0,-0.5,0,0),nrow=1)

comparison75 <-matrix(c(-0.5,0,0,0.5,0,0,0,
                       0.5,0,0,-0.5,0,0,0),nrow=1)

comparison76 <-matrix(c(-0.5,0,0.5,0,0,0,0,
                       0.5,0,-0.5,0,0,0,0),nrow=1)

comparison77 <-matrix(c(-0.5,0.5,0,0,0,0,0,
                       0.5,-0.5,0,0,0,0,0),nrow=1)

comparison78 <-matrix(c(0,0,0,0,0,0.5,-0.5,
                       0,0,0,0,0,-0.5,0.5),nrow=1)

comparison79 <-matrix(c(0,0,0,0,0.5,0,-0.5,
                       0,0,0,0,-0.5,0,0.5),nrow=1)

comparison80 <-matrix(c(0,0,0,0.5,0,0,-0.5,
                        0,0,0,-0.5,0,0,0.5),nrow=1)

comparison81 <-matrix(c(0,0,0.5,0,0,0,-0.5,
                        0,0,-0.5,0,0,0,0.5),nrow=1)

comparison82 <-matrix(c(0,0.5,0,0,0,0,-0.5,
                        0,-0.5,0,0,0,0,0.5),nrow=1)

comparison83 <-matrix(c(0,0,0,0,0.5,-0.5,0,
                        0,0,0,0,-0.5,0.5,0),nrow=1)

comparison84 <-matrix(c(0,0,0,0.5,0,-0.5,0,
                        0,0,0,-0.5,0,0.5,0),nrow=1)

comparison85 <-matrix(c(0,0,0.5,0,0,-0.5,0,
                        0,0,-0.5,0,0,0.5,0),nrow=1)

comparison86 <-matrix(c(0,0.5,0,0,0,-0.5,0,
                        0,-0.5,0,0,0,0.5,0),nrow=1)

comparison87 <-matrix(c(0,0,0,0.5,-0.5,0,0,
                        0,0,0,-0.5,0.5,0,0),nrow=1)

comparison88 <-matrix(c(0,0,0.5,0,-0.5,0,0,
                        0,0,-0.5,0,0.5,0,0),nrow=1)

comparison89 <-matrix(c(0,0.5,0,0,-0.5,0,0,
                        0,-0.5,0,0,0.5,0,0),nrow=1)

comparison90 <-matrix(c(0,0,0.5,-0.5,0,0,0,
                        0,0,-0.5,0.5,0,0,0),nrow=1)

comparison91 <-matrix(c(0,0.5,0,-0.5,0,0,0,
                        0,-0.5,0,0.5,0,0,0),nrow=1)

comparison92 <-matrix(c(0,0.5,-0.5,0,0,0,0,
                        0,-0.5,0.5,0,0,0,0),nrow=1)

comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6, comparison7,
                    comparison8, comparison9, comparison10, comparison11, comparison12, comparison13, comparison14,
                    comparison15, comparison16, comparison17, comparison18, comparison19, comparison20, comparison21,
                    comparison22, comparison23, comparison24, comparison25, comparison26, comparison27, comparison28,
                    comparison29, comparison30, comparison31, comparison32, comparison33, comparison34, comparison35,
                    comparison36, comparison37, comparison38, comparison39, comparison40, comparison41, comparison42, 
                    comparison43, comparison44, comparison45, comparison46, comparison47, comparison48, comparison49, 
                    comparison50, comparison51, comparison52, comparison53, comparison54, comparison55, comparison56, 
                    comparison57, comparison58, comparison59, comparison60, comparison61, comparison62, comparison63, 
                    comparison64, comparison65, comparison66, comparison67, comparison68, comparison69, comparison70, 
                    comparison71, comparison72, comparison73, comparison74, comparison75, comparison76, comparison77, 
                    comparison78, comparison79, comparison80, comparison81, comparison82, comparison83, comparison84, 
                    comparison85, comparison86, comparison87, comparison88, comparison89, comparison90, comparison91, 
                    comparison92)
# Set the column names
colnames(comparison)<- levels(as.factor(data$Condition))
# Set the names of each row
row.names(comparison)<-c("Mutant-Control", 
                         "F10-F9",
                         "F10-F8",
                         "F10-F7",
                         "F10-F6",
                         "F10-F5",
                         "F10-F4",
                         "F9-F8",
                         "F9-F7",
                         "F9-F6",
                         "F9-F5",
                         "F9-F4",
                         "F8-F7",
                         "F8-F6",
                         "F8-F5",
                         "F8-F4",
                         "F7-F6",
                         "F7-F5",
                         "F7-F4",
                         "F6-F5",
                         "F6-F4",
                         "F5-F4",
                         "Mutant.F10-Mutant.F9",
                         "Mutant.F10-Mutant.F8",
                         "Mutant.F10-Mutant.F7",
                         "Mutant.F10-Mutant.F6",
                         "Mutant.F10-Mutant.F5",
                         "Mutant.F10-Mutant.F4",
                         "Mutant.F9-Mutant.F8",
                         "Mutant.F9-Mutant.F7",
                         "Mutant.F9-Mutant.F6",
                         "Mutant.F9-Mutant.F5",
                         "Mutant.F9-Mutant.F4",
                         "Mutant.F8-Mutant.F7",
                         "Mutant.F8-Mutant.F6",
                         "Mutant.F8-Mutant.F5",
                         "Mutant.F8-Mutant.F4",
                         "Mutant.F7-Mutant.F6",
                         "Mutant.F7-Mutant.F5",
                         "Mutant.F7-Mutant.F4",
                         "Mutant.F6-Mutant.F5",
                         "Mutant.F6-Mutant.F4",
                         "Mutant.F5-Mutant.F4",
                         "Control.F10-Control.F9",
                         "Control.F10-Control.F8",
                         "Control.F10-Control.F7",
                         "Control.F10-Control.F6",
                         "Control.F10-Control.F5",
                         "Control.F10-Control.F4",
                         "Control.F9-Control.F8",
                         "Control.F9-Control.F7",
                         "Control.F9-Control.F6",
                         "Control.F9-Control.F5",
                         "Control.F9-Control.F4",
                         "Control.F8-Control.F7",
                         "Control.F8-Control.F6",
                         "Control.F8-Control.F5",
                         "Control.F8-Control.F4",
                         "Control.F7-Control.F6",
                         "Control.F7-Control.F5",
                         "Control.F7-Control.F4",
                         "Control.F6-Control.F5",
                         "Control.F6-Control.F4",
                         "Control.F5-Control.F4",
                         "Mutant.F10-Control.F10",
                         "Mutant.F4-Control.F4",
                         "Mutant.F5-Control.F5",
                         "Mutant.F6-Control.F6",
                         "Mutant.F7-Control.F7",
                         "Mutant.F8-Control.F8",
                         "Mutant.F9-Control.F9",
                         "Mutant-ControlF10-F9",
                         "Mutant-ControlF10-F8",
                         "Mutant-ControlF10-F7",
                         "Mutant-ControlF10-F6",
                         "Mutant-ControlF10-F5",
                         "Mutant-ControlF10-F4",
                         "Mutant-ControlF9-F8",
                         "Mutant-ControlF9-F7",
                         "Mutant-ControlF9-F6",
                         "Mutant-ControlF9-F5",
                         "Mutant-ControlF9-F4",
                         "Mutant-ControlF8-F7",
                         "Mutant-ControlF8-F6",
                         "Mutant-ControlF8-F5",
                         "Mutant-ControlF8-F4",
                         "Mutant-ControlF7-F6",
                         "Mutant-ControlF7-F5",
                         "Mutant-ControlF7-F4",
                         "Mutant-ControlF6-F5",
                         "Mutant-ControlF6-F4",
                         "Mutant-ControlF5-F4")

comparison

# add simulated subject variance
fx1 <- formula("Abundance ~ 1 + Group + (1|Mixture) + (1|Mixture:Subject)")
fx2 <- formula("Abundance ~ 1 + Group + (1|Subject)")
fx3 <- formula("Abundance ~ 1 + Group + (1|Mixture)")
fx4 <- formula("Abundance ~ 1 + Group")
model_list <- c(fx1, fx2, fx3, fx4)
methods <- c("CMS", "CS", "CM", "C")

for(k in seq_along(model_list)){
  # without moderation
  data.res <- groupComparisonTMT.v2(data = data,
                                    formula = model_list[[k]],
                                    contrast.matrix = comparison,
                                    moderated = FALSE)
  
  filename1 <- paste0('TMT-Swip-PD-by-moderated', methods[k], ".testing.rda")
  save(data.res, file=filename1)
  rm(data.res)
  rm(filename1)
  
  # with moderation
  data.res <- groupComparisonTMT.v2(data = data,
                                    formula = model_list[[k]],
                                    contrast.matrix = comparison,
                                    moderated = TRUE)
  
  filename2 <- paste0('TMT-Swip-PD-by-moderated', methods[k], ".testing-EB.rda")
  save(data.res, file=filename2)
  rm(data.res)
  rm(filename2)
}

data.res <- groupComparison(data = data,
                            method = "limma+timecourse",
                            contrast.matrix = comparison,
                            moderated = TRUE,
                            adj.method = "BH")
filename4 <- paste0('TMT-Swip-PD-by-limmaCS.testing.rda')
save(data.res, file=filename4)
rm(data.res)
rm(filename4)

data.res <- groupComparison(data = data,
                            method = "limma+timecourse+twoway",
                            contrast.matrix = comparison,
                            moderated = TRUE,
                            adj.method = "BH")
filename5 <- paste0('TMT-Swip-PD-by-limmaCMS.testing.rda')
save(data.res, file=filename5)
rm(data.res)
rm(filename5)
