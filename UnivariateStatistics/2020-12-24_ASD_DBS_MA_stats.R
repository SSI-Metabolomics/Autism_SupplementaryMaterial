###############################################################
#                                                             #
#     Script to                                               #
#         1- format your feature table from the GNPS output   #
#            into MetaboAnalyst format EXCLUDING Pools        #
#         2- run MetaboAnalyst missing values imputation,     #
#             normalization/transf./scaling                   #
#         3- run PCA                                          #
#                                                             #
#                                                             #  
###############################################################


#Making sure to get error messages in English
Sys.setenv(LANG = "en")
#Setting the working directory where you have your files and want to save the output
setwd("C:/ASD_DBS") 
library(tidyverse)
library(amap)

##################################
# 1- formatting of table for MA  #
##################################

ft <- read.table ("PEGPPG_cur_68Samples_Pools.csv", sep=";", header=TRUE, dec=".", check.names = FALSE)
##creating the feature ID unique column
ft2 <- ft %>% unite(feature_ID, "row ID", "row m/z", "row retention time", sep = "_", remove = TRUE)  
##making the feature ID the row names and removing the feature ID column
rownames(ft2) <- ft2$feature_ID
ft2 <- select(ft2, -feature_ID)
##transposing the dataframe to have features in columns and samples in row. 
ft2 <- as.data.frame(t(ft2))
class(ft2[1,1]) #should be "numeric"
##creating the "type" column and placing it first
ft2$type <- rownames(ft2)
which(colnames(ft2) == "type")
ft2 <- ft2[ ,c(which(colnames(ft2) == "type"),1:(which(colnames(ft2) == "type")-1))]

ft2$type <- gsub("\\_case.*","case",ft2$type)
ft2$type <- gsub(".*case","case",ft2$type)
ft2$type <- gsub("\\_control.*","control",ft2$type)
ft2$type <- gsub(".*control","control",ft2$type)
ft2$type <- gsub("\\_Pool.*","Pool",ft2$type)
ft2$type <- gsub(".*Pool","Pool",ft2$type)

#making sure zero values are replaced by NA
ft2[ft2 == 0] <- NA

#removing pools
ft_nopools <- ft2[-which(ft2$type == "Pool"),]

##creating the CSV file
write.table(ft_nopools, 'PEGPPG_cur_68Samples_MA.csv', sep=',',col.names = NA, row.names = TRUE, quote = TRUE, na="NA")

######################################################################################################
# 2- run missing values imputation, filtering, transformation, scaling on 72 samples including Pools #
######################################################################################################

# Load the MetaboAnalystR package
library("MetaboAnalystR")
library(pls)

# First step is to create the mSet Object, specifying that the data to be uploaded
# is a peak table ("pktable") and that statistical analysis will be performed ("stat").
mSet<-InitDataObjects("pktable", "stat", FALSE)

# Second step is to read in the filtered peak list, please set the path right first
mSet<-Read.TextData(mSet, "PEGPPG_cur_68Samples_MA.csv", "rowu", "disc")

# The third step is to perform data processing using MetaboAnalystR (filtering/normalization)
# Performing data processing - Data checking
mSet<-SanityCheckData(mSet)

# Performing data processing - Minimum Value Replacing
mSet<-ReplaceMin(mSet)


mSet<-PreparePrenormData(mSet)
#"Normalization" here includes "row-wise" normalization, which is sample-based (sum/mean/median of the samples)
#then transformation (generalized log), then column-wise scaling (feature-based) including Pareto
#Pareto scaling is mean-centered and divided by the square root of the SD.
mSet<-Normalization(mSet, "NULL", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)

#visualizing (saves the figures in your wd)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

#saving your tables in your working directory
mSet<-SaveTransformedData(mSet)

#######################################################################
# 3- run fold-change analysis, t-tests and Wilcoxon rank sum tests    #
#######################################################################

#perform fold-change analysis
mSet <- FC.Anal.unpaired(mSet, 2.0, 0)
mSet <- PlotFC(mSet, "fc_0_", "png", 300, width=NA)

#perform t-test analysis (with FDR correction)
mSet <- Ttests.Anal(mSet, FALSE, 0.1, FALSE, TRUE, TRUE)
mSet <- PlotTT(mSet, "tt_0_", "png", 300, width=NA)

#perform Wilcoxon rank sum tests (with FDR correction)
mSet <- Ttests.Anal(mSet, TRUE, 0.1, FALSE, TRUE, TRUE)
mSet <- PlotTT(mSet, "wilc_0_", "png", 300, width=NA)

#Check PLS-DA feasibility
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);
mSet<-PlotPLS3DLoading(mSet, "pls_loading3d_0_", "json", 1,2,3)
mSet<-PLSDA.CV(mSet, "L",5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
