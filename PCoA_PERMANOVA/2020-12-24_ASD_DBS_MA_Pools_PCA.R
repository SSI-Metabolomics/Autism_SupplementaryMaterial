###############################################################
#                                                             #
#     Script to                                               #
#         1- format your feature table from the GNPS output   #
#            into MetaboAnalyst format including Pools        #
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

##creating the CSV file
write.table(ft2, 'PEGPPG_cur_68Samples_Pools_MA.csv', sep=',',col.names = NA, row.names = TRUE, quote = TRUE, na="NA")

######################################################################################################
# 2- run missing values imputation, filtering, transformation, scaling on 72 samples including Pools #
######################################################################################################

# Load the MetaboAnalystR package
library("MetaboAnalystR")

# First step is to create the mSet Object, specifying that the data to be uploaded
# is a peak table ("pktable") and that statistical analysis will be performed ("stat").
mSet<-InitDataObjects("pktable", "stat", FALSE)

# Second step is to read in the filtered peak list, please set the path right first
mSet<-Read.TextData(mSet, "PEGPPG_cur_68Samples_Pools_MA.csv", "rowu", "disc")

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

###############
# 3- run PCA  #
###############

#The fourth step is to perform Principal Component Analysis (PCA)
mSet<-PCA.Anal(mSet)

#Editing colours
colVec<-c("#a80347","#ffcd82", "#808080")
mSet<-UpdateGraphSettings(mSet, colVec, shapeVec)

#saving plots
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 300, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 300, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 300, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 300, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 300, width=NA, 1,2)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
