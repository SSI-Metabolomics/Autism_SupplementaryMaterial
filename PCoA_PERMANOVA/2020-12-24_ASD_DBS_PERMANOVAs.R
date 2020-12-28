###############################################################
#                                                             #
#     Script to run PERMANOVAs                                #
#                                                             #  
###############################################################


#Make sure to get error messages in English
Sys.setenv(LANG = "en")
#set the working directory where you have your files and want to save the output
setwd("C:/ASD_DBS") 

#load libraries
library(ggfortify)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(rlang)
library(vegan)
library(viridis)

#load feature table and format it
ft <- read.table ("PEGPPG_cur_68Samples_Pools.csv", sep=";", header=TRUE, dec=".", check.names = FALSE)
colnames(ft) <- gsub(' filtered Peak area','',colnames(ft))
if (length(which(is.na(colSums(ft)))) != 0){
  ft <- ft[,-which(is.na(colSums(ft)))]
}

ft[1:5,1:5]

ft <- ft[,-c(2:3)]
ft <- t(ft)
colnames(ft) <- ft[1,]
ft <- ft[-1,]
class(ft) <-"numeric"

ft[1:5,1:5]

ft[ft==0] <- NA

#Assess sparsity of featuretable
length(which(is.na(ft_nopools)))/length(ft_nopools)

#load metadata table
md <- read.table('PEGPPG_cur_68Samples_Pools_metadata.txt', sep = '\t', header = T)
head(md)

colnames(ft)[-which(colnames(ft) %in% md$filename)]

md <- md[match(rownames(ft),md$filename),]
identical(as.character(md$filename),as.character(rownames(ft)))

md <- md[which(md$filename%in%rownames(ft)),]
md <- md[match(rownames(ft),md$filename),]

sels <- c('ATTRIBUTE_injOrder',
          'ATTRIBUTE_type',
          'ATTRIBUTE_ACAT',
          'ATTRIBUTE_XX',
          'ATTRIBUTE_BW',
          'ATTRIBUTE_GA',
          'ATTRIBUTE_AAS',
          'ATTRIBUTE_MoB')

#clean metadata 
for (i in 1:length(sels)){
  
  if('na' %in% names(table(md[,which(colnames(md) == sels[i])], exclude = NULL))){
    
    md[,which(colnames(md) == sels[i])][which(md[,which(colnames(md) == sels[i])] == 'na')] <- NaN
    md[,which(colnames(md) == sels[i])] <- droplevels(md[,which(colnames(md) == sels[i])])
  }
}

md$ATTRIBUTE_BW <- as.numeric(as.character(md$ATTRIBUTE_BW))
md$ATTRIBUTE_GA <- as.numeric(as.character(md$ATTRIBUTE_GA))
md$ATTRIBUTE_AAS <- as.numeric(as.character(md$ATTRIBUTE_AAS))

#run PERMANOVAs
distmetrics <- c('euclidean')

#remove pools
ft_nopools <- ft[-which(rownames(ft) %in% as.character(md$filename[which(md$ATTRIBUTE_type == 'Pool')])),]
md_nopools <- md[match(rownames(ft_nopools),md$filename),]
identical(as.character(md_nopools$filename), rownames(ft_nopools)) 

Rsqs <- c()
pvals <- c()
#samprem <- c()
#samptype <- c()

for (j in 1:length(distmetrics)){
  
  for (i in 1:length(sels)){
    
    catcor <- sels[i]
    message(catcor)
    
    distm <- vegdist(ft_nopools,method = distmetrics[j], na.rm = T)
    message(distmetrics[j])
    
    if (length(which(is.na(md_nopools[,colnames(md_nopools) == catcor]))) !=0){
      
      red <- as.dist(as.matrix(distm)[-which(is.na(md_nopools[,colnames(md_nopools) == catcor])),-which(is.na(md_nopools[,colnames(md_nopools) == catcor]))])
      md_red <- md_nopools[-which(is.na(md_nopools[,colnames(md_nopools) == catcor])),]
      Rsq <- adonis(red ~ na.omit(md_red[,colnames(md_red) == catcor]))$aov.tab$R2[1]
      pval <- adonis(red ~ na.omit(md_red[,colnames(md_red) == catcor]))$aov.tab$'Pr(>F)'[1]
      
    } else {
      
      Rsq <- adonis(distm ~ md_nopools[,colnames(md_nopools) == catcor])$aov.tab$R2[1]
      pval <- adonis(distm ~ md_nopools[,colnames(md_nopools) == catcor])$aov.tab$'Pr(>F)'[1]
      
    }
    
    Rsqs <- c(Rsqs,Rsq)
    pvals <- c(pvals,pval)
    
  }
  
}

out <- cbind(sels,Rsqs,pvals,rep(distmetrics,each = length(sels)))
out <- as.data.frame(out, stringsAsFactors = F)
colnames(out) <- c("metadata","AdonisR2","pvalue","metric")
head(out)

out$AdonisR2 <- as.numeric(out$AdonisR2)
out$pvalue <- as.numeric(out$pvalue)

out$sig <- out$pvalue
out$sig[out$sig>=0.05] <- "p>0.05"
out$sig[out$sig<=0.05] <- "p<0.05"
out$sig <- as.factor(out$sig)

out$stars <- as.character(out$sig)
out$stars[which(out$stars == 'p>0.05')] <- ''
out$stars[which(out$stars == 'p<0.05')] <- '*'


keys <- cbind(sels, c('Injection order','ASD (yes/no)','ASD subtype', 'Gender','Birthweight',
                      'Gestational age', 'Age at sampling','Month of birth'))
keys <- keys[order(keys[,1]),]
out$metadatanice <- keys[as.numeric(as.factor(out$metadata)),2]
out

#generate plot for figure 4
ggplot(out, aes(x=metadatanice, y=AdonisR2, fill=metric)) + 
  geom_bar(colour="black", stat="identity", position='dodge') +
  geom_text(aes(label=stars), position=position_dodge(width=0.9), vjust=-0.25) +
  ggpubr::rotate_x_text() +
  labs(y= "Adonis R2") +
  theme(axis.title.x=element_blank(), legend.position = "none")


pdf(file="Figure4.pdf", width=4, height=5)
ggplot(out, aes(x=metadatanice, y=AdonisR2, fill=metric)) + 
  geom_bar(colour="black", stat="identity", position='dodge') +
  geom_text(aes(label=stars), position=position_dodge(width=0.9), vjust=-0.25) +
  ggpubr::rotate_x_text() +
  labs(y= "Adonis R2") +
  theme(axis.title.x=element_blank(), legend.position = "none")
dev.off()

write.table(out, 'PERMANOVA_Summary.txt' ,sep = '\t', quote = F, row.names = F)

