#R-codes for Figures and statistical analyses

###############################################################################
#PCoA statistical analyses and biplot generation
###############################################################################
#Random forest method to measure variable importance
#Extract top10 most discriminating ARG-MRG types across all genomes for biplot
install.packages('randomForest')
install.packages('ggplot2')
install.packages('vegan')
install.packages('BiodiversityR')
install.packages('scales')
library(randomForest)
library(ggplot2)
library(vegan)
library(BiodiversityR)
library(scales)
options(stringsAsFactors = F)

ARGMRG_tab<-read.table(file="file.name",header=T,sep="\t") #the file to read is abundance data frame with column names of genome ID,all ARG-MRG types
ARGMRG_tab<-ARGMRG_pathogen_tab[,c(TRUE,colSums(ARGMRG_pathogen_tab[,2:272])>=10)] #filter the ARG-MRG types with less than 10 across all geomes
ARGMRG_pathogen_tab<-merge(ARGMRG_tab,Pathogen_list,by="genomeID") #Pathogen_list is pathogen status information of all genomes with column names of genome ID and pathogen status
ARGMRG_pathogen_tab$pathogen<-as.factor(ARGMRG_pathogen_tab$pathogen) #convert to factor for random forest analysis
ARGMRG_pathogen_rf<-randomForest(pathogen~., data=ARGMRG_pathogen_tab,importance=TRUE)
ARGMRG_pathogen_imp<-importance(ARGMRG_pathogen_rf) #extract importance values from random forest analysis
ARGMRG_pathogen_imp<-ARGMRG_pathogen_imp[order(ARGMRG_pathogen_imp[,c("MeanDecreaseAccuracy")],decreasing=T),] #order importance value 
Top10_ARGMRG_pathogen<-ARGMRG_pathogen_imp[1:10,1] #extract top 10 most discriminating ARG-MRG

# PCoA statistical analyses
rownames(ARGMRG_tab)<-ARGMRG_tab[,1]
ARGMRG_tab<-ARGMRG_tab[,-1]
ARGMRG_tab<-ARGMRG_tab[!rowSums(ARGMRG_tab)==0,]
dist_ARGMRG_tab<-vegdist(ARGMRG_tab,method="bray") #Bray-Curtis distance calculation 
cmdscale_ARGMRG_tab<-cmdscale(dist_ARGMRG_tab, k=nrow(ARGMRG_tab)-1,eig=T, add=F) #PCoA calculation
cmdscale_ARGMRG_tab<- add.spec.scores(cmdscale_ARGMRG_tab,ARGMRG_tab,method="wa.scores", Rscale=T, scaling=1, multi=1) #Calculate scores (coordinates) for ARG-MRG types

#prepare for plotting
pcoa_ARGMRG_tab<-as.data.frame(cmdscale_ARGMRG_tab[[1]][,1:2]) #extract data frame for PCoA plotting
colnames(pcoa_ARGMRG_tab)<-c("PC1","PC2")
pcoa_ARGMRG_tab$genomeID<-rownames(pcoa_ARGMRG_tab)
pcoa_ARGMRG_pathogen_tab<-merge(pcoa_ARGMRG_tab,Pathogen_list,by="genomeID") #get pathogen status of genomes
pcoa_ARGMRG_pathogen_tab$pathogen<-gsub("pathogen","red",pcoa_ARGMRG_pathogen_tab$pathogen,fixed=T) #color for pathogen
pcoa_ARGMRG_pathogen_tab$pathogen<-gsub("nonred","gold",pcoa_ARGMRG_pathogen_tab$pathogen,fixed=T) #color for non-pathogen

pcoa_ARGMRG_pathogen_10rf<-as.data.frame(pcoa_ARGMRG_pathogen_tab[[6]][rownames(pcoa_ARGMRG_pathogen_tab[[6]])%in%Top10_ARGMRG_pathogen,1:2]) #extract scores (coordinates) for top 10 most discriminating ARG-MRG types
pcoa_ARGMRG_pathogen_10rf$name<-rownames(pcoa_ARGMRG_pathogen_10rf)
Top10_ARGMRG_abundsize<-data.frame(colSums(ARGMRG_tab[,Top10_ARGMRG_pathogen])) #total abundance of the top 10 ARG-MRG types across genomes, used as biplot point size  
Top10_ARGMRG_abundsize$name<-rownames(Top10_ARGMRG_abundsize)
colnames(Top10_ARGMRG_abundsize)[1]<-c("abundancescale")
pcoa_ARGMRG_pathogen_10rf<-merge(pcoa_ARGMRG_pathogen_10rf,Top10_ARGMRG_abundsize,by="name") #merge with coordinate table of the top10 ARG-MRG types

#plotting
png("MRGARG_pathogen_biplot.png", res = 300,width = 11, height = 9, units = 'in')
ggplot()+geom_point(data=pcoa_ARGMRG_pathogen_tab,aes(x=PC1, y=PC2),size=3,colour=alpha(pcoa_ARGMRG_pathogen_tab$pathogen,0.6))+
  theme(text=element_text(size=15))+
  xlab("PC1 - Percentage variation explained (--%)")+
  ylab("PC2 - Percentage variation explained (--%)")+
  geom_point(data=pcoa_ARGMRG_pathogen_10rf,aes(x=Dim1, y=Dim2,size=abundancescale),colour=alpha("blue",0.7))
  geom_text(data=pcoa_ARGMRG_pathogen_10rf,aes(x=Dim1,y=Dim2,label=name))
dev.off()





