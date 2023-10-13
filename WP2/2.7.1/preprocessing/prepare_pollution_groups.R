# In this script there is the code for detecting the two pollution groups (HI and LO)
# in the control samples for which pollution information is available

# Descriptive analysis
pheno_original<-read.table("/home/cmengoni/INC_BOTH_AD_CL/clinical_encoded.txt",header = 1)
pollutants<-read.csv("/home/cmengoni/pollutants/e5.dat",sep='\t')
pollutants_e4<-read.csv("/home/cmengoni/pollutants/e4.dat",sep='\t')
pollutants_e4$ergoid<-paste0('ID_',pollutants_e4$ergoid)
pollutants$ergoid<-paste0('ID_',pollutants$ergoid)

pollutants<-pollutants[!is.na(pollutants$pm10tselect),]
pollutants_e4<-pollutants_e4[!is.na(pollutants_e4$pm10tselect),]

length(pheno_original$ID) #766
length(intersect(pheno_original$ID,pollutants$ergoid)) #400

only_e5 <-intersect(pheno_original$ID,pollutants$ergoid)
only_e4<-intersect(setdiff(pollutants_e4$ergoid,pollutants$ergoid),pheno_original$ID) #268

pollutants_e4_ft <- pollutants_e4[pollutants_e4$ergoid %in% only_e4,c('ergoid','pm10tselect',
                                                                      'pm25tselect','absorbancet',
                                                                      'noxtselecti','no2tselecti')]
pollutants<-pollutants[pollutants$ergoid %in% only_e5,c('ergoid','pm10tselect',
              'pm25tselect','absorbancet',
              'noxtselecti','no2tselecti')]

all_pollutants<-rbind(pollutants_e4_ft,pollutants) 
 

# check AD status
MERGED_DATASET<-merge(pheno_original,all_pollutants, by.x='ID',by.y='ergoid') #668
write.table(MERGED_DATASET,'/home/cmengoni/pollutants/clinical_encoded_with_pollution.tsv') #clinical data for pollutants encoded
MERGED_DATASET_AD<-MERGED_DATASET[MERGED_DATASET$STATUS=='DE',] #342
MERGED_DATASET_CON<-MERGED_DATASET[MERGED_DATASET$STATUS!='DE',] #326


######## PREPARE DATASET FOR HI vs LO IN CONTROL################
pca<-prcomp(MERGED_DATASET_CON[,c("pm10tselect","pm25tselect","absorbancet",
                              "noxtselecti","no2tselecti")], center=TRUE, scale. = TRUE)
summary(pca)

library(devtools)
install_github('vqv/ggbiplot')
library(ggbiplot)
png("/home/cmengoni/pollutants/pca.png")
ggbiplot(pca, groups=MERGED_DATASET_CON$STATUS,ellipse=TRUE)
dev.off()

MERGED_DATASET_CON$PC1<-pca$x[,1]
table(MERGED_DATASET_CON[MERGED_DATASET_CON$PC1<=median(MERGED_DATASET_CON$PC1),'STATUS'])
table(MERGED_DATASET_CON[MERGED_DATASET_CON$PC1>median(MERGED_DATASET_CON$PC1),'STATUS'])

pollutants_cols<-MERGED_DATASET_CON[,c("pm10tselect","pm25tselect","absorbancet",
                                   "noxtselecti","no2tselecti")]
pollutants_cols<-scale(pollutants_cols)
k<-kmeans(pollutants_cols, 2,nstart = 10)
MERGED_DATASET_CON$Kmeans<-k$cluster
table(MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==1,"STATUS"])
table(MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==2,"STATUS"]) #high
mean(MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==2,]$pm10tselect)
mean(MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==1,]$pm10tselect)

png("/home/cmengoni/pollutants/pca_kmeansgroups.png")
ggbiplot(pca, groups=MERGED_DATASET_CON$Kmeans,ellipse=FALSE)
dev.off()

#MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==1,c(1:14)]
MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==1,'Kmeans']<-'HI'
MERGED_DATASET_CON[MERGED_DATASET_CON$Kmeans==2,'Kmeans']<-'LO'

MERGED_DATASET_CON[,'ID'] = unlist(lapply(MERGED_DATASET_CON$ID, function(x){str_split(x,'_')[[1]][2]}))
MERGED_DATASET_CON[,'STATUS']<-MERGED_DATASET_CON[,'Kmeans']
write.table(MERGED_DATASET_CON[,c(1:14)],'/home/cmengoni/pollutants/CN_all_dataset.tsv',sep='\t')
write.table(MERGED_DATASET_CON[,c(1,3)],'/home/cmengoni/pollutants/CN_RSID_all_dataset.tsv',sep='\t')


training_ids_LO_subset<-sample(MERGED_DATASET_CON[MERGED_DATASET_CON['STATUS']=='LO','ID'],size=length(MERGED_DATASET_CON[MERGED_DATASET_CON['STATUS']=='HI','ID']))
MERGED_DATASET_CON_BALANCED <- MERGED_DATASET_CON[MERGED_DATASET_CON$ID%in%c(training_ids_LO_subset,MERGED_DATASET_CON[MERGED_DATASET_CON['STATUS']=='HI','ID']),]
write.table(MERGED_DATASET_CON_BALANCED[,c(1:14)],'/home/cmengoni/pollutants/CN_all_dataset_BALANCED.tsv',sep='\t')
write.table(MERGED_DATASET_CON_BALANCED[,c(1,3)],'/home/cmengoni/pollutants/CN_RSID_all_dataset_BALANCED.tsv',sep='\t')

