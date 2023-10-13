# This script contains code to preprocess metabolomics and make it ready for netDx.

# read metabolomics
m <-read.table("/home/cmengoni/INC_BOTH_AD_CL/metabolomics.txt",sep=',',header = 1)
# recode ids
rownames(m)<-paste0('ID_',m$RS_ID)
# retain metabolomics columns
m<-m[,c(5:ncol(m))]
# recode metabolomics names
colnames(m)<-sapply(colnames(m), function(x){substr(x,7,nchar(x))})

##PREPROCESSING
# Remove metabolites having more than 20% NA values (none)
colnames(m[colSums(is.na(m) | m==0)>0.2*dim(m)[1]])

# Impute missing values to half of the minimum for that metabolite
minimums<-apply(m,2,function(x){min(x[!is.na(x) & x != 0])})
minimums<-minimums/2
names(minimums)<-colnames(m)

for (c in colnames(m)){
  idxs<-which(is.na(m[,c]))
  m[idxs,c]<-minimums[c]
}

m[is.na(m)] # 0

write.table(t(m),"/home/cmengoni/INC_BOTH_AD_CL/metabolomics_preprocessed.txt",sep='\t')
