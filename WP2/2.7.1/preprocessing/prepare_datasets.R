# Prepare datasets (split training and tests) for AD vs CN
########################### DATASET SPLITTING ##################################
set.seed(1)

# load data
all_ids<-read.table('/home/cmengoni/DiGAS/pheno/RS_pheno.txt')
colnames(all_ids)<-c('ID','PHENO')

# save length of dataset and number of CN/AD
n<-dim(all_ids)[[1]]
n_pheno<-n/2

# create testing set (20% -> 152 samples (76+76)) used only for netDx validation
test_CN<-sample(all_ids[all_ids$PHENO == 'CN','ID'],n_pheno*0.2)
test_AD<-sample(all_ids[all_ids$PHENO != 'CN','ID'],n_pheno*0.2)

test_set<-all_ids[all_ids$ID %in% c(test_CN,test_AD),]

# training set (80% -> 614 samples (307+307))
training_set<-all_ids[!(all_ids$ID %in% test_set$ID),]
write.table(training_set, '/home/cmengoni/DiGAS/pheno/test_missingness/original.txt', quote=F,row.names = F,col.names=F, sep='\t')
 