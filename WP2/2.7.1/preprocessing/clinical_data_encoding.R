# This script contains clinical data preprocessing

pheno<-read.csv("/home/cmengoni/INC_BOTH_AD_CL/clinical.txt",header = 1)
names(pheno)<-c('ID','AGE','STATUS','EDUCATION','SMOKING','WAIST','DIABETES','SYSTOLIC','DIASTOLIC','PULSE','LIPID_MED',
                'GLUCOSE_MED','HYPER','APOE4')

pheno$ID<-paste0('ID_',pheno$ID)
typeof(pheno$ID)

typeof(pheno$AGE)
pheno$AGE<-round(pheno$AGE,4)

typeof(pheno$STATUS)


pheno[pheno$SMOKING=='Never',"SMOKING"]<-0#0
pheno[pheno$SMOKING=='Former',"SMOKING"]<-1#0.8
pheno[pheno$SMOKING=='Current',"SMOKING"]<-1#1
pheno$SMOKING<-as.numeric(pheno$SMOKING)
typeof(pheno$SMOKING)pathFile

table(pheno$EDUCATION)
pheno[pheno$EDUCATION=='primary education',"EDUCATION"]<-1
pheno[pheno$EDUCATION=='lower/intermediate general education OR lower vocational education',"EDUCATION"]<-2
pheno[pheno$EDUCATION=='intermediate vocational education OR higher general education',"EDUCATION"]<- 2
pheno[pheno$EDUCATION=='higher vocational education OR university',"EDUCATION"]<-3
pheno$EDUCATION<-as.numeric(pheno$EDUCATION)
typeof(pheno$EDUCATION)

typeof(pheno$WAIST)
typeof(pheno$SYSTOLIC)
typeof(pheno$DIASTOLIC)
typeof(pheno$PULSE)

typeof(pheno$DIABETES)

pheno[pheno$LIPID_MED=='yes','LIPID_MED']<-1
pheno[pheno$LIPID_MED=='no','LIPID_MED']<-0


pheno[pheno$HYPER=='yes','HYPER']<-1
pheno[pheno$HYPER=='no','HYPER']<-0

pheno[pheno$GLUCOSE_MED=='yes','GLUCOSE_MED']<-1
pheno[pheno$GLUCOSE_MED=='no','GLUCOSE_MED']<-0

pheno[pheno$APOE4==2,"APOE4"]<-1
typeof(pheno$APOE4)

write.table(pheno, "/home/cmengoni/INC_BOTH_AD_CL/clinical_encoded.txt",sep='\t',quote=F,row.names = F)
