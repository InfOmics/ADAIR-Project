# In this script genotype data is preprocessed and made ready to be run in DiGAS.
# Below, there is also the code to preprocess DiGAS results and made the matrix ready for netDx

snpstogenes<-read.csv2("/home/cmengoni/digas_out/digastmp/SnpToGene_RSID.txt",sep='\t') #28,418 genes 

library(data.table)
mat<-data.table::fread("/home/cmengoni/digas_out/digastmp/merged.raw")
mat<-data.frame(mat)

#remove final underscores
colnames(mat)<-substr(colnames(mat),start=1,stop=nchar(colnames(mat))-2)
length(intersect(colnames(mat),unlist(strsplit(snpstogenes$SNPsNames_Converter,','))))  
mat_only<-(setdiff(colnames(mat),unlist(strsplit(snpstogenes$SNPsNames_Converter,',')))) 
snp_only<-(setdiff(unlist(strsplit(snpstogenes$SNPsNames_Converter,',')),colnames(mat))) 

#remove final underscores in genes
library(EnsDb.Hsapiens.v79)
fix_ensembl<- strsplit(snpstogenes$Genes,'.',fixed=T)
fix_ensembl<-unlist(lapply(fix_ensembl,function(x){x[1]}))

#transform names ens to id
snpstogenes$Genes<-fix_ensembl
symbols <- ensembldb::select(EnsDb.Hsapiens.v79, keys= fix_ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
trans_snpstogenes<-merge(symbols,snpstogenes,by.x='GENEID',by.y='Genes') #18,653

############### RUN DIGAS (external script) #######################

# Once digas has been executed these are the steps necessary to prepare the matrix for netDx:

digas_genes<-read.table("/home/cmengoni/digas_new_runs/digas_original.txt") 
digas_genes<-digas_genes$V3[-c(1,2)] 
digas_genes<- strsplit(digas_genes,'.',fixed=T)
digas_genes<-data.frame(unlist(lapply(digas_genes,function(x){x[1]})))
colnames(digas_genes)<-'Genes'
  
trans_snpstogenes<-trans_snpstogenes[trans_snpstogenes$GENEID%in%digas_genes$Genes,]

GENETOPATIENT<-c()
for (i in unique(trans_snpstogenes$SYMBOL)){
  SNPS<-trans_snpstogenes[trans_snpstogenes$SYMBOL==i,]$SNPsNames_Converter

  SNPS<-unique(unlist(strsplit(SNPS,','))) 
  SNPS<-intersect(SNPS,colnames(mat))

  if(length(SNPS)>0){
    values<-rowSums(mat[SNPS])
    GENETOPATIENT[[i]]<-values/length(SNPS)  
  }
}

GENETOPATIENT_df<-data.frame(GENETOPATIENT) 
GENETOPATIENT_df_binary<-GENETOPATIENT_df
GENETOPATIENT_df_binary[GENETOPATIENT_df_binary > 0] = 1


hist(rowSums(GENETOPATIENT_df_binary)) # patients has x snp
hist(colSums(GENETOPATIENT_df_binary)) # snp has x patients


pheno<-read.table("/home/cmengoni/DiGAS/pheno/RS_pheno.txt")
de_IDX<-pheno[pheno$V2=='DE',]$V1
cn_IDX<-pheno[pheno$V2=='CN',]$V1
IDX<-c(de_IDX,cn_IDX)
IDX<-paste0('ID_',IDX)
rownames(GENETOPATIENT_df_binary)<-IDX
GENETOPATIENT_df_binary<-t(GENETOPATIENT_df_binary)
saveRDS(GENETOPATIENT_df_binary,"/home/cmengoni/digas_new_runs/out/genotype_matrix_all.rds")


