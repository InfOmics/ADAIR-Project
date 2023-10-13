# In this script the pathways are prepared for grouping metabolites and SNPs

pathways_path<-"/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways.rds"
#Prepare pathways 
pathList=readRDS(pathways_path)
excl=-grep("KEGG",names(pathList))
pathList=pathList[excl]
min=30;max=180;
dim_sets=sapply(pathList,length)
keep_bol=(dim_sets>=min & dim_sets<=max);cat(sum(keep_bol),"\n")
pathList=pathList[keep_bol]
#Simplify the name of the pathways
isi_names=paste("pathwayK",seq(1,length(pathList)),sep="")
#Create map and replace
map_pathways=cbind(names(pathList),isi_names)
colnames(map_pathways)=c("pathway_ids","isi_names")
names(pathList)=isi_names
saveRDS(pathList,'/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_mock.rds')
saveRDS(map_pathways,'/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_mapper.rds')

#################################################################################
metab_pathways<-read.table(file = '/home/cmengoni/INC_BOTH_AD_CL/meta_pathways.txt')
metab_pathways_l<-c()
for (m in unique(metab_pathways$V2)){
  metab_pathways_l[[m]]<- metab_pathways[metab_pathways$V2==m,]$V1
}
saveRDS(metab_pathways_l,'/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_meta.rds')

#################################################################################
pathList <- readPathways(fetchPathwayDefinitions("April",2022)) #FILTER: sets with num genes in [10, 200]
saveRDS(pathList,'/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_NETDX.rds')

