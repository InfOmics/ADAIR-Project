# This script contains the commands and functions used to run NetDx on the Rotterdam Study data.
# This takes as input the genotype matrix preprocessed with DiGAS,  metabolomics data and pollution data.
# In this script netDx is trained on 80% of samples (AD vs CN)
 
#!/usr/bin/Rscript

assign('.lib.loc','/usr/local/lib/R/site-library/',envir=environment(.libPaths))
.libPaths(.libPaths()[2:4])
.libPaths()

# Packages
library(netDx)
library(SummarizedExperiment)
library(parallel)
library(doParallel)
library(MultiAssayExperiment)

# Parameters
nco<-5
run_name<-'pollutants'

# 1) Prepare data
## list with assays and phenotype. Should be named with assay names and  "pheno"
pheno<-read.table("/home/cmengoni/pollutants/clinical_encoded_with_pollution.tsv",header = 1)
rownames(pheno)<-pheno$ID

geno<-readRDS("/home/cmengoni/digas_new_runs/out/genotype_matrix_all.rds")
geno<-geno[,rownames(pheno)] #put them in the same order as meta and pheno
#head(geno)

meta<-read.table("/home/cmengoni/INC_BOTH_AD_CL/metabolomics_preprocessed.txt",header = 1)  
meta<-meta[,rownames(pheno)]
meta<-as.matrix(meta+0.0000000001)
meta<-log10(meta)
meta<-meta-min(meta)
head(meta)

## split clinical data from phenotype
clinical<-pheno[,!colnames(pheno)%in%c('ID','STATUS',   "pm10tselect", "pm25tselect", "absorbancet", "noxtselecti", "no2tselecti")]
#clinical<-pheno[,!colnames(pheno)%in%c('ID','STATUS')]

pheno<-pheno[,c('ID','STATUS')]

## objList
clinicaldata<-SummarizedExperiment(t(clinical),colData=pheno)  
snpdata<-SummarizedExperiment(geno,colData=pheno)  
metabdata<-SummarizedExperiment(meta,colData=pheno) 
objList<-list(genetic=snpdata, clin=clinicaldata,meta=metabdata)

## create groupList 
groupList<-list()
### create GroupList for snps
pathFile_g<-'/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_NETDX.rds'
pathList_g<-readRDS(pathFile_g)
 
i=0
COUNTER=0
keep<-c()
for (p in pathList_g){
  i=i+1
  inte<-(intersect(p,rownames(geno)))
  if (length(inte>0)){
    COUNTER=COUNTER+1
    keep<-c(keep,names(pathList_g[i]))
  }
}
 
sprintf("PATHWAYS WITH AT LEAST ONE GENE: %i",COUNTER)
pathList_g<-pathList_g[keep]
groupList[['genetic']]=pathList_g


# create GroupList for clinical
groupList[['clin']]<-list(AGE='AGE',EDUCATION='EDUCATION',SMOKING='SMOKING',
                          WAIST='WAIST',DIABETES='DIABETES',SYSTOLIC='SYSTOLIC',
                          DIASTOLIC='DIASTOLIC',PULSE='PULSE',LIPID_MED='LIPID_MED',
                          GLUCOSE_MED='GLUCOSE_MED',HYPER='HYPER',APOE4='APOE4')
                          #PM10='pm10tselect',PM25='pm25tselect',ABSORBANCE='absorbancet',NOX='noxtselecti',NO2='no2tselecti')

pathFile_m<-"/home/cmengoni/INC_BOTH_AD_CL/ext_data/Pathways_meta.rds"
pathList_m<-readRDS(pathFile_m)
groupList[['meta']]=pathList_m

## create multiassayexpriment  
dataList<-MultiAssayExperiment(objList,pheno)   
summary(groupList)

## 2) define similarities 
makeNets<-function(dataList,groupList,netDir,numCores,...){
  # make genotype nets
  netList1<-c()
  if (!is.null(groupList[['genetic']])){
    netList1<-makeMutNets(assay(dataList[['genetic']]),
                          groupList[['genetic']], netDir,
                          numC = nco)
  }
  
  # # make metabolite net
  # netList2<-c() 
  # if(!is.null(groupList[['meta']])){ 
  #   netList2<-makePSN_NamedMatrix_meta(assay(dataList[['meta']]),
  #                                 rownames(dataList[['meta']]),
  #                                 groupList[['meta']],netDir, 
  #                                 simMetric = 'custom',customFunc=avgNormDiff,
  #                                 verbose=T,writeProfiles = F, sparsify = T)
  # }
  netList2<-c()
  if(!is.null(groupList[['meta']])){
    netList2<-makePSN_NamedMatrix_meta(assay(dataList[['meta']]),
                                       rownames(dataList[['meta']]),
                                       groupList[['meta']],netDir, simMetric = 'custom',customFunc=avgNormDiff,
                                       verbose=T,sparsify=T, writeProfiles = F)
  }
  
  
  # make clinical nets
  netList3<-c() 
  if(!is.null(groupList[['clin']])){
    netList3<-makePSN_NamedMatrix_clinical(assay(dataList[['clin']]),
                                           rownames(dataList[['clin']]),
                                           groupList[['clin']],outDir=netDir,
                                           simMetric = 'custom',customFunc=normDiff,
                                           writeProfiles = F, sparsify=T, verbose=T)
  }
  netList<-c(unlist(netList1),unlist(netList2),unlist(netList3))
}


makeMutNets<- function(g,pList,oDir, numC){ 
  g<-t(g) # transpose to have genes as columns
  cl<-makeCluster(numC)
  registerDoParallel(cl)
  
  numPat<-c()
  netList<-foreach(k=1:length(pList)) %do% { #for each pathway
    idx<-which(colnames(g)%in%pList[[k]])
    cat(sprintf('%s \t %i \n',names(pList)[k],length(idx)))
    
    if(length(idx)>0){ # if at least one gene of the pathway exists in the dataset
      
      has_mut<-rowSums(g[,idx,drop=F]) #count how many patients have anygene mutated in the pathway
      has_mutp<-names(has_mut)[which(has_mut>0)] #patients with mutation
      
      if(length(has_mutp)>=6){ #at least 6 patients
        
        cat(sprintf("%s: %i patients\n",names(pList)[k],length(has_mutp))) # this pathway has x patients with genes mutated
        pat_pairs<-t(combinat::combn(has_mutp,2)) # assign similarity of each couple of mutated patients to 1
        pat_pairs<-cbind(pat_pairs,1)
        outFile<-sprintf("%s/%s_cont.txt",oDir,names(pList)[k])
        write.table(pat_pairs,file=outFile,sep='\t',col=F,row=F,quote=F)
        basename(outFile)
      } else NULL
    } else NULL
  }
  stopCluster(cl)
  unlist(netList)
}

makePSN_NamedMatrix_clinical <- function (xpr, nm, namedSets, outDir = tempdir(), simMetric = "pearson", 
                                          verbose = TRUE, numCores = 1L, writeProfiles = TRUE, sparsify = FALSE, 
                                          useSparsify2 = FALSE, cutoff = 0.3, sparsify_edgeMax = Inf, 
                                          sparsify_maxInt = 50, minMembers = 1L, runSerially = FALSE, 
                                          ...) 
{ 
  
  if ((!simMetric %in% c("pearson", "MI")) & writeProfiles == 
      TRUE) {
    print(simMetric)
    stop(paste("writeProfiles must only be TRUE with simMetric", 
               " set to pearson or MI. For all other metrics, ", 
               "set writeProfiles=FALSE", sep = ""))
  }
  cl <- makeCluster(numCores, outfile = paste(outDir, "makePSN_log.txt", 
                                              sep = getFileSep()))
  if (!runSerially) {
    registerDoParallel(cl)
  }
  else {
    message("running serially")
  }
  if (simMetric == "pearson") {
    message(paste("Pearson similarity chosen - ", "enforcing min. 5 patients per net.", 
                  sep = ""))
    minMembers <- 5
  }
  #`%myinfix%` <- ifelse(runSerially, `%do%`, `%dopar%`)
  outFiles<-c()
  i = 0
  for(curSet in names(namedSets))
  { library(netDx, lib.loc = "/usr/local/lib/R/site-library")
    i =i+1
    if (verbose) 
      message(sprintf("%s: ", curSet))
    idx <- which(nm %in% namedSets[[curSet]])
    if (verbose) 
      message(sprintf("%i members", length(idx)))
    oFile <- NULL
    if (length(idx) >= minMembers) {
      if (writeProfiles) {
        outFile <- paste(outDir, sprintf("%s.profile", 
                                         curSet), sep = getFileSep())
        write.table(t(xpr[idx, , drop = FALSE]), file = outFile, 
                    sep = "\t", dec = ".", col.names = FALSE, 
                    row.names = TRUE, quote = FALSE)
      }
      else {
        outFile <- paste(outDir, sprintf("%s_cont.txt", 
                                         curSet), sep = getFileSep())
        message(sprintf("computing sim for %s", curSet))
        sim <-normDiff(xpr[idx, , drop=FALSE]) 
        if (is.null(sim)) {
          stop(sprintf(paste("makePSN_NamedMatrix:%s: ", 
                             "similarity matrix is empty (NULL).\n", 
                             "Check that there isn't a mistake in the ", 
                             "input data or similarity method of choice.\n", 
                             sep = ""), curSet))
        }
        pat_pairs <- sim
        if (sparsify) {
          if (useSparsify2) {
            tryCatch({
              spmat <- sparsify2(pat_pairs, cutoff = cutoff, 
                                 EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                                 maxInt = sparsify_maxInt)
            }, error = function(ex) {
              stop("sparsify2 caught error\n")
            })
          }
          else {
            message("sparsify3")
            tryCatch({
              sp_t0 <- Sys.time()
              spmat <- sparsify3(pat_pairs, cutoff = cutoff, 
                                 EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                                 maxInt = sparsify_maxInt, verbose = FALSE)
              print(Sys.time() - sp_t0)
            }, error = function(ex) {
              stop("sparsify3 caught error\n")
            })
          }
        }
        else {
          write.table(pat_pairs, file = outFile, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE)
          print(basename(outFile))
          message("done")
        }
      }
      oFile <- basename(outFile)
    }
    outFiles[[i]]<-oFile
  }
  stopCluster(cl)
  outFiles
}


makePSN_NamedMatrix_meta <- function (xpr, nm, namedSets, outDir = tempdir(), simMetric = "pearson", 
                                      verbose = TRUE, numCores = 1L, writeProfiles = TRUE, sparsify = FALSE, 
                                      useSparsify2 = FALSE, cutoff = 0.3, sparsify_edgeMax = Inf, 
                                      sparsify_maxInt = 50, minMembers = 1L, runSerially = FALSE, 
                                      ...) 
{ 
  
  if ((!simMetric %in% c("pearson", "MI")) & writeProfiles == 
      TRUE) {
    print(simMetric)
    stop(paste("writeProfiles must only be TRUE with simMetric", 
               " set to pearson or MI. For all other metrics, ", 
               "set writeProfiles=FALSE", sep = ""))
  }
  cl <- makeCluster(numCores, outfile = paste(outDir, "makePSN_log.txt", 
                                              sep = getFileSep()))
  if (!runSerially) {
    registerDoParallel(cl)
  }
  else {
    message("running serially")
  }
  if (simMetric == "pearson") {
    message(paste("Pearson similarity chosen - ", "enforcing min. 5 patients per net.", 
                  sep = ""))
    minMembers <- 5
  }
  #`%myinfix%` <- ifelse(runSerially, `%do%`, `%dopar%`)
  outFiles<-c()
  i = 0
  for(curSet in names(namedSets))
  { library(netDx, lib.loc = "/usr/local/lib/R/site-library")
    i =i+1
    if (verbose) 
      message(sprintf("%s: ", curSet))
    idx <- which(nm %in% namedSets[[curSet]])
    if (verbose) 
      message(sprintf("%i members", length(idx)))
    oFile <- NULL
    if (length(idx) >= minMembers) {
      if (writeProfiles) {
        outFile <- paste(outDir, sprintf("%s.profile", 
                                         curSet), sep = getFileSep())
        write.table(t(xpr[idx, , drop = FALSE]), file = outFile, 
                    sep = "\t", dec = ".", col.names = FALSE, 
                    row.names = TRUE, quote = FALSE)
      }
      else {
        outFile <- paste(outDir, sprintf("%s_cont.txt", 
                                         curSet), sep = getFileSep())
        message(sprintf("computing sim for %s", curSet))
        sim <-avgNormDiff(xpr[idx, , drop=FALSE]) 
        #sim <- round(cor(na.omit(xpr[idx, , drop = FALSE]), method = "pearson"), digits = 3) 
        #getSimilarity(xpr[idx, , drop = FALSE],  #
        #type = simMetric, ...)
        if (is.null(sim)) {
          stop(sprintf(paste("makePSN_NamedMatrix:%s: ", 
                             "similarity matrix is empty (NULL).\n", 
                             "Check that there isn't a mistake in the ", 
                             "input data or similarity method of choice.\n", 
                             sep = ""), curSet))
        }
        pat_pairs <- sim
        if (sparsify) {
          if (useSparsify2) {
            tryCatch({
              spmat <- sparsify2(pat_pairs, cutoff = cutoff, 
                                 EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                                 maxInt = sparsify_maxInt)
            }, error = function(ex) {
              stop("sparsify2 caught error\n")
            })
          }
          else {
            message("sparsify3")
            tryCatch({
              sp_t0 <- Sys.time()
              spmat <- sparsify3(pat_pairs, cutoff = cutoff, 
                                 EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                                 maxInt = sparsify_maxInt, verbose = FALSE)
              print(Sys.time() - sp_t0)
            }, error = function(ex) {
              stop("sparsify3 caught error\n")
            })
          }
        }
        else {
          write.table(pat_pairs, file = outFile, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE)
          print(basename(outFile))
          message("done")
        }
      }
      oFile <- basename(outFile)
    }
    outFiles[[i]]<-oFile
  }
  stopCluster(cl)
  outFiles
}

# 3) Create Holdout
create_holdout_from_pheno<-function (dataMAE,training_IDs)
{
  pheno <- colData(dataMAE)
  idx_holdout <- which(pheno$ID %in% training_IDs)
  
  holdout <- dataMAE[, rownames(pheno)[idx_holdout]]
  colData(holdout)$ID <- as.character(colData(holdout)$ID)
  tokeep <- setdiff(1:nrow(pheno), idx_holdout)
  dataMAE <- dataMAE[, rownames(pheno)[tokeep]]
  return(list(trainMAE = dataMAE, validationMAE = holdout))
}

pheno_original<-read.table("/home/cmengoni/INC_BOTH_AD_CL/clinical_encoded.txt",header = 1)
pheno_2<-read.table('/home/cmengoni/pollutants/digas_original_pheno.txt')
training_IDs<-pheno_original[!(pheno_original$ID %in% paste0('ID_',pheno_2$V1)),]$ID
dsets<-create_holdout_from_pheno(dataList,training_IDs)
dataList<-dsets$trainMAE
holdout<-dsets$validationMAE
dataList



# 4) Build predictor
outDir<-paste('/home/cmengoni',run_name,'pred_output',sep=getFileSep())
#if(!file.exists(outDir)) unlink(outDir,recursive=T)

out<-buildPredictor(dataList=dataList,
                      groupList=groupList,
                      makeNetFunc = makeNets,
                      outDir=outDir,
                      debugMode=FALSE,
                      numCores=nco, logging='all',
                      numSplits = 100, JavaMemory = 20L,
                      featScoreMax = 10L,
                      trainProp = 0.8,
                      featSelCutoff = 9L
)

save.image(paste0("/home/cmengoni/netDx_out_",run_name))
