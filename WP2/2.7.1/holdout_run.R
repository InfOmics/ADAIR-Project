# This script is used to run the classifcation on the test set (holdout). 
# By changing run name and classes it was ran on the different tested datasets
.libPaths(.libPaths()[2:4])
.libPaths()

# Packages
library(netDx)
library(SummarizedExperiment)
library(parallel)
library(doParallel)
library(MultiAssayExperiment)

# Parameters
run_name<-'ADvsCN'
classes<-c('DE','CN')
load((paste0("/home/cmengoni/netDx_out_",run_name)))
outDir = paste0('/home/cmengoni/holdout_',run_name)
dir.create(outDir)
  
# define similarities 
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
                                         groupList[['meta']],netDir,
                                         verbose=T,sparsify=T)
   }
   
   
   # make clinical nets
   netList3<-c() 
   if(!is.null(groupList[['clinical']])){
      netList3<-makePSN_NamedMatrix_clinical(assay(dataList[['clinical']]),
                                             rownames(dataList[['clinical']]),
                                             groupList[['clinical']],netDir,
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
            sim <- round(cor(na.omit(xpr[idx, , drop = FALSE]), method = "pearson"), digits = 3) 
            #getSimilarity(xpr[idx, , drop = FALSE],  #avgNormDiff(xpr[idx, , drop=FALSE]) #
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
 
# Run classification
for ( cutoff in c(0.5,0.6,0.7,0.8,0.9,1)){

   results <- getResults(
      out,
      classes,
      featureSelCutoff=9L,
      featureSelPct=cutoff
    ) 
   
   outDir_cutoff = paste0(outDir,paste0("tmp_",cutoff))
   dir.create(outDir_cutoff)
   
   predModel<-predict(trainMAE=dataList,testMAE=holdout,
                       groupList = groupList,
                       selectedFeatures = results$selectedFeatures, 
                       outDir=outDir_cutoff,makeNetFunc = makeNets)
   
   saveRDS(predModel, paste0(outDir,paste0("predModel_",cutoff)))
} 
