# Preparatory script for methylkit_analysis.R 

library(bsseq)
files_dir='/home/cmengoni/methyl/bsseq_data'
files=list.files(files_dir)
x=list.files(files_dir)
lapply(files,function(x){
    print(paste(files_dir,x,sep='/'))
    name <- strsplit(x,'.',fixed=TRUE)[[1]][1]
    
    load(paste(files_dir,x,sep='/'))
    seqs<-seqnames(bismarkBSseq@rowRanges)
    base<-ranges(bismarkBSseq@rowRanges)
    
    CHR<-c()
    
    counter<-0
    for (i in seqs@values){
        counter<-counter+1
        print(i)
        CHR<-c(CHR,rep(i,seqs@lengths[counter]))
    }
    
    CHR <- unlist(lapply(CHR, function(x){paste0('chr',x)}))
    hea.con.M<-getCoverage(bismarkBSseq, type="M")
    hea.con.cov<-getCoverage(bismarkBSseq, type="Cov")
    
    for ( i in c(1:dim(hea.con.cov)[2])){
        cov_value <-as.numeric(hea.con.cov[,i])
        freqC_value <- as.numeric(hea.con.M[,i]/hea.con.cov[,i])
        freqT_value<- 1-freqC_value
        
        tmp_df <- data.frame('chrBase'=paste(CHR, base, sep='.'),
                             'chr'=CHR,
                             'base'=base@start,
                             'strand'=rep('*',dim(hea.con.cov)[1]),
                             'coverage'=cov_value,
                             'freqC'=freqC_value*100,
                             'freqT'=freqT_value*100)
        tmp_df[is.na(tmp_df)] <- 0
        
        write.table(tmp_df,
                    paste0(paste(paste0('/home/cmengoni/methyl/methylkit_data/',name),colnames(hea.con.cov)[i],sep='_'),'.txt'), sep='\t', quote=FALSE,row.names = FALSE)
    }
    
})
