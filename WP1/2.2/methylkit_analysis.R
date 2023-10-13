# DIFFERENTIAL METHYLATION ANALYSIS OF hOM CELLS (Claudia Mengoni)
# Data: methylation of healthy (HEA) and AD (DIS) hOM cells
# exposed to clean air (CON) and Diesel Exhaust (EXP) through ALI system.
# 32 samples of which 8 HEA.CON, 8 HEA.EXP, 8 DIS.CON, 8 DIS.EXP

# Aims:
# - detect biomarkers for acute air pollutant exposure in healthy and disease.
# - detect differences between healthy and disease

# Note: before this we run the mthylkit_preparation.R script to prepare data for the analysis

library(methylKit)
library(bsseq)
library(BRGenomics) 
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


## READ IN DATA (preprocessed into tab sep file)
load('/home/cmengoni/methyl/data.rds')
# Read in data 
files_dir='/home/cmengoni/methyl/methylkit_data'
file_list=list.files(files_dir)
ids=lapply(file_list,function(x){ strsplit(strsplit(x,'.',fixed=TRUE)[[1]][1], '_', fixed=TRUE)[[1]][4] })
treat=lapply(file_list,function(x){ paste(strsplit(x, '_', fixed=TRUE)[[1]][2],strsplit(x, '_', fixed=TRUE)[[1]][3], sep='_') })
file_list=lapply(file_list,function(x){paste0('/home/cmengoni/methyl/methylkit_data/',x)})

# read the files to a methylRawList. Saved in: '/home/cmengoni/methyl/data.rds'
ad.expVSad.con=methRead(file_list[1:16],
                sample.id=ids[1:16],
                assembly="hg18",
                treatment=c(rep(0,8),rep(1,8)),
                context="CpG",
                mincov = 5, dbtype='tabix',
                dbdir='/home/cmengoni/methyl/methylDB/gene-body')
he.expVShe.con=methRead(file_list[17:32],
                        sample.id=ids[17:32],
                        assembly="hg18",
                        treatment=c(rep(0,8),rep(1,8)), #exp=1
                        context="CpG",
                        mincov = 5, dbtype='tabix',
                        dbdir='/home/cmengoni/methyl/methylDB/gene-body')
ad.conVShe.con=methRead(file_list[c(1:8,17:24)],
                        sample.id=ids[c(1:8,17:24)],
                        assembly="hg18",
                        treatment=c(rep(1,8),rep(0,8)), #ad=1
                        context="CpG",
                        mincov = 5, dbtype='tabix',
                        dbdir='/home/cmengoni/methyl/methylDB/gene-body')


# Descriptive stats
getMethylationStats(ad.expVSad.con[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(ad.expVSad.con[[1]],plot=TRUE,both.strands=FALSE)


# Unite
AD_meth=unite(ad.expVSad.con, destrand=TRUE)
HE_meth=unite(he.expVShe.con, destrand=TRUE)
CON_meth=unite(ad.conVShe.con, destrand=TRUE)

# Differential analysis on BP
AD_s_diff=calculateDiffMeth(AD_meth,mc.cores=1)
HE_s_diff=calculateDiffMeth(HE_meth,mc.cores=1)
CON_s_diff=calculateDiffMeth(CON_meth,mc.cores=1)

myDiff_s_AD=getMethylDiff(AD_s_diff,difference=0,qvalue=0.05)
myDiff_s_HE=getMethylDiff(HE_s_diff,difference=0,qvalue=0.05)
myDiff_s_CON=getMethylDiff(CON_s_diff,difference=0,qvalue=0.05)

# Annotation
aCR_s_AD<-assignChromosomeRegion(as(myDiff_s_AD,'GRanges'), nucleotideLevel = TRUE,
                                  precedence=c("Promoters", "immediateDownstream",
                                               "fiveUTRs", "threeUTRs",
                                               "Exons", "Introns"), proximal.promoter.cutoff = c(upstream = 1500, downstream = 500),
                                  TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)

tiff("/home/cmengoni/methyl/plots/aCR_s_AD.tiff")
par(mar = c(15, 4, 4, 2) + 0.1)
bp<-barplot(aCR_s_AD$percentage,main='DMP in AD EXP vs AD CON',las=2)
dev.off()

aCR_s_HE<-assignChromosomeRegion(as(myDiff_s_HE,'GRanges'), nucleotideLevel = TRUE,
                                  precedence=c("Promoters", "immediateDownstream",
                                               "fiveUTRs", "threeUTRs",
                                               "Exons", "Introns"), proximal.promoter.cutoff = c(upstream = 1500, downstream = 500),
                                  TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
tiff("/home/cmengoni/methyl/plots/aCR_s_HE.tiff")
par(mar = c(15, 4, 4, 2) + 0.1)
bp<-barplot(aCR_s_HE$percentage,main='DMP in HE EXP vs HE CON',las=2)
dev.off()

aCR_s_CON<-assignChromosomeRegion(as(myDiff_s_CON,'GRanges'), nucleotideLevel = TRUE,
                              precedence=c("Promoters", "immediateDownstream",
                                           "fiveUTRs", "threeUTRs",
                                           "Exons", "Introns"), proximal.promoter.cutoff = c(upstream = 1500, downstream = 500),
                              TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
tiff("/home/cmengoni/methyl/plots/aCR_s_CON.tiff")
par(mar = c(15, 4, 4, 2) + 0.1)
bp<-barplot(aCR_s_CON$percentage,main='DMP in AD CON vs HE CON',las=2)
dev.off()

load('/home/cmengoni/methyl/singlebase.rds')

################################################################################

## REGION ANALYSIS
library(rtracklayer)

regions =c('gene-body','promoters','mirna')#'gene-body' promoters both mirna

lapply(regions, function(x){

genes<-readGFFAsGRanges("/home/cmengoni/methyl/Homo_sapiens.GRCh38.104.gff3")
genes<-genes[!is.na(genes$biotype), ]
genes<-genes[ genes$biotype == "protein_coding" & genes$type=='gene', ]
genes<-renameSeqlevels(genes,paste0('chr',levels(genes@seqnames)))
genes<-tidyChromosomes(genes,keep.X=FALSE,keep.Y=FALSE) #18920

if (x=='gene-body'){
  all_genes<-genes
  genes<-reduce(genes)
} else if (x=='promoters'){
  genes<-promoters(genes, upstream=1500,downstream=500)
  all_genes<-genes
  genes<-reduce(genes) #make sure no overlapping sequences are present
#
# } else if (region=='both'){
#   genes<-resize(genes,width(genes) + 1500L, fix = "end")
#   all_genes<-genes
#   genes<-reduce(genes) #make sure no overlapping sequences are present

} else if (x=='mirna'){
  genes<-readGFF('/home/cmengoni/methyl/hsa.gff3')
  genes<-GRanges(genes)
  genes<-tidyChromosomes(genes,keep.X=FALSE,keep.Y=FALSE)
  #keep<-unlist(lapply(genes$Name,function(x){strsplit(x,split = '-')[[1]][2]})=='miR')
  genes<-genes[genes$type=='miRNA',]
  all_genes<-genes
}

# Map bases to regions
AD_regions=regionCounts(ad.expVSad.con,genes,strand.aware = F)#F)#gene.obj$promoters,strand.aware = F)#
HE_regions=regionCounts(he.expVShe.con,genes,strand.aware = F)#gene.obj$promoters,strand.aware = F)#
CON_regions=regionCounts(ad.conVShe.con,genes,strand.aware = F)#gene.obj$promoters,strand.aware = F)#

# Unite
AD_regions=unite(AD_regions, destrand=TRUE)
HE_regions=unite(HE_regions, destrand=TRUE)
CON_regions=unite(CON_regions, destrand=TRUE)

# # Differential analysis on regions
AD_meth_diff=calculateDiffMeth(AD_regions,mc.cores=1)
HE_meth_diff=calculateDiffMeth(HE_regions,mc.cores=1)
CON_meth_diff=calculateDiffMeth(CON_regions,mc.cores=1)

# # Significative
myDiff_AD=getMethylDiff(AD_meth_diff,difference=0,qvalue=0.05)
myDiff_HE=getMethylDiff(HE_meth_diff,difference=0,qvalue=0.05)
myDiff_CON=getMethylDiff(CON_meth_diff,difference=0,qvalue=0.05)

#save.image(paste0(paste0('/home/cmengoni/methyl/',x),'.rds'))
})

################################################################################
## ANNOTATE REGIONS

region="mirna"
#load(paste0(paste0('/home/cmengoni/methyl/',region),'.rds'))
genes<-readGFFAsGRanges("/home/cmengoni/methyl/Homo_sapiens.GRCh38.104.gff3")
genes<-genes[!is.na(genes$biotype), ]
genes<-genes[ genes$biotype == "protein_coding" & genes$type=='gene', ]
genes<-renameSeqlevels(genes,paste0('chr',levels(genes@seqnames)))
genes<-tidyChromosomes(genes,keep.X=FALSE,keep.Y=FALSE) 

if (region=='gene-body'){
    all_genes<-genes
    genes<-reduce(genes)
    myDiff_CON<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d845d2dfc6_all.txt.bgz")
    myDiff_HE<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d856d610a0_all.txt.bgz")
    myDiff_AD<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d82c26a3d3_all.txt.bgz")
} else if (region=='promoters'){
    genes<-promoters(genes, upstream=1500,downstream=500)
    all_genes<-genes
    genes<-reduce(genes) #make sure no overlapping sequences are present
    myDiff_AD<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d818d60232_all.txt.bgz")
    myDiff_HE<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d878e8bfa2_all.txt.bgz")
    myDiff_CON<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d828057900_all.txt.bgz")
    # 
    # } else if (region=='both'){
    #   genes<-resize(genes,width(genes) + 1500L, fix = "end")
    #   all_genes<-genes
    #   genes<-reduce(genes) #make sure no overlapping sequences are present
    
} else if (region=='mirna'){
    genes<-readGFF('/home/cmengoni/methyl/hsa.gff3')
    genes<-GRanges(genes)
    genes<-tidyChromosomes(genes,keep.X=FALSE,keep.Y=FALSE)
    #keep<-unlist(lapply(genes$Name,function(x){strsplit(x,split = '-')[[1]][2]})=='miR')
    genes<-genes[genes$type=='miRNA',]
    all_genes<-genes
    
    myDiff_AD<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d861d75534_all.txt.bgz")
    myDiff_HE<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d86cfae5b6_all.txt.bgz")
    myDiff_CON<-readMethylDB("/home/cmengoni/methyl/methylDB/gene-body/methylDiff_670d8767a24bc_all.txt.bgz")
}


library(splicejam)
gr_AD<-GRanges(myDiff_AD)
gr_AD<-annotateGRfromGR(reduce(gr_AD),gr_AD)
tmp_genes<-keepSeqlevels(all_genes,unique(seqnames(gr_AD)),pruning.mode='coarse')
tmp_genes<-tmp_genes[,c('Name','gene_id')]
gr_AD <- annotateGRfromGR(GR1=gr_AD, GR2=tmp_genes,type = ifelse(region=='mirna', 'end', 'any'))##any')#type='end') #'any
gr_AD

gr_HE<-GRanges(myDiff_HE)
gr_HE<-annotateGRfromGR(reduce(gr_HE),gr_HE)
tmp_genes<-keepSeqlevels(all_genes,unique(seqnames(gr_HE)),pruning.mode='coarse')
tmp_genes<-tmp_genes[,c('Name','gene_id')]
gr_HE <- annotateGRfromGR(gr_HE, tmp_genes,type = ifelse(region=='mirna', 'equal', 'any'))#'equal')
gr_HE

gr_CON<-GRanges(myDiff_CON)
gr_CON<-annotateGRfromGR(reduce(gr_CON),gr_CON)
tmp_genes<-keepSeqlevels(all_genes,unique(seqnames(gr_CON)),pruning.mode='coarse')
tmp_genes<-tmp_genes[,c('Name','gene_id')]
gr_CON <- annotateGRfromGR(gr_CON, tmp_genes,type = ifelse(region=='mirna', 'equal', 'any'))#'equal')
gr_CON

# #save.image('/home/cmengoni/methyl/methylkit_data/gene-body.RData')
# ################################################################################
## Compare with RNASeq

library(openxlsx)

if(region=='mirna'){
    mrna_he_hc<-openxlsx::read.xlsx('methyl/DE_mirna.xlsx',sheet = 3)
    mrna_de_dc<-openxlsx::read.xlsx('methyl/DE_mirna.xlsx',sheet = 4)
    mrna_dc_hc<-openxlsx::read.xlsx('methyl/DE_mirna.xlsx',sheet = 1)
    
    mrna_dc_hc$SYMBOL<-mrna_dc_hc[,1]
    mrna_de_dc$SYMBOL<-mrna_de_dc[,1]
    mrna_he_hc$SYMBOL<-mrna_he_hc[,1]
} else {
    mrna_he_hc<-openxlsx::read.xlsx('methyl/de_mrna.xlsx',sheet=2)
    mrna_de_dc<-openxlsx::read.xlsx('methyl/de_mrna.xlsx',sheet=3)
    mrna_dc_hc<-openxlsx::read.xlsx('methyl/de_mrna.xlsx',sheet=5)
}

#png(paste0(paste('methyl/plots/de_dc',region,sep='_'),'.png'),width = 400, height = 400, units='mm', res = 300)
de_dc_common<-intersect(gr_AD$gene_id,mrna_de_dc$ENSEMBL) #29 #0 #25
mrna_de_dc<-mrna_de_dc[order(mrna_de_dc$ENSEMBL),]
df_AD<-data.frame(gr_AD)
rownames(df_AD)<-df_AD$gene_id
# plot(df_AD[de_dc_common[order(de_dc_common)],'meth.diff'],mrna_de_dc[mrna_de_dc$ENSEMBL%in%de_dc_common,'logFC'],
#      ylab='mRNA',xlab='methylation',main='AD EXP vs AD CON')
# text(df_AD[de_dc_common[order(de_dc_common)],'meth.diff'],
#      mrna_de_dc[mrna_de_dc$SYMBOL%in%de_dc_common,'logFC'],
#      labels=de_dc_common[order(de_dc_common)],cex= 0.7, pos=4)
# abline(v = 0, col="red", lwd=1, lty=1)
# abline(h = 0, col="red", lwd=1, lty=1)
# dev.off()


# png(paste0(paste('methyl/plots/he_hc',region,sep='_'),'.png'),width = 400, height = 400, units='mm', res = 300)
he_hc_common<-intersect(gr_HE$gene_id,mrna_he_hc$ENSEMBL) #29 #0 #25
mrna_he_hc<-mrna_he_hc[order(mrna_he_hc$ENSEMBL),]
df_HE<-data.frame(gr_HE)
rownames(df_HE)<-df_HE$gene_id
plot(df_HE[he_hc_common[order(he_hc_common)],'meth.diff'],mrna_he_hc[mrna_he_hc$ENSEMBL%in%he_hc_common,'logFC'],
     ylab='mRNA',xlab='methylation',main='HE EXP vs HE CON')
text(df_HE[he_hc_common[order(he_hc_common)],'meth.diff'],
     mrna_he_hc[mrna_he_hc$ENSEMBL%in%he_hc_common,'logFC'],
     labels=df_HE[he_hc_common[order(he_hc_common)],'Name'],cex= 0.7, pos=4)
abline(v = 0, col="red", lwd=1, lty=1)
abline(h = 0, col="red", lwd=1, lty=1)
# dev.off()

#png(paste0(paste('methyl/plots/dc_hc',region,sep='_'),'.png'),width = 400, height = 400, units='mm', res = 300)
dc_hc_common<-intersect(gr_CON$gene_id,mrna_dc_hc$ENSEMBL) #29 #0 #25
mrna_dc_hc<-mrna_dc_hc[order(mrna_dc_hc$ENSEMBL),]
df_CON<-data.frame(gr_CON)
rownames(df_CON)<-df_CON$gene_id
plot(df_CON[dc_hc_common[order(dc_hc_common)],'meth.diff'],mrna_dc_hc[mrna_dc_hc$ENSEMBL%in%dc_hc_common,'logFC'],
     ylab='mRNA',xlab='methylation',main='AD CON vs HE CON')
text(df_CON[dc_hc_common[order(dc_hc_common)],'meth.diff'],
     mrna_dc_hc[mrna_dc_hc$SYMBOL%in%dc_hc_common,'logFC'],
     labels=dc_hc_common[order(dc_hc_common)],cex= 0.7, pos=4)
abline(v = 0, col="red", lwd=1, lty=1)
abline(h = 0, col="red", lwd=1, lty=1)
# dev.off()


#################################################################################
## SAVE DE GENES

de_list<- list(merge(df_CON[,c('gene_id','Name','qvalue','meth.diff')],mrna_dc_hc[,c('ENSEMBL','logFC')],by.x='gene_id',by.y='ENSEMBL',all.x=T,all.y=F),
               merge(df_HE[,c('gene_id','Name','qvalue','meth.diff')],mrna_he_hc[,c('ENSEMBL','logFC')],by.x='gene_id',by.y='ENSEMBL',all.x=T,all.y=F),
               merge(df_AD[,c('gene_id','Name','qvalue','meth.diff')],mrna_de_dc[,c('ENSEMBL','logFC')],by.x='gene_id',by.y='ENSEMBL',all.x=T,all.y=F))
names(de_list) <- c("DisvsHeaForCon","ExpvsConForHea","ExpvsConForDis")

wb <- createWorkbook()
lapply(seq_along(de_list), function(i){
    addWorksheet(wb=wb, sheetName = names(de_list[i]))
    writeData(wb, sheet = i, de_list[[i]])
})
saveWorkbook(wb, paste0(paste0('/home/cmengoni/methyl/results/',region),'.xlsx'), overwrite = TRUE )
# #View(merge(df_HE[,c('Name','qvalue','meth.diff','type')],mrna_de_dc[,c('SYMBOL','logFC')],by.x='Name',by.y='SYMBOL',all.x=T,all.y=F))
