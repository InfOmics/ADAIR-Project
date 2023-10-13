# DIFFERENTIAL EXPRESSION ANALYSIS OF hOM CELLS (Claudia Mengoni)
# Data: (KK group) mRNA of healthy (HEA) and AD (DIS) hOM cells
# exposed to clean air (CON) and Diesel Exhaust (EXP) through ALI system.
# 48 samples of which 12 HEA.CON, 12 HEA.EXP, 12 DIS.CON, 12 DIS.EXP

# Aims:
# - detect biomarkers for acute air pollutant exposure in healthy and disease.
# - detect differences between healthy and disease

# 1) IMPORT DATA AND PREPARE OBJECTS

## load packages
library(edgeR)
library(RColorBrewer)
library(Homo.sapiens)
#library(proBatch)
#library(sva)
library(gplots)
#install.packages('/home/claudia/Downloads/SPARQL_1.16/SPARQL/',type='source',repos=NULL)
#install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)


setwd("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/")

## Read in metadata
meta<-read.delim("metadata.txt.csv",strip.white=TRUE,sep=',')

## Factorize
meta$REP <- as.factor(meta$REP)
meta$EXPOSURE <- as.factor(meta$EXPOSURE)
meta$CELL_LINE <- as.factor(meta$CELL_LINE)
meta$GENOTYPE <- as.factor(meta$GENOTYPE)
meta$GENOTYPE<-relevel(meta$GENOTYPE, 'HEA')
meta$GROUP<- as.factor(paste0(meta$GENOTYPE,".",meta$EXPOSURE))
meta$GROUP <- relevel(meta$GROUP,'HEA.CON')
meta$NAME<- paste0(meta$GROUP,".",meta$CELL_LINE,".",meta$REP)
meta$COLOR<-brewer.pal(9,'Set1')[meta$GROUP]
head(meta)

## Read in raw counts obtained from STAR + HTSeq-count
counts <- read.delim("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/featurecounts/featurecounts.txt", skip = 1)
rownames(counts)<-counts[,1]
counts <- counts[,7:length(counts)]
colnames(counts) <-meta$NAME

counts <- counts[ , colSums(is.na(counts)) < nrow(counts)]  # Remove rows with NA only


## Create DGE object
d0<-DGEList(counts)
d0$samples$REP <- meta$REP
d0$samples$group <- meta$GROUP
d0$samples$CELL_LINE <- meta$CELL_LINE

## Create an Ensembl to symbol to entrezID mapping dataframe and save it in the DGE object
geneid <- rownames(d0)
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),] #remove dups
genes<- genes[!is.na(genes$ENTREZID),] #remove NAs
genes<- genes[!is.na(genes$SYMBOL),]

dim(genes)
head(genes)
#d0$genes<-genes

# 2) FILTER AND NORMALIZE

## Show count per milion (cpm) before across samples normalization
lcpm <- cpm(d0, log=T)
boxplot(lcpm, las=2, col=meta$COLOR, main="Unnormalised data", ylab="Log-cpm")

## Filtering
dim(d0) #60664
keep.exprs <- rowSums(counts>15)>=12 # Empirical: Retain rows having at least 15 counts across 12 samples (cardinality of groups)
d0 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0) #13679

## TMM normalization
norm_counts <- calcNormFactors(d0, method = "TMM")
lcpm <- cpm(norm_counts, log=TRUE)
boxplot(lcpm, las=2, col=meta$COLOR, main="Normalised data", ylab="Log-cpm",cex.axis = 0.7)


# 3) DATA VISUALIZATION

plotMDS(d0,col =brewer.pal(9,'Set1')[meta$CELL_LINE],main='CELL-LINE') # strong separation!
plotMDS(d0,col =brewer.pal(9,'Set1')[meta$REP],main='REPLICATE') # some separation within cell-line
plotMDS(d0,col =brewer.pal(9,'Set1')[meta$GROUP],main='GROUP') # weak separation => use samples paired by cell-line and replicate

nb.cols <- 24
meta$PAIRED<-as.factor(paste0(meta$CELL_LINE,meta$REP))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)[meta$PAIRED]
plotMDS(d0,col =mycolors, main='PAIRED')

## Visualize other possible effects (e.g. smoking, hyposmia)
patients_metadata<-read.csv('/home/claudia/Desktop/ADAIR/ADAIR_exp/extra_metadata.csv',sep=',')
patients_metadata[patients_metadata$Smoking!='never','Smoking']<-'Y'
patients_metadata[patients_metadata$Smoking=='never','Smoking']<-'N'
patients_metadata[patients_metadata$Sense.of.smell=='hyposmia','Sense.of.smell']<-'HYP'
patients_metadata[patients_metadata$Sense.of.smell!='HYP','Sense.of.smell']<-'NOR'
patients_metadata<-patients_metadata[,c(1,7,8)]

new_meta<-merge(meta,patients_metadata,by.x='CELL_LINE',by.y='cell.line.code.UEM',)
new_meta<-new_meta[order(new_meta$SAMPLE,decreasing=F),]
new_meta$SAMPLE == meta$SAMPLE

new_meta$Sense.of.smell<-as.factor(new_meta$Sense.of.smell)
new_meta$Sense.of.smell<-relevel(new_meta$Sense.of.smell,ref = 'NOR')
new_meta$Smoking<-as.factor(new_meta$Smoking)
new_meta$Smoking<-relevel(new_meta$Smoking,ref = 'N')

mycols<-colorRampPalette(brewer.pal(3, "Set1"))(2)[new_meta$Smoking]
plotMDS(d0, col=mycols,labels = new_meta$NAME, main='SMOKERS') # => Smoking shows some level of separation, tested in the model
mycols<-colorRampPalette(brewer.pal(3, "Set1"))(2)[new_meta$Sense.of.smell]
plotMDS(d0, col=mycols, labels = new_meta$NAME)
# => Hyposmia is embedded within one condition (disease) so it is not useful to use as confounder between health and disease
# In disease EXP vs CON it is less powerful than pairing by cell-line and rep, so it is not used

# 4) MODEL
# 4.1) Paired comparison for within health and within disease comparisons

## design matrix
design<-model.matrix(~0+meta$PAIRED)
h.e <- meta$GENOTYPE =='HEA' & meta$EXPOSURE == 'EXP'
d.e <- meta$GENOTYPE =='DIS' & meta$EXPOSURE == 'EXP'
design <- cbind(design, h.e,d.e)

## HEA.EXP vs HEA.CON and DIS.EXP vs DIS.CON
v <- voomWithQualityWeights(norm_counts, design, normalize.method = 'cyclicloess')
fit <- lmFit(v,design)
fit_paired<-contrasts.fit(fit,coefficients = c(25,26))
fit_paired <- eBayes(fit_paired)
summary(decideTests(fit_paired,p.value=0.05))
hist(fit_paired$p.value)

# 4.2) Unpaired comparison for healthy vs disease in control condition (with correction for smoking)
design<-model.matrix(~new_meta$Smoking+meta$GROUP)
colnames(design)[3:5] <-  substr(colnames(design)[3:5], start = 11, stop=17)

## DIS.CON vs HEA.CON
v <- voomWithQualityWeights(norm_counts, design,normalize.method = 'cyclicloess')
fit <- lmFit(v,design)
fit_unpaired <- contrasts.fit(fit, coefficients=3)
fit_unpaired <- eBayes(fit_unpaired)
summary(decideTests(fit_unpaired,p.value=0.05))
hist(fit_unpaired$p.value)


# 5) PLOT RESULTS
h.e_h.c  = rownames(topTreat(fit_paired, coef=1, n=Inf,p.value=0.05,sort.by = 'logFC'))
d.e_d.c  = rownames(topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC'))
d.c_h.c  = rownames(topTreat(fit_unpaired, coef=1, n=Inf,p.value=0.05,sort.by = 'logFC'))

dt_paired <- decideTests(fit_paired, p.value=0.05)
plotMD(fit_paired, column=1, status=dt_paired[,1], main=colnames(fit_paired)[1], xlim=c(-8,13))
plotMD(fit_paired, column=2, status=dt_paired[,2], main=colnames(fit_paired)[2], xlim=c(-8,13))
vennDiagram(dt_paired[,1:2], include=c("up", "down"),
            counts.col=c("red", "blue"),names=c('HEA.EXP vs HEA.CON','DIS.EXP vs DIS.CON'))

dt_unpaired <- decideTests(fit_unpaired, p.value=0.05)
vennDiagram(dt_unpaired[,1], include=c("up", "down"),
            counts.col=c("red", "blue"), names = 'DIS.CON vs HEA.CON')

# Plot
he_hc_for_plot<-topTreat(fit_paired, coef=1, n=Inf,p.value=0.05,sort.by = 'logFC')
he_hc_for_plot$adj.P.Val<--log10(he_hc_for_plot$adj.P.Val)
he_hc_for_plot$names<-rownames(he_hc_for_plot)

de_dc_for_plot<-topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC')
de_dc_for_plot$adj.P.Val<--log10(de_dc_for_plot$adj.P.Val)
de_dc_for_plot$names<-rownames(de_dc_for_plot)

merge_for_plot<- merge(de_dc_for_plot[,c('names','adj.P.Val')],he_hc_for_plot[,c('names','adj.P.Val')],
                       by.x='names',by.y='names', all.x=FALSE,all.y=FALSE)

hea_col="#4DAF4A"
dis_col="#984EA3"


svg('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/common_DE.svg',width=20, height=26)
merge_for_plot<-merge_for_plot[order(merge_for_plot$adj.P.Val.x,decreasing=F),]
NAMES<-merge_for_plot$names
NAMES<-genes[genes$ENSEMBL%in%NAMES,'SYMBOL']
merge_for_plot$names<-NA
par(oma = c(8, 16, 0, 0))

barplot(t(merge_for_plot), names.arg= NAMES,
        horiz = T,xlab = '-log(p-value)',
        las=1,col=c(dis_col,hea_col),beside=T,cex.names = 1.5,cex.axis = 1.5)
abline(v=1.3,lty='dashed')
dev.off()
### HEATMAP OF DE GENES
library(Seurat)
ipa_palette<-CustomPalette(low = '#0000ea', high = '#ea6901', mid = 'white', k = 100)

#png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_hehc_top20_p.png",units='mm',width=250, height=200,res = 300)
#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_hehc_top20_p.svg",width=10, height=8)#units='mm',,res = 300)
h.e_h.c_10<-rownames(topTreat(fit_paired, coef=1, n=Inf,p.value=0.05,sort.by = 'P')[1:20,])
lcpm_de <- lcpm[h.e_h.c_10,grep('HEA',colnames(lcpm))]
#heatmap(lcpm_de,ColSideColors = )
exp_con <-brewer.pal(9,'Set1')[meta[grep('HEA',meta$NAME),'EXPOSURE']]
heatmap.2(lcpm_de, scale="row", xlab = 'Samples', ylab='Genes', main='Healthy',
          labCol='',labRow = '',#colnames(lcpm_de),
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con, key.xlab = 'LogCPM',
          col = ipa_palette)#hcl.colors(200,palette='PuOr',rev=T))
legend("bottomleft",legend=c('Exposed','Controls'),fill=unique(exp_con))
#dev.off()

# merge replicates
HEA.EXP.A<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.EXP.A.2','HEA.EXP.A.3','HEA.EXP.A.1')])
HEA.EXP.B<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.EXP.B.2','HEA.EXP.B.3','HEA.EXP.B.1')])
HEA.EXP.C<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.EXP.C.2','HEA.EXP.C.3','HEA.EXP.C.1')])
HEA.EXP.D<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.EXP.D.2','HEA.EXP.D.3','HEA.EXP.D.1')])

HEA.CON.A<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.A.2','HEA.CON.A.3','HEA.CON.A.1')])
HEA.CON.B<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.B.2','HEA.CON.B.3','HEA.CON.B.1')])
HEA.CON.C<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.C.2','HEA.CON.C.3','HEA.CON.C.1')])
HEA.CON.D<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.D.2','HEA.CON.D.3','HEA.CON.D.1')])

meta_sub<-meta[meta$REP==1,]
exp_con <-brewer.pal(9,'Set1')[meta_sub[grep('HEA',meta_sub$NAME),'EXPOSURE']]
png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_hehc_top20_p_REP.png",units='mm',width=250, height=200,res = 300)
lcpm_de_avg<-data.frame(HEA.EXP.A,HEA.EXP.B,HEA.EXP.C,HEA.EXP.D,HEA.CON.A,HEA.CON.B,HEA.CON.C,HEA.CON.D)
heatmap.2(as.matrix(lcpm_de_avg), scale="row", xlab = 'Samples', ylab='Genes', main='Healthy',
          labCol='',labRow = '',#colnames(lcpm_de),
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con, key.xlab = 'LogCPM',
          col = ipa_palette)
legend("bottomleft",legend=c('Exposed','Controls'),fill=unique(exp_con))
dev.off()


#png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_dedc_top20_p.png",units='mm',width=250, height=200,res = 300)
#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_dedc_top20_p.svg",width=10, height=8)#units='mm',,res = 300)#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_dedc_top20_logFC.svg",width=16, height=12)#,units="in",res=300)
d.e_d.c_10<-rownames(topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'P')[1:20,])
lcpm_de <- lcpm[d.e_d.c_10,grep('DIS',colnames(lcpm))]
#heatmap(lcpm_de,ColSideColors = )
exp_con <-brewer.pal(9,'Set1')[meta[grep('DIS',meta$NAME),'EXPOSURE']]
heatmap.2(lcpm_de, scale="row", xlab = 'Samples', ylab='Genes', main='AD',
          labCol='',labRow = '',
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con,key.xlab = 'LogCPM',
          col = ipa_palette)#hcl.colors(200,palette='PuOr',rev=T))
legend("bottomleft",legend=c('Exposed','Controls'),fill=unique(exp_con))
#dev.off()

#SUB BY REP
DIS.EXP.E<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.EXP.E.2','DIS.EXP.E.3','DIS.EXP.E.1')])
DIS.EXP.F<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.EXP.F.2','DIS.EXP.F.3','DIS.EXP.F.1')])
DIS.EXP.G<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.EXP.G.2','DIS.EXP.G.3','DIS.EXP.G.1')])
DIS.EXP.H<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.EXP.H.2','DIS.EXP.H.3','DIS.EXP.H.1')])

DIS.CON.E<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.E.2','DIS.CON.E.3','DIS.CON.E.1')])
DIS.CON.F<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.F.2','DIS.CON.F.3','DIS.CON.F.1')])
DIS.CON.G<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.G.2','DIS.CON.G.3','DIS.CON.G.1')])
DIS.CON.H<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.H.2','DIS.CON.H.3','DIS.CON.H.1')])

meta_sub<-meta[meta$REP==1,]
exp_con <-brewer.pal(9,'Set1')[meta_sub[grep('DIS',meta_sub$NAME),'EXPOSURE']]
png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_DeDc_top20_p_REP.png",units='mm',width=250, height=200,res = 300)
lcpm_de_avg<-data.frame(DIS.EXP.E,DIS.EXP.F,DIS.EXP.G,DIS.EXP.H,DIS.CON.E,DIS.CON.F,DIS.CON.G,DIS.CON.H)
heatmap.2(as.matrix(lcpm_de_avg), scale="row", xlab = 'Samples', ylab='Genes', main='AD',
          labCol='',labRow = '',#colnames(lcpm_de),
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con, key.xlab = 'LogCPM',
          col = ipa_palette)
legend("bottomleft",legend=c('Exposed','Controls'),fill=unique(exp_con))
dev.off()


#png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_dchc_top20_p.png",units='mm',width=250, height=200,res = 300)
#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_dchc_top20_p.svg",width=10, height=8)#,units="in",res=300)
d.c_h.c_10<-rownames(topTreat(fit_unpaired, coef=1, n=Inf,p.value=0.05,sort.by = 'p')[1:20,])
lcpm_de <- lcpm[d.c_h.c_10,grep('CON',colnames(lcpm))]
#heatmap(lcpm_de,ColSideColors = )
exp_con <-brewer.pal(9,'Set1')[3:9][meta[grep('CON',meta$NAME),'GENOTYPE']]
heatmap.2(lcpm_de, scale="row", xlab = 'Samples', ylab='Genes', main='Controls',
          labCol='',labRow = '',
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con,key.xlab = 'LogCPM',
          col = ipa_palette)
legend("bottomleft",legend=c('Healthy','AD'),fill=unique(exp_con))
#dev.off()

#SUB BY REP
HEA.CON.A<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.A.1','HEA.CON.A.2','HEA.CON.A.3')])
HEA.CON.B<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.B.1','HEA.CON.B.2','HEA.CON.B.3')])
HEA.CON.C<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.C.1','HEA.CON.C.2','HEA.CON.C.3')])
HEA.CON.D<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('HEA.CON.D.1','HEA.CON.D.2','HEA.CON.D.3')])

DIS.CON.E<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.E.2','DIS.CON.E.3','DIS.CON.E.1')])
DIS.CON.F<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.F.2','DIS.CON.F.3','DIS.CON.F.1')])
DIS.CON.G<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.G.2','DIS.CON.G.3','DIS.CON.G.1')])
DIS.CON.H<-rowMeans(lcpm_de[,colnames(lcpm_de)%in%c('DIS.CON.H.2','DIS.CON.H.3','DIS.CON.H.1')])

meta_sub<-meta[meta$REP==1,]
exp_con <-brewer.pal(9,'Set1')[3:9][meta_sub[grep('CON',meta_sub$NAME),'GENOTYPE']]
png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/heatmap_HcDc_top20_p_REP.png",units='mm',width=250, height=200,res = 300)
lcpm_de_avg<-data.frame(HEA.CON.A,HEA.CON.B,HEA.CON.C,HEA.CON.D,
                        DIS.CON.E,DIS.CON.F,DIS.CON.G,DIS.CON.H)
heatmap.2(as.matrix(lcpm_de_avg), scale="row", xlab = 'Samples', ylab='Genes', main='Controls',
          labCol='',labRow = '',#colnames(lcpm_de),
          margin=c(8,6), lhei=c(2,10),
          trace="none", density.info="none",
          dendrogram="column", ColSideColors = exp_con, key.xlab = 'LogCPM',
          col = ipa_palette)
legend("bottomleft",legend=c('Healthy','AD'),fill=unique(exp_con))
dev.off()

##########
# SAVE DE GENES
# de_list<- list(topTreat(fit_paired, coef=1, n=Inf, p.value=0.05, sort.by = 'logFC'),
#                topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC'),
#                topTreat(fit_unpaired, coef=1, n=Inf,p.value=0.05,sort.by = 'logFC'))
# names(de_list) <- c("ExpvsConForHea","ExpvsConForDis","DisvsHeaForCon")
#
# wb <- createWorkbook()
# lapply(seq_along(de_list), function(i){
#   addWorksheet(wb=wb, sheetName = names(de_list[i]))
#   writeData(wb, sheet = i, de_list[[i]])
# })
#saveWorkbook(wb, "out/DE.xlsx", overwrite = TRUE )

# 6) ENRICHMENT ANALYSIS

# ## GSEA
# foldchanges <- topTreat(fit_paired, coef=2, n=Inf)$logFC
# names(foldchanges) <- topTreat(fit_paired, coef=1, n=Inf)$ENTREZ
# comp_foldchanges <- sort(foldchanges, decreasing = TRUE)
#
# GSEA using gene sets from KEGG pathways
# set.seed(1)
# gseaKEGG <- gseKEGG(geneList =  comp_foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
#                     organism = "hsa", # supported organisms listed below
#                     minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
#                     pvalueCutoff = 0.05, # padj cutoff value
#                     verbose = FALSE)
#
# gseaKEGG@result


## Create background dataset for hypergeometric testing using all genes tested for significance in the results
allOE_genes <- read.csv('/home/claudia/Desktop/ADAIR/ADAIR_exp/genes.csv',sep='\t',col.names = F)
allOE_genes<-allOE_genes$FALSE.

## Run enrichment analysis on HEA.EXP vs HEA.CON
sig_OE_genes<-h.e_h.c

### GO
go_h.e_h.c <- enrichGO(gene = sig_OE_genes,
                       universe = allOE_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)

dotplot(go_h.e_h.c,font.size=12,showCategory=15)
go_h.e_h.c_df<-data.frame(go_h.e_h.c)
# write.csv(go_d.c_h.c_df,'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/mrna.csv')

### KEGG
kegg_h.e_h.c <- enrichKEGG(gene = genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'],
                        universe = genes[genes$ENSEMBL %in% allOE_genes,'ENTREZID'],
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

dotplot(kegg_h.e_h.c, showCategory=30,font.size=12)
kegg_h.e_h.c_df<-data.frame(kegg_h.e_h.c)
# write.csv(kegg_h.e_h.c_df,'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/08022022/kegg_mrna.csv')

### WIKIPATH
wikipath_h.e_h.c<-enrichWP(genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], organism = "Homo sapiens")
dotplot(wikipath_h.e_h.c, showCategory=30,font.size=12)

### REACTOME
reactome_h.e_h.c<-enrichPathway(gene=genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], pvalueCutoff = 0.05, readable=TRUE)

## Run enrichment analysis on DIS.EXP vs DIS.CON
sig_OE_genes<-d.e_d.c

### GO
go_d.e_d.c <- enrichGO(gene = sig_OE_genes,
                    universe = allOE_genes,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)

dotplot(go_d.e_d.c, showCategory=15, font.size=10)
go_d.e_d.c_df<-data.frame(go_d.e_d.c)
#write.csv(go_d.e_d.c_df,'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/08022022/p.csv')

### KEGG
ekegg_d.e_d.c <- enrichKEGG(gene = genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'],
                        universe = genes[genes$ENSEMBL %in% allOE_genes,'ENTREZID'],
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

dotplot(ekegg_d.e_d.c, showCategory=30)
ekegg_d.e_d.c_df<-data.frame(ekegg_d.e_d.c)
#write.csv(ekegg_d.e_d.c_df,'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/08022022/p.csv')

### WikiPath
wikipath_d.e_d.c<-enrichWP(genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], organism = "Homo sapiens")
dotplot(wikipath_d.e_d.c, showCategory=30,font.size=12)

### Reactome
reactome_d.e_d.c<-enrichPathway(gene=genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], pvalueCutoff = 0.05, readable=TRUE)

## Run enrichment analysis on DIS.EXP vs DIS.CON - HEA.EXP vs HEA.CON
# sig_OE_genes<-setdiff(d.e_d.c,h.e_h.c)
# #up
# #sig_OE_genes<-rownames(topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC')[topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC')$ENSEMBL %in% sig_OE_genes & topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC')$logFC<0,])
# ego_diff <- enrichGO(gene = sig_OE_genes,
#                      universe = allOE_genes,
#                      keyType = "ENSEMBL",
#                      OrgDb = org.Hs.eg.db,
#                      ont = "BP",
#                      pAdjustMethod = "BH",
#                      qvalueCutoff = 0.05,
#                      readable = TRUE)
#
# dotplot(ego_diff, showCategory=30)
# ego_diff_df<-data.frame(ego_diff)
#
# ekegg_diff <- enrichKEGG(gene = genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'],
#                            universe = genes[genes$ENSEMBL %in% allOE_genes,'ENTREZID'],
#                            pAdjustMethod = "BH",
#                            qvalueCutoff = 0.05)
#
# dotplot(ekegg_diff, showCategory=30)
# ekegg_diff_df<-data.frame(ekegg_diff)

## Run enrichment analysis on DIS.CON vs HEA.CON
sig_OE_genes<-d.c_h.c

### GO
go_d.c_h.c <- enrichGO(gene = sig_OE_genes,
                       universe = allOE_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
go_d.c_h.c_df<-data.frame(go_d.c_h.c)

dotplot(go_d.c_h.c, showCategory=25, font.size=10)
#write.csv(data.frame(go_d.c_h.c),'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/mrna.csv')

### KEGG
kegg_d.c_h.c <- enrichKEGG(gene = genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'],
                           universe = genes[genes$ENSEMBL %in% allOE_genes,'ENTREZID'],
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)

dotplot(kegg_d.c_h.c, showCategory=30)

### WikiPath
wikipath_d.c_h.c<-enrichWP(genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], organism = "Homo sapiens")
dotplot(wikipath_d.c_h.c, showCategory=30,font.size=12)

### reactome
reactome_d.c_h.c<-enrichPathway(gene=genes[genes$ENSEMBL %in% sig_OE_genes,'ENTREZID'], pvalueCutoff = 0.05, readable=TRUE)


# dis_enrich_group4<-disease_enrichment(entities =genes[genes$ENSEMBL %in% d.c_h.c,'ENTREZID'], vocabulary = "ENTREZ", database = "CURATED" )
# ext.dis_enrich_group4 <- extract(dis_enrich_group4)
# ext.dis_enrich_group4 <- ext.dis_enrich_group4[ext.dis_enrich_group4$FDR<0.05,]
# plot(dis_enrich_group4, class = "Enrichment", count =3,  cutoff= 0.05, nchars=70)


# # SAVE GO AND KEGG
# wb <- createWorkbook()
# lapply(seq_along(go_list), function(i){
#   addWorksheet(wb=wb, sheetName = names(go_list[i]))
#   writeData(wb, sheet = i, go_list[[i]],rowNames = TRUE)
# })
# saveWorkbook(wb, "out/GO.xlsx", overwrite = TRUE )
#
# wb <- createWorkbook()
# lapply(seq_along(kegg_list), function(i){
#   addWorksheet(wb=wb, sheetName = names(kegg_list[i]))
#   writeData(wb, sheet = i, kegg_list[[i]],rowNames = TRUE)
# })
# saveWorkbook(wb, "out/KEGG.xlsx", overwrite = TRUE )


## Retrieve disease associated to target genes
disgenet_api_key <- get_disgenet_api_key(
  email = "cla.mengoni1@gmail.com",
  password = "DisGeNet1" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

key='Alzheimer' #smell
dis_d.c_h.c <- gene2disease( gene = genes[genes$ENSEMBL %in% d.c_h.c,'SYMBOL'], verbose = TRUE,database = 'ALL')
dis_d.c_h.c<-extract(dis_d.c_h.c)
dis_d.c_h.c<-dis_d.c_h.c[grep('Alzheimer',dis_d.c_h.c$disease_name),]
print(unique(dis_d.c_h.c$gene_symbol))
#"Alzheimer" 70 and NONE (ALL, CURATED)
#"smell" 3 and NONE (ALL, CURATED)   sense of smell impaired SEMA3A,SPRY4,SEMA7A

dis_h.e_h.c <- gene2disease( gene = genes[genes$ENSEMBL %in% h.e_h.c,'SYMBOL'], verbose = TRUE, database = 'ALL')
dis_h.e_h.c<-extract(dis_h.e_h.c)
dis_h.e_h.c<-dis_h.e_h.c[grep('Alzheimer',dis_h.e_h.c$disease_name),]
print(unique(dis_h.e_h.c$gene_symbol))#42 + 3
#"Alzheimer" 42 and 3 (ALL, CURATED)
#"smell,hyposmia,anosmia" 1 and NONE  (ALL, CURATED)   sense of smell impaired UCHL1

dis_d.e_d.c <- gene2disease( gene =  genes[genes$ENSEMBL %in% d.e_d.c,'SYMBOL'], verbose = TRUE, database = 'ALL')
dis_d.e_d.c <-extract(dis_d.e_d.c)
dis_d.e_d.c<-dis_d.e_d.c[grep('Alzheimer',dis_d.e_d.c$disease_name),]
print(unique(dis_d.e_d.c$gene_symbol)) #123 + 6
#"Alzheimer" 123 and 6 (ALL, CURATED)
#"smell,hyposmia,anosmia" none and NONE (ALL, CURATED)

### enrichment
dis_h.e_h.c_ENRICH <-disease_enrichment( entities = genes[genes$ENSEMBL %in% h.e_h.c,'ENTREZID'] , vocabulary = 'ENTREZ',database = "CURATED")
dis_h.e_h.c_ENRICH <-extract(dis_h.e_h.c_ENRICH)
dis_h.e_h.c_ENRICH<-dis_h.e_h.c_ENRICH[grep('Alzheimer',dis_h.e_h.c_ENRICH$Description),]
print(dis_h.e_h.c_ENRICH$shared_symbol) # 3 (PLAU;VEGFA;HMOX1)

dis_d.e_d.c_ENRICH <-disease_enrichment( entities = genes[genes$ENSEMBL %in% d.e_d.c,'ENTREZID'] , vocabulary = 'ENTREZ',database = "CURATED")
dis_d.e_d.c_ENRICH <-extract(dis_d.e_d.c_ENRICH)
dis_d.e_d.c_ENRICH<-dis_d.e_d.c_ENRICH[grep('Alzheimer',dis_d.e_d.c_ENRICH$Description),]
print(dis_d.e_d.c_ENRICH$shared_symbol) # 6 (HMOX1;CST3;IGF2R;PRNP;NECTIN2;TOMM40)

dis_d.c_h.c_ENRICH <-disease_enrichment( entities = genes[genes$ENSEMBL %in% d.c_h.c,'ENTREZID'] , vocabulary = 'ENTREZ',database = "CURATED")
dis_d.c_h.c_ENRICH <-extract(dis_d.c_h.c_ENRICH)
dis_d.c_h.c_ENRICH<-dis_d.c_h.c_ENRICH[grep('Alzheimer',dis_d.c_h.c_ENRICH$Description),]
print(dis_d.c_h.c_ENRICH$shared_symbol) # 0


## GO OX genes heatmap
sig_OE_genes = topTreat(fit_paired, coef=1, n=Inf,p.value=0.05,sort.by = 'logFC')
sig_OE_genes<-rownames(sig_OE_genes[sig_OE_genes$logFC>0,])
go_h.e_h.c <- enrichGO(gene = sig_OE_genes,
                       universe = allOE_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE, )

go_h.e_h.c_df<-data.frame(go_h.e_h.c)
go_h.e_h.c_df<-go_h.e_h.c_df[grep('oxidative',go_h.e_h.c_df$Description),]
go_OX_h.e_h.c<-unique(c(strsplit(go_h.e_h.c_df$geneID[[1]],split='/')[[1]],strsplit(go_h.e_h.c_df$geneID[[2]],split='/')[[1]]))
print(go_OX_h.e_h.c)

sig_OE_genes  = topTreat(fit_paired, coef=2, n=Inf,p.value=0.05,sort.by = 'logFC')
sig_OE_genes<-rownames(sig_OE_genes[sig_OE_genes$logFC>0,])
go_d.e_d.c <- enrichGO(gene = sig_OE_genes,
                       universe = allOE_genes,
                       keyType = "ENSEMBL",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE, )

go_d.e_d.c_df<-data.frame(go_d.e_d.c)
go_d.e_d.c_df<-go_d.e_d.c_df[grep('oxidative',go_d.e_d.c_df$Description),]
go_OX_d.e_d.c<-unique(c(strsplit(go_d.e_d.c_df$geneID[[1]],split='/')[[1]],strsplit(go_d.e_d.c_df$geneID[[2]],split='/')[[1]]))
print(unique(go_OX_d.e_d.c))


all_OX<-c(go_OX_d.e_d.c,go_OX_h.e_h.c)
all_OX<-unique(all_OX)

#h.e_h.c_df<-h.e_h.c_df[,-c(length(h.e_h.c_df),length(h.e_h.c_df)-1)]
h.e_h.c_df<-topTreat(fit_paired, coef=1, n=Inf,p.value=1,sort.by = 'logFC')
h.e_h.c_df$ENSEMBL<-rownames(h.e_h.c_df)
h.e_h.c_df<-merge(h.e_h.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
h.e_h.c_df_plot_OX <- h.e_h.c_df[h.e_h.c_df$SYMBOL %in% all_OX,]#$logFC
h.e_h.c_df_plot_OX$sign <-  ifelse(h.e_h.c_df_plot_OX$adj.P.Val < 0.05, T,F)

 d.e_d.c_df<-topTreat(fit_paired, coef=2, n=Inf,p.value=1,sort.by = 'logFC')
d.e_d.c_df$ENSEMBL<-rownames(d.e_d.c_df)
d.e_d.c_df$id  <- 1:nrow(d.e_d.c_df)
d.e_d.c_df<-merge(d.e_d.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
d.e_d.c_df_plot_OX <- d.e_d.c_df[d.e_d.c_df$SYMBOL %in% all_OX,]#$logFC
d.e_d.c_df_plot_OX$sign <-  ifelse(d.e_d.c_df_plot_OX$adj.P.Val < 0.05, T,F)


plot_OX<-merge(d.e_d.c_df_plot_OX,h.e_h.c_df_plot_OX,by='SYMBOL')
nms<-plot_OX$SYMBOL
rownames(plot_OX) <- nms
plot_OX<-plot_OX[order(plot_OX$id), ]
plot_OX<-plot_OX[,c("logFC.x","logFC.y")]
colnames(plot_OX)<- c('AD Exposed','Healthy Exposed')



library(gplots)
library(RColorBrewer)
library(ggplot2)
#cust_col<-CustomPalette(low = '#0000FF', high = '#FF6600', mid = 'white', k = 200) #library(Seurat)

#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/ox_genes_up.svg",width=4, height=13)#,units="in",res=300) #height 16 all
#png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/ox_genes_up.png",units='mm',width=50, height=250,res = 300)
par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))
#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/ox_genes_all.svg",width=10, height=8)#,units="in",res=300)
heatmap.2(as.matrix(plot_OX[all_OX,]), Rowv = F, Colv = F, trace = "none", na.color = "Grey",density.info="none",
          scale="none", dendrogram = "none",col = ipa_palette, key = F, cexCol = 1,cexRow=1)#,main='Oxidative stress related genes')
#hcl.colors(200,palette='PuOr',rev=T))
#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)#
#dev.off()


##################################################################
# enrichment of overlap
library(xopen)
length(intersect(d.e_d.c,h.e_h.c))
sig_OE_genes  = intersect(d.e_d.c,h.e_h.c)
#write.table(sig_OE_genes_g,'/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/common_genes.tsv',quote=F, row.names = F,sep='\t')
common_pathways<-readxl::read_xls('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/IPA_common/Claudia_pathways.xls',col_names = T,skip = 1)
common_pathways <- common_pathways[common_pathways$`-log(p-value)`>1.3,]
common_pathways<-common_pathways[order(common_pathways$`-log(p-value)`,decreasing=F),]
svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/pathways_common_scale.svg",width=10, height=8)#,units="in",res=300)
par(oma = c(4, 32, 0, 4) , mar=c(0,0,1,4))
barplot(common_pathways$`-log(p-value)`, names.arg= common_pathways$`Ingenuity Canonical Pathways`,
        horiz = T,xlab = '-log(p-value)',ylab='Canonical Pathways',las=1)#, main='Shared AD, Control genes in Exposed')
abline(v=1.3,lty='dashed')
dev.off()

common_upst<-readxl::read_xls('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/IPA_common/Claudia_upstream.xls',col_names = T,skip = 1)
common_upst <- common_upst[common_upst$`p-value of overlap`<0.05,]
common_upst$`-log(p-value)` <- -log(common_upst$`p-value of overlap`,10)
common_upst<-common_upst[order(common_upst$`-log(p-value)`,decreasing=F),]
# common_upst<-common_upst[-(grep('mir-',common_upst$`Upstream Regulator`)),]
# common_upst<-common_upst[-(grep('hemoglobin',common_upst$`Upstream Regulator`)),]
common_upst<-common_upst[common_upst$`Molecule Type`=='transcription regulator',]
svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/upst_common.svg",width=10, height=8)#,units="in",res=300)
par(oma = c(2, 10, 0, 4) , mar=c(0,1,1,4))
barplot(common_upst$`-log(p-value)`, names.arg= common_upst$`Upstream Regulator`,
        horiz = T,xlab = '-log(p-value)',ylab='Transcription Factors',las=1,  xaxp = c(0,8, 8))#,cex.names = 0.6)#, main='Shared AD, Control genes in Exposed')
abline(v=1.3,lty='dashed')
dev.off()
##################################################################
# VOLCANO PLOTS
library(ggplot2)
library(ggrepel)
##H.E_H.C
h.e_h.c_ALZ<- topTreat(fit_paired, coef=1, n=Inf,p.value=1,sort.by = 'logFC')
h.e_h.c_ALZ$ENSEMBL<-rownames(h.e_h.c_ALZ)
SYMBOLS<-genes[genes$ENSEMBL %in% rownames(h.e_h.c_ALZ),c('ENSEMBL','SYMBOL')]
h.e_h.c_ALZ<-merge(h.e_h.c_ALZ,SYMBOLS,by='ENSEMBL')
#write.table(h.e_h.c_ALZ,"/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/14022022/res/h.e_h.c_unfiltered.tsv")

h.e_h.c_ALZ$diffexpressed <- "Not significant"
h.e_h.c_ALZ$diffexpressed[h.e_h.c_ALZ$logFC > 0 & h.e_h.c_ALZ$adj.P.Val  < 0.05] <- "Up-regulated"
h.e_h.c_ALZ$diffexpressed[h.e_h.c_ALZ$logFC < 0 & h.e_h.c_ALZ$adj.P.Val  < 0.05] <- "Down-regulated"

lit_search<-c('HSPA6','HMOX1')
AD_enr<-unique(unlist(strsplit(dis_h.e_h.c_ENRICH$shared_symbol,';')))#c('PLAU','VEGFA','HMOX1')
GO_enr<-unique(go_OX_h.e_h.c)


words<-AD_enr

h.e_h.c_ALZ[!(h.e_h.c_ALZ$SYMBOL %in% words), 'SYMBOL'] <- F

#png("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/14022022/res/LIT_genes_H.E_H.C.png",width=1024, height=800,res = 100)
#tiff("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/AD_genes_H.E_H.C.tif",width=10, height=8, units="in",res=300)
ggplot(data=h.e_h.c_ALZ, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=SYMBOL))+
  geom_point(size = 1.5) +
  theme_minimal() +
  geom_text_repel(aes(logFC, -log10(adj.P.Val)),
                  label = ifelse(h.e_h.c_ALZ$SYMBOL != F,
                                 as.character(h.e_h.c_ALZ$SYMBOL),""),
                  box.padding = unit(.7, "lines"),hjust= 0.30,col='black', max.overlaps = 1000) +

  xlim(min(h.e_h.c_ALZ$logFC)-0.5,max(h.e_h.c_ALZ$logFC)+0.5)+
  ylim(min(-log10(h.e_h.c_ALZ$adj.P.Val))-0.5,max(-log10(h.e_h.c_ALZ$adj.P.Val))+0.5)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_hline(yintercept=-log10(0.05), col="red",linetype='dotted',) +
  labs(title='healthy exp vs con', y='-log10(p-adj)', x='log2FC')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + theme(legend.title = element_blank())
#dev.off()

# D.E_D.C
d.e_d.c_ALZ<- topTreat(fit_paired, coef=2, n=Inf,p.value=1,sort.by = 'logFC')
d.e_d.c_ALZ$ENSEMBL<-rownames(d.e_d.c_ALZ)
SYMBOLS<-genes[genes$ENSEMBL %in% rownames(d.e_d.c_ALZ),c('ENSEMBL','SYMBOL')]
d.e_d.c_ALZ<-merge(d.e_d.c_ALZ,SYMBOLS,by='ENSEMBL')
#write.table(d.e_d.c_ALZ,"/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/14022022/res/d.e_d.c_unfiltered.tsv")

d.e_d.c_ALZ$diffexpressed <- "Not significant"
d.e_d.c_ALZ$diffexpressed[d.e_d.c_ALZ$logFC > 0 & d.e_d.c_ALZ$adj.P.Val  < 0.05] <- "Up-regulated"
d.e_d.c_ALZ$diffexpressed[d.e_d.c_ALZ$logFC < 0 & d.e_d.c_ALZ$adj.P.Val  < 0.05] <- "Down-regulated"


lit_search<-c('HSPA6','HMOX1','FOSB','NR4A3')
AD_enr<-unique(unlist(strsplit(dis_d.e_d.c_ENRICH$shared_symbol,';')))
GO_enr<-unique(go_OX_d.e_d.c)

words<-AD_enr

d.e_d.c_ALZ[!(d.e_d.c_ALZ$SYMBOL %in% words), 'SYMBOL'] <- F

#svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/AD_genes_D.E_D.C.svg",width=10, height=8)#, units="in", res=300)

ggplot(data=d.e_d.c_ALZ, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=SYMBOL))+
  geom_point(size = 1.5) +
  theme_minimal() +
  geom_text_repel(aes(logFC, -log10(adj.P.Val)),
                  label = ifelse(d.e_d.c_ALZ$SYMBOL != F,
                                 as.character(d.e_d.c_ALZ$SYMBOL),""),
                  box.padding = unit(.7, "lines"),hjust= 0.30,col='black',max.overlaps = 9000) +
  xlim(min(d.e_d.c_ALZ$logFC)-0.5,max(d.e_d.c_ALZ$logFC)+0.5)+
  ylim(min(-log10(d.e_d.c_ALZ$adj.P.Val))-0.5,max(-log10(d.e_d.c_ALZ$adj.P.Val))+0.5)+
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_hline(yintercept=-log10(0.05), col="red",linetype='dotted',) +
  labs(title='AD exp vs con', y='-log10(p-adj)', x='log2FC')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + theme(legend.title = element_blank())
#dev.off()


##############
# NFE2L2 Upstream
h.e_h.c_NFE2L2<-c('FTH1','FTL','GCLC','GCLM','HMOX1','NQO1','OSGIN1','PRDX1')
d.e_d.c_NFE2L2<-c('AKR1B10','CHGB','GCLM','HMOX1','NQO1','NUCB2','SCARB1','TGFB1')
all_NFE2L2<-c(d.e_d.c_NFE2L2,h.e_h.c_NFE2L2)
length(all_NFE2L2)
all_NFE2L2<-unique(all_NFE2L2)
length(all_NFE2L2)

h.e_h.c_df<-topTreat(fit_paired, coef=1, n=Inf,p.value=1,sort.by = 'logFC')
h.e_h.c_df$ENSEMBL<-rownames(h.e_h.c_df)
h.e_h.c_df<-merge(h.e_h.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
h.e_h.c_df_plot_NFE2L2 <- h.e_h.c_df[h.e_h.c_df$SYMBOL %in% all_NFE2L2,]#$logFC

d.e_d.c_df<-topTreat(fit_paired, coef=2, n=Inf,p.value=1,sort.by = 'logFC')
d.e_d.c_df$ENSEMBL<-rownames(d.e_d.c_df)
d.e_d.c_df$id  <- 1:nrow(d.e_d.c_df)
d.e_d.c_df<-merge(d.e_d.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
d.e_d.c_df_plot_NFE2L2 <- d.e_d.c_df[d.e_d.c_df$SYMBOL %in% all_NFE2L2,]#$logFC

plot_NFE2L2<-merge(d.e_d.c_df_plot_NFE2L2,h.e_h.c_df_plot_NFE2L2,by='SYMBOL')
nms<-plot_NFE2L2$SYMBOL
rownames(plot_NFE2L2) <- nms
plot_NFE2L2<-plot_NFE2L2[,c("logFC.x","logFC.y")]
colnames(plot_NFE2L2)<- c('DIS.EXP','HEA.EXP')
plot_NFE2L2<-plot_NFE2L2[order(abs(plot_NFE2L2$DIS.EXP),decreasing = T), ]
colnames(plot_NFE2L2)<- c('AD Exposed','Healthy Exposed')


svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/NFE2L2_genes.svg",width=4, height=8)
par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))

heatmap.2(as.matrix(plot_NFE2L2), Rowv = F, Colv = F, trace = "none", na.color = "Grey", cexCol = 1,cexRow=1,density.info="none",
          scale="none", dendrogram = "none",col =ipa_palette,key = F)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)
dev.off()


##############
# NRF2 CANONICAL
h.e_h.c_NRF2<-c('BACH1','DNAJB1','DNAJB4','DNAJB9','FTH1', 'FTL','GCLC','GCLM','HMOX1','MAFG','MGST1','NQO1','PRDX1','SQSTM1','TXN','TXNRD1')
d.e_d.c_NRF2<-c('CLPP','DNAJA4','DNAJB11','DNAJC3','FTL','GCLM','HMOX1','HSP90AA1','HSP90B1','JUND','MAFF','MAP2K2','NQO1','PIK3R6','SCARB1','TXNRD1')
all_NRF2<-c(d.e_d.c_NRF2,h.e_h.c_NRF2)
length(all_NRF2)
all_NRF2<-unique(all_NRF2)
length(all_NRF2)

#h.e_h.c_df<-h.e_h.c_df[,-c(length(h.e_h.c_df),length(h.e_h.c_df)-1)]
h.e_h.c_df<-topTreat(fit_paired, coef=1, n=Inf,p.value=1,sort.by = 'P')
h.e_h.c_df$ENSEMBL<-rownames(h.e_h.c_df)
h.e_h.c_df<-merge(h.e_h.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
h.e_h.c_df_plot_NRF2 <- h.e_h.c_df[h.e_h.c_df$SYMBOL %in% all_NRF2,]#$logFC


#d.e_d.c_df<-d.e_d.c_df[,-c(length(d.e_d.c_df),length(d.e_d.c_df)-1)] #questo si puo togliere quando tolgo il volcano plot
d.e_d.c_df<-topTreat(fit_paired, coef=2, n=Inf,p.value=1,sort.by = 'P')
d.e_d.c_df$ENSEMBL<-rownames(d.e_d.c_df)
d.e_d.c_df<-merge(d.e_d.c_df,genes,by.x="ENSEMBL",by.y="ENSEMBL")
d.e_d.c_df_plot_NRF2 <- d.e_d.c_df[d.e_d.c_df$SYMBOL %in% all_NRF2,]#$logFC

plot_NRF2<-merge(d.e_d.c_df_plot_NRF2,h.e_h.c_df_plot_NRF2,by='SYMBOL')
nms<-plot_NRF2$SYMBOL
rownames(plot_NRF2) <- nms
plot_NRF2<-plot_NRF2[,c("logFC.x","logFC.y")]
colnames(plot_NRF2)<- c('DIS.EXP','HEA.EXP')
plot_NRF2<-plot_NRF2[order(abs(plot_NRF2$DIS.EXP),decreasing = T),]
colnames(plot_NRF2)<- c('AD Exposed','Healthy Exposed')

svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/NRF2_canonical.svg",width=3, height=8)
par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))
heatmap.2(as.matrix(plot_NRF2), Rowv = all_NRF2, Colv = F, trace = "none", na.color = "Grey", cexCol = 0.7,cexRow=0.7,density.info="none",
          scale="none", dendrogram = "none",col = ipa_palette,key = F)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)

dev.off()

# test simone's method
blues <- ipa_palette[1:47]#50]
oranges <- ipa_palette[54:100]
my.breaks1 <- c(seq(-2, -0.75, by=0.25),  seq(0.75, 2, by=0.25))
my.breaks2 <- c(seq(-0.75, 0, by=0.1), seq(0, 0.75, by=0.1))
#my.breaks3 <- seq(1, 2, by=1)
my.colors <- c(colorRampPalette(colors = c("#0000EA", "#0000EA"))((length(my.breaks1)/2)),
               colorRampPalette(blues)(length(my.breaks1)),
               colorRampPalette(colors=c('#FFFFFF'))(1),
               colorRampPalette(oranges)(length(my.breaks2)-4),
               colorRampPalette(colors = c("#EA6901", "#EA6901" ))(length(my.breaks1)-5)) #aumenta numero x diminuire arancioni
              # colorRampPalette(colors = c("#EA6901"))(length(my.breaks3)))
par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))
heatmap.2(as.matrix(plot_NRF2), Rowv = all_NRF2, Colv = F, trace = "none", na.color = "Grey", cexCol = 0.7,cexRow=0.7,density.info="none",
          scale="none", dendrogram = "none",col = my.colors,key = T)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)

#subset by logfc
# plot_NRF2<-merge(d.e_d.c_df_plot_NRF2,h.e_h.c_df_plot_NRF2,by='SYMBOL')
# nms<-plot_NRF2$SYMBOL
# rownames(plot_NRF2) <- nms
# plot_NRF2<-plot_NRF2[,c("logFC.x","logFC.y")]
# colnames(plot_NRF2)<- c('DIS.EXP','HEA.EXP')
# plot_NRF2<-plot_NRF2[order(abs(plot_NRF2$DIS.EXP),decreasing = T),]
# colnames(plot_NRF2)<- c('AD Exposed','Healthy Exposed')
# #plot_NRF2<-plot_NRF2[plot_NRF2$`AD Exposed`>0.2 | plot_NRF2$`Healthy Exposed`>0.2,]
# svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/NRF2_canonical_LOGFC0.2.svg",width=3, height=8)
# ipa_palette_ORANGE<-CustomPalette(low = 'white', high = '#ea6901',mid='white',k = 200)
#
# par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))
# heatmap.2(as.matrix(plot_NRF2), Rowv = all_NRF2, Colv = F, trace = "none", na.color = "Grey", cexCol = 0.7,cexRow=0.7,density.info="none",
#           scale="none", dendrogram = "none",col = ipa_palette_ORANGE,key = F)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)
#
# dev.off()

############
# Recreate IPA
path_data<-readxl::read_xlsx("/home/claudia/Desktop/ADAIR/ADAIR_exp/IPA/pathway_heatmap.xlsx",sheet = 2)
path_data<-data.frame(path_data)
rownames(path_data)<-path_data$Canonical.Pathways
path_data<-path_data[2:3]
colnames(path_data)<-c('AD Exposed','Healthy Exposed')
path_data<-path_data[path_data[,1]+path_data[,2] != 0,]
svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/ipa_pathways_R.svg",width=3, height=12)
par(oma = c(4, 0, 0, 11) , mar=c(0,0,0,0))
heatmap.2(as.matrix(path_data), Rowv = rownames(path_data),Colv = F, trace = "none", na.color = "Grey", cexCol = 0.7,cexRow=0.7,density.info="none",
          scale="none", dendrogram = "none",col = ipa_palette,key = F)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)
dev.off()

upst_data<-readxl::read_xlsx("/home/claudia/Desktop/ADAIR/ADAIR_exp/IPA/upstream_analysis_heatmap.xlsx",sheet = 2)
upst_data<-data.frame(upst_data)
rownames(upst_data)<-upst_data$Transcription.Factors
upst_data<-upst_data[2:3]
colnames(upst_data)<-c('AD Exposed','Healthy Exposed')
svg("/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/figures/ipa_upstream_R.svg",width=3, height=8)
par(oma = c(4, 0, 0, 0) , mar=c(0,0,0,0))
heatmap.2(as.matrix(upst_data), Rowv = rownames(upst_data),Colv = F, trace = "none", na.color = "Grey", cexCol = 0.7,cexRow=0.7,density.info="none",
          scale="none", dendrogram = "none",col = ipa_palette,key = F)#hcl.colors(200,palette='RdBU',rev=T))#direction=-1#,breaks=colors,col=my_palette)
dev.off()
