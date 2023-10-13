library(VennDiagram)
library(openxlsx)
library(stringr)


# Rotterdam
AD<-read.table('/home/claudia/Desktop/ADAIR/Rotterdam/excel_res/all_AD.txt',header = 1)
HI<-read.table('/home/claudia/Desktop/ADAIR/Rotterdam/excel_res/CN_HI.txt',header=1)

intersect(AD$features,HI$features) 
# 
# venn.diagram(list("HI"=HI[HI$genes!='L',]$features, "LO"=LO[LO$genes!='L',]$features, "AD"=AD[AD$genes!='L',]$features), cat.pos=2,
#              filename="/home/claudia/Desktop/ADAIR/Rotterdam/overlap/rotterdam.tiff")


# hOM
h.e_h.c<-readxl::read_xlsx('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/res/res/de_mrna.xlsx',sheet=2)
d.e_d.c<-readxl::read_xlsx('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/res/res/de_mrna.xlsx',sheet=3)
d.c_h.c<-readxl::read_xlsx('/home/claudia/Desktop/ADAIR/ADAIR_exp/DE_analysis/res/res/de_mrna.xlsx',sheet=5)


# venn.diagram(list("h.e_h.c"=h.e_h.c$ENSEMBL, "d.e_d.c"=d.e_d.c$ENSEMBL, "d.c_h.c"=d.c_h.c$ENSEMBL), cat.pos=2,
#              filename="/home/claudia/Desktop/ADAIR/Rotterdam/overlap/hOM.tiff")
# 

# check if h.e vs h.c DE genes are in any altered pathway of HI
common_h.e_h.c_HI<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
h.e_h.c_genes<-h.e_h.c$SYMBOL
for (g in h.e_h.c_genes){
    h.e_h.c_genes<-h.e_h.c_genes[h.e_h.c_genes!='NA']
    p_number<-0
    for (p in HI$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(HI[p_number,"features"])
            print(strsplit(HI[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(HI[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_h.e_h.c_HI<-rbind(common_h.e_h.c_HI,c(inter,HI[p_number,"features"],HI[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_h.e_h.c_HI<-common_h.e_h.c_HI[2:dim(common_h.e_h.c_HI)[1],]
common_h.e_h.c_HI$long <- unlist(lapply(common_h.e_h.c_HI$long, function(x){str_replace_all(x,'_',' ')}))

# svg("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_h.e_h.c_HI.svg")#,width = 800)
# #png("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_h.e_h.c_HI.png",width = 800)
# #par(mar = c(3, 35, 4, 2) + 0.2, mgp = c(25, 1, 0))    # Modify plotting options
# par(mar = c(3, 25, 1, 3) + 0.2, mgp = c(24, 1, 0))    # Modify plotting options
# barplot(sort(table(common_h.e_h.c_HI$long),decreasing=F), cex.names =0.6,
#         xlab="Control exposed DE gene frequency",ylab="Pathways altered in class HI",las=1,horiz=T, col='#4cae4aff')
# dev.off()

df_h.e_h.c_HI<-data.frame('Features'=NA,'Ctrl DE VS Ctrl Air'=NA)
for (path in unique((common_h.e_h.c_HI$long))){
    acutes<-common_h.e_h.c_HI[common_h.e_h.c_HI$long==path,'acute']
    df_h.e_h.c_HI<-rbind(df_h.e_h.c_HI, c(path, toString(acutes)))
}

# check if D.e vs D.c DE genes are in any altered pathway of HI
common_d.e_d.c_HI<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
d.e_d.c_genes<-d.e_d.c$SYMBOL
for (g in d.e_d.c_genes){
    d.e_d.c_genes<-d.e_d.c_genes[d.e_d.c_genes!='NA']
    p_number<-0
    for (p in HI$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(HI[p_number,"features"])
            print(strsplit(HI[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(HI[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_d.e_d.c_HI<-rbind(common_d.e_d.c_HI,c(inter,HI[p_number,"features"],HI[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_d.e_d.c_HI<-common_d.e_d.c_HI[2:dim(common_d.e_d.c_HI)[1],]

common_d.e_d.c_HI$long <- unlist(lapply(common_d.e_d.c_HI$long, function(x){str_replace_all(x,'_',' ')}))

#common_d.e_d.c_HI[120,'long']<-'BIOSYNTHESIS OF THE N-GLYCAN PRECURSOR AND TRANSFER TO A NASCENT PROTEIN'
#png("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.e_d.c_HI.png",width = 800)
# svg("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.e_d.c_HI.svg")#,width = 800)
# par(mar = c(3, 25, 1, 3) + 0.2, mgp = c(24, 1, 0))    # Modify plotting options
# barplot(sort(table(common_d.e_d.c_HI$long),decreasing=F), 
#         cex.names =0.6, col='#974ca1ff',
#         xlab="AD exposed DE gene frequency",ylab="Pathways altered in class HI",las=1,horiz=T)
# dev.off()


df_d.e_d.c_HI<-data.frame('Features'=NA,'AD DE VS AD Air'=NA)
for (path in unique((common_d.e_d.c_HI$long))){
    acutes<-common_d.e_d.c_HI[common_d.e_d.c_HI$long==path,'acute']
    df_d.e_d.c_HI<-rbind(df_d.e_d.c_HI, c(path, toString(acutes)))
}


# check if d.c vs h.c DE genes are in any altered pathway of AD
common_d.c_h.c_AD<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
d.c_h.c_genes<-d.c_h.c$SYMBOL
for (g in d.c_h.c_genes){
    d.c_h.c_genes<-d.c_h.c_genes[d.c_h.c_genes!='NA']
    p_number<-0
    for (p in AD$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(AD[p_number,"features"])
            print(strsplit(AD[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(AD[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_d.c_h.c_AD<-rbind(common_d.c_h.c_AD,c(inter,AD[p_number,"features"],AD[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_d.c_h.c_AD<-common_d.c_h.c_AD[2:dim(common_d.c_h.c_AD)[1],]
common_d.c_h.c_AD$long <- unlist(lapply(common_d.c_h.c_AD$long, function(x){str_replace_all(x,'_',' ')}))

# svg("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.c_h.c_AD.svg")#,width = 800)
# par(mar = c(3, 25, 1, 3) + 0.2, mgp = c(24, 1, 0))    # Modify plotting options
# barplot(sort(table(common_d.c_h.c_AD$long),decreasing=F), cex.names =0.3,
#         xlab="AD DE gene frequency",ylab="Pathways altered in class AD",las=1,horiz=T, col='#974ca1ff')
# dev.off()

df_d.c_h.c_AD<-data.frame('Features'=NA,'AD Air VS Ctrl Air'=NA)
for (path in unique((common_d.c_h.c_AD$long))){
    acutes<-common_d.c_h.c_AD[common_d.c_h.c_AD$long==path,'acute']
    df_d.c_h.c_AD<-rbind(df_d.c_h.c_AD, c(path, toString(acutes)))
}


# intersect pathways of h.e VS h.c in HI and d.e VS d.c in HI
# venn.diagram(list("h.e VS h.c in HI"=common_h.e_h.c_HI$long, "d.e VS d.c in HI"=common_d.e_d.c_HI$long),
#              filename="/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_HI.tiff",fill = c("darkred", "darkgreen"))
# 
# venn.diagram(list("h.e VS h.c in LO"=common_h.e_h.c_LO$long, "d.e VS d.c in LO"=common_d.e_d.c_LO$long),
#              filename="/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_LO.tiff",fill = c("darkred", "darkgreen"))



### DIFF MARKER

# check if h.c vs h.e DE genes are in any altered pathway of AD
common_h.e_h.c_AD<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
h.e_h.c_genes<-h.e_h.c$SYMBOL
for (g in h.e_h.c_genes){
    h.e_h.c_genes<-h.e_h.c_genes[h.e_h.c_genes!='NA']
    p_number<-0
    for (p in AD$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(AD[p_number,"features"])
            print(strsplit(AD[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(AD[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_h.e_h.c_AD<-rbind(common_h.e_h.c_AD,c(inter,AD[p_number,"features"],AD[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_h.e_h.c_AD<-common_h.e_h.c_AD[2:dim(common_h.e_h.c_AD)[1],]
common_h.e_h.c_AD$long <- unlist(lapply(common_h.e_h.c_AD$long, function(x){str_replace_all(x,'_',' ')}))

# svg("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_h.e_h.c_AD.svg")#,width = 800)
# #png("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_h.e_h.c_AD.png",width = 800,height=1200)
# par(mar = c(3, 25, 1, 3) + 0.2, mgp = c(24, 1, 0))    # Modify plotting options
# barplot(sort(table(common_h.e_h.c_AD$long),decreasing=F), cex.names =0.6,
#         xlab="AD exposed DE gene frequency",ylab="Pathways altered in class AD",las=1,horiz=T, col='#4cae4aff')
# dev.off()


df_h.e_h.c_AD<-data.frame('Features'=NA,'Ctrl DE VS Ctrl Air'=NA)
for (path in unique((common_h.e_h.c_AD$long))){
    acutes<-common_h.e_h.c_AD[common_h.e_h.c_AD$long==path,'acute']
    df_h.e_h.c_AD<-rbind(df_h.e_h.c_AD, c(path, toString(acutes)))
}

# check if d.c vs d.e DE genes are in any altered pathway of AD
common_d.e_d.c_AD<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
d.e_d.c_genes<-d.e_d.c$SYMBOL
for (g in d.e_d.c_genes){
    d.e_d.c_genes<-d.e_d.c_genes[d.e_d.c_genes!='NA']
    p_number<-0
    for (p in AD$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(AD[p_number,"features"])
            print(strsplit(AD[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(AD[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_d.e_d.c_AD<-rbind(common_d.e_d.c_AD,c(inter,AD[p_number,"features"],AD[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_d.e_d.c_AD<-common_d.e_d.c_AD[2:dim(common_d.e_d.c_AD)[1],]
common_d.e_d.c_AD$long <- unlist(lapply(common_d.e_d.c_AD$long, function(x){str_replace_all(x,'_',' ')}))

#png("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.e_d.c_AD.png",width = 800,height=1200, )
# svg("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.e_d.c_AD.svg")
# par(mar = c(3, 25, 1, 3) + 0.2, mgp = c(24, 1, 0))    # Modify plotting options
# barplot(sort(table(common_d.e_d.c_AD$long),decreasing=F), cex.names =0.3,
#         xlab="AD exposed DE gene frequency",ylab="Pathways altered in class AD",las=1,horiz=T, col='#974ca1ff')
# dev.off()

df_d.e_d.c_AD<-data.frame('Features'=NA,'AD DE VS AD Air'=NA)
for (path in unique((common_d.e_d.c_AD$long))){
    acutes<-common_d.e_d.c_AD[common_d.e_d.c_AD$long==path,'acute']
    df_d.e_d.c_AD<-rbind(df_d.e_d.c_AD, c(path, toString(acutes)))
}

# venn.diagram(list("h.e VS h.c in AD"=common_h.e_h.c_AD$long, "d.e VS d.c in AD"=common_d.e_d.c_AD$long),
#              filename="/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_AD.tiff",fill = c("darkred", "darkgreen"))



# check if d.c vs d.e DE genes are in any altered pathway of HI
common_d.c_h.c_HI<-data.frame(acute = NA, long = NA, sing_genes = NA, genes_overlap=NA)
d.c_h.c_genes<-d.c_h.c$SYMBOL
for (g in d.c_h.c_genes){
    d.c_h.c_genes<-d.c_h.c_genes[d.c_h.c_genes!='NA']
    p_number<-0
    for (p in HI$genes){
        p_number<-p_number+1
        p_genes<-strsplit(p,',')[[1]]
        p_genes<-p_genes[p_genes!='NA']
        inter<-intersect(g,p_genes)
        if (length(inter)!=0){
            print(inter)
            #print(HI[p_number,"features"])
            print(strsplit(HI[p_number,"sign_genes"],',')[[1]])
            sign_genes_inter<-intersect(g,strsplit(HI[p_number,"sign_genes"],',')[[1]])
            if (length(sign_genes_inter)==0){
                sign_genes_inter=NA
            }
            common_d.c_h.c_HI<-rbind(common_d.c_h.c_HI,c(inter,HI[p_number,"features"],HI[p_number,"sign_genes"],sign_genes_inter))
        }
        
    }
}
common_d.c_h.c_HI<-common_d.c_h.c_HI[2:dim(common_d.c_h.c_HI)[1],]

# png("/home/claudia/Desktop/ADAIR/Rotterdam/overlap/common_d.c_h.c_HI.png",width = 800,height=1200)
# par(mar = c(7, 45, 4, 6) + 0.2, mgp = c(43, 1, 0))    # Modify pHItting options
# barplot(sort(table(common_d.c_h.c_HI$long),decreasing=F), 
#         main = "N° of genes\n in pathways altered\n in HI",
#         xlab="Acute genes frequency",ylab="Pathways altered in cohort (HI)",las=1,horiz=T, cex.names = 0.75)
# dev.off()

df_d.c_h.c_HI<-data.frame('Features'=NA,'AD Air VS Ctrl Air'=NA)
for (path in unique((common_d.c_h.c_HI$long))){
    acutes<-common_d.c_h.c_HI[common_d.c_h.c_HI$long==path,'acute']
    df_d.c_h.c_HI<-rbind(df_d.c_h.c_HI, c(path, toString(acutes)))
}
 
#############
## SAVE FILES

# HI
HI['features'] <- unlist(lapply(HI['features'], function(x){str_replace_all(x,'_',' ')}))

HI<-merge(HI['features'],df_h.e_h.c_HI, by.x='features',by.y='Features',all.x = T)
HI<-merge(HI, df_d.e_d.c_HI, by.x='features',by.y='Features',all.x = T)
HI<-merge(HI, df_d.c_h.c_HI, by.x='features',by.y='Features',all.x = T)

# AD
AD['features'] <- unlist(lapply(AD['features'], function(x){str_replace_all(x,'_',' ')}))

AD<-merge(AD['features'],df_h.e_h.c_AD, by.x='features',by.y='Features',all.x = T)
AD<-merge(AD, df_d.e_d.c_AD, by.x='features',by.y='Features',all.x = T)
AD<-merge(AD, df_d.c_h.c_AD, by.x='features',by.y='Features',all.x = T)

library(openxlsx)

# iniate workbook object
wb <- createWorkbook()

# add 2 worksheets to the workbook object, arbitrary names
addWorksheet(wb, "AD Significant Features")
addWorksheet(wb, "HI Significant Features")

# write data to worksheets in the workbook object
writeData(wb, 1, AD)
writeData(wb, 2, HI)

# save the workbook to a file
saveWorkbook(wb, "/home/claudia/Desktop/ADAIR/Rotterdam/overlap/overlap.xlsx")
