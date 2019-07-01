# title: DE analysis roots cold stress data
# author: Alexander Vergara
# Date: Jul 10 2017

#' ---
#' title: "DE analysis roots cold stress data"
#' author: "Alexander Vergara"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("~/cold-stress-roots")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/cold-stress-roots")
#' ```


### ============= 1. Start and set-up of the data =================
## 1.1 load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(knitr))
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/rmd.R")  
library(limma)

### ==============================
## 1.2 set the working directory
### ==============================
#' Create a result dir
system("mkdir DE_roots")
outdir <- file.path("DE_roots")
dir.create(outdir,showWarnings=FALSE)

### ==============================
## 1.3 read the samples details
### ==============================
samples <- read.csv("~/samples2.csv")

### ==============================
## 1.4 read the HTSeq files in a matrix
## names are set according to the sample.csv!
### ==============================
res <- mclapply(dir("htseq_technical_replicates",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=9)
names(res) <- sub(".*_P1406_","",sub("_sortmerna_trimmomatic_STAR.txt","",dir("htseq_technical_replicates",pattern="*.txt")))
#names(res) <- sub(".*_P1406_","P1406",sub("_sortmerna_trimmomatic_STAR.txt","",dir("htseq_technical_replicates",pattern="*.txt")))

names(res) <- samples$SampleName[match(names(res),samples$ID)]

### ==============================
## 1.5 get the count table 
### ==============================
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]


count.table <- count.table[,order(colnames(count.table))]
head(count.table)
#ACA VOY y count table esta Ok
dim(count.table)
### ===========2. Creat the frame of work ===================
## 2.1 Create filters
### ==============================

names <- sub("P1406_","",colnames(count.table))
samples <- substr(x=colnames(count.table),1,10)
treatments <- substr(x=colnames(count.table),1,8)

### ==============================
## 2.2 create the design matrix - this is the design for my biological question
### ==============================
df <- data.frame (
  name=names,
  sample=samples,

  treatment=treatments)
df
head(df)

### ==============================
## 2.2 create dds object
### ==============================
dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = df,
                              design = ~treatment)
                  
## check the size factors

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## estimate the disperison (takes a long time!)

dds <- estimateDispersions(dds)
plotDispEsts(dds)

### ===========3. Negative binomial Wald Test===================
## run the test - negative binomial Wald test (takes a long time!)
### ==============================
dds <- nbinomWaldTest(dds)

### ==============================
## create ddsC (for comparison by condition)
### ==============================
ddsC <- dds
### ===========4. Methods and Functions===================
## 
### ==============================
setGeneric(name="VolcanoPlotMA",def=function(object,alpha=0.01){
  standardGeneric("VolcanoPlotMA")})

setMethod(f="VolcanoPlotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha
            
            ## plot
            heatscatter(object$log2FoldChange[sel],
                        -log10(object$padj[sel]),
                        main="Volcano",xlab="Log2 Fold Change",
                        ylab="- log(10) adj. p-value")
            
            ## legend
            legend("topleft",bty="n",paste("cutoff @",alpha),lty=2,col="gray")
            
            ## points
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="lightblue",pch=19)
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="dodgerblue3",pch=19,cex=0.5)
            
            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,col="gray")
          })

setMethod(f="plotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha
            
            ## graphic params
            orig.par <- par(no.readonly=TRUE)
            par(mfrow=c(2,1))
            
            ## plots
            kde2dplot(log10(object$baseMean[sel]),
                      object$log2FoldChange[sel],
                      grid=250,ncol=30,nlevels=10,
                      main="MA density estimation"
            )
            
            heatscatter(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        add.contour=TRUE,main="MA")
            
            mtext(paste(sum(sel2),"sig. feats. @",alpha,"cutoff"),
                  side=1,line=2)
            
            points(log10(object$baseMean[sel][sel2]),
                   object$log2FoldChange[sel][sel2],
                   col="darkred",pch=19,cex=.5)
            
            par(orig.par,no.readonly=TRUE)
            invisible(TRUE)
          })

##########################################################
####============5. Comparison of groups ==================
##########################################################

ddsG <- DESeq(dds)

### ==============================
## Comparing control_ vs 6h 5°C__
### ==============================
df
res0 <- results(ddsG,contrast=c("treatment","6h 5°C__","control_"))
head(res0) 

VolcanoPlotMA(res0,alpha=0.01)
sum(res0$padj <= 0.01,na.rm=TRUE)
plot(density(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj)]))
table(sign(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj)]))
table(sign(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj) & abs(res0$log2FoldChange)>=2]))

write.table(res0,file.path(outdir,"control_vs_6h_5_°C.txt"))

######################
### ==============================
## Comparing control vs 24h_5_°C
### ==============================

res1 <- results(ddsG,contrast=c("treatment","24h 5°C_","control_"))
head(res1)

VolcanoPlotMA(res1,alpha=0.01)
sum(res1$padj <= 0.01,na.rm=TRUE)
plot(density(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj)]))
table(sign(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj)]))
table(sign(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj) & abs(res1$log2FoldChange)>=2]))

write.table(res1,file.path(outdir,"control_vs_24h_5_°C.txt"))

### ==============================
## Comparing control vs 3d_5_°C
### ==============================

res2 <- results(ddsG,contrast=c("treatment","3d 5°C__","control_"))
head(res2)

VolcanoPlotMA(res2,alpha=0.01)
sum(res2$padj <= 0.01,na.rm=TRUE)
plot(density(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj)]))
table(sign(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj)]))
table(sign(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj) & abs(res2$log2FoldChange)>=2]))

write.table(res2,file.path(outdir,"control_vs_3d_5_°C.txt"))

####
### ==============================
## Comparing control vs 10d_5_°C
### ==============================
res3 <- results(ddsG,contrast=c("treatment", "10d 5°C_","control_"))
head(res3)

VolcanoPlotMA(res3,alpha=0.01)
sum(res3$padj <= 0.01,na.rm=TRUE)
plot(density(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj)]))
table(sign(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj)]))
table(sign(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj) & abs(res3$log2FoldChange)>=2]))

write.table(res3,file.path(outdir,"control_vs_10d_5_°C.txt"))
####
### ==============================
## Comparing control vs 6h_-5_°C
### ==============================
####
res4 <- results(ddsC,contrast=c("treatment","6h -5°C_","control_"))
head(res4)

VolcanoPlotMA(res4,alpha=0.01)
sum(res4$padj <= 0.01,na.rm=TRUE)
plot(density(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj)]))
table(sign(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj)]))
table(sign(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj) & abs(res4$log2FoldChange)>=2]))

write.table(res4,file.path(outdir,"control_vs_6h_-5_°C.txt"))

####
### ==============================
## Comparing control vs 24h_-5_°C
### ==============================
####
res5 <- results(ddsC,contrast=c("treatment","24h -5°C", "control_"))
head(res5)

VolcanoPlotMA(res5,alpha=0.01)
sum(res5$padj <= 0.01,na.rm=TRUE)
plot(density(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj)]))
table(sign(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj)]))
table(sign(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj) & abs(res5$log2FoldChange)>=2]))

write.table(res5,file.path(outdir,"control_vs_24h_-5_°C.txt"))
########
### ==============================
## Comparing control vs 3d_-5_°C
### ==============================

res6 <- results(ddsC,contrast=c("treatment","3d -5°C_","control_"))
head(res6)

VolcanoPlotMA(res6,alpha=0.01)
sum(res6$padj <= 0.01,na.rm=TRUE)
plot(density(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj)]))
table(sign(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj)]))
table(sign(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj) & abs(res6$log2FoldChange)>=2]))


write.table(res6,file.path(outdir,"control_vs_3d_-5_°C.txt"))
###
### ==============================
## Comparing control vs 10d_-5_°C
### ==============================
####


res7 <- results(ddsC,contrast=c("treatment","10d -5°C","control_"))
head(res7)

VolcanoPlotMA(res7,alpha=0.01)
sum(res7$padj <= 0.01,na.rm=TRUE)
plot(density(res7$log2FoldChange[res7$padj <= .01 & !is.na(res7$padj)]))
table(sign(res7$log2FoldChange[res7$padj <= .01 & !is.na(res7$padj)]))
table(sign(res7$log2FoldChange[res7$padj <= .01 & !is.na(res7$padj) & abs(res7$log2FoldChange)>=2]))


write.table(res7,file.path(outdir,"control_vs_10d_-5_°C.txt"))

# Session Info
sessionInfo()


###### End
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.2 LTS
##
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
##
## locale:
## [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C
## [3] LC_TIME=en_US.UTF-8 LC_COLLATE=en_US.UTF-8
## [5] LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8
## [7] LC_PAPER=en_US.UTF-8 LC_NAME=C
## [9] LC_ADDRESS=C LC_TELEPHONE=C
23
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
##
## attached base packages:
## [1] parallel stats4 stats graphics grDevices utils datasets
## [8] methods base
##
## other attached packages:
## [1] limma_3.32.2 knitr_1.16
## [3] LSD_3.0 scatterplot3d_0.3-40
## [5] vsn_3.44.0 RColorBrewer_1.1-2
## [7] DESeq2_1.16.1 SummarizedExperiment_1.6.3
## [9] DelayedArray_0.2.7 matrixStats_0.52.2
## [11] Biobase_2.36.2 GenomicRanges_1.28.3
## [13] GenomeInfoDb_1.12.2 IRanges_2.10.2
## [15] S4Vectors_0.14.3 BiocGenerics_0.22.0
##
## loaded via a namespace (and not attached):
## [1] Rcpp_0.12.11 locfit_1.5-9.1
## [3] lattice_0.20-35 rprojroot_1.2
## [5] digest_0.6.12 plyr_1.8.4
## [7] backports_1.1.0 acepack_1.4.1
## [9] RSQLite_2.0 evaluate_0.10
## [11] BiocInstaller_1.26.0 ggplot2_2.2.1
## [13] zlibbioc_1.22.0 rlang_0.1.1
## [15] lazyeval_0.2.0 data.table_1.10.4
## [17] annotate_1.54.0 blob_1.1.0
## [19] rpart_4.1-11 Matrix_1.2-10
## [21] preprocessCore_1.38.1 checkmate_1.8.2
## [23] rmarkdown_1.6 splines_3.4.0
## [25] BiocParallel_1.10.1 geneplotter_1.54.0
## [27] stringr_1.2.0 foreign_0.8-68
## [29] htmlwidgets_0.8 RCurl_1.95-4.8
## [31] bit_1.1-12 munsell_0.4.3
## [33] compiler_3.4.0 base64enc_0.1-3
## [35] htmltools_0.3.6 nnet_7.3-12
## [37] tibble_1.3.3 gridExtra_2.2.1
## [39] htmlTable_1.9 GenomeInfoDbData_0.99.0
## [41] Hmisc_4.0-3 XML_3.98-1.9
## [43] bitops_1.0-6 grid_3.4.0
## [45] xtable_1.8-2 gtable_0.2.0
## [47] affy_1.54.0 DBI_0.7
## [49] magrittr_1.5 scales_0.4.1
## [51] stringi_1.1.5 XVector_0.16.0
## [53] genefilter_1.58.1 affyio_1.46.0
## [55] latticeExtra_0.6-28 Formula_1.2-1
## [57] tools_3.4.0 bit64_0.9-7
## [59] survival_2.41-3 AnnotationDbi_1.38.1
## [61] colorspace_1.3-2 cluster_2.0.6
## [63] memoise_1.1.0
###### End
