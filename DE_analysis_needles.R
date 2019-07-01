##This script is doing DE analysis for needles
##alexander Vergara
## July 10 2017
#' ---
#' title: "DE analysis needles cold stress data"
#' author: "Alexander Vergara"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/vhurry/COLD_STRESS/cold-stress-needles")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/vhurry/COLD_STRESS/cold-stress-needles")
#' ```
#' 
#' 
### ============= 1. Start and set-up of the data =================
## 1.1 load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(stringr))
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/rmd.R")  
library(limma)

### ==============================
## 1.2 set the working directory
### ==============================
setwd("/mnt/picea/projects/spruce/vhurry/COLD_STRESS/cold-stress-needles")

#' Create a result dir
system("mkdir DE_Needles_DATA_2017_manuscript_REMOVING_6_OUTLIERS")
outdir <- file.path("DE_Needles_DATA_2017_manuscript_REMOVING_6_OUTLIERS")
dir.create(outdir,showWarnings=FALSE)


### ==============================
## 1.3 read the samples details
### ==============================/
samples <- read.csv("~/Git/UPSCb/projects/spruce_cold_stress/spruce-needles-cold-stress/samples2.csv", head=T, sep=",")
samples
### ==============================
## 1.4 read the HTSeq files in a matrix
## names are set according to the sample.csv!
### ==============================
res <- mclapply(dir("htseq_without_6_outliers",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=9)
names(res) <- sub(".*_P1554_","",sub("_sortmerna_trimmomatic_STAR.txt","",dir("htseq_without_6_outliers",pattern="*.txt")))

names(res) <- samples$SampleName[match(names(res),samples$ID)]

### ==============================
## 1.5 get the count table 
### ==============================
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
#addInfo <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]  #antes era 3 en vez de 2 !!!
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

count.table <- count.table[,order(colnames(count.table))]
head(count.table)
#ACA VOY y count table esta Ok
dim(count.table)
### ===========2. Creat the frame of work ===================
## 2.1 Create filters
### ==============================
names <- colnames(count.table)
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
###### Simple way
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

dds <- DESeq(dds)

### ==============================
## Comparing control_ vs 6h 5°C__
### ==============================
df
res0 <- results(dds,contrast=c("treatment","6h 5°C__","control_"))
head(res0)

VolcanoPlotMA(res0,alpha=0.01)
sum(res0$padj <= 0.01,na.rm=TRUE)
plot(density(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj)]))
table(sign(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj)]))
table(sign(res0$log2FoldChange[res0$padj <= .01 & !is.na(res0$padj) & abs(res0$log2FoldChange)>=2]))

write.table(res0,file.path(outdir,"control_vs_6h_5_°C.txt"))


### ==============================
## Comparing control vs 24h_5_°C
### ==============================

res1 <- results(dds,contrast=c("treatment","24h 5°C_","control_"))
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

res2 <- results(dds,contrast=c("treatment","3d 5°C__","control_"))
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
res3 <- results(dds,contrast=c("treatment", "10d 5°C_","control_"))
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
res4 <- results(dds,contrast=c("treatment","6h -5°C_","control_"))
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
res5 <- results(dds,contrast=c("treatment","24h -5°C", "control_"))
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

res6 <- results(dds,contrast=c("treatment","3d -5°C_","control_"))
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


res7 <- results(dds,contrast=c("treatment","10d -5°C","control_"))
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

