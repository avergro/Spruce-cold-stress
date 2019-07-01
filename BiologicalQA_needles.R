#' ---
#' title: "Needles cold stress data from Spruce"
#' author: "Nicolas Delhomme-Alexander Vergara"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("~/cold-stress-needles")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/cold-stress-needles")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
#' ```{r detach, echo=FALSE, eval=FALSE}
#' detach("package:limma")
#' ```
#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
samples <- read.csv("~/samples_needles.csv")
head(samples)
samples$ID
#' Create a result dir
 system("mkdir ANALYSIS_NEEDLES")
outdir <- file.path("ANALYSIS_NEEDLES","manuscript_version")
dir.create(outdir,showWarnings=FALSE)


#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq_raw",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub(".*_P1554_","",sub("_sortmerna.*\\.txt","",dir("htseq_raw",pattern="*.txt")))

#' Reorder the sample data.frame according to the way
#' the results were read 
samples_ok <- samples[match(names(res),samples$ID),]
samples_ok
#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
count.table <- count.table[,match(samples$ID,colnames(count.table))]

write.csv(count.table,file.path(outdir,"raw-unormalised-data_45_samples_Needles.csv"))

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' Plot the stats
#' 
#' 
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 75%
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected, around 100X
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#'
#' 
plot.multidensity(log10(count.table),
                  col=pal[as.integer(samples$SampleID)],
                  legend.x="topright",
                  legend=levels(samples$SampleID),
                  legend.col=pal[1:8],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
conditions <- colnames(count.table)
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)

write.csv(vst,file.path(outdir,"VST_library-size-normalized_data_45_samples_Needles.csv"))


#' Validate the VST 
#' 
#' Visualize the corrected mean - sd
#' relationship. It is fairly linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is
#' due to genes having few counts in
#' a few subset of the samples and hence 
#' having a higher variability. This is
#' expected.
meanSdPlot(vst[rowSums(count.table)>0,], ylim = c(0,2.5))

#' # QC on the normalised data
#' 
#' #################################################################################
#' ## PCA
#' 
#' First perform a Principal Component
#' Analysis (PCA) of the data
#'  to do a quick quality assessment; 
#' i.e. replicate should cluster
#' and the first 2-3 dimensions should 
#' be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA 3 first dimensions

mar=c(5.1,4.1,4.1,8.1, xpd=TRUE)
#mar – A numeric vector of length 4, which sets the margin sizes in the 

s3d <- scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$SampleID)],
              pch=19)

s3d 

legend(5, 5, pch=19,
       col=pal[1:8],
       legend=levels(samples$SampleID))

###########################################
# Now we try using a new column od samples object called temperature like sub-population
mar=c(5.1,4.1,4.1,2.1)
popIDs <- unique(samples$Temperature)

# make plot
s3d <- scatterplot3d(pc$x[,1],
                     pc$x[,2],
                     pc$x[,3],
                     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                     zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                     color=as.numeric(samples$Temperature),
                     pch=as.numeric(samples$Temperature),
                     y.margin.add=1)
s3d
# add legend 
legend("top", legend = popIDs,
       pch = as.numeric(popIDs), 
       col=as.numeric(popIDs),
       inset = -0.25, xpd = TRUE, horiz = TRUE)
##DONE
# 
###########################################
# Now we try using all the labels of sub-populations
head(samples)
popIDs <- unique(samples$SampleID)

s3d <- scatterplot3d(pc$x[,1],
                     pc$x[,2],
                     pc$x[,3],
                     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                     zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                     color=as.numeric(samples$SampleID),
                     pch=as.numeric(samples$SampleID),
                     y.margin.add=1)
s3d
# 
# add legend
legend(6,5, legend = popIDs,
       pch = as.numeric(popIDs), 
       col=as.numeric(popIDs),
       inset = -0.25, xpd = TRUE, horiz = FALSE)

###########Done

##########################################################################################
#Now PCA using colours to distigish how deep the dots are in 3D and their names
library(plot3D)
pc3 <- cbind(pc$x[,1],pc$x[,2],pc$x[,3])
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
#Saving coordinates to add text
pcOK<-as.data.frame((pc3))
x <- pcOK$V1
y <- pcOK$V2
z <- pcOK$V3
#ploting 3D and names
scatter3D(pc$x[,1],
          pc$x[,2],
          pc$x[,3],
          xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
          ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
          zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                    pch=19)
text3D(x,y,z,labels = row.names(pc$x), add= TRUE, colkey= FALSE, cex= 0.5)

#######################################################################
#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]

#' First with their samples numbers IDs
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#######################################################################
#Now we explore the data to filter low expressed genes
#In total we have 70736 gene models
#######################################################################
#The function takes 4 parameters:
#(i) the expression matrix; 
#(ii) a vector of "conditions" matching the columns in the expression matrix;
#(iii) an expression cutoff X and 
#(iv) a number of replicates cutoff Y.

#' ## Raw data tissues
#' analyizing raw counts from phase1 (5°C) and and phase2 (-5°C)
#' Keep genes with 2 or more reads in at least 3 replicates
head(samples)
#' In each condition we have 5 reps
#' We considered "outliers" for DE analysis the samples: 105(control), 116(3d 5°C), 119(24h 5°C),
#' 130(6h -5°C),135(24h -5°C), 142(10d -5°C)
#' So we set 3 like the minimum number of biological replicates because 3 = (4-1) 
#' So althought in some cases we removed 1 replicate we have almost 4 replicates anyway in these cases. 
#' So we used that minus 1 like criterion.
#' Function geneSelect to remove low expressed genes:
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

head(samples)
head(count.table)
dim(count.table)
#sum(sel)
####
#' get the boolean selection
plot(density(log10(as.matrix(count.table)+1))) 
plot(density(log10(as.matrix(count.table))))
sel <- geneSelect(count.table,samples$SampleID,0,0)
#' how many genes are selected
sum(sel)  #no genes removed
class(sel)
head(sel)
plot(density(log10(as.matrix(count.table[sel,]))))

sel <- geneSelect(count.table,samples$SampleID,2,3)
sum(sel)
count.table[sel,]
survivors <- row.names(count.table[sel,])
dim(count.table[sel,])
plot(density(log10(as.matrix(count.table[sel,]))))    #,main="survivors genes counts density"))
#saving survivors genes from needles using 2,3 filter (num min reads, num min of replicates with 2 or more reads)
write.table(survivors, file="Survivors_genes_needles_38316.txt", quote = F, row.names = F, col.names = F)

#pdf("HeatMap_Counts_Needles_38316_survivors_genes.pdf")
hpal <- colorRampPalette(c("blue","white","red"))(100)
heatmap.2(as.matrix(count.table[sel,]),scale="row",
          labRow = FALSE,labCol = samples$SampleName,
          trace="none",col=hpal)
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#'  sessionInfo()
#' ```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.2 LTS
## 
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] plot3D_1.1                 hexbin_1.27.1             
##  [3] vsn_3.44.0                 scatterplot3d_0.3-40      
##  [5] RColorBrewer_1.1-2         pander_0.6.0              
##  [7] gplots_3.0.1               DESeq2_1.16.1             
##  [9] SummarizedExperiment_1.6.3 DelayedArray_0.2.7        
## [11] matrixStats_0.52.2         Biobase_2.36.2            
## [13] GenomicRanges_1.28.3       GenomeInfoDb_1.12.2       
## [15] IRanges_2.10.2             S4Vectors_0.14.3          
## [17] BiocGenerics_0.22.0       
## 
## loaded via a namespace (and not attached):
##  [1] bit64_0.9-7             splines_3.4.0          
##  [3] gtools_3.5.0            Formula_1.2-1          
##  [5] affy_1.54.0             latticeExtra_0.6-28    
##  [7] blob_1.1.0              GenomeInfoDbData_0.99.0
##  [9] yaml_2.1.14             RSQLite_2.0            
## [11] backports_1.1.0         lattice_0.20-35        
## [13] limma_3.32.2            digest_0.6.12          
## [15] XVector_0.16.0          checkmate_1.8.2        
## [17] colorspace_1.3-2        preprocessCore_1.38.1  
## [19] htmltools_0.3.6         Matrix_1.2-10          
## [21] plyr_1.8.4              XML_3.98-1.9           
## [23] misc3d_0.8-4            genefilter_1.58.1      
## [25] zlibbioc_1.22.0         xtable_1.8-2           
## [27] scales_0.4.1            gdata_2.18.0           
## [29] affyio_1.46.0           BiocParallel_1.10.1    
## [31] htmlTable_1.9           tibble_1.3.3           
## [33] annotate_1.54.0         ggplot2_2.2.1          
## [35] nnet_7.3-12             lazyeval_0.2.0         
## [37] survival_2.41-3         magrittr_1.5           
## [39] memoise_1.1.0           evaluate_0.10          
## [41] foreign_0.8-68          BiocInstaller_1.26.0   
## [43] tools_3.4.0             data.table_1.10.4      
## [45] stringr_1.2.0           munsell_0.4.3          
## [47] locfit_1.5-9.1          cluster_2.0.6          
## [49] AnnotationDbi_1.38.1    compiler_3.4.0         
## [51] caTools_1.17.1          rlang_0.1.1            
## [53] grid_3.4.0              RCurl_1.95-4.8         
## [55] htmlwidgets_0.8         labeling_0.3           
## [57] bitops_1.0-6            base64enc_0.1-3        
## [59] rmarkdown_1.6           gtable_0.2.0           
## [61] DBI_0.7                 gridExtra_2.2.1        
## [63] knitr_1.16              bit_1.1-12             
## [65] Hmisc_4.0-3             rprojroot_1.2          
## [67] KernSmooth_2.23-15      stringi_1.1.5          
## [69] Rcpp_0.12.11            geneplotter_1.54.0     
## [71] rpart_4.1-11            acepack_1.4.1
