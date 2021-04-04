########################################################################################################
# Response to Technical Comment by Abruzzi et al. (2021). Code for:
# Figure 1 & 2
########################################################################################################

# set the working directory
setwd("~/Downloads/")

# ######################################################################################################
# load all the needed packages
# NB: you may need to install some packages before running this script
# ######################################################################################################

library(R.utils)
library(tools)
library(tidyr)
library(dplyr)
library(reshape2)
library(rain)
library(scales)
library(gplots)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(heatmap3)
library(eulerr)
library(RColorBrewer)
library(plotrix)
library(UpSetR)
library(GenomicFeatures)
library(GeneOverlap)


options(stringsAsFactors=FALSE)

# ######################################################################################################
# Download the data from GEO and unpack it to data directory
# ######################################################################################################

dir.create("data", showWarnings=FALSE)
setwd("data")

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111696/suppl/GSE111696%5FBmal1%5FexLiver%5Fgenes%2Efpkm%5Ftracking%2Etxt%2Egz", "GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt.gz")
gunzip("GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt.gz", overwrite = T)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111696&format=file&file=GSE111696%5FBmal1%5FMSF%5Fgenes%2Efpkm%5Ftracking%2Etxt%2Egz", "GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt.gz")
gunzip("GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt.gz", overwrite = T)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134333&format=file&file=GSE134333%5FDataset%5FAM%5FPM%5FTempComp%5FFPKM%5Fvalues%2Etxt%2Egz", "GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz")
gunzip("GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz", overwrite = T)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE114943&format=file", "GSE114943_RAW.tar")
untar("GSE114943_RAW.tar", exdir="./GSE114943_RAW")

setwd("..")

# ######################################################################################################
# reading the unpacked data 
# ######################################################################################################

exLiver_fpkm <- read.delim("data/GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt", quote = "")
MSF_fpkm <- read.delim("data/GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt", quote = "")
AM_PM_TempComp_FPKM <- read.delim("data/GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt", quote = "")

fileDir="./data/GSE114943_RAW"
fileID = "*._counts.csv.gz$"
seqIDcol=1

    ## Make sure that fileDir has trailing slash
    if(length(grep("/$",fileDir))==0){
      fileDir <- paste(fileDir,"/",sep="")
    }
    
    ## Get Files
    files <- list.files(fileDir,pattern=fileID)

    indFiles <- list()
    for(k in 1:length(files)){
      indFiles[[files[k]]] <- read.table(paste(fileDir,files[k],sep=""), col.names=c("GeneID", files[k]), header=TRUE, sep="\t")
    }
    
    ## Confirm that rows are in the same order:
    namesTemp <- lapply(indFiles,function(x) {as.character(x[,1])})
    names <- do.call(cbind,namesTemp)
    
    testName <- apply(names,1,function(x) {sum(x == x[1]) })

    nNotMatchAll <- sum(testName != length(indFiles))
    nNotMatchAll
    
    ## Limit down to just columns to keep 
    colsToKeep=2
    colsInt <- lapply(indFiles,function(x) { 
      out <- x[,colsToKeep]
      ## works only when only one col kept
      if(length(colsToKeep == 1)){
        out <- as.data.frame(out)
        if(class(colsToKeep) == "numeric"){
          names(out) <- names(x)[colsToKeep]
        } else {
          names(out) <- colsToKeep
        }
      }
      row.names(out) <- x[,1]
      return(out)
      })

GSE114943_Counts <- as.data.frame(colsInt)
head(GSE114943_Counts)
GSE114943_Counts <- GSE114943_Counts[ order(rownames(GSE114943_Counts)) , ]
head(GSE114943_Counts)

rm(fileDir, fileID, seqIDcol, files, indFiles, names, testName, colsToKeep, colsInt, k)

# ######################################################################################################
# subset the data
# ######################################################################################################

exLiver_WT_data <- exLiver_fpkm[, grep("WT_.*_FPKM", colnames(exLiver_fpkm))]
exLiver_KO_data <- exLiver_fpkm[, grep("KO_.*_FPKM", colnames(exLiver_fpkm))]
rownames(exLiver_WT_data) <- exLiver_fpkm$gene_id
rownames(exLiver_KO_data) <- exLiver_fpkm$gene_id

Bmal_msf_WT_data <- MSF_fpkm[, grep("WT_.*_FPKM", colnames(MSF_fpkm))]
Bmal_msf_KO_data <- MSF_fpkm[, grep("KO_.*_FPKM", colnames(MSF_fpkm))]
rownames(Bmal_msf_WT_data) <- MSF_fpkm$gene_id
rownames(Bmal_msf_KO_data) <- MSF_fpkm$gene_id

Bmal_KO_37C <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.37C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_37C <- AM_PM_TempComp_FPKM[, grep("WT.*.37C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_KO_32C <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.32C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_32C <- AM_PM_TempComp_FPKM[, grep("WT.*.32C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_KO_27C <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.27C.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_27C <- AM_PM_TempComp_FPKM[, grep("WT.*.27C.FPKM", colnames(AM_PM_TempComp_FPKM))]
rownames(Bmal_KO_37C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_37C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_KO_32C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_32C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_KO_27C) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_27C) <- AM_PM_TempComp_FPKM$gene_id

Bmal_KO_AM <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.AM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_AM <- AM_PM_TempComp_FPKM[, grep("WT.*.AM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_KO_PM <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.PM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_PM <- AM_PM_TempComp_FPKM[, grep("WT.*.PM.FPKM", colnames(AM_PM_TempComp_FPKM))]
rownames(Bmal_KO_AM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_AM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_KO_PM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_PM) <- AM_PM_TempComp_FPKM$gene_id

epidermal_WT_data <- GSE114943_Counts[, grep(".*WT._counts.csv.gz", colnames(GSE114943_Counts))]
epidermal_KO_data <- GSE114943_Counts[, grep(".*KO._counts.csv.gz", colnames(GSE114943_Counts))]
epidermal_RE_data <- GSE114943_Counts[, grep(".*RE._counts.csv.gz", colnames(GSE114943_Counts))]

download.file("ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz", "Mus_musculus.GRCm38.81.gtf.gz")
gunzip("Mus_musculus.GRCm38.81.gtf.gz", overwrite = T)

gtf_txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.81.gtf")
gene_list <- genes(gtf_txdb)
gene_list <- as.data.frame(gene_list)
head(gene_list)
rm(gtf_txdb)

    ## Confirm that gene_list rows and counts rows are in the same order:
    indFiles <- list()
    indFiles[["gene_list"]] <- gene_list
    indFiles[["epidermal_WT_data"]] <- epidermal_WT_data
    indFiles[["epidermal_KO_data"]] <- epidermal_KO_data
    indFiles[["epidermal_RE_data"]] <- epidermal_RE_data
    
    namesTemp <- lapply(indFiles,function(x) {as.character(row.names(x))})
    names <- do.call(cbind,namesTemp)
    
    testName <- apply(names,1,function(x) {sum(x == x[1]) })

    nNotMatchAll <- sum(testName != length(indFiles))
    nNotMatchAll


rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}


epidermal_WT_rpkm <- apply(epidermal_WT_data, 2, function(x) rpkm(x, gene_list$width))
epidermal_KO_rpkm <- apply(epidermal_KO_data, 2, function(x) rpkm(x, gene_list$width))
epidermal_RE_rpkm <- apply(epidermal_RE_data, 2, function(x) rpkm(x, gene_list$width))

rm(rpkm, gene_list, indFiles, namesTemp, names, testName, nNotMatchAll)

# ######################################################################################################
# check/view the data structure
# ######################################################################################################

colnames(exLiver_WT_data); head(exLiver_WT_data)
colnames(exLiver_KO_data); head(exLiver_KO_data)
colnames(Bmal_msf_WT_data); head(Bmal_msf_WT_data)
colnames(Bmal_msf_KO_data); head(Bmal_msf_KO_data)
colnames(Bmal_WT_37C); head(Bmal_WT_37C)
colnames(Bmal_WT_32C); head(Bmal_WT_32C)
colnames(Bmal_WT_27C); head(Bmal_WT_27C)
colnames(Bmal_KO_37C); head(Bmal_KO_37C)
colnames(Bmal_KO_32C); head(Bmal_KO_32C)
colnames(Bmal_KO_27C); head(Bmal_KO_27C)
colnames(Bmal_WT_AM); head(Bmal_WT_AM)
colnames(Bmal_WT_PM); head(Bmal_WT_PM)
colnames(Bmal_KO_AM); head(Bmal_KO_AM)
colnames(Bmal_KO_PM); head(Bmal_KO_PM)
colnames(epidermal_WT_rpkm); head(epidermal_WT_rpkm)
colnames(epidermal_KO_rpkm); head(epidermal_KO_rpkm)
colnames(epidermal_RE_rpkm); head(epidermal_RE_rpkm)

# ######################################################################################################
# data cleanup, remove rows with any zero values, Log Transformation
# ######################################################################################################

# data cleanup, remove rows with any zero values
exLiver_WT_data <- exLiver_WT_data[apply(exLiver_WT_data,1,function(y) !any(y==0)),]
exLiver_KO_data <- exLiver_KO_data[apply(exLiver_KO_data,1,function(y) !any(y==0)),]

Bmal_msf_WT_data <- Bmal_msf_WT_data[apply(Bmal_msf_WT_data,1,function(y) !any(y==0)),]
Bmal_msf_KO_data <- Bmal_msf_KO_data[apply(Bmal_msf_KO_data,1,function(y) !any(y==0)),]

Bmal_KO_37C <- Bmal_KO_37C[apply(Bmal_KO_37C,1,function(y) !any(y==0)),]
Bmal_WT_37C <- Bmal_WT_37C[apply(Bmal_WT_37C,1,function(y) !any(y==0)),]
Bmal_KO_32C <- Bmal_KO_32C[apply(Bmal_KO_32C,1,function(y) !any(y==0)),]
Bmal_WT_32C <- Bmal_WT_32C[apply(Bmal_WT_32C,1,function(y) !any(y==0)),]
Bmal_KO_27C <- Bmal_KO_27C[apply(Bmal_KO_27C,1,function(y) !any(y==0)),]
Bmal_WT_27C <- Bmal_WT_27C[apply(Bmal_WT_27C,1,function(y) !any(y==0)),]

Bmal_KO_AM <- Bmal_KO_AM[apply(Bmal_KO_AM,1,function(y) !any(y==0)),]
Bmal_WT_AM <- Bmal_WT_AM[apply(Bmal_WT_AM,1,function(y) !any(y==0)),]
Bmal_KO_PM <- Bmal_KO_PM[apply(Bmal_KO_PM,1,function(y) !any(y==0)),]
Bmal_WT_PM <- Bmal_WT_PM[apply(Bmal_WT_PM,1,function(y) !any(y==0)),]

epidermal_WT_rpkm <- epidermal_WT_rpkm[apply(epidermal_WT_rpkm,1,function(y) !any(y==0)),]
epidermal_KO_rpkm <- epidermal_KO_rpkm[apply(epidermal_KO_rpkm,1,function(y) !any(y==0)),]
epidermal_RE_rpkm <- epidermal_RE_rpkm[apply(epidermal_RE_rpkm,1,function(y) !any(y==0)),]


# Log2 Transformation
exLiver_WT_data_log <- log2(exLiver_WT_data)
exLiver_KO_data_log <- log2(exLiver_KO_data)

Bmal_msf_WT_data_log <- log2(Bmal_msf_WT_data)
Bmal_msf_KO_data_log <- log2(Bmal_msf_KO_data)

Bmal_KO_37C_log <- log2(Bmal_KO_37C)
Bmal_WT_37C_log <- log2(Bmal_WT_37C)
Bmal_KO_32C_log <- log2(Bmal_KO_32C)
Bmal_WT_32C_log <- log2(Bmal_WT_32C)
Bmal_KO_27C_log <- log2(Bmal_KO_27C)
Bmal_WT_27C_log <- log2(Bmal_WT_27C)

Bmal_KO_AM_log <- log2(Bmal_KO_AM)
Bmal_WT_AM_log <- log2(Bmal_WT_AM)
Bmal_KO_PM_log <- log2(Bmal_KO_PM)
Bmal_WT_PM_log <- log2(Bmal_WT_PM)

epidermal_WT_log <- log2(epidermal_WT_rpkm)
epidermal_KO_log <- log2(epidermal_KO_rpkm)
epidermal_RE_log <- log2(epidermal_RE_rpkm)

# ######################################################################################################
# Running RAIN
# ######################################################################################################
dir.create("output", showWarnings=FALSE)

# wrapper function for running RAIN with defaults settings
runRAIN <- function(dat, deltat, period.delta = 0) {
  rain(x = t(dat), deltat = deltat, period = 24, period.delta = period.delta, adjp.method = "ABH", peak.border = c(0.1,0.9)) 
}

# ex-vivo liver
WTexLiver <- runRAIN(dat = exLiver_WT_data_log, deltat = 3)
WTexLiver$p.adj <- round(p.adjust(WTexLiver$pVal,method='BH'), digit = 6)
KOexLiver <- runRAIN(dat = exLiver_KO_data_log, deltat = 3)
KOexLiver$p.adj <- round(p.adjust(KOexLiver$pVal,method='BH'), digit = 6)

# msf
WTmsf <- runRAIN(dat = Bmal_msf_WT_data_log, deltat = 3)
WTmsf$p.adj <- round(p.adjust(WTmsf$pVal,method='BH'), digit = 6)
KOmsf <- runRAIN(dat = Bmal_msf_KO_data_log, deltat = 3)
KOmsf$p.adj <- round(p.adjust(KOmsf$pVal,method='BH'), digit = 6)

# temperature compensation msf
WT27 <- runRAIN(dat = Bmal_WT_27C_log, deltat = 2)
WT32 <- runRAIN(dat = Bmal_WT_32C_log, deltat = 2)
WT37 <- runRAIN(dat = Bmal_WT_37C_log, deltat = 2)
KO27 <- runRAIN(dat = Bmal_KO_27C_log, deltat = 2)
KO32 <- runRAIN(dat = Bmal_KO_32C_log, deltat = 2)
KO37 <- runRAIN(dat = Bmal_KO_37C_log, deltat = 2)

WT27$p.adj <- round(p.adjust(WT27$pVal,method='BH'), digit = 6)
WT32$p.adj <- round(p.adjust(WT32$pVal,method='BH'), digit = 6)
WT37$p.adj <- round(p.adjust(WT37$pVal,method='BH'), digit = 6)
KO27$p.adj <- round(p.adjust(KO27$pVal,method='BH'), digit = 6)
KO32$p.adj <- round(p.adjust(KO32$pVal,method='BH'), digit = 6)
KO37$p.adj <- round(p.adjust(KO37$pVal,method='BH'), digit = 6)
      
# opposite entrainment msf
WTAM <- runRAIN(dat = Bmal_WT_AM_log, deltat = 3)
WTPM <- runRAIN(dat = Bmal_WT_PM_log, deltat = 3)
KOAM <- runRAIN(dat = Bmal_KO_AM_log, deltat = 3)
KOPM <- runRAIN(dat = Bmal_KO_PM_log, deltat = 3)

WTAM$p.adj <- round(p.adjust(WTAM$pVal,method='BH'), digit = 6)
WTPM$p.adj <- round(p.adjust(WTPM$pVal,method='BH'), digit = 6)
KOAM$p.adj <- round(p.adjust(KOAM$pVal,method='BH'), digit = 6)
KOPM$p.adj <- round(p.adjust(KOPM$pVal,method='BH'), digit = 6)


# wrapper function for running RAIN with defaults settings
runRAIN <- function(dat, deltat, period.delta = 0, replicates) {
  rain(x = t(dat), deltat = deltat, period = 24, period.delta = period.delta, nr.series = replicates, adjp.method = "ABH", peak.border = c(0.1,0.9)) 
}

# ex-vivo liver
epidermal_WT_log_rain_results <- runRAIN(dat = epidermal_WT_log, deltat = 4, replicates=4)
epidermal_WT_log_rain_results$p.adj <- round(p.adjust(epidermal_WT_log_rain_results$pVal,method='BH'), digit = 6)

epidermal_KO_log_rain_results <- runRAIN(dat = epidermal_KO_log, deltat = 4, replicates=4)
epidermal_KO_log_rain_results$p.adj <- round(p.adjust(epidermal_KO_log_rain_results$pVal,method='BH'), digit = 6)

epidermal_RE_log_rain_results <- runRAIN(dat = epidermal_RE_log, deltat = 4, replicates=4)
epidermal_RE_log_rain_results$p.adj <- round(p.adjust(epidermal_RE_log_rain_results$pVal,method='BH'), digit = 6)

# ################################################################################################################
# save all the output dataframes
# ################################################################################################################

current_ession <- sessionInfo()
save(current_ession, WTexLiver, KOexLiver, Bmal_msf_WT_data, Bmal_msf_KO_data,
    Bmal_WT_37C, Bmal_WT_32C, Bmal_WT_27C, Bmal_KO_37C, Bmal_KO_32C, Bmal_KO_27C,
    Bmal_WT_AM, Bmal_WT_PM, Bmal_KO_AM, Bmal_KO_PM,
    WTmsf, KOmsf, WT37, WT32, WT27, KO37, KO32, KO27, WTAM, WTPM, KOAM, KOPM,
    epidermal_WT_log_rain_results, epidermal_KO_log_rain_results, epidermal_RE_log_rain_results,
    file = "output/RAIN_results_Abruzzi_1&2.RData" )

# load('~/Downloads/output/RAIN_results_Abruzzi_1&2.RData')

# ########################################################
# Plots for Abruzzi Response et al.
# ########################################################

common <- unique(c(rownames(WTmsf),rownames(WT37),rownames(WT32),rownames(WT27),rownames(WTAM),rownames(WTPM),rownames(KOmsf),rownames(KO37),rownames(KO32),rownames(KO27),rownames(KOAM),rownames(KOPM)))

allP <- cbind(
  WT.exp1 = WTmsf[common, ]$pVal,
  WT.37C = WT37[common, ]$pVal,
  WT.32C = WT32[common, ]$pVal,
  WT.27C = WT27[common, ]$pVal,
  WT.AM = WTAM[common, ]$pVal,
  WT.PM = WTPM[common, ]$pVal,
  KO.exp1 = KOmsf[common, ]$pVal,
  KO.37C = KO37[common, ]$pVal,
  KO.32C = KO32[common, ]$pVal,
  KO.27C = KO27[common, ]$pVal,
  KO.AM = KOAM[common, ]$pVal,
  KO.PM = KOPM[common, ]$pVal
)

allFDR <- apply(allP, 2, p.adjust, "fdr")

# columns with the control condition 
contCols <- c( "WT.exp1", "WT.37C", "WT.AM", "WT.PM",
                "KO.exp1", "KO.37C", "KO.AM", "KO.PM" )

allFDRcont <- allFDR[, contCols]

# make list of sets to input into upsetR
FDRlist <- apply(allFDR, 2, function(q) { which(q < 0.1) })

# metadata for upsetR
metaData <- data.frame( sets = names(FDRlist),
                        experiment = rep(c(1, 2, 2, 2, 3, 3), 2),
                        WTKO = substr(names(FDRlist), 1, 2),
                        cond = rep(c("cont", "cont", "exp", "exp", "cont", "cont"), 2)   )

metaData$WTKO.cond <- paste(metaData$WTKO, metaData$cond, sep = ".")
rownames(metaData) <- metaData$sets

# ########################################################
# Figure 2A
# ########################################################

cont_KO_Cols <- c( "KO.exp1", "KO.37C", "KO.AM", "KO.PM"  )

allFDRcont_KO <- allFDR[, cont_KO_Cols]

pdf("output/Figure_2A_MSF_KO.pdf", height = 4, width = 4)
upset(fromList(FDRlist[cont_KO_Cols]), keep.order = T, nsets = 4, nintersects = 50, 
      sets = rev(cont_KO_Cols), order.by = "freq", sets.bar.color = c("red", "red", "red", "red"),
      set.metadata = list(data = metaData[contCols, ], plots = list(list(type = "matrix_rows", column = "WTKO.cond", alpha = 0.5))) )

dev.off()

# ########################################################
# Figure 2B
# ########################################################

temp_KO_Cols <- c( "KO.27C", "KO.32C", "KO.37C" )

allFDRtemp_KO <- allFDR[, temp_KO_Cols]

pdf("output/Figure_2B_MSF_KO.pdf", height = 4, width = 4)
upset(fromList(FDRlist[temp_KO_Cols]), keep.order = T, nsets = 3, nintersects = 50, 
      sets = rev(temp_KO_Cols), order.by = "freq", sets.bar.color = c("red", "red", "red"),
      set.metadata = list(data = metaData[contCols, ], plots = list(list(type = "matrix_rows", column = "WTKO.cond", alpha = 0.5)))  )

dev.off()

# ################################################################################################################
# Figure 1E    ### Fold Change 
# ################################################################################################################

source("https://raw.githubusercontent.com/ejanalysis/analyze.stuff/master/R/rowMaxs.R")
source("https://raw.githubusercontent.com/ejanalysis/analyze.stuff/master/R/rowMins.R")

plotFoldChange <- function(raw_data, subset_genes_list, condition_name="Fold Change", my_color="blue", bin_width=0.4) {
  my_subset <- data.frame(matrix(0, ncol=ncol(raw_data), nrow=length(subset_genes_list)))
  colnames(my_subset) <- colnames(raw_data)

  for(k in 1:length(subset_genes_list)){
    my_subset[k,] <- raw_data[which(rownames(raw_data) == subset_genes_list[k]),]
  }

  fold_change <- ( rowMaxs(my_subset)/rowMins(my_subset) )

  my_median <- median(fold_change)
  my_n_total <- length(fold_change)

  fold_change <- fold_change[fold_change<10]
  my_n_lower <- length(fold_change)
  my_n_upper <- my_n_total - my_n_lower

  bin_breakpoints <- seq((min((fold_change[fold_change>0])) ), (max(fold_change[fold_change>0])+bin_width), bin_width)

  hist(as.matrix(fold_change[fold_change>0]), breaks=bin_breakpoints, col= my_color, include.lowest=TRUE, las=1, xlim=c(0,(max(bin_breakpoints)+bin_width)), ylim=c(0,(max(hist(as.matrix(fold_change[fold_change>0]),breaks=bin_breakpoints,plot=FALSE)$counts)+(0.2*(max(hist(as.matrix(fold_change[fold_change>0]),breaks=bin_breakpoints,plot=FALSE)$counts))) ) ), xlab=("Fold Change Distribution\nRaw Data"), ylab=("Frequency"), main= paste(condition_name,"\nMedian Value: ",my_median,"\nTotal Number: ",my_n_total,"\nLess than ten: " ,my_n_lower,sep="") )

}


pdf(file=paste("output/Figure1E_Histogram_for_FoldChange","_exLiver",".pdf",sep=""),width = 10.69, height = 3.76* 1.5 )
par( mfrow = c(1, 2),oma = c(0,0,4,0), mar = c(4,4,6,3) )

plotFoldChange(raw_data=exLiver_WT_data, subset_genes_list=rownames(WTexLiver[which(WTexLiver$p.adj<0.1),]), condition_name="exLiver WT", my_color="blue", bin_width=0.4)

plotFoldChange(raw_data=exLiver_KO_data, subset_genes_list=rownames(KOexLiver[which(KOexLiver$p.adj<0.1),]), condition_name="exLiver KO", my_color="red", bin_width=0.4)

dev.off()



# ################################################################################################################
# Figure 2F&G    #### Comparing our data with GSE114943 DD Mouse skin data
# ################################################################################################################


WTexLiver_epidermal_rain_cyclic_genes_fdr10 <- list( "WTexLiver FDR<0.1" = rownames(WTexLiver[which(WTexLiver$p.adj<0.1),]), 
  "WT epidermal FDR<0.1" = rownames(epidermal_WT_log_rain_results[which(epidermal_WT_log_rain_results$p.adj<0.1),]) )
summary(WTexLiver_epidermal_rain_cyclic_genes_fdr10)

KOexLiver_epidermal_rain_cyclic_genes_fdr10 <- list( "KOexLiver FDR<0.1" = rownames(KOexLiver[which(KOexLiver$p.adj<0.1),]), 
  "KO epidermal FDR<0.1" = rownames(epidermal_KO_log_rain_results[which(epidermal_KO_log_rain_results$p.adj<0.1),]) )
summary(KOexLiver_epidermal_rain_cyclic_genes_fdr10)

v2 <- plot(venn(WTexLiver_epidermal_rain_cyclic_genes_fdr10))
v5 <- plot(venn(KOexLiver_epidermal_rain_cyclic_genes_fdr10))

WTmsf_epidermal_rain_cyclic_genes_fdr10 <- list( "WTmsf FDR<0.1" = rownames(WTmsf[which(WTmsf$p.adj<0.1),]), 
  "WT37 FDR<0.1" = rownames(WT37[which(WT37$p.adj<0.1),]), 
  "WTAM FDR<0.1" = rownames(WTAM[which(WTAM$p.adj<0.1),]),  
  "WTPM FDR<0.1" = rownames(WTPM[which(WTPM$p.adj<0.1),]),  
  "WT epidermal FDR<0.1" = rownames(epidermal_WT_log_rain_results[which(epidermal_WT_log_rain_results$p.adj<0.1),]) )
WTmsf_ours_fdr10 <- WTmsf_epidermal_rain_cyclic_genes_fdr10[-c(5)]
WTmsf_ours_fdr10_cyclic_genes <- unique(unlist(WTmsf_ours_fdr10[grep(".*.",WTmsf_ours_fdr10)]))

KOmsf_epidermal_rain_cyclic_genes_fdr10 <- list( "KOmsf FDR<0.1" = rownames(KOmsf[which(KOmsf$p.adj<0.1),]), 
  "KO37 FDR<0.1" = rownames(KO37[which(KO37$p.adj<0.1),]), 
  "KOAM FDR<0.1" = rownames(KOAM[which(KOAM$p.adj<0.1),]),  
  "KOPM FDR<0.1" = rownames(KOPM[which(KOPM$p.adj<0.1),]),  
  "KO epidermal FDR<0.1" = rownames(epidermal_KO_log_rain_results[which(epidermal_KO_log_rain_results$p.adj<0.1),]) )
KOmsf_ours_fdr10 <- KOmsf_epidermal_rain_cyclic_genes_fdr10[-c(5)]
KOmsf_ours_fdr10_cyclic_genes <- unique(unlist(KOmsf_ours_fdr10[grep(".*.",KOmsf_ours_fdr10)]))


WTmsfUnion_epidermal_rain_cyclic_genes_fdr10 <- list( "Union WT MSF FDR<0.1" = WTmsf_ours_fdr10_cyclic_genes, 
  "WT epidermal FDR<0.1" = rownames(epidermal_WT_log_rain_results[which(epidermal_WT_log_rain_results$p.adj<0.1),]) )

KOmsfUnion_epidermal_rain_cyclic_genes_fdr10 <- list( "Union KO MSF FDR<0.1" = KOmsf_ours_fdr10_cyclic_genes, 
  "KO epidermal FDR<0.1" = rownames(epidermal_KO_log_rain_results[which(epidermal_KO_log_rain_results$p.adj<0.1),]) )


# ### GeneOverlap calculations for Fisher's test

WT.genome.size = length(unique(c(rownames(WTmsf), rownames(WT37), rownames(WTAM), rownames(WTPM), rownames(epidermal_WT_log_rain_results) ) ))
WT.go.obj.fdr10 <- newGeneOverlap(unlist(WTmsfUnion_epidermal_rain_cyclic_genes_fdr10[1]), 
                         unlist(WTmsfUnion_epidermal_rain_cyclic_genes_fdr10[2]), 
                         genome.size=WT.genome.size )
WT.go.obj.fdr10 <- testGeneOverlap(WT.go.obj.fdr10)


KO.genome.size = length(unique(c(rownames(KOmsf), rownames(KO37), rownames(KOAM), rownames(KOPM), rownames(epidermal_KO_log_rain_results) ) ))
KO.go.obj.fdr10 <- newGeneOverlap(unlist(KOmsfUnion_epidermal_rain_cyclic_genes_fdr10[1]), 
                         unlist(KOmsfUnion_epidermal_rain_cyclic_genes_fdr10[2]), 
                         genome.size=KO.genome.size )
KO.go.obj.fdr10 <- testGeneOverlap(KO.go.obj.fdr10)


v20 <- plot(venn(WTmsfUnion_epidermal_rain_cyclic_genes_fdr10), main=paste("Overlapping p-value:",WT.go.obj.fdr10@pval,"\nUniverse:",WT.genome.size,sep="") )
v23 <- plot(venn(KOmsfUnion_epidermal_rain_cyclic_genes_fdr10), main=paste("Overlapping p-value:",KO.go.obj.fdr10@pval,"\nUniverse:",KO.genome.size,sep="") )


pdf("output/Figure2F+G_VennDiagrams.pdf", width=15, height=7 )
grid.arrange(v23,v5,ncol=2)
dev.off()


# ################################################################################################################



