# ##########################################################################################################################
# Ray S#, Valekunja UK#, Stangherlin A, Howell SA, Snijders AP, Damodaran G, Reddy AB. (2020). 
# Circadian rhythms in the absence of the clock gene Bmal1. Science. 2020 Feb 14;367(6479):800-806. 
# doi: 10.1126/science.aaw7365.
# ##########################################################################################################################
rm(list = ls(all = TRUE))
# ######################################################################################################
# set the working directory
setwd("~/Downloads/")
# ######################################################################################################

# ######################################################################################################
# load all the needed packages
# ######################################################################################################

library(R.utils)
library(tools)
library(tidyr)
library(dplyr)
library(reshape2)
library(gnm)
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

source("~/Downloads/JTK_CYCLEv3.1.R") # downloaded from https://openwetware.org/wiki/HughesLab:JTK_Cycle
source("https://raw.githubusercontent.com/gangwug/MetaCycle/master/R/ARS.R")


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
setwd("..")

# ######################################################################################################
# reading the unpacked data 
# ######################################################################################################

exLiver_fpkm <- read.delim("data/GSE111696_Bmal1_exLiver_genes.fpkm_tracking.txt", quote = "")
MSF_fpkm <- read.delim("data/GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt", quote = "")
AM_PM_TempComp_FPKM <- read.delim("data/GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt", quote = "")


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

# ######################################################################################################
# Running JTK
# ######################################################################################################
# wrapper function for running JTK 
runJTK <- function(dat, sampling.time.sequence, sampling.resolution, replicates, period.low, period.high){

  jtkdist(timepoints=length(sampling.time.sequence),reps=replicates)
  periods <- (period.low/sampling.resolution):(period.high/sampling.resolution)
  jtk.init(periods, sampling.resolution)
  
  res <- apply(dat,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  
  res = as.data.frame(t(res))
  bhq = p.adjust(unlist(res[,1]),"BH")
  res = cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- res
} 

exLiver_WT_jtk_results <- runJTK(dat=exLiver_WT_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
exLiver_KO_jtk_results <- runJTK(dat=exLiver_KO_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)

WTmsf_jtk_results <- runJTK(dat = Bmal_msf_WT_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
KOmsf_jtk_results <- runJTK(dat = Bmal_msf_KO_data_log, sampling.time.sequence=seq(0,69,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)

WT27_jtk_results <- runJTK(dat = Bmal_WT_27C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)
WT32_jtk_results <- runJTK(dat = Bmal_WT_32C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)
WT37_jtk_results <- runJTK(dat = Bmal_WT_37C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)
KO27_jtk_results <- runJTK(dat = Bmal_KO_27C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)
KO32_jtk_results <- runJTK(dat = Bmal_KO_32C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)
KO37_jtk_results <- runJTK(dat = Bmal_KO_37C_log, sampling.time.sequence=seq(0,46,by=2), sampling.resolution=2, replicates=1, period.low=24, period.high=24)

WTAM_jtk_results <- runJTK(dat = Bmal_WT_AM_log, sampling.time.sequence=seq(0,48,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
WTPM_jtk_results <- runJTK(dat = Bmal_WT_PM_log, sampling.time.sequence=seq(0,48,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
KOAM_jtk_results <- runJTK(dat = Bmal_KO_AM_log, sampling.time.sequence=seq(0,48,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)
KOPM_jtk_results <- runJTK(dat = Bmal_KO_PM_log, sampling.time.sequence=seq(0,48,by=3), sampling.resolution=3, replicates=1, period.low=24, period.high=24)


# ######################################################################################################
# Running ARS
# ######################################################################################################


exLiver_WT_ars_results <- runARS(cbind(row.names(exLiver_WT_data_log),data.frame(exLiver_WT_data_log)), ARStime=rep(seq(0,69,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
exLiver_KO_ars_results <- runARS(cbind(row.names(exLiver_KO_data_log),data.frame(exLiver_KO_data_log)), ARStime=rep(seq(0,69,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)

WTmsf_ars_results <- runARS(cbind(row.names(Bmal_msf_WT_data_log),data.frame(Bmal_msf_WT_data_log)), ARStime=rep(seq(0,69,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KOmsf_ars_results <- runARS(cbind(row.names(Bmal_msf_KO_data_log),data.frame(Bmal_msf_KO_data_log)), ARStime=rep(seq(0,69,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)

WT27_ars_results <- runARS(cbind(row.names(Bmal_WT_27C_log),data.frame(Bmal_WT_27C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
WT32_ars_results <- runARS(cbind(row.names(Bmal_WT_32C_log),data.frame(Bmal_WT_32C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
WT37_ars_results <- runARS(cbind(row.names(Bmal_WT_37C_log),data.frame(Bmal_WT_37C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KO27_ars_results <- runARS(cbind(row.names(Bmal_KO_27C_log),data.frame(Bmal_KO_27C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KO32_ars_results <- runARS(cbind(row.names(Bmal_KO_32C_log),data.frame(Bmal_KO_32C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KO37_ars_results <- runARS(cbind(row.names(Bmal_KO_37C_log),data.frame(Bmal_KO_37C_log)), ARStime=rep(seq(0,46,by=2), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)

WTAM_ars_results <- runARS(cbind(row.names(Bmal_WT_AM_log),data.frame(Bmal_WT_AM_log)), ARStime=rep(seq(0,48,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
WTPM_ars_results <- runARS(cbind(row.names(Bmal_WT_PM_log),data.frame(Bmal_WT_PM_log)), ARStime=rep(seq(0,48,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KOAM_ars_results <- runARS(cbind(row.names(Bmal_KO_AM_log),data.frame(Bmal_KO_AM_log)), ARStime=rep(seq(0,48,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)
KOPM_ars_results <- runARS(cbind(row.names(Bmal_KO_PM_log),data.frame(Bmal_KO_PM_log)), ARStime=rep(seq(0,48,by=3), each=1), minper=24, maxper=24, arsper=24, arsmet="auto", releaseNote=TRUE)


# ######################################################################################################
# save all the output dataframes
# ######################################################################################################

current_session <- sessionInfo()
save(current_session, WTexLiver, KOexLiver, Bmal_msf_WT_data, Bmal_msf_KO_data,
    Bmal_WT_37C, Bmal_WT_32C, Bmal_WT_27C, Bmal_KO_37C, Bmal_KO_32C, Bmal_KO_27C,
    Bmal_WT_AM, Bmal_WT_PM, Bmal_KO_AM, Bmal_KO_PM,
    WTmsf, KOmsf, WT37, WT32, WT27, KO37, KO32, KO27, WTAM, WTPM, KOAM, KOPM,
    exLiver_WT_jtk_results, exLiver_KO_jtk_results, WTmsf_jtk_results, KOmsf_jtk_results, 
    WT27_jtk_results, WT32_jtk_results, WT37_jtk_results, 
    KO27_jtk_results, KO32_jtk_results, KO37_jtk_results,
    WTAM_jtk_results, WTPM_jtk_results, KOAM_jtk_results, KOPM_jtk_results,
    exLiver_WT_ars_results, exLiver_KO_ars_results, WTmsf_ars_results, KOmsf_ars_results, 
    WT27_ars_results, WT32_ars_results, WT37_ars_results, 
    KO27_ars_results, KO32_ars_results, KO37_ars_results,
    WTAM_ars_results, WTPM_ars_results, KOAM_ars_results, KOPM_ars_results,
    file = "output/out_bmal1_results.RData" )

# load('~/Downloads/output/out_bmal1_results.RData')
# ######################################################################################################
# list of cyclic genes in each condition
# ######################################################################################################

WTexLiver_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WTexLiver[which(WTexLiver$p.adj<0.05),]),
  "FDR<0.1" = rownames(WTexLiver[which(WTexLiver$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WTexLiver[which(WTexLiver$pVal<0.05),]) )
summary(WTexLiver_cyclic_genes_list)

KOexLiver_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KOexLiver[which(KOexLiver$p.adj<0.05),]),
  "FDR<0.1" = rownames(KOexLiver[which(KOexLiver$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KOexLiver[which(KOexLiver$pVal<0.05),]) )
summary(KOexLiver_cyclic_genes_list)


WTmsf_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WTmsf[which(WTmsf$p.adj<0.05),]),
  "FDR<0.1" = rownames(WTmsf[which(WTmsf$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WTmsf[which(WTmsf$pVal<0.05),]) )
summary(WTmsf_cyclic_genes_list)

KOmsf_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KOmsf[which(KOmsf$p.adj<0.05),]),
  "FDR<0.1" = rownames(KOmsf[which(KOmsf$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KOmsf[which(KOmsf$pVal<0.05),]) )
summary(KOmsf_cyclic_genes_list)


WT27_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WT27[which(WT27$p.adj<0.05),]),
  "FDR<0.1" = rownames(WT27[which(WT27$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WT27[which(WT27$pVal<0.05),]) )
summary(WT27_cyclic_genes_list)

KO27_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KO27[which(KO27$p.adj<0.05),]),
  "FDR<0.1" = rownames(KO27[which(KO27$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KO27[which(KO27$pVal<0.05),]) )
summary(KO27_cyclic_genes_list)

WT32_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WT32[which(WT32$p.adj<0.05),]),
  "FDR<0.1" = rownames(WT32[which(WT32$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WT32[which(WT32$pVal<0.05),]) )
summary(WT32_cyclic_genes_list)

KO32_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KO32[which(KO32$p.adj<0.05),]),
  "FDR<0.1" = rownames(KO32[which(KO32$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KO32[which(KO32$pVal<0.05),]) )
summary(KO32_cyclic_genes_list)

WT37_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WT37[which(WT37$p.adj<0.05),]),
  "FDR<0.1" = rownames(WT37[which(WT37$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WT37[which(WT37$pVal<0.05),]) )
summary(WT37_cyclic_genes_list)

KO37_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KO37[which(KO37$p.adj<0.05),]),
  "FDR<0.1" = rownames(KO37[which(KO37$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KO37[which(KO37$pVal<0.05),]) )
summary(KO37_cyclic_genes_list)


WTAM_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WTAM[which(WTAM$p.adj<0.05),]),
  "FDR<0.1" = rownames(WTAM[which(WTAM$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WTAM[which(WTAM$pVal<0.05),]) )
summary(WTAM_cyclic_genes_list)

KOAM_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KOAM[which(KOAM$p.adj<0.05),]),
  "FDR<0.1" = rownames(KOAM[which(KOAM$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KOAM[which(KOAM$pVal<0.05),]) )
summary(KOAM_cyclic_genes_list)

WTPM_cyclic_genes_list <- list(  "FDR<0.05" = rownames(WTPM[which(WTPM$p.adj<0.05),]),
  "FDR<0.1" = rownames(WTPM[which(WTPM$p.adj<0.1),]), 
  "pVal<0.05" = rownames(WTPM[which(WTPM$pVal<0.05),]) )
summary(WTPM_cyclic_genes_list)

KOPM_cyclic_genes_list <- list(  "FDR<0.05" = rownames(KOPM[which(KOPM$p.adj<0.05),]),
  "FDR<0.1" = rownames(KOPM[which(KOPM$p.adj<0.1),]), 
  "pVal<0.05" = rownames(KOPM[which(KOPM$pVal<0.05),]) )
summary(KOPM_cyclic_genes_list)


# ######################################################################################################
# plots
# ######################################################################################################




pdf(file=paste("output/Figure1C_a_WTliver_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_b_KOliver_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_c_KOliver_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_d_WTliver_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)




pdf(file=paste("output/Figure1F_a_WTmsf_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_b_KOmsf_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_c_KOmsf_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_d_WTmsf_heatmap_pVal_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$pVal<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)





pdf(file=paste("output/Figure1C_a_WTliver_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_b_KOliver_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_c_KOliver_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_d_WTliver_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)




pdf(file=paste("output/Figure1F_a_WTmsf_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_b_KOmsf_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_c_KOmsf_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_d_WTmsf_heatmap_FDR_0.05.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$p.adj<0.05),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)





pdf(file=paste("output/Figure1C_a_WTliver_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_b_KOliver_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTexLiver[which(WTexLiver$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_c_KOliver_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_KO_data_log[current_cyclic_genes,]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1C_d_WTliver_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOexLiver[which(KOexLiver$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- exLiver_WT_data_log[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO exLiver  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)




pdf(file=paste("output/Figure1F_a_WTmsf_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_b_KOmsf_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- WTmsf[which(WTmsf$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in WT MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_c_KOmsf_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_KO_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  -/-",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)

pdf(file=paste("output/Figure1F_d_WTmsf_heatmap_FDR_0.10.pdf",sep=""), width=6, height=12 )
current_cyclic_genes <- KOmsf[which(KOmsf$p.adj<0.1),]
current_cyclic_genes <- rownames(current_cyclic_genes[order(current_cyclic_genes$phase,current_cyclic_genes$pVal), ])
current_cyclic <- Bmal_msf_WT_data[current_cyclic_genes,]
current_cyclic <- current_cyclic[complete.cases(current_cyclic), ]
heatmap3(as.matrix(current_cyclic ), useRaster=TRUE, Rowv=NA, Colv=NA, scale="row", hclustfun=function(d) hclust(d, method="ward"), col=colorRampPalette(c("blue","black","red"))(21),
    margins = c(6, 12), labCol= colnames(current_cyclic), labRow=rownames(current_cyclic), xlab=paste("Rhythmic in KO MSF  +/+",sep=""), ylab="", balanceColor=TRUE, showRowDendro=FALSE, keep.dendro=FALSE, symm=F)
dev.off()
rm(current_cyclic_genes)
rm(current_cyclic)




# Chart
pdf(file=paste("output/Figure1B_1_venn.pdf",sep="") )
plot(venn(WTexLiver_cyclic_genes_list))
dev.off()







