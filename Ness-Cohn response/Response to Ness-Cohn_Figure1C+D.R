# ##########################################################################################################################
# Response to Technical Comment by Ness-Cohn et al. (2021). Code for:
# Figure 1C-D
# ##########################################################################################################################

# assumes we're standing in the R subdirectory
setwd("~/Downloads/")
dir.create("data", showWarnings=FALSE)
dir.create("output", showWarnings=FALSE)


rm(list = ls())

# should data be re-pulled from GEO? set F to use data from git repo (faster)
reGEO <- TRUE

# should RAIN be re-run? set F to use pre-computed results
reRAIN <- TRUE

# =====================================================================
# Dependencies
# NB: you may need to install some packages before running this script
# =====================================================================

library(R.utils)
library(rain)
library(gplots)
library(scales)
library(ggplot2)
library(UpSetR)
# library(cosinor) # not used in this script
library(plotrix)
library(tidyverse)

# =====================================================================
# Downloading & formatting the data
# =====================================================================

if (reGEO) { # re-download the data?
  setwd("data")
  # Experiment 1
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111696&format=file&file=GSE111696%5FBmal1%5FMSF%5Fgenes%2Efpkm%5Ftracking%2Etxt%2Egz", "GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt.gz")
  gunzip("GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt.gz", overwrite = T)
  # Experiments 2 & 3
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134333&format=file&file=GSE134333%5FDataset%5FAM%5FPM%5FTempComp%5FFPKM%5Fvalues%2Etxt%2Egz", "GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz")
  gunzip("GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt.gz", overwrite = T)
  setwd("..")
}

# reading the data from disk:
MSF_fpkm <- read.delim("data/GSE111696_Bmal1_MSF_genes.fpkm_tracking.txt", quote = "")
AM_PM_TempComp_FPKM <- read.delim("data/GSE134333_Dataset_AM_PM_TempComp_FPKM_values.txt", quote = "")

# reformatting data from Expt 1:
Bmal_msf_WT_data <- MSF_fpkm[, grep("WT_.*_FPKM", colnames(MSF_fpkm))]
Bmal_msf_KO_data <- MSF_fpkm[, grep("KO_.*_FPKM", colnames(MSF_fpkm))]
rownames(Bmal_msf_WT_data) <- MSF_fpkm$gene_id
rownames(Bmal_msf_KO_data) <- MSF_fpkm$gene_id
# making sure times are in order
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_msf_WT_data))) == seq(0, 69, by = 3))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_msf_KO_data))) == seq(0, 69, by = 3))

# reformatting data from Expt 2:
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
# making sure times are in order
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_WT_37C))) == seq(0, 46, by = 2))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_WT_32C))) == seq(0, 46, by = 2))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_WT_27C))) == seq(0, 46, by = 2))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_KO_37C))) == seq(0, 46, by = 2))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_KO_32C))) == seq(0, 46, by = 2))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_KO_27C))) == seq(0, 46, by = 2))

# reformatting data from Expt 3:
Bmal_KO_AM <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.AM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_AM <- AM_PM_TempComp_FPKM[, grep("WT.*.AM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_KO_PM <- AM_PM_TempComp_FPKM[, grep("Bmal1.*.PM.FPKM", colnames(AM_PM_TempComp_FPKM))]
Bmal_WT_PM <- AM_PM_TempComp_FPKM[, grep("WT.*.PM.FPKM", colnames(AM_PM_TempComp_FPKM))]
rownames(Bmal_KO_AM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_AM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_KO_PM) <- AM_PM_TempComp_FPKM$gene_id
rownames(Bmal_WT_PM) <- AM_PM_TempComp_FPKM$gene_id
# making sure times are in order
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_WT_AM))) == seq(0, 48, by = 3))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_WT_PM))) == seq(0, 48, by = 3))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_KO_AM))) == seq(0, 48, by = 3))
all(as.numeric(gsub(".*CT(..).*", "\\1", colnames(Bmal_KO_PM))) == seq(0, 48, by = 3))

# =====================================================================================
# data cleanup, remove rows with any zero values and Log Transformation
# =====================================================================================

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


Bmal_msf_WT_data <- log2(Bmal_msf_WT_data)
Bmal_msf_KO_data <- log2(Bmal_msf_KO_data)

Bmal_KO_37C <- log2(Bmal_KO_37C)
Bmal_WT_37C <- log2(Bmal_WT_37C)
Bmal_KO_32C <- log2(Bmal_KO_32C)
Bmal_WT_32C <- log2(Bmal_WT_32C)
Bmal_KO_27C <- log2(Bmal_KO_27C)
Bmal_WT_27C <- log2(Bmal_WT_27C)

Bmal_KO_AM <- log2(Bmal_KO_AM)
Bmal_WT_AM <- log2(Bmal_WT_AM)
Bmal_KO_PM <- log2(Bmal_KO_PM)
Bmal_WT_PM <- log2(Bmal_WT_PM)

# =====================================================================
# Running RAIN
# =====================================================================


if (reRAIN) { # NB: this takes some time.  Canned results are in output/

  # wrapper function for running RAIN with defaults as in paper
  runRAIN <- function(x, deltat, period.delta = 0) {
    rain(x = t(x), deltat = deltat, period = 24, period.delta = period.delta, adjp.method = "ABH", peak.border = c(0.1,0.9))
  }

  # Exp 1
  WT1 <- runRAIN(Bmal_msf_WT_data, deltat = 3)
  KO1 <- runRAIN(Bmal_msf_KO_data, deltat = 3)

  # Exp 2
  WT27 <- runRAIN(Bmal_WT_27C, deltat = 2)
  WT32 <- runRAIN(Bmal_WT_32C, deltat = 2)
  WT37 <- runRAIN(Bmal_WT_37C, deltat = 2)
  KO27 <- runRAIN(Bmal_KO_27C, deltat = 2)
  KO32 <- runRAIN(Bmal_KO_32C, deltat = 2)
  KO37 <- runRAIN(Bmal_KO_37C, deltat = 2)

  # Exp 3
  WTAM <- runRAIN(Bmal_WT_AM, deltat = 3)
  WTPM <- runRAIN(Bmal_WT_PM, deltat = 3)
  KOAM <- runRAIN(Bmal_KO_AM, deltat = 3)
  KOPM <- runRAIN(Bmal_KO_PM, deltat = 3)

  # save it all
  runSession <- sessionInfo()
  save(runSession, Bmal_msf_WT_data, Bmal_msf_KO_data,
    Bmal_WT_37C, Bmal_WT_32C, Bmal_WT_27C, Bmal_WT_AM, Bmal_WT_PM,
    Bmal_KO_37C, Bmal_KO_32C, Bmal_KO_27C, Bmal_KO_AM, Bmal_KO_PM,
    WT1, WT37, WT32, WT27, WTAM, WTPM,
    KO1, KO37, KO32, KO27, KOAM, KOPM,
    file = "output/RAINres.RData"
  )
} else {
  load("output/RAINres.RData")
}

# =====================================================================
# analysis functions
# =====================================================================

timeDiff <- function(trueTimes, predTimes, ...) {
  predTimes <- predTimes %% 24
  trueTimes <- trueTimes %% 24
  timeDiffs <- abs(predTimes - trueTimes)
  timeDiffs <- pmin(timeDiffs, 24 - timeDiffs)
  return(timeDiffs)
}

# =====================================================================
# Figures and table data for tech comment response
# =====================================================================

# Get common genes among the datasets

# common <- intersect(rownames(WT1), rownames(WT37))            # Ness-Cohen common method

# get union of all gene names in all samples (removing duplicates), instead of using intersect of WT1+WT37 (Ness-Cohen analysis)

common <- unique(c(rownames(WT1),rownames(WT37),rownames(WT32),rownames(WT27),rownames(WTAM),rownames(WTPM),rownames(KO1),rownames(KO37),rownames(KO32),rownames(KO27),rownames(KOAM),rownames(KOPM)))

allP <- cbind(
  WT1 = WT1[common, ]$pVal,
  WT37 = WT37[common, ]$pVal,
  WT32 = WT32[common, ]$pVal,
  WT27 = WT27[common, ]$pVal,
  WTAM = WTAM[common, ]$pVal,
  WTPM = WTPM[common, ]$pVal,
  KO1 = KO1[common, ]$pVal,
  KO37 = KO37[common, ]$pVal,
  KO32 = KO32[common, ]$pVal,
  KO27 = KO27[common, ]$pVal,
  KOAM = KOAM[common, ]$pVal,
  KOPM = KOPM[common, ]$pVal
)

allFDR <- apply(allP, 2, p.adjust, "fdr")
colnames(allFDR) <- c(
  "WT.exp1", "WT.37C", "WT.32C", "WT.27C", "WT.AM", "WT.PM",
  "KO.exp1", "KO.37C", "KO.32C", "KO.27C", "KO.AM", "KO.PM"
)
# columns with the control condition [37C, 48h post-DEX] (* added WT PM set in Response)
contCols <- c(
  "WT.exp1", "WT.37C", "WT.AM", "WT.PM",
  "KO.exp1", "KO.37C", "KO.AM", "KO.PM"
)
allFDRcont <- allFDR[, contCols]

# make list of sets to input into upsetR

FDRlist <- apply(allFDR, 2, function(q) {
  which(q < 0.1)
})

# metadata for upsetR
metaData <- data.frame(
  sets = names(FDRlist),
  experiment = rep(c(1, 2, 2, 2, 3, 3), 2),
  WTKO = substr(names(FDRlist), 1, 2),
  cond = rep(c("cont", "cont", "exp", "exp", "cont", "cont"), 2)
)
metaData$WTKO.cond <- paste(metaData$WTKO, metaData$cond, sep = ".")
rownames(metaData) <- metaData$sets

# Figure 1C: WT sets only   ---------------------------------------------------------
cont_WT_Cols <- c(
  "WT.exp1", "WT.37C", "WT.AM", "WT.PM"
)
allFDRcont_WT <- allFDR[, cont_WT_Cols]

pdf("output/Figure1C_MSF_WT_only.pdf", height = 4, width = 4)
upset(fromList(FDRlist[cont_WT_Cols]),
      keep.order = T, nsets = 4, nintersects = 50, sets = rev(cont_WT_Cols), order.by = "freq",
      sets.bar.color = c("blue", "blue", "blue", "blue"),
      set.metadata = list(data = metaData[contCols, ], plots = list(list(type = "matrix_rows", column = "WTKO.cond", alpha = 0.5)))
)
dev.off()

# Figure 1C: KO sets only    ---------------------------------------------------------
cont_KO_Cols <- c(
  "KO.exp1", "KO.37C", "KO.AM", "KO.PM"
)
allFDRcont_KO <- allFDR[, cont_KO_Cols]

pdf("output/Figure1C_MSF_KO_only.pdf", height = 4, width = 4)
upset(fromList(FDRlist[cont_KO_Cols]),
      keep.order = T, nsets = 4, nintersects = 50, sets = rev(cont_KO_Cols), order.by = "freq",
      sets.bar.color = c("red", "red", "red", "red"),
      set.metadata = list(data = metaData[contCols, ], plots = list(list(type = "matrix_rows", column = "WTKO.cond", alpha = 0.5)))
)
dev.off()


# Figure 1D -------------------------------------------------------------

# RAIN phase estimates
APp.RAIN <- cbind(
  WT_AM = WTAM$phase,
  WT_PM = WTPM$phase,
  KO_AM = KOAM$phase,
  KO_PM = KOPM$phase
)

# RAIN fdr values
APf <- apply(cbind(
  WT_AM = WTAM$pVal,
  WT_PM = WTPM$pVal,
  KO_AM = KOAM$pVal,
  KO_PM = KOPM$pVal
), 2, p.adjust, method = "fdr")

# RAIN ranks
APr <- apply(cbind(
  WT_AM = WTAM$pVal,
  WT_PM = WTPM$pVal,
  KO_AM = KOAM$pVal,
  KO_PM = KOPM$pVal
), 2, rank, ties = "min")

APp <- APp.RAIN

# genes passing FDR < 0.1
WTf01 <- (APf[, 1] < 0.1) & (APf[, 2] < 0.1)
KOf01 <- (APf[, 3] < 0.1) & (APf[, 4] < 0.1)
table(round(timeDiff(APp[WTf01, 1], APp[WTf01, 2])))
table(round(timeDiff(APp[KOf01, 3], APp[KOf01, 4])))

